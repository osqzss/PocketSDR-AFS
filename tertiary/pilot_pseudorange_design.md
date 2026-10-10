# AFS-Q Pilot Channel Pseudorange - Design Study

Status: design only (not implemented yet).

This document studies the source code changes needed to generate pseudoranges
from the data-less AFS pilot channel (AFSP) by using the tertiary code sync
implemented in `tertiary_code_changes.md`.

## 1. Principle

A pseudorange needs the signal transmission time `t_tx` at the receiver epoch:

```
P = c * (t_rx - t_tx)
```

For AFSP, there is no navigation data, but the tertiary code gives the
transmission time modulo 12 s:

- The tertiary code (1500 chips x 8 ms = 12 s) starts at the start of the AFS-I
  data frame (LSIS-220).
- AFS frames are aligned to the AFS time of week: TOW = ITOW * 1200 + TOI * 12,
  and 604800 / 12 = 50400 is an integer, so a frame start is always a multiple
  of 12 s in the week.
- Therefore, after the tertiary code sync, `t_tx mod 12 s` is known with the
  same resolution as the AFSD time tag.

The remaining 12 s ambiguity can be resolved with the receiver time, which is
already known by the PVT epoch (`pvt->time`), because the signal travel time
plus the receiver clock error is far below 6 s for lunar users (LNSP ranges are
about 5000 - 10000 km, i.e. 17 - 35 ms; even the Earth-Moon distance is 1.3 s).
This is the same idea as the existing 100 ms ambiguity resolution for
L1CP/E1C (`tow_v == 2` in `gen_prng()`), with a 12 s period instead.

Benefits:

- Pseudoranges are available about 2 s after the secondary code sync (tertiary
  sync), without waiting for the AFSD frame sync (2 frames = 24 s).
- AFSP is BPSK(5) with 5 times the chip rate of AFSD, so its code tracking
  noise and multipath error are smaller than those of AFSD.
- AFSP can be tracked and used at a lower C/N0 than the AFSD data demodulation
  needs.

## 2. Current Implementation (branch tertiary)

### Time tag of a channel

- `ch->tow` (int, ms): transmission time at the end of the current epoch
  (primary code period). `ch->week`: week number (0: unknown).
- `ch->tow_v`: 0: invalid, 1: valid, 2: ambiguity unresolved (100 ms).
- `sdr_ch.c` `update_tow()` advances `ch->tow` by `T` every epoch (and by
  -T/+T on code offset wrap), modulo one week.
- AFSD (`sdr_nav.c` `decode_AFSD_frame()`): sets `week` and
  `tow = ITOW * 1200 + TOI * 12 + 0.136` (68 sync symbols x 2 ms after the
  frame start) when SB2 is decoded.
- AFSP (`sdr_nav.c` `decode_AFSP()`): empty, so `ch->tow = -1` and AFSP never
  produces observations.

### Pseudorange generation (src/sdr_pvt.c)

- `sdr_pvt_udobs()`: initializes the PVT epoch by `init_epoch()` (needs
  `ch->week > 0`), then calls `update_obs()` if
  `state == LOCK && tow >= 0 && tow_v > 0 && (fsync > 0 || sec_sync > 0)`.
- `gen_prng()`:
  - `week > 0`: `tau = (week diff) + tow_rx - tow * 1e-3 + coff`
  - `tow_v == 2`: same without week, then the 100 ms ambiguity is resolved
    (`0.05 <= tau < 0.15`).
- `update_obs()`: `data_idx()` puts AFSD (`CODE_L9D`) at index 0 and AFSP
  (`CODE_L9P`) at index 1. LLI bit 2 (half-cycle ambiguity) is cleared when
  `fsync > 0 || sec_sync > 0`.
- `sdr_pvt_udsol()`: uses `update_sol_afs()` only if some observation has
  `code[0] == CODE_L9D`.

### AFS PVT (src/sdr_pvt_afs.c)

- `update_sol_afs()`: for each satellite, uses the first non-empty code in
  `data->code[]`, i.e. AFSD if available, otherwise AFSP. Only satellites with
  a valid ephemeris (`eph[prn-1].vflg`, decoded from AFSD SB2) are used.

### Baseline test result

In the 30 s simulation (afs_sim, PRN 2-8), AFSD SB2 was first decoded at
about 24 s, so no `$OBS` / `$POS` was output. A run of 60 - 90 s is needed for
the verification.

## 3. Design Options

### Option A: tertiary code time + 12 s ambiguity resolution (recommended)

AFSP sets `ch->tow` = transmission time modulo 12 s from the tertiary code
sync, with a new flag `tow_v = 3` (12 s ambiguity unresolved). `gen_prng()`
resolves the 12 s ambiguity with the receiver epoch time.

- Pros: independent of other channels, no cross-thread access, same pattern
  as L1CP/E1C, small changes.
- Cons: PVT epoch initialization still needs at least one AFSD channel with a
  decoded week/TOW (also needed for the ephemeris anyway).

### Option B: copy week/TOW from the AFSD channel of the same PRN

AFSP gets `week` and `tow` from the AFSD channel of the same PRN (both are
coherent within 0.5 ns, LSIS-180).

- Pros: full time tag without ambiguity.
- Cons: AFSD and AFSP run in different threads with different `lock` and
  `time`, so mapping the epochs needs locking and care. It also needs the AFSD
  frame sync of the same PRN, which removes most of the benefits in section 1.

Option B is better used only as a consistency check (see 4.7).

## 4. Required Changes (Option A)

### 4.1 src/pocket_sdr.h

- Update the comment of `tow_v`:
  `// TOW flag (0:invalid,1:valid,2:amb-unresolved,3:12s-amb-unresolved)`
- Optional: `#define SDR_TOW_AMB12 3` or similar, to avoid a magic number.

### 4.2 src/sdr_nav.c: decode_AFSP()

Set the TOW modulo 12 s every epoch from the tertiary code sync:

```c
// decode AFSP nav data --------------------------------------------------------
static void decode_AFSP(sdr_ch_t *ch)
{
    if (ch->trk->ter_sync <= 0) {
        ch->tow = -1;
        ch->tow_v = 0;
        return;
    }
    // transmission time within the 12 s tertiary code period at the end of
    // the current epoch (ms) (LSIS-220: tertiary code start = frame start)
    int L = ch->len_sec_code * ch->len_ter_code; // 6000 epochs
    int n = ((ch->lock - ch->trk->ter_sync) % L + L) % L;
    ch->tow = (int)(n * ch->T / 1e-3 + 0.5); // 0 - 11998 ms
    ch->tow_v = 3; // 12 s ambiguity unresolved
}
```

Notes:

- `ter_sync` is the lock count at which a tertiary code period starts, i.e.
  the end of the epoch at `lock == ter_sync` is a frame start (`tow = 0`).
  `ter_sync` may be ahead of `lock` just after the sync, so a positive modulo
  is used.
- Recomputing every epoch (instead of only at the frame start like L1CP)
  keeps `tow` consistent after `lock--` / `lock++` on code offset wrap, and
  avoids a negative `tow` from `update_tow(ch, -T)` at `tow = 0`.
- `decode_AFSP()` is called only when `lock * T >= T_NPULLIN`, after
  `sync_ter_code()` in the same epoch, so `ter_sync` is up to date.
- When the tertiary sync is lost (`ter_sync = 0`), `tow` becomes invalid and
  the observation stops automatically.

### 4.3 src/sdr_pvt.c: gen_prng()

Add the 12 s ambiguity resolution:

```c
    else if (ch->tow_v == 3) { // resolve 12 s ambiguity (-6 <= tau < 6)
        tau = tow - ch->tow * 1e-3 + ch->coff;
        tau -= floor(tau / 12.0 + 0.5) * 12.0;
    }
```

- `tow` is the GPST TOW of the receiver epoch (`pvt->time`). AFS time is
  aligned with GPST in the simulator (`gpst2afst()` only regroups the TOW).
- A negative `tau` gives `P <= 0`, which is already rejected by
  `update_obs()` (`if (P > 0.0 && ...)`).
- The window `[-6, 6)` s is valid as long as |travel time + receiver clock
  error| < 6 s. The PVT epoch is initialized from AFSD (`init_epoch()`), so
  the clock error is small.
- Alternatively, use a window starting slightly below 0 (e.g. `[-1, 11)` s)
  if negative clock biases larger than the travel time are expected. To be
  decided at implementation.

### 4.4 src/sdr_pvt.c: sdr_pvt_udobs()

The current condition `(ch->nav->fsync > 0 || ch->trk->sec_sync > 0)` is
already true for AFSP after the secondary code sync. With 4.2, `tow >= 0 &&
tow_v > 0` holds only after the tertiary sync, so no change is strictly
needed. For clarity and robustness, the condition can be extended:

```c
            (ch->nav->fsync > 0 || ch->trk->sec_sync > 0) &&
            (ch->len_ter_code < 2 || ch->trk->ter_sync > 0)) {
```

### 4.5 src/sdr_pvt.c: update_obs() (LLI half-cycle flag)

For AFSP, `sec_pol` contains the sign of the tertiary chip at the moment of
the secondary code sync, so the half-cycle ambiguity of the carrier phase is
resolved only after the tertiary sync. Change the LLI bit 2 condition:

```c
        int sync = ch->nav->fsync > 0 ||
            (ch->trk->sec_sync > 0 &&
             (ch->len_ter_code < 2 || ch->trk->ter_sync > 0));
        if (!sync) {
            obs->data[i].LLI[j] |= 2; // half-cyc-amb unresolved
        }
```

Optional (separate step): correct the carrier phase by the resolved polarity
`sec_pol * ter_pol` (e.g. `L = -adr + (pol < 0 ? 0.5 : 0.0)`). The sign
convention (`CHIP[]` maps bit 1 to +1 in both secondary and tertiary codes)
must be checked against the simulator before enabling it.

### 4.6 src/sdr_pvt.c: sdr_pvt_udsol() (AFS PVT selection)

Use the AFS PVT also when only AFSP observations exist:

```c
            if (pvt->obs->data[i].code[0] == CODE_L9D ||
                pvt->obs->data[i].code[1] == CODE_L9P) {
```

### 4.7 src/sdr_pvt_afs.c: update_sol_afs() (observation selection)

Currently the first code (AFSD) is used. A selection policy is needed:

| Policy             | Description                                        |
|--------------------|----------------------------------------------------|
| AFSD preferred     | Current behavior. AFSP only if AFSD is missing     |
| AFSP preferred     | Lower noise (BPSK(5)). AFSD only if AFSP is missing |
| Selectable         | Global variable (e.g. `sdr_afs_obs_sel`) and a `pocket_trk` option |

Recommendation: selectable, default AFSP preferred. To be decided.

Optional consistency check: if both AFSD and AFSP exist for the same
satellite, log the difference `P(9P) - P(9D)`, and reject AFSP if it differs by
more than e.g. 1 km (wrong ambiguity resolution or wrong tertiary sync).

### 4.8 Not changed

- `init_epoch()`: still needs a channel with `week > 0` (AFSD). AFSP alone
  cannot initialize the PVT epoch.
- Ephemeris: still decoded from AFSD SB2 only (`decode_afs_frame()`).
- `out_log_obs()`: already outputs AFSP as signal code `9P`.
- `sdr_ch.c`: no change. `ter_sync` already provides the needed timing.

## 5. Constraints and Assumptions

- At least one AFSD channel must be frame-synced and decoded to give the week,
  the receiver time (PVT epoch) and the ephemeris.
- |travel time + receiver clock error| < 6 s for the 12 s ambiguity
  resolution.
- AFS time is treated as GPST in the receiver, as for AFSD.
- The interim code assignment (LSIS Table 11) with tertiary phase offset 0 is
  assumed. If a non-zero tertiary phase offset is assigned in the future, it
  must be added to `ch->tow` in `decode_AFSP()`.

## 6. Verification Plan

1. Generate a 90 s 2-bit IQ file with afs_sim (PRN 2-8):
   `./afs_sim -t 90 -b 2 afssim_iq2.bin`
2. Run `pocket_trk -sig AFSD -prn 2-8 -sig AFSP -prn 2-8 -f 12 -IQ 2 -log log.txt`.
3. Check the `$OBS` log:
   - `9P` observations appear for all PRNs once the PVT epoch is initialized.
   - `P(9P) - P(9D)` for the same satellite and epoch is within a few meters
     (both components are coherent in the simulator).
   - LLI bit 2 of `9P` is cleared only after the tertiary sync.
4. Check `$POS`: positions computed with AFSP preferred and with AFSD
   preferred are consistent, and the AFSP solution shows smaller scatter.
5. Check the 12 s ambiguity: `ch->tow` of AFSP modulo 12000 ms equals that of
   AFSD of the same PRN (debug trace in `gen_prng()` already prints
   `ch->tow`).
6. Optional robustness test: start the PVT with AFSD channels of some PRNs
   only, and confirm that AFSP-only satellites are used in the solution.

## 7. Decisions Needed Before Implementation

1. Observation selection policy in `update_sol_afs()` (4.7).
2. Ambiguity window in `gen_prng()`: `[-6, 6)` s or `[-1, 11)` s (4.3).
3. Whether to include the carrier phase half-cycle correction (4.5, optional).
4. Whether to add the AFSD/AFSP consistency check (4.7, optional).
