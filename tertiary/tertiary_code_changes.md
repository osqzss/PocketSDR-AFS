# AFS-Q Tertiary Code Synchronization - Source Code Changes

This document summarizes the source code changes made to support the tertiary
code synchronization of the LANS AFS pilot channel (AFS-Q) in PocketSDR-AFS.

Reference: LSIS V1.0 (doc/LSIS_AFS_v1.pdf), LSIS-171, LSIS-200, LSIS-210,
LSIS-220, LSIS-224 and Appendix E.

## 1. AFS-Q Tiered Code Structure

| Code      | Length [chips] | Chip period | Code period |
|-----------|----------------|-------------|-------------|
| Primary   | 10230          | -           | 2 ms        |
| Secondary | 4              | 2 ms        | 8 ms        |
| Tertiary  | 1500           | 8 ms        | 12 s        |

- One secondary code period (4 primary code periods) corresponds to one
  tertiary code chip (LSIS-171).
- The start of the tertiary code is synchronized with the start of the AFS-I
  data frame (LSIS-220).
- Interim code assignment (Table 11): tertiary PRN = node ID, phase offset = 0.

## 2. Changed Files

| File                  | Change                                               |
|-----------------------|------------------------------------------------------|
| `src/sdr_code_afs.c`  | New. Tertiary code HEX strings for PRN 1-210         |
| `src/sdr_code.c`      | Tertiary code generation, new API `sdr_ter_code()`   |
| `src/pocket_sdr.h`    | New fields in `sdr_trk_t` / `sdr_ch_t`, API prototype |
| `src/sdr_ch.c`        | Tertiary code sync and removal, sec-code index fix   |
| `src/sdr_rcv.c`       | Show "T" in the sync status for tertiary code sync   |
| `lib/build/libsdr.mk` | Add `sdr_code_afs.o` to the library                  |
| `test/utest/makefile` | Add `sdr_code_afs.o` (needed to link `sdr_code.o`)   |

## 3. Tertiary Code Generation

### src/sdr_code_afs.c (new)

- `const char *code_afs_ter[210]`: the 1500-chip tertiary codes as HEX strings
  (375 hex digits per PRN), in the same format as `sdr_code_gal.c`.
- Generated from `tertiary/008_Weil1500hex210prns.txt`. All bits were verified
  against `tertiary/tcode.txt`, and the first/last 24 chips of PRN 1
  (`E4D09E` / `F2164C`) match Table E-2 of the LSIS.

### src/sdr_code.c

- `extern const char *code_afs_ter[];` added.
- `ter_code_AFSP()`: completed. Valid for PRN 1-210. The code is converted by
  `read_code_hex()` and cached in `AFSPT[]`. The chip mapping is the same as
  the secondary code (`CHIP[] = {-1, 1}`, i.e. bit 1 -> +1).
- `ter_code()` (static): returns `ter_code_AFSP()` for `AFSP`, otherwise a
  dummy code `{1}` with length 1 (no tertiary code).
- `sdr_ter_code()` (new API, thread-safe with a mutex like `sdr_sec_code()`):

```c
int8_t *sdr_ter_code(const char *sig, int prn, int *N);
```

## 4. Data Structures (src/pocket_sdr.h)

New fields in `sdr_trk_t`:

| Field      | Description                                                    |
|------------|----------------------------------------------------------------|
| `ter_sync` | Tertiary code sync status (0: no-sync, >0: lock count at the start of a tertiary code period) |
| `ter_pol`  | Tertiary code polarity (+1 / -1)                               |
| `ter_err`  | Tertiary symbol error count (for sync lost detection)          |
| `ter_nsym` | Number of tertiary symbols accumulated for sync                |
| `ter_syms` | Tertiary symbol buffer (`float`, length 1500)                  |

New fields in `sdr_ch_t`: `ter_code` (tertiary code) and `len_ter_code`
(length of tertiary code, 1: no tertiary code).

## 5. Tertiary Code Sync (src/sdr_ch.c)

### Channel setup

- `sdr_ch_new()`: gets the tertiary code by `sdr_ter_code()` and allocates
  `trk->ter_syms` only if `len_ter_code >= 2`.
- `trk_new()`, `trk_init()`, `trk_free()` and the signal-lost branch in
  `track_sig()` initialize or free the new fields.

### Processing flow in track_sig()

```
sdr_corr_std()        correlation (P history updated, lock++)
sync_sec_code()       sync and remove secondary code
sync_ter_code()       sync and remove tertiary code   <-- new
FLL/PLL, DLL, CN0
sdr_nav_decode()
```

`sync_ter_code(ch, N, M)` is called when `len_ter_code >= 2` and
`lock * T >= T_NPULLIN` (N = 4: secondary code length, M = 1500: tertiary code
length).

### Algorithm

1. Requires secondary code sync (`sec_sync > 0`). Otherwise the tertiary sync
   state and symbol buffer are reset (`reset_ter_code()`).
2. Tertiary chip boundaries are given by the secondary code sync:
   `n = lock - sec_sync`, and a tertiary chip ends when `n % N == 0`.
   The tertiary chip count is `k = n / N - 1`.
3. Tertiary symbol: `S = sum of the last N P correlations` (I component) with
   the secondary code removed (8 ms coherent sum).
4. Before sync, `S` is stored in `ter_syms[k % M]`. Since the buffer is indexed
   by `k % M`, duplicated or skipped epochs caused by code offset wrapping
   (`lock--` / `lock++`) do not break the chip alignment.
5. When at least `N_TSYNC` (250 symbols = 2 s) are accumulated, the phase is
   searched every `N_TCHK` (125 symbols = 1 s) by `search_ter_code()`:
   `P(p) = sum_j ter_syms[j] * ter_code[(j + p) % M]` for all p = 0, ..., M-1.
   Slots not yet received are zero and do not contribute.
6. Sync condition:
   `|P_max| >= THRES_TSYNC * R` and `R / m >= THRES_SYNC * N`,
   where `R = sum |ter_syms[j]|` and `m = min(ter_nsym, M)`.
7. On sync, chip count k corresponds to tertiary chip `(k + p) % M`, so the
   tertiary code starts at chip count `k0 = (M - p) % M`:
   `ter_sync = sec_sync + N * k0`, `ter_pol = sign(P_max)`.
8. After sync, the tertiary code is removed from the correlations `C[]` and the
   latest P history in every epoch:
   `i = ((lock - ter_sync - 1) mod (N * M)) / N`,
   `C = ter_code[i] * ter_pol`.
   After removal, the tertiary symbols should always be positive.
9. Sync lost: negative symbols are counted in `ter_err`. If the count exceeds
   `THRES_TLOST * N_TCHK` within `N_TCHK` symbols (1 s), the sync is reset and
   re-searched. A half-cycle slip of the PLL is detected in the same way.

### Parameters

| Macro         | Value | Description                                      |
|---------------|-------|--------------------------------------------------|
| `THRES_TSYNC` | 0.8   | Threshold for sync (normalized correlation)      |
| `THRES_TLOST` | 0.25  | Threshold for sync lost (symbol error rate)      |
| `N_TSYNC`     | 250   | Min number of tertiary symbols for sync (2 s)    |
| `N_TCHK`      | 125   | Interval of sync search / lost check (1 s)       |

### Log messages (log level 3)

```
$LOG,<time>,AFSP,<prn>,TERTIARY CODE SYNC (<chip index>,<polarity>,<corr>)
$LOG,<time>,AFSP,<prn>,TERTIARY CODE LOST (<error count>)
```

### Sync status display (src/sdr_rcv.c)

`sync_stat()` builds the 4-character sync status (SYNC column) shown by
`pocket_trk`. The first character now shows the tertiary code sync:

| 1st char | Condition                                              |
|----------|--------------------------------------------------------|
| `T`      | Tertiary code synced (`ter_sync > 0`), implies secondary code sync |
| `S`      | Secondary code synced only (`sec_sync > 0`)            |
| `-`      | No secondary code sync                                 |

The other characters (`B`: symbol sync, `F`: frame sync, `R`: reversed
polarity) are unchanged.

### Computational load

The phase search costs 1500 x 1500 multiply-adds per attempt. It runs at most
once per second per AFSP channel, and only while the tertiary code is not
synced. The inner loop is split in two to avoid the modulo operation.

## 6. Fix of Existing Code: sync_sec_code()

When `lock == sec_sync` (the epoch at which the secondary code sync is
established), the index `(lock - sec_sync - 1) % N` was -1 in C, which read
`sec_code[-1]` (out of bounds).

```c
// before
int8_t C = ch->sec_code[(ch->lock - ch->trk->sec_sync - 1) % N] *
// after
int8_t C = ch->sec_code[(ch->lock - ch->trk->sec_sync - 1 + N) % N] *
```

This affects all signals with a secondary code, but only that one epoch. The
results of later epochs are unchanged.

## 7. Verification

The sync logic was tested by a simulation driving `sync_sec_code()` and
`sync_ter_code()` with synthetic P correlations (PRN 1-12, random tertiary
phases, noise levels equivalent to about 33, 35 and 38 dB-Hz, and noise-free).

- Every sync was declared at the correct tertiary chip index (no false sync).
- Noise-free, 35 and 38 dB-Hz: sync about 2 s after the secondary code sync.
- Half-cycle slip: detected as sync lost, then re-synced correctly.
- Tertiary phase 0 and runs longer than one tertiary code period: OK.
- 33 dB-Hz: in 6 of 12 runs the existing secondary code lost check
  (`THRES_LOST`) reset the secondary sync near the end of the run, and the
  tertiary sync was reset with it.

Not yet tested with real IF data or with a full link of `pocket_trk`.

## 8. Known Limitations and Future Work

- A secondary code re-sync resets the tertiary sync, so the tertiary code needs
  another 2 s to re-sync. At low C/N0, frequent secondary code re-sync may
  delay the tertiary sync. The tertiary sync could be kept when the new
  `sec_sync` keeps the same chip alignment.
- `decode_AFSP()` in `src/sdr_nav.c` is still empty, and `ch->tow` is not set
  for AFSP. The tertiary sync gives only the time within the 12 s frame, and
  `gen_prng()` with `tow_v = 2` assumes a 100 ms ambiguity, so setting the TOW
  here would give wrong pseudoranges. The week number and TOW should be taken
  from the AFSD channel of the same PRN.
- The combined polarity `sec_pol * ter_pol` resolves the half-cycle ambiguity
  of the carrier phase, but it is not yet used for the observation data.
- The Python version (`python_afs/`) is not updated.
