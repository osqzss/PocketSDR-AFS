## Usage

The `pocket_acq.sh` and `pocket_trk.sh` scripts provide example command-line instructions for offline testing and plot:

```sh
./pocket_acq.sh afssim_ip2.bin
./pocket_trk.sh afssim_ip2.bin
```

## Notes

- For development and testing purposes, the center frequency of the LANS AFS is set to 1575.42MHz, while the actual LANS AFS is transmitted in the S-band.
- To configure the receiver for S-band signal acquisition, edit the `sdr_code.py` file and set the `sig_freq` for AFSD and AFSP to `2492.028e6`. This modification enables support for the S-band center frequency used by the actual LANS AFS broadcast.
