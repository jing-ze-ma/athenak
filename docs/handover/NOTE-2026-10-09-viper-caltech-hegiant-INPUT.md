# NOTE viper -> Caltech: He giant fresh N897 input REVISED (one key) before link 1

Inner-boundary IC check (10-09): the IC column's eint came from the python IC tool's EOS; with the code EOS the
T round trip was off by up to 1.7e-3 (61.6 Rsun). Fix = an existing key, no new binary:
- In the input written by SETUP.sh (hegiant_fresh897.athinput), add in `<problem>`:
  `he_ic_eint_from_t = true`
- Everything else unchanged. (viper md5 with the key: 9034dac91ab256e697dd6d01602cf312, viper paths.)
- Smoke with the revised input before link 1 (as in the TASK). The remaining 5e-3 in the 4 wall cells and ~1e-4
  inside are grid truncation of the HSE march (2nd order, converges with resolution) — expected, not a bug.
