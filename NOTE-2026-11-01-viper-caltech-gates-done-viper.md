# NOTE 2026-11-01 viper -> Caltech (USER): G1 (POL) and G2 (CSV) from TASK 4bdc9e767 already PASSED on viper; skip them

viper apudev finally ran two of the gates we moved to you.

## G1. POL (fork/rkgpolar-test-1101 @ 33a4d519f) vs fin51-r-fix 65658289c: PASS
- Job 12554957, MI300A, W121 cs MHD restarts at rot 320, key-free.
- Final rst and dhj.mhd.hst IDENTICAL for all three:
  - b100 (explicit);
  - s10 (RKG);
  - s10 + rkg_lean + rkg_dt_safety 0.5.
- The polar fix goes to fin51-r-fix once its CPU suite is done. We will NOTE the hash.

## G2. CSV defaults (fork/csvertex-defaults-1101 @ 6029a2aa7): PASS
- **Old restart without the key vs b935e145:** GPU job 12553770 (b100, s10) gives final rst and mhd.hst IDENTICAL; CPU gives 3/3 rst identical. The log says "restart input has no key: the old default, key not added".
- **Fresh start, no key vs explicit true:** rst data identical (4/4). Explicit false = old (4/4).
- **cs MHD blast:** the default run equals the explicit-on run (max p 19.822, |B|max 0.52981, div B 4.9e-14).
- **mhd *_cpu suite:** 31 passed.
- The default still waits for your csvertex A/B (TASK 76efbd1d5).

## Still wanted from you, in priority order
1. max_eta R15/R14/R13 + X13;
2. csvertex A/B;
3. ambipolar T15/TA15;
4. tag-only floor arm;
5. G3 TRC defaults;
6. G4/rkghalo.

Raven is also queued on rkghalo and TRC; the first result wins.
