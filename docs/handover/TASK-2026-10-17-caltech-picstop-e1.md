# TASK raven -> Caltech: Picard stop experiment 2 on the REAL agcE1r restart (picstop-1017), GPU test, NOT production

From the Raven worker PICSTOP (USER approved experiment 2 on 10-17). Tests only: do NOT touch, continue, or edit agcE1r (copy its
last restart to a new test dir). Raven runs the same arms on the S3/P268 half-res restarts; this TASK adds the real E1 536 grid.

## What is tested (branch picstop-1017 @ 4775edf5 on fork jing-ze-ma/athenak; all keys OPT-IN, default = bitwise b16b66db/bc5d24c1)
- `rad_m1/implicit_final_resid_log = true`: one extra assembly per solve, logs the scaled nonlinear residual of the accepted iterate
  (PICRES lines; diagnostic, trajectory unchanged).
- (a) `rad_m1/implicit_conv_est_nofloor = true`: the conv_est contraction q ignores cells at e_floor.
- (b) `rad_m1/implicit_final_resid_stop = 1e-4`: accept a solve only if max(rc, rT) over non-floor cells <= 1e-4; else the
  check assembly becomes the next Picard pass (up to implicit_maxit; then NON-CONVERGED -> usual dt_retry).

## Steps (one GPU job, same GPU count/layout as agcE1r links, <= 2 h total)
1. Build `athena_gpu_he_star_m1_4775edf5_nofma` with exactly the options of athena_gpu_he_star_m1_b16b66db_nofma (CUDA 12.9,
   HOPPER90, MPI, host -ffp-contract=off). Do not use it for production.
2. New dir `picstop_1017/`; R = copy of agcE1r `agcar3d.00302.rst` (t 1.48e8). Common keys for every arm:
   `time/ndiag=1 output2/dt=DT output3/dt=DT output5/dt=1e30 time/tlim=1.48e8+DT` with DT chosen so the OFF arm takes ~8-10 min
   (start with DT = 2e5; outputs 2/3 must land exactly at tlim so the arms are compared at the same time).
   The production keys come from the restart (ROL2U: implicit_tol_rel 0.1, opac_lag pred/2, gas_newton_u4) -- do not repeat them.
3. Arms (each from R, own run dir, `-t` = 2x the OFF wall):
   - `off`  : `rad_m1/implicit_final_resid_log=true`
   - `aN`   : off + `rad_m1/implicit_conv_est_nofloor=true`
   - `b`    : off + `rad_m1/implicit_final_resid_stop=1e-4`
   - `ab`   : off + both
   - `ref`  : off + `rad_m1/implicit_tol=1e-10 rad_m1/implicit_tol_rel=0.0 rad_m1/implicit_conv_est=false`
   - `off_nl`, `b_nl`: as off / b WITHOUT implicit_final_resid_log (clean wall cost)
   Write `rc=$? wall=<seconds>` into `<arm>/rc.txt` after each run (stat.py reads it).
4. Analysis (docs/handover/picstop-1017, check MD5SUMS; python3 + numpy):
   `python3 stat.py picstop_1017/{off,aN,b,ab,ref,off_nl,b_nl}` and
   `python3 cmpbin.py picstop_1017/ref picstop_1017/{off,aN,b,ab} --ib 380:480` (band = the cells above the photosphere on the
   536 grid; adjust if you know the photosphere index, and say which you used).

## Report back (NOTE-2026-10-17-caltech-picstop-e1.md)
Both tables verbatim, the `<rad_m1> picstop:` and `implicit transport:` summary lines of every arm, DT/tlim, GPU count, binary md5,
and any FATAL / NON-CONVERGED / DIVERGED / NaN. Nothing else; no production change (the USER decides).
