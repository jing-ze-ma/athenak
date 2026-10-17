# TASK 2026-10-18 raven -> caltech: longer agcE1r comparison of implicit_conv_est_nofloor (follow-up to NOTE-2026-10-17-caltech-picstop-e1)

USER (10-18): hold `implicit_conv_est_nofloor` out of production for now; first check whether its E1 effect on the
solution is real or noise. Your 10-17 run (1 DT = 7e5 s from rst 302) showed: tail >1e-2 off 12.4 % -> nofloor 4.8 %
(ref 4.3 %), +20 % wall, but mean density difference vs ref 2.6e-3 (nofloor) vs 2.3e-3 (off), max larger.

GPU test only; never touch agcE1r (read-only COPY of the restart in a scratch dir, as last time).

## Do
Same binary (athena_gpu_he_star_m1_4775edf5_nofma, md5 b375ccf1...), same 1 node x 2 H200 launch, from the same
copy of agcar3d.00302.rst. Three arms only:
- `off`  : production keys (as your `off`);
- `aN`   : + `rad_m1/implicit_conv_est_nofloor=true`;
- `ref`  : tol 1e-10, tol_rel 0, conv_est off (as your `ref`).
Run each to **t_rst + 5 DT = t_rst + 3.5e6 s** (5x longer than last time), outputs (hydro_w + m1 bins) every DT, so
the comparison is available at 1, 2, 3, 4, 5 DT. Walltime as needed (OFF took 623 s per DT).

## Report (NOTE on bsg-files-1009)
- per arm: stat.py line (solves, passes, NC, retries, tail >1e-4 / >1e-2, wall);
- cmpbin.py of off and aN vs ref **at each of the 5 output times** (dens/eint/m1_e/m1_f1 max and mean, full and
  band --ib 380:480), plus the same for off vs aN;
- whether the off-vs-ref and aN-vs-ref differences grow with time and which is closer to ref at 5 DT.
