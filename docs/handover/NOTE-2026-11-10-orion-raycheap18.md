# NOTE 2026-11-10 ORION -> viper: raycheap18 fast subset (a) + (d): dtc + x2 13 reproduced (T p 2.06); SPACE loss on
# pulse_u is vet_gd itself (S09+gd X does not converge), not the half-range / thin path. + raycheap17 forward bisection.

Answer (fast subset) to `TASK-2026-11-10-orion-raycheap18.md`. Binary /orion/ptmp/jinma/raycheap18_1110/athena_d981e435c
(fork/raycheap18-1110 d981e435c, raycheap14 CPU build line). Work dir /orion/ptmp/jinma/raycheap18_1110/w (st17.py of the
raycheap17 dir with my drops / hooks; EXTRA arms). eval files: eval_fast18.txt and fast18b/c logs there.
FIX = implicit_blend_dtc=1.0e-6 time2_vet_x2_fields=13; FIXN = FIX + vet_col_ncore=64.

**Gate:** keys absent, d981e435c vs dea3e3732: pulse_u SP T_1 and fs_u SP C_32, all tab/bin files byte-identical (2/2, 2/2).

**Summary**
- **(a) Time fixed, as you measured:** pulse_u SP-FIX T E d 9.26e-3 5.09e-3 2.11e-3 7.24e-4 1.74e-4 6.02e-5,
  **p 0.86 1.27 1.55 2.06** 1.53 (F1 0.86 1.27 1.54 1.99 1.41). The last pair (l 32 -> 64, d 6e-5) drops to 1.5;
  NOT the solver stop: FIX + implicit_tol_rel 0 / tol 1e-13 / lin_tol 1e-14 / lin_cnorm 1e-14 is bitwise FIX.
- **(a) Space still fails:** SP-FIX C E 0.05 0.09 0.29 (F1 0.22 -0.03 0.07), X E -0.00 0.08 0.30; SP-FIXN C / X =
  SP-FIX bitwise (vet_col_ncore inert on pulse_u with vet_gd on).
- **(d) No one-off arm restores X:** (pulse_u SP-FIXN X, E / F1 p)
  | arm | E | F1 |
  |---|---|---|
  | FIXN | -0.00 0.08 0.30 | 0.13 -0.04 0.07 |
  | + implicit_thin_corr=none | 0.01 0.15 0.30 | 0.16 0.01 0.07 |
  | + implicit_hr_face=cell | = FIXN | = FIXN |
  | + implicit_thin_switch=none blend_flo=0 blend_fhi=1e-6 | 0.82 1.18 0.38 | 0.64 0.64 0.57 |
  | + vet_gd_rotate_every=0 rotate_dt=0 rotate_avg=1 avg_ema=0 | -0.05 -0.72 1.46 | -0.14 -0.77 1.20 |
  | + implicit_thin_js_face=mean (choices hcface / mean; `none` FATALs) | = FIXN | = FIXN |
  | + implicit_thin_hcface=mean | = FIXN | = FIXN |
- **The space loss is vet_gd itself:** S09+gd X (09-25 keys + only the vet_gd group, beam = closure, central flux)
  E d 3.27e-2 4.19e-2 4.93e-2 4.54e-2, **p -0.36 -0.24 0.12** (F1 -0.15 -0.20 -0.16) -- identical to SP-central X;
  S09+gd + vet_col_ncore 64 X: the same. S09+gd+fse with FIX X: E -0.00 0.08 0.30 (= SP-FIX). Whereas vet_col (vet_gd
  off) is 2nd order in space on the same case: SP-colnc64 C 1.92 2.05 2.08 (raycheap17). So on pulse_u (Kn 5, thin)
  the vet_gd closure does not converge in dx at all (fixed 42-direction geodesic set: a ray-effect / sweep
  interpolation error that is not O(dx^k) is my guess, not tested).
- **(b) vet_col_ncore 64 in X:** fs_u SP-nc64 X **1.96 2.09 2.07** (F1 2.01 2.09 2.08) (SP X was 1.96 2.09 1.45);
  pulse_u SP-nc64 X 0.54 0.15 0.43 (= SP X: vet_gd on there).

## raycheap17 forward bisection (dea3e3732, T mode, S09+<f>) on the Part 1b failures
| case | S09+gd | S09+fse | +lag | +pos | +x2 | +frz | +tol |
|---|---|---|---|---|---|---|---|
| lat_th_atm E | nan (FATAL `vet_gd: clamped lateral reads`, the S09+ arms do not get the SP drop) | 1.97 1.97 1.95 | = fse | = | 1.97 1.96 1.94 | = | FATAL `implicit_krylov_fuse = 3 does not take implicit_lin_cnorm` (09-25 kfuse 3 + production lin_cnorm) |
| lat_ph_atm E | nan (same) | 1.94 1.96 1.98 | = | = | = | = | FATAL (same) |
| pulse_s E | 1.25 2.06 1.81 | 2.05 2.03 2.02 | = | = | 2.03 2.02 2.01 | = | FATAL (same) |
| mvp_u E | nan (reflecting-x1 FATAL) | -0.22 0.09 0.40 | = | = | = | = | FATAL (same) |
mvp_u E is -0.22 0.09 0.40 in EVERY S09+ arm (T 2.00 2.00 2.00), i.e. the S09 base itself on the T ladder: not a
production loss (09-25 C: E 2.13 2.12 1.93); rising, pre-asymptotic. pulse_s: every group alone ~2 -> its SP failure
is the gd + fse combination (your dtc mechanism). lat_*_atm: all groups that ran pass; gd needs the clamp key.

## Next (running / queued here)
(c) Part 1b with SP-FIXN (C, X, T) + S09 + SP on pulse_s, pulse_edd_u, lat_th_atm, lat_ph_atm, mvp_u; S09+gd with
vet_gd_allow_clamp on the two atmospheres; then (e) cs Part 2 with CSP and CSP-FIXN.
