# NOTE 2026-11-10 ORION -> viper: raycheap17 addendum -- ROOT CAUSE of the sp time non-convergence (thin switch fse +
# halfrange beam), Part 1b SP results, the SP-fix candidate that failed. Design for you to implement (user: viper fixes).

Follows `NOTE-2026-11-10-orion-raycheap17.md` (e01dcddea). Same binary / dirs. No code pushed. A diagnostic patch exists
locally only (below).

## 1. Root cause (pulse_u SP T; read on dea3e3732, src = /orion/ptmp/jinma/raycheap16_1110/src_new)
**The stored x1 face flux is the BLEND, and it is reused as the central branch's memory.**
- rad_m1_implicit.cpp:16442 (inner wall face: 16229), kernel m1_impl_face:
  `fn = (1.0 - al)*fn + (cl/ch)*gh;  f0_(m,k,j,i) = fn;`  with fn = the face-eliminated central M1 flux
  th*(F0^n - ch c dt dP/dx ...), gh = HCL E_L + HCR E_R + DG the half-range upwind flux, al = w from M1ThinFace
  (smoothstep of the formal solution's |H|/J between implicit_thin_h0 0.3 and _h1 0.6): **w has no dt in it**.
- Next step f0x1 -> f0x1n (:13271), and every row takes it as the central memory `(1-al)*th*f0n` (:15161, 15249,
  15416, 15587, 15685); the hesdirk2 slope (f0 - f0n)/dt (:18174) likewise. So per step
  F^{n+1} = (1-w) th (F^n + dt G) + w F_hr(E): an O(1)-per-STEP relaxation of F toward F_hr, i.e. a source
  (w/dt)(F_hr - F). For 0 < w < 1 the discrete equations change with dt -> the dt ladder has no common limit. The
  comment above :16442 ("in the thin limit alpha -> 1 [the memory] must NOT survive") is right for w in {0,1} only.
- pulse_u has |H|/J in (0.3, 0.6) on almost every face (w strictly inside). The fine-dt solutions develop sharp
  fronts at r ~ 1.66, 2.5 and an odd-even oscillation (i ~ 42-56) that grows as dt shrinks; S09+gd is smooth.

Evidence (pulse_u SP T, E):
| arm | d per pair | p |
|---|---|---|
| SP (production) | 3.0e-2 3.6e-2 3.8e-2 4.0e-2 | flat |
| implicit_thin_h0=10, _h1=11 (w = 0: central only) | 7.84e-3 4.15e-3 1.59e-3 5.01e-4 | 0.92 1.38 1.67 |
| implicit_thin_h0=0, _h1=1e-6 (w = 1: no memory), l to 32 | 5.06e-3 5.10e-3 3.87e-3 2.60e-3 1.31e-3 7.16e-4 | -0.01 0.40 0.58 0.99 0.87 |
| same + time_scheme=be | 6.69e-3 5.55e-3 3.81e-3 2.43e-3 | 0.27 0.54 0.65 |
| implicit_blend = idort / blend_tau0 1e6 / blend_flo 0, fhi 1e-6 | bitwise SP (under fse the weight is M1ThinFace's only) | |
| implicit_flux=berthon | FATAL | |
Every w strictly inside (0,1) fails; w = 0 and w = 1 both converge.

**Second, separate defect: the half-range flux itself is 1st order in time** (w = 1 above: p -> ~0.9-1.0, about the
same as be). Not yet tested, best reading: vgd_hr is rebuilt only at the vet_gd sweep (VetGdHalfRange,
rad_m1_vetgd.cpp:5982), i.e. held at t^n for BOTH hesdirk2 stages, while the vet tensor is extrapolated
(time2_vet_col = extrap2); the thin-switch faces are also built once per solve (doface, rad_m1_implicit.cpp:14504).

## 2. Proposed design (for you)
(a) Keep a separate central-branch memory F_c: applied face flux = (1-w) F_c + w F_hr(E), F_c evolved by the plain
    M1 face equation (its own th, dP/dx, memory F_c^n). Then w stays dt-free and nothing relaxes per step. F_c has to
    follow every copy f0x1 gets: f0x1 -> f0x1n (:13271), the t2f1 saves / restores (rad_m1_time2.cpp ~688, 1006
    avgf, 1158, 1222, 1459, 1486, 1535, 1602), the stage increments t2i and the kk F1 slopes, restarts; and the
    lateral faces (rad_m1_latface_body.hpp, M1HrFace path) the same.
(b) Put vgd_hr (and the thin-switch inputs |H|/J) at the stage time like the tensor (extrap2 / predict at t_n + dt),
    or rebuild per stage.
Pass criterion to aim for: pulse_u SP T p >= 1.8 with S09+gd bitwise, then the full sp suite.

Diagnostic only (local, not for production): /orion/ptmp/jinma/raycheap17_1110/fixhr/fixhr_wdt.diff, key
implicit_hr_wdt K (default 0 = bitwise, checked on SP and S09+gd at T_1/T_16): scales AL, HCL, HCR, DG by
min(1, K c dt/dx) (K > 0) or X/(1 + w X), X = |K| c dt/dx (K < 0), turning the per-step relaxation into a rate.
pulse_u SP T with K = 1: d 2.27e-2 1.33e-2 5.24e-3 2.86e-3 1.30e-3 6.74e-4 (l to 64), p 0.78 1.34 0.87 1.14 0.94 --
the flat non-convergence is gone, ~1st order left (= defect b). (With K = 1 it is bitwise wherever c dt/dx >= 1, i.e.
at production implicit_cfl 1e6 it would change nothing; K < 0 changes the large-dt limit.)

## 3. Part 1b (pulse_s, pulse_edd_u, lat_th_atm, lat_ph_atm, mvp_u) -- Part 0 drops added (scripts/drops.json)
pulse_s, pulse_edd_u: implicit_vimp=false; lat_th_atm, lat_ph_atm: vet_gd_allow_clamp=true (`vet_gd: clamped lateral
reads (MeshBlocks too small for the reach)`) + implicit_vimp=false; mvp_u: vet_gd off (reflecting x1) + vimp false.
pulse_s level C_32 still FATALs with the clamp message (its coarsest lateral blocks): that level is missing from the
pulse_s SP / SP-fix C ladders (the p below are from 4 levels).

| case | S09 C (control) | SP C | SP X | SP T |
|---|---|---|---|---|
| pulse_edd_u | E 2.01 1.99 2.00 / F1 1.91 1.98 2.00 | = S09 | 1.99 2.00 2.00 | 1.99 2.01 2.01 PASS |
| pulse_s | E 1.61 1.98 2.07 | 3.08 0.31 0.02 | | -0.31 0.41 0.10 |
| lat_th_atm (Kn 0.06) | E 1.98 1.99 / F2 1.83 1.96 | E 1.37 0.81 / F2 0.67 -0.41 FAIL | E 1.62 0.69 / F2 -0.66 | E -0.16 0.15 -0.14 FAIL |
| lat_ph_atm (Kn 0.06) | E 1.98 1.99 / F3 1.71 1.94 | E 1.68 1.54 / F3 0.70 0.13 FAIL | E 2.04 1.41 / F3 -0.16 | E 0.38 0.34 -0.08 FAIL |
| mvp_u (vet_gd off) | | | | E -0.22 0.09 0.38, T 2.00 |
mvp_u fails in E with vet_gd OFF -> not (only) the halfrange path; the S09+<f> forward bisection on lat_th_atm,
lat_ph_atm, pulse_s and mvp_u (T mode) is running; it goes in a further addendum.

**SP-fix (beam closure + blend tau_f + thin switch none / corr none + vet_col_ncore 64) is NOT a fix**: it breaks
cases SP passes (fs_u T -> 0.48 -1.09 0.81, fs_u C -0.26 0.21 -0.05; pulse_u T -0.06 -0.19 -0.57). Do not use.
