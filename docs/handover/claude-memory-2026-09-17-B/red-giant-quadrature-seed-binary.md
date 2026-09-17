---
name: red-giant-quadrature-seed-binary
description: "Three findings 2026-09-09 04:00: (1) ck_nquad=1 (hemispheric mean, D=1.66) gives 0.904x the diffusion flux -- use ck_nquad=2 (2-pt Gauss, sum w mu^2 = 1/3 exactly) for grey; (2) the RCB sup dipole is NOT a flux deficit (R0 == MLT v3 to 4 digits) -- suspect the vpert SEED advecting the entropy step across the RCB; (3) rebuilding a binary under a running MPI job segfaults it -- pin binaries per run."
metadata:
  type: project
---

**Quadrature.** `two_stream_rt.hpp` ~1066: `ck_nquad = 1` -> mu = 1/1.66, weight pi (F = pi I);
`ck_nquad = 2` -> Gauss-Legendre on [0,1]: mu = 0.2113, 0.7887, w = 0.5, 0.5, so
sum w mu^2 = 1/3 exactly and the thick limit is F = -(4 pi/3) dB/dtau, i.e. diffusion.
The hemispheric mean gives (1/1.66)(0.5)/(1/3) = 0.904 of the diffusion flux -- a fixed 10 %
mismatch exactly where the two-stream hands over to the conduction operator (tau 10-100).
This is the "handover consistent to 10 %" of [[radiative-conduction-deep-interior]].  The
production input had ck_nquad = 1.  Use 2 for grey (user's call; no code change needed).

**RCB dipole is not a flux deficit.** Closure v3 (Schwarzschild-gated on the noise-free
grad_rad, deficit D = F_req - F_rad_cond - F_conv_res, no x gate) closes F_used + F_rad =
F_req at t = 0 to 1e-4 at faces 7-12, yet sup(5)/sup(6) in R1g equal R0f to FOUR digits at
every time.  Dump budget at t = 4e5, faces 4-12, identical in both: F_res/F_req = -0.16,
-0.10, +0.31, +1.64, +2.28 at i = 4, 5, 6, 8, 10.  Whatever drives the dipole is in the
resolved dynamics.  Hypothesis under test (R2a/R2c): the SEED -- vpert imposes
v1 = sin(3*2pi*x) over the whole block including across the RCB entropy step, so a standing
radial flow advects entropy across it at a constant rate: linear from t = 0, closure-blind.
New `problem/vpert_rmin` seeds only r > rmin (1.6e12 = 0.5 R in the tests).

**v3's own defect:** D did not subtract the two-stream flux; above the blend F_rad_cond -> 0
while the two-stream carries L, so MLT was asked for the full L in thin gas -> Mach-1
outflow at the top -> dt collapse at 4e5.  Fix: D also subtracts the shell-mean
`two_stream_rt::rt_face_flux()`.  With v3, L_rad,out/L = 0.72-0.88 (R0f: 0.30-0.60) and
T(tau=2/3) = 3802 K (R0f 3564; target 4000) -- the direction is right.

**Operational:** `make` into build_rg_prod while R0f (srun, 6 ranks) was running SEGFAULTED
it at t = 2.1e6 (the executable's pages were replaced).  Every run must copy the binary
into its own directory and point sub.sh at the copy.  R0f's from-scratch data to 2.1e6 is
intact.  Also: one RT column is NOT L_out -- the single-column emergent flux read 0.20 L
where the all-column diagnostic (`L_rad,out/L` in the face-budget print) gives 0.30-0.60.
