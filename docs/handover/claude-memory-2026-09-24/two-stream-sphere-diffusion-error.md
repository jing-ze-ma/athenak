---
name: two-stream-sphere-diffusion-error
description: The spherical two-stream (J = A I, commit 1159a8f3) is exact in the transparent limit but in the DIFFUSION limit gives F/F_exact = 1 - H_T/(2r) (verified 1-4 % on the He4 star; wrong SIGN if H_T > 2r). Production: hand tau > 20..300 to radiative diffusion; long-term fix = variable Eddington factor moment scheme
metadata:
  type: project
---

Derivation (09-17, tests_3d/handover, wt_he4): transporting J = A I with source A B solves
mu dI/dr = -kappa rho (I - B) - 2 mu I/r, i.e. the curvature term (1-mu^2)/r dI/dmu replaced by a
sink. Transparent limit: I ~ 1/r^2, L conserved (exact). Diffusion limit: D = -(mu/kappa rho)
d(AB)/dr -> F_2s = -(mu/kappa rho)(dB/dr + 2B/r) -> F_2s/F_exact = 1 - H_T/(2r), H_T = T/|dT/dr|
(B ~ T^4). The exact first moment has no such term (f = 1/3 kills the sphericity term).
Measured on the He4 smoke dump: predicted vs measured agree 1-4 % at every face (min 0.10-0.14
at 0.65 R where H_T/r = 1.7). ck_nquad 4 and rt_impl_mixed 0 change nothing.
HANDOVER SCAN (one-cycle mlt_dump budgets): rad_tau_lo/hi 20/300 (H1) closes the budget < 1 % in
the stable layer 0.50-0.67 R that killed arms A/B/D; residual 0.74-0.91 over 0.72-0.90 R is the
convective share (Rosseland itself is 0.85 there) -> MLT closure. 5/50 overshoots 31 % near the
photosphere (Rosseland invalid); 50/1000 leaves the interior on the two-stream. Above tau 20 the
two-stream error is 3-4 % (H_T/r < 0.1).
**How to apply:** on any sphere spanning a sizeable fraction of R: rad_tau_lo/hi ~ 20/300, L injected
at the conduction wall (rad_flux_inner), rt_bottom_flux false, mlt_alpha 1.5 shell-mean closure.
Never run the whole-column two-stream on a sphere with H_T ~ r. Long-term: rewrite the spherical
two-stream as moment equations with a variable Eddington factor (sphericity only in the zeroth
moment), gated by tests_r2/rg1d (transparent) + tests_r2/thick (diffusion, being built).
See [[he4-presn-global-plan]], [[rt-shallow-ramp-missplit]].

UNIT TEST (tests_r2/thick, inputs/tests/two_stream_sph_thick.athinput, 76e21c7d): 8 cases
(rho ~ r^-n, n 0/1/3/7; tau 30/100/300; r_out/r_in 1.1/1.5/2), exact radiative-equilibrium T(r)
as ic, one cycle, mode 3 pcr. F_2s/F_req = 1 - H_T/(2r) face by face, residual 2e-5..7e-4; tau- and
thickness-independent except through H_T/r; Rosseland control 0.9997-1.0002. 20x deficit at the
bottom of n=0. Gotcha: srun swallows a piped stdin case list -> give srun < /dev/null.

FIXED 09-17 evening (8bca3dfa on he4-presn-global, tests_r4/README.md): the layer solve stays
plane-parallel per unit area with source B; ALL geometry is at the faces: S = (I+ + I-)/2 and
A D = A (I+ - I-)/2 passed continuously (u_above = u_below + c, d_below = d_above + c,
c = beta (d_above - u_below), beta = (A_above - A_below)/(A_above + A_below)). Phi = A F is
single-valued at faces -> deposit -(F_hi - F_lo)/dx telescopes to -(A_top F_top - A_bot F_bot)/V.
Mode 3: exact (new RTCol3::FormalSph Riccati pair; PCR path calls it in Kokkos::single on
curvilinear meshes); mode 0: exact up-sweep, down-sweep lagged to O(beta^2). Legacy whole-cell
layers refused under rt_grey on curvilinear meshes. GATES: box bitwise (mode 0 and 3), rg1d L
const 1.0002, thick test median 1e-4 / worst 1.2 % (was 16-97 %), He star t=0 F_2s/F_raddiff
0.983-1.016. TRAP: mlt_dump col 17 F_2s is the EXPLICIT sweep's flux unless rt_col3_skip_sweep.
MLT closure double count fixed (F_cond + (1-w) F_2s). Arm B' now survives 1 turnover (old 0.37);
remaining: upper FeCZ drains 54 %/turnover, efloor -> EOS-garbage detonation with vpert 1e-2.
