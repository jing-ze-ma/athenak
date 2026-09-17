---
name: he4-presn-global-plan
description: Plan (09-17) for the GLOBAL cubed-sphere 4 Msun presupernova He-star FeCZ model; base = red_giant.cpp; two code blockers (ADI transverse is Cartesian-only fatal on cs; two-stream deposit has no r^2 dilution); full analysis in bench/hestar_presn/GAP_ANALYSIS.md
metadata:
  type: project
---

User decided 09-17: build the 4 Msun presupernova He star (Woosley 2019; M 3.15, log L 4.78,
Teff 49 kK, R 2.37e11) as a full-sphere cubed-sphere model, 2 cells per surface Hp horizontally,
two-stream all along the column, Strang, tapers clear of the FeCZ (FeCZ tau 12..300, base 0.48 R).
GAP ANALYSIS (bench/hestar_presn/GAP_ANALYSIS.md, verified 09-17):
- base = src/pgen/red_giant.cpp (spherical IC march with point mass, open inner/outer BCs, cs
  seeding, WB+point mass on cs); port the box physics (~1000-1200 lines) rather than teach
  box_convection.cpp spheres. Envelope mass 1e-7 Msun -> mstar point mass exact.
- BLOCKER 1: rad_implicit_ang (ADI lod2 transverse diffusion) is a startup FATAL on cs/sp
  (conduction.cpp:339-347, "Cartesian-only in this version"); branch implicit-transverse-raddiff
  has NO commits beyond rt-integration, so the cs generalization does not exist. Interim:
  explicit rad_cap_ang.
- BLOCKER 2: the two-stream split path deposits -(F_top-F_bot)/dx1 with no r^2 dilution
  (two_stream_rt.hpp:331-334, kernel ~:3158): factor 6.4 over r_out/r_in 2.54 -> must become
  (A F)_top-(A F)_bot / V, + mode-3 Jacobian rows; bitwise-inert under rt_plane_parallel.
- first grid: nx1 192 stretched (0.4 R .. tau 1e-2 = 2.4057e11), 6 x 320^2, blocks 192x80x80
  (96 blocks / 8 GPUs), dt ~3 s, ~2000 cycles/turnover (6060 s), ~1 h wall / 8 GPU-h per
  turnover (calibrated 1.6e7 zc/s/GPU on prod_w7); smoke 96 x 6x32^2.
- build order G.0-G.7 in the doc: decide transverse op -> baseline rg -> r^2 fix (bitwise box
  test + 1-D spherical RE test) -> mode-3 plumbing -> IC -> BCs/sponges -> seed + smoke -> prod.

PROGRESS 09-17 07:30: (1) r^2 fix DONE 1159a8f3 on he4-presn-global (area-weighted intensities
J=A I through the sweep + volume deposit + mode-3 rows; box production config BITWISE; 1-D
spherical test L const to 0.03 % vs x6.25 before; ck kernel + mode-1/2 Jacobian NOT converted).
(2) plumbing DONE (he4-rg-plumbing, 6 commits, merged into he4-presn-global): 56 box switches
in red_giant.cpp, rt_surface (panel,x2v,x3v,F_top rows) / rt_profile (8 shell means), input
inputs/hydro/he4_presn_cs.athinput; red giant bitwise. FOUND: red_giant's built-in 1-D march
cannot build a Prad-dominated envelope (fatal aT^4/3 > 0.9 ptop; ignores taper/force; 1e6 rho
error at r_in) -> ic_profile route. (3) SPHERICAL 1-D structure DONE: bench/hestar_presn/
column_sph.py, column_he4_presn_sph.txt, make_ic_sph.py, ic_he4_presn_sph.txt (r rho eint,
taper applied, covers 0.35-1.026 R). On the sphere: FeCZ 0.635-0.968 R (tau 11-208), v_MLT
1.45e7 (Ma 0.42), Fconv/F 0.16, turnover 4705 s, thermal 5.7 d, DENSITY INVERSION 0.72-0.95 R
contrast 3.9 (physical, in both columns), envelope mass 4e-5 M. Plane-parallel column is only
good above 0.98 R. (4) cs ADI transverse: implemented on cs-implicit-transverse (355015e6),
tests pending. (5) NEXT: ic_profile reader + 1-D He gate (agent running on wt_he4).
08:30: inner boundary moved to 0.5 R (user; 0.4-0.5 R = 6.6 Hp of stable 0.24-1 MK gas = half the
radial cells; 1 Hp / 2 Hrho of buffer remain below the base). ic_profile reader DONE 526791b5
(+ input 1512b882): read-back 6e-9 rho, photosphere R/Teff to 0.02/0.45 %. 1-D He gate: mode 0
RUNS (L_out/L 0.95), MODE 3 BLOWS UP on the sphere (L_out/L 8.6e5 in a few cycles; box bitwise)
-> bug in the spherical arithmetic of the mode-3 rows (Av/Ac/Vc); agent fixing on the
constant-kappa 1-D test. rt_force_verbose normalises by g(r_in) not g(r) (cosmetic, being fixed).
09:30: cs IMPLICIT TRANSVERSE DONE on cs-implicit-transverse (2e7b30a4, 355015e6, 19fbff67,
418f0bb3; docs/dev/cs_implicit_transverse.md; test tst/test_suite/rad/test_rad_cs_implicit_ang_cpu.py):
diagonal metric implicit, cross term g^{23} explicit (|cos a| <= 1/2 -> unconditionally stable),
seam faces as isolated two-cell backward-Euler pairs (rad_adi_seam_w 0.5), new rad_adi_cross_iter.
Cartesian BITWISE (incl. box_w8 production config). cs rate error 0.3 % (l=2) / 1.3 % (l=6) at 32,
conservation 1e-6 (seam resample), stiff z=7e5 ADI monotone; RKL1 (sts) BLOWS UP when its
substage cap clamps -> on the cs use rad_ang_solver = adi. LESSON: dividing by an exact 1.0 or
hoisting into a branch changes FMA contraction under hipcc and breaks bitwise tests; append
curvilinear forms as overwrites. Not yet merged into he4-presn-global (waiting for the mode-3 fix).
10:00: mode-3 spherical "bug" was the ABSOLUTE-area scaling (J = A I with A ~ 1e22 next to unit
entries -> single-precision pivot test failed on every block, mixed precision silently disabled);
fixed bbcb06a0 (areas relative to the column top). Merge of cs-implicit-transverse into
he4-presn-global CLEAN (2c59e35f); input on the ADI (e24f9932). rt_rad_force transverse bug
fixed 4bcdc855 (see rt-rad-force-transverse-cs-bug). 1-D He GATE IS IMPOSSIBLE: a super-Eddington
Prad-dominated FeCZ has no hydrostatic radiative-equilibrium state (needs 15-40 % convective flux);
even with the base flux injected the column runs away radially inside the FeCZ at 0.47 turnover.
DECISION (mine, user asleep): skip the 1-D gate, run the 3-D smoke arms from the MLT ic
(A: inner_bc open; B: wall + bottom flux; C: survivor with vpert 1e-2). tests_3d/mk_ic_from_profile.py
builds an ic from rt_profile.bin if a relaxed state is ever wanted. box_w9 deleted (user).
11:00: 3-D SMOKE ARMS A (inner open) and B (wall + base flux) BOTH DIE (0.12 / 0.37 turnover,
tests_3d/arms/README.md): the upper FeCZ DRAINS from t=0 (L_rad,out/L 0.94 -> 0.66 by 0.06
turnover, rho at 0.97 R -43 %) because the two-stream sees div F_rad != 0 where the MLT column had
15-40 % of L convective; efloor cascade in the emptied top; convection never starts (KEh/KE1 0.3 %).
Not resolution-limited (eos_fail 0, fofc 0, vertices silent). Open inner BC as configured = mass
source + energy sink -> use WALL + rt_bottom_flux. cs-seam-4x4 merged (86aff7f9). NEXT (running):
carry the convective flux from t=0 with red_giant's MLT SUB-GRID FLUX (mlt_alpha/mlt_tau/ramp) and
ramp it off over turnovers 1-4 while vpert 1e-2 convection grows (arms D/E). NO production launch
without the user.
11:30: MLT SUB-GRID FLUX arm D FAILS the same way (tests_3d/mlt/README.md; new switches
mlt_ramp_down_time/mlt_hold_time f4daa23c; NOTE vpert_mlt defaults true, gated on mlt_alpha_ic).
DECISIVE MEASUREMENT (mlt_dump t=0 face budget, D/mltfaces_D.txt): Rosseland diffusion on the
shell-mean state carries 0.85-1.07 F_req at every face; the grey TWO-STREAM carries only
0.11-0.74 over 0.55-0.90 R (0.97 at 0.968 R) -> 2-9x short in the deep FeCZ = where every arm
collapses. The J = A I spherical scheme (1159a8f3) is exact in the transparent limit but in the
DIFFUSION limit adds a spurious ~2B/(r dB/dr) term (O(1) where H_T ~ 0.3-0.6 R). So "two-stream
all along" (my 09:00 answer) is WRONG on a sphere spanning half the star: hand the deep interior
to radiative diffusion (rad_tau_lo/hi ~ 20/300 as in dhj/rg), two-stream only in the thin layers.
Agent running: analytic thick-shell unit test of the two-stream on a sphere + handover budget
scan + arm B/E with the handover.
12:30: handover budget scan committed (6d918992) then REVERTED (fe429a58) because arm BH
(20/300 + wall + rad_flux_inner, mlt_alpha 0) died at 0.26 turnover at the FeCZ BASE (0.633 R),
worse than arm B: the live conduction feeds the base while the FeCZ interior (0.72-0.90 R, the
unstable shells, Rosseland 0.85 of F_req) has no carrier. Meter caveats: mlt_dump fires after the
Strang dt/2 pre-step (F_2s and T not simultaneous); deep T drifts 0.4 % between configs.
PENDING: arm DH (20/300 + mlt_alpha 1.5 + vpert 1e-3, job 11760307) = the only coherent config;
if it clears 1 turnover -> restore 6d918992 values and continue smoke; if it dies -> the
variable-Eddington moment rewrite of the spherical two-stream is LOAD-BEARING. Thick-shell
analytic test (tests_r2/thick, job 11760190) also pending. Input is at the arm-A/B values.
