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
21:00 (tests_r5, e247210d): THE DRAIN IS PHYSICAL INFLATION. Gamma = kappa L/(4 pi c G M) = 1.46 kappa
> 1 over 0.66-0.94 R (max 1.20 at 0.82 R) on radiation alone -> with mlt_alpha 0 the layer expands,
sweeps a dense shell at 0.92-0.96 R (rho x3), the top atmosphere falls onto it; mass through the top
only 8e-4/turnover; velocities 0.08 c_s. With mlt_alpha 1.5 (arm D'): Gamma_rad 1.00-1.05, FeCZ mass
-9 % (was -36 %), no dt collapse, L_out/L 0.94 at 1 turnover -> a SETTLING inflated state (residual
Gamma_rad - 1 ~ 0.03, ~30-turnover timescale). RETRACTION: the "EOS garbage T = 5e11" was wtemp in
code units (p/rho at mu=1): T = 6310 K = eos_logt_min clamp, a CONSISTENT neutral-He state. Real
defects: tfloor 5e3 K BELOW the table edge (unreachable; now raised to the table's lowest T at
startup) and dfloor 1e-13 = 4e-5 of the FeCZ rho (cells empty then v = m/dfloor -> 212 c). New
switch efloor_as_tfloor (default off; floors set e = e(rho, max(T, tfloor)) on the tabulated EOS;
counter eos_tset, appended last). Box bitwise. NOTE the pgen face_budget g in/out numbers are ~1e11
too large (unit bug) - use hst mass. NEXT: dfloor_keep_velocity true + vceil ~1e8, then D' 3-5
turnovers (gate Gamma_rad -> 1, Mdot(0.9 R) -> 0), then resolved convection takes over the closure.
22:00 (tests_r6, 5b71b144): dfloor_keep_velocity + vceil 1e8 in the input: null on the structure,
eos_vceil 0 in the healthy phase, converts the 212 c cell into one pinned at vceil - BUT the star
still dies at 1.04-1.08 turnover in EVERY arm (seed 1e-3/1e-2, vceil 1e8/5e7, floors on/off): a
DETERMINISTIC THERMAL RUNAWAY of one cell at r/R 0.79 (i=36) at t ~ 5.0e3 s, T 1.3-2.3e6 K
(ambient 1.7e5), dt set by its c_s. Not settling: Mdot(0.9 R) still rising (2.8e21 g/s), Gamma at
the kappa peak 1.15 -> 1.25 (mass piles onto the rising side of the Fe bump), F_res/F_req = 0.
NEW BLOCKER 2: the MLT closure's relaxed profile fmlt1d (+ seeded flag) is pgen state, not restart
state: a chained restart resets F_mlt/F_req 0.10 -> 0.71 and kicks the star (KE1 x2.5-4.6). Must be
written to the restart (or mlt_relax_time 0). NEXT: reproduce the 0.79 R heating in a 1-D/thin
wedge (subgrid flux divergence at that face? two-stream/conduction blend? opacity table on the
rising Fe-bump side?) and add fmlt1d to the restart.
23:30 (tests_r7, 52d16d6b + 65cd5bd1): THE 0.79 R KILLER IS THE RADIATIVE FORCE. Reproduces in
1-D (1.18 turnover); shell means flat; the collapsing cell is one cell 4 orders above its shell,
ALWAYS on a panel-edge row (11/11); rank count picks the cell. Bisection: only rt_rad_force=false
removes it (alive > 2.4 turnovers); mlt_alpha 0 / mlt_mean false (WORSE) / strang / uform / fofc /
dc / nx1 64 only shift it. Opacity feedback is STABILISING (cell sits on the Fe-bump maximum).
Mechanism: the EOS/force taper weight w = w(rho) only (eos_rad_rho_hi/lo 1.93e-9/4.97e-10 sized at
tau 3/0.3 on the IC); the FeCZ density inversion min (2.3e-9) sits at the window, so a cell thinned
by a seam-row flux error gets the full thin-region force (1-w) rho kappa F/c at kappa F/(cg) 1.23
inside the CZ and evacuates. FIX TO DO: gate w on tau or radius as well (consistently in the EOS
and the force), or cap the force at g as a probe. RESTART: new optional marked pgen state block
("PGENST01", src/pgen/pgen.hpp, restart.cpp STEP 3); red_giant stores fmlt1d + seeded; chained
1-KE continuity 1.65e-1 -> 1e-15. Pre-existing: restart.cpp writes ~21 uninitialised bytes of
RegionIndcs (cnx1..cke unset on uniform grids) -> zero-init. Do NOT pursue per-column MLT.
00:30 09-18 (tests_r8, a6d66c3c): TAPER GATED ON TEMPERATURE: w_eff = max(w_rho, w_T) with
eos_rad_t_hi/lo = T at tau 3/0.3 (6.33e4/4.52e4 K); T is an argument at all ~150 EOS sites so it
is automatically consistent (a radius/tau gate would need an EOS interface change: follow-up).
rad_taper::WeightGated; dw/dlnT terms in chi_T, c_v, the root-find derivative and the force
(P_rad dw/dlnT grad T/T). ic rebuilt (nearly a no-op: |dw| 6e-4). Box bitwise (except 2 new
default lines in the parameter dump). RegionIndcs zero-init -> restart files byte-identical.
RESULT: 1-D x_gate alive 2.26+ turnovers (x_gate5 long arm running, 2.5 turnovers at dt 0.075);
3-D Gnr8 1.53 turnovers, no collapse, Gamma at the kappa peak flat (1.154 -> 1.168), Gamma_rad
at 0.8 R 1.007 -> 1.000, Mdot(0.9 R) turned over, restart kick gone. Convection NOT started
(F_res = 0, vr_rms/v_MLT 0.02). THIRD MECHANISM: a 0.97 R shell piles up (+47 %/turnover, rho
+134 % by 1.5 turnovers, local Gamma 1.05, dt -> 0.3 s), force-independent (the noforce control
drains/heats the same shell); it is where two-stream, MLT closure and taper all hand over.
NEXT (1-D, cheap): mlt_alpha 0 and rt_bottom_flux false vs the gated baseline; rt_profile eint
slot is NOT eint (1.02x at 0.5 R, 80x at 1.01 R) - an6 energy columns untrustworthy.
01:00 09-18 (tests_r9, e73c6add; handover docs/handover/HANDOVER-2026-09-18.md 4ff8bce8): the
0.97 R pile-up = EMERGENT-LUMINOSITY SHORTFALL: L_out/L relaxes to 0.956-0.981 while the domain
thermal time is 0.46 turnover -> envelope banks 4-9 %/turnover (+20 % by 2.5), lifts 0.94-0.99 R,
dt -> 0.08 s; six input arms (mlt 0, wall flux, top wall, no sponge, taper tau 1/0.1, rad_tau
20/300, nx1 144) change nothing. The shortfall is in the last 5 thin cells + the TOP FACE
(F_2s/F_req 1.28 at t=0) -> suspect the spherical sweep's top boundary (rt_top_re off: ghost
mirrors the top cell), not the interior form (fixed). Also the ic vs solver disagree 28 % at the
top face. rt_profile slots 5/6 now T[K] and eint. NEXT: 1-D grey atmosphere test incl. the top
face gating L_out/L = 1.000; fix the top-face treatment; rebuild the ic against the solver; then
3-D 3-5 turnovers. OPUS QUOTA EXHAUSTED until 2026-09-19 23:00 (Brussels).
02:00 09-18 (tests_r10, 1ebcbc3a, Sonnet, measurement only): the spherical sweep's TOP FACE is
wrong in the thin layers: +20-35 % excess at t ~ 0.1 t_therm(top cell) in every thick-shell case
(thin or thick top); run long, cases with dtau_top < 0.1 RUN AWAY in mode 3 (L_out/L 5-77x),
thicker tops decay to a 0.4-0.7x deficit; the He4 star decays to 0.978. Plane-parallel box:
F_top/F_imposed 1.004 from cycle 0 -> the defect is the spherical top boundary / thin-layer
treatment. rt_top_re=true fatals under mode 3; rt_top_vacuum is not read by red_giant. NEXT (Opus,
after 09-19 23:00): fix the top-face condition of the spherical sweep (mode 3 FormalSph + mode 0)
so that a thin outer layer emits L/(4 pi r^2); gate = tests_r10 matrix L_out/L = 1.000 at t = 0
and after relaxation, box bitwise.
04:00 09-18 (tests_r11, Fable): rt_top_vacuum exposed in red_giant: a NULL on the He4 star
(dtau above ~6e-3). mlt_alpha 3.0 is WORSE (L_out 0.36 at 1.6 turnovers, outer shells fall in and
heat to 2e5 K). Reading: the 2-3 % L shortfall = expansion/stored-radiation power of a layer at
Gamma_rad ~ 1; only resolved convection can hold it. SUBMITTED 3-D smoke run to 5 turnovers in
one piece (job 11793509, tests_r11/g3d, 4 h apu, rst every half turnover): check alive, L_out/L,
lnKEh vs lnKE1, F_res/F_req, v_r/v_MLT. The tests_r10 thick-top/thin-top numbers are the
unresolved-column model with gas pressure (separate item, not the star).
09:30 09-18: 3-D g3d run (tests_r11/g3d, one piece) reached 2.2 turnovers, no collapse, horizontal
KE growing x30/turnover (convection starting), L_out/L 0.97; dt eroded to 0.07 s. CAUSE (dump 3 at
1.5 turnovers): ~160 cells (1 % of columns, NOT seam rows) at 0.96-0.99 R just above the swept
dense shell are VOIDS: rho 1e-12..1e-11 (median 1e-9), e/rho up to 4e17 (T ~ 1e7 K), v1 up to
6e7 outward -> c_s sets dt. Mechanism: the shell pushes up, the low-density gas above is expelled
(super-Eddington), the near-empty cell then heats under the implicit RT deposit (rt_de_max caps
only EXPLICIT updates). ARMS SUBMITTED (2 turnovers, 3.5 h apu, rst every 0.2 turnover, script
tests_r11/r11_3d_arm.sh): d1 dfloor 1e-11 (job 11798237), d2 dfloor 1e-10 (11798238); IC top rho
is 5.5e-11 so 1e-10 is above it (check the top cells). Gate: dt not eroding below ~1 s, no hot
voids, KEh growth, L_out/L.
