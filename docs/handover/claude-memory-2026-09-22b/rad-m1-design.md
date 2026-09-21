---
name: rad-m1-design
description: "Grey photon M1 module plan (2026-09-21): design note docs/dev/rad_m1_design.md on rt-integration (uncommitted at writing); key decisions and open risks"
metadata: 
  node_type: memory
  type: project
  originSessionId: e67ee046-aa0b-46f4-881b-84099bfda99a
  modified: 2026-09-20T23:28:58.033Z
---

User wants to replace the two-stream + transverse conduction stitch with an evolved M1
radiation field; agreed order: explicit Cartesian M1 first, implicit later. Design note written
09-21: docs/dev/rad_m1_design.md (no code yet).

Decisions in the note: mixed frame, sources evaluated in the COMOVING frame and transformed
(G0 = g0 + beta.g, G = g + beta g0, with E0 and F0/c = F/c - beta E - beta.P) = AREPO-IDORT
(user's own paper, arXiv:2503.16627) eqs 56-57 = Krumholz 2007 coefficient (kappa_R - 2 kappa_P);
the beta^2 E terms are leading order at beta*tau > 1 (O(v/c) truncations of Skinner-Ostriker,
Rosdahl-Teyssier, AREPO-RT, QUOKKA I give spurious heating). Implicit solve = scalar root find in
T with the table c_v, T^4 not linearised, algebraic energy conservation; F linear backward Euler.
Thick-limit flux (revised cde8940e after reading Berthon-Turpault 2011, Bloch 2021, Jiang 2021, QUOKKA 2404.08247): default ap_hll = alpha on the E-FLUX ONLY (Bloch variant; alpha on c^2 P would make the discrete radiation force alpha*grad P_rad), compact stencil, ARITHMETIC face opacity; thick_flux=scaled = Jiang/IDORT modified HLLE in moment form (both speeds scaled, fixed 4-20 % excess diffusion, not AP); uncorrected HLL (QUOKKA) rejected: tau_cell x excess flux at limiter-clipped extrema (own assessment, test T3c). Moving fluid: upwinded enthalpy-flux split vE+v.P in the E equation only (own construction; IDORT found Jiang's split unstable with local time-stepping). Time integration: IMEX PD-ARS with the source solve INSIDE both stages; transport-then-source splitting gives spurious c^2 dt/3 diffusion.
RSLA only valid for chat >> v*tau_max -> useless in our deep cells; explicit M1 is a validation
vehicle, production needs implicit transport (column + ADI line solves as preconditioner).
EOS must be gas-only (eos_radiation=false, no taper; rt_rad_force path goes away). No Planck-mean
opacity exists in the tree. Upstream origin/project/ccsn src/radiation_m1: reuse tasks/bvals
skeleton only; do not copy its root finders (GSL/GPL-looking); it has a dx1 bug in x2/x3 A-factor.

**Why:** porosity/force-consistency problems of the He4 envelope come from the stitched operators.
**How to apply:** follow the note's gates T1-T9 in order; open risks: R1 admissibility of the E-only alpha variant, R2 the advective split (gates T4/T4b), R4 Planck-mean table source; R3 (PD-ARS) closed. Committed 089a1e7c + cde8940e on rt-integration (unpushed). See [[use-rt-integration-branch]], [[he4-porosity-luminosity-excess]].

09-21 later: milestone 1a committed 6b9a14f9 (skeleton + free-streaming transport; box G1
bitwise, beam FWHM 11.4 vs exact 11.3 cells, pulse order 1.63-1.76, 4 ranks bitwise); gate
scripts tests_m1/ b8d8112d (Su-Olson reference is the script's own S_N solve, published table
not sourced). R4 CLOSED: Planck-mean tables in bench/m1_opac/ (TOPS/ATOMIC, scriptable via
fetch_tops.py; planck_gs98_x0.7_z0.014.txt, planck_he_x0.0_z0.02.txt in the repo table_rho
format; Ferguson 2005 low-T Planck (user's pointer) spliced below log T 4.2 in *_ferg+tops.txt;
TOPS Rosseland agrees with our tables to median 0.005 dex). kappa_P/kappa_R ~ 100-1000 at
tau <= 1 (108 at tau=1 on the He column), 40-70 at the Fe bump; TOPS Planck is absorption-only
(no scattering). Published practice: Jiang 2015/2017 and Moens 2022 set kappa_P = kappa_E =
kappa_R - kappa_es; Jiang 2018 (Nature) and Goldberg 2022 SAY they use an OPAL Planck table, but
OPAL ships Rosseland only. User 09-21: they probably also used OPLIB/TOPS for the Planck mean; not crucial, do not pursue.

SIZING 09-21 (bench/m1_sizing/sizing.py, reproduces measured production dt): explicit M1 + RSLA
IS affordable on the He BOX (box_w8: tau_bot 322, v_max 8.6e4, FeCZ = tau 4.4-64): first
comparison = He box cut at tau_bot 100 (84 cells), chat/c 2.9e-3 (K=10 in chat = K v tau),
N_sub 8.8, ~13 node-h per 10 turnovers = same cost as the two-stream; full He box chat/c 9e-3,
N_sub 28, 58 node-h. K=3 is marginal (t_diff/t_flow = 3/K); prefer K >= 10-30. B star NOT
affordable explicitly: FeCZ at tau 579-3268 needs tau_bot 1e4 where chat = c, ~2900 node-h ->
B star (and the global He4 model) need implicit transport. My earlier "RSLA useless" statement
holds for the global model and the B star, not for the He box.

STAGE 2 PLAN 09-21: bench/m1_stage2/PLAN.md (+check_ic.py). Key facts: (1) WB scheme walks a
hydrostatic background whose p includes w aT^4/3; with gas-only EOS an operator-split M1 force
leaves a per-step excursion up to 18 v_MLT (|grad P_rad|/(rho g) up to 2.85 at the Fe bump) ->
need a second effective potential Phi_eff = int(g - kappa_R F/c)dz for the WB walk/source only
(~145 lines, default off). (2) box_w8 column is hydrostatic under the TAPERED EOS: w=0 at the
box top, 0.81 at tau 2/3 -> the top ~15 cells have no radiation pressure today; M1 IC must be
rebuilt (ic_hse_retune=2). (3) x1 BCs are `user` -> BoxConvBC needs an M1 branch. (4) opacity
set: kappa_F+kappa_s := kappa_R(total), kappa_s=0, kappa_P=kappa_E=Planck table (15-527x kappa_R
in this box); table lookup wants KELVIN, code T unit = 1.202724e-8 K; arad(code)=1.583122e-46.
(5) keep the Conduction object alive (tables, rad_tauf, rad_flux_inner) with blend weight 0;
box_convection has NO MLT closure. Cost 12.9 node-h/10 turnovers vs 11.9 two-stream.
1b committed f271d57d (ap_hll rate 1.000015 at tau_cell 1e3; findings in design note sect. 10).

RESULTS PAGE: https://claude.ai/artifact/6b5dXN5p9ESoNZkNh6rxs8 ("M1 Gate Ledger"; source
scratchpad m1_gates.html of session 09-21; numbers only, no profile plots yet because dumps are
deleted for the inode quota; user asked for test/convergence plots -> add profile figures after
1c-A). MULTIGRID PORT: branch mg-port (bench/wt_mgport, cb39c07c) = upstream 8a6a8efa cherry-picked
onto rt-integration, all gates passed (box G1 10/10 bitwise, poisson3d/jeans3d/MPI-4 pass, rad_m1
beam ok); grav_phi output index moved to 160 (m1_* keep 155-159; new m1 outputs must stay < 160 or
bump the guard); to be merged into rt-integration after 1c-A is committed. Other branches in
flight: wb-phieff (bench/wt_wbphi), IC builder bench/m1_stage2/ic.

M1 HE-BOX IC 09-21 (bench/m1_stage2/ic/build_ic.py; use ic_m1_V2_pgen.txt + arad_V2.txt, with
ic_hse_retune = 0): V2 = fully self-consistent radiative-equilibrium + hydrostatic column on the
production grid, cut at face j=50 (tau 99.0, 84 cells). Findings: (1) CORRECTS the stage-2 PLAN:
a_rad/g peaks at 0.75 (tau~30), NO super-Eddington layer and no density inversion; the "2.85" was
the gradient of the TAPER weight, not a force. (2) production column is already in radiative
equilibrium (F_rad/F 0.97-1.03 in the FeCZ), FeCZ extent unchanged (tau 4.5-64.0 for V2 vs 4.65-64.0 production; the agent's 10.8-54.4 and its 'gas-only grad_ad = 0.5' came from a factor-2 bug in build_ic.py's log derivatives, fixed by me 09-21: gas-only grad_ad = 0.397-0.400, mixture 0.250-0.268; IC files unchanged by the fix); top differs:
P_rad(top) 7.7e4 vs 0 today, T_top 68400 vs 65975 K. (3) box EOS table is built at run time from
eos_composition.hpp (analytic Saha) + w aT^4. (4) the box column reader takes `z rho eint` only.
(5) vacuum BC is a plain copy for outgoing flux: the top closure f(0)=1/2 is set by the IC.
(6) my own caution: with S&O-type RSLA the conserved energy is E_gas + (c/chat)E, i.e. radiation's
effective heat capacity is inflated by c/chat unless dE/dt is quasi-static (same condition
chat >> v tau, error ~3/K) -> in the radiation-dominated FeCZ (beta 0.09-0.6) run K = 10/30/100.
MILESTONES: 1c-A 446918e0 (split works: lag of half the advection distance without it; source form
confirmed, O(v/c) control runs away; alpha2 stays default; T3b 1.4 % flux deficit at tau_cell~1
open); multigrid port merged 7861549c.

WB EFFECTIVE POTENTIAL 09-21: branch wb-phieff (bench/wt_wbphi, 2 commits on 65a0569f, NOT yet
merged: merge after 1c-B is committed). <problem>/wb_phi_eff + wb_arad_file + wb_arad_force;
Hydro::phicc_wb / phi_wb_x1f default to shallow copies of the true potential (bitwise inert, G1
10/10). 1-D adiabatic He column, 35 sound crossings: fix 0.46 v_MLT, unfixed split 10.6 v_MLT
(rho drifts 40 %). KEY for the stage-2 wiring: once the WB pair uses Phi_eff it already delivers
rho*a_rad_ref, so the M1 coupling must apply only the RESIDUAL rho*(kappa_R F/c - a_rad_ref) to the
momentum (else double counting); its work term goes explicitly into the total energy (a_rad is not
in the etotgrav potential). Trap: Kokkos::realloc keeps the allocation when extents match (aliased
phicc0). Restart of this 1-D config is not bitwise even in the control (4.7e-13): pre-existing.
Untested: spherical path, non-zero residual with a live M1.

STATE 09-21 end of stage 1: rt-integration HEAD 2d5c6ba1 (unpushed beyond cde8940e): 1c-B aea14fc0
(restart bitwise, reconstruct dc|plm|ppm4|ppmx|wenoz with plm kept default -- wenoz order 2.0 free
streaming but slower in the diffusion limit; split_vel=recon default; 1/2/4 blocks and 1/4 ranks
bitwise; conservation 2e-13/200 cycles), multigrid port merged 7861549c, wb-phieff merged 2d5c6ba1
(both worktrees removed). Results page v3 has profile plots (source copy:
bench/m1_stage2/m1_gates_page.html; data tests_m1/plots/*.json). Open: T3b 1.4 % flux deficit at a
1e3 jump at tau_cell~1; ParameterInput records defaulted Reals with 6 digits -> pgen state derived
from a defaulted Real is not bitwise restartable (code-wide); stage 1 gates T2 (shadow), T7, T9 not
run. NEXT: wire <rad_m1> into box_convection (BC branch, table opacities with kelvin conversion,
residual force vs a_rad_ref, diagnostics), 1-D He column relaxation with K = 10/30/100; stage 3
design note; handover doc; RKL2 port (#779).

STAGE 2a 09-21 late: committed 44c8cf55 (M1 in box_convection: table opacities, BC branch, residual
force vs a_rad_ref, M1 IC file, history, RSLA check; fixed: etotgrav rho*Phi read as internal
energy, Conduction injecting rad_flux_inner twice, bottom ghost reading a zero opacity array).
1-D He column (84 cells, K=10, sponges off, closed walls) is NOT quiet: a growing OSCILLATION,
period ~72 s, e-folding ~195 s (I checked the hst: oscillatory, so NOT convection as the agent
guessed, impossible in 1-D, and not a monotone drain). Suspects: physical kappa/strange mode
(L/L_Edd 0.75), RSLA phase lag (c/chat), split lag (cf. fmode-dt-taper-test, fmode-wb-cache-culprit),
frozen a_rad_ref, closed lid. Diagnostic agent running (K scan, dt scan, frozen opacities, work
integral, sponges, two-stream twin) -> tests_m1/runs_2a_diag/. The module's `vacuum` fill (copy
when flux points out) lets a diffuse surface float (emergent flux 0.36 F_in): use dark ghost;
being changed to always-dark. Other branches in flight: pin-real-precision (bench/wt_pinreal),
m1-open-gates (bench/wt_m1gates: shadow, radiative shocks, T3b remedy f_source=wb).

PULSATION DIAGNOSED 09-21 (commit after 13cf9772; tests_m1/runs_2a_diag/RESULTS.txt): the 1-D He
column's growing ~65-75 s mode is a KAPPA-MECHANISM AT THE Fe BUMP ENABLED BY THE REDUCED SPEED OF
LIGHT: gamma = 5.7e-3/8.5e-3/4.2e-3/6.3e-4/5e-5 per s at K = 10/30/100/300/1000, -> 0 as chat -> c
(true-c diffusion time 2.2 s = 0.03 P); independent of hydro dt (not split lag); frozen opacities
kill it; needs force AND heat exchange; sponges remove only 4 %; WB Phi_eff only reduces the kick.
Work integral: driving tau 53->1.6, peak tau~21; dF lags drho by 117 deg. LESSON: RSLA validity in a
radiation-dominated zone is t_diff(chat) << pulsation period, K >~ 300-1000 (K=300: ~390 node-h per
10 turnovers in 3-D), NOT chat >> v*tau with K~10 -> explicit+RSLA is not the way for the He box
either; implicit true-c (stage 3a, branch m1-implicit-x1) is the path. Vacuum BC now always dark.

RESTART FIXES merged 09-21: pin-real-precision 13cf9772 (GetOrAddReal/SetReal wrote 6 digits) and
box-restart-bitwise 6ae9ccbe: (1) tabulated-EOS T inversion is NOT idempotent + Driver::Initialize
runs an extra ConsToPrim -> restart file now carries wder behind wtemp and the first ConToPrim
reads the caches (ConsToPrimFrozen); (2) rt_surface/rt_profile dump clocks now in the pgen state
block; (3) NOT fixed: rt_impl_warm=1 (PRODUCTION setting) keeps per-cell Newton warm-start history
outside the restart file -> mode 3 restarts bitwise only with rt_impl_warm=0. Same wder loss exists
for general-EOS MHD. INODES: each git worktree of this repo costs ~8.5k inodes (tests_* dirs);
the GPFS hard limit was hit once with 4 worktrees: keep <= 2-3 worktrees, remove when merged.

STAGE 3a MERGED 09-21 (f2d4bf47; design docs/dev/rad_m1_implicit_design.md incl. findings;
tests_m1/runs_3a/RESULTS.txt): transport = implicit_x1 = face-eliminated scalar E tridiagonal
(M-matrix) column solve at TRUE c, backward Euler, Picard. Thick-pulse rate within 0.1 % for CFL
0.4..1e4; T3b opacity-jump flux error 1.31 % -> 1.9e-3 %; equilibration 1e-12; restart + 1v4
blocks bitwise. **1-D He column at true c does NOT pulsate** (gamma -2.9e-5/s; F1/F_in
0.9968-1.0012; residual 10 v_MLT sits in the bottom cell at the imposed-flux BC, period ~50 s,
constant amplitude) at HALF the cost of explicit K=10 (4.7e-3 vs 9.5e-3 s/step; Picard 13.5 mean).
Limits: free streaming heavily damped even at CFL<1 (peak 0.37 after one crossing; explicit 0.95)
-> not for thin regions; Marshak 1.4 % vs 0.19 %; one MeshBlock per x1 column; F2=F3=0; MPI/GPU
untested; needs efix BC for pure-scattering problems (singular otherwise). Bugs fixed on the way:
Sherman-Morrison corner entries swapped; gas energy must come from the row's applied source.
NEXT: 3b (3-D Krylov + this line solve as preconditioner; partitioned line across MeshBlocks),
then He box 2-D/3-D at true c vs the two-stream twin.

09-21 night: m1-open-gates merged 7a2c9382 (post-merge gate tests_gate_merge/postmerge.sh: box G1
bitwise + implicit He column): T2 shadow front 1.005c pass, shadow depth 7.7 % (gate 1e-3 NOT
met; Eddington 89 %); T7 Mach 2/5 L1 0.25-1.0 % pass (Eddington closure 8-100x better there;
steady states independent of chat); T3b: most of the 1.4 % was the TEST's fixed-E ghost
(exact_ghost), the cell-F mechanism is real and <rad_m1>/f_source = wb fixes it (0.0018/0.049 %;
recommended default for stage 2, currently off). IN FLIGHT: m1-implicit-2 (bench/wt_m1impl2):
implicit ap_hll blend with lagged f (thin regions/Marshak front), steady atmosphere gate, bottom
imposed-flux BC oscillation, partitioned line solve + MPI, GPU run; rt-warm-restart
(bench/wt_warmrst). After that: 3b transverse transport (user agreed to address the limits;
my proposal = column solve + implicit ADI transverse lines on E first, Krylov later).
rt-warm-restart merged (after 7a2c9382): warm-start history now in the restart file (RTWARM01
block); He box_w8 and B-star prod_w7 configurations restart BITWISE with production settings
(warm 0/1/2, 1/4 blocks, 4 ranks). Gates: tests_gate_merge/postmerge.sh + postmerge_rst.sh.

3a2 MERGED a766ef98 (09-21 night): implicit_bmom_half (boundary-face force shared correctly):
He column max|v| 10.3 -> 0.99 v_MLT (verified by me at 300 s), fluxes within 1.7e-4; steady M1
atmosphere gate (t9_atmosphere.py): CENTRAL face-flux form within 1.35e-3 of the closed form
incl. top cells at every CFL, all upwind forms ~40x worse, EXPLICIT scheme 49 % off there ->
central stays default for quasi-static problems; implicit_flux=berthon only for propagating
fronts (pulse peak 0.82/0.70 at CFL 0.4/1); the literal alpha F_HLL+(1-alpha)F_diff blend double
counts (my spec error, second time); Marshak 1.0 % best (explicit 0.19 %). NOT done there:
partitioned x1 line solve, GPU. IN FLIGHT: stage 3b agent in the MAIN checkout (uncommitted):
phase 0 bmom default + He column restart, phase 1 gather-partitioned line solve + MPI, phase 2
transverse transport (line_jacobi in the Picard loop, then BiCGStab with x1-line preconditioner),
phase 3 He box 2-D/3-D smoke test with the entropy seed, phase 4 GPU. Logs tests_m1/runs_3b/.

USER 09-21 late night (going to bed, left me to continue autonomously; nothing to be pushed, no
production jobs without approval): (1) fix the remaining issues of the implicit scheme; (2) idea:
blend smoothly from the central form to Berthon based on local tau like Jiang 2021
(w = 1 - exp(-tau^2)); my note: the atmosphere gate shows central wins even at tau_cell < 1 when
the field is quasi-isotropic (f ~ 0.5), Berthon wins only for beam-like fronts (f -> 1), so the
blend weight should depend on the reduced flux f (and tau), not tau alone: test both. To do after
stage 3b returns (agent in the main checkout): gates = atmosphere (keep 1.35e-3), free-streaming
pulse, Marshak (target 0.5 %), T3b, thick pulse, He column. Results page v4 published (He column,
growth rate vs K, implicit gates, shadow, shocks).

09-22 early: 3b phases 0-1 committed 441762eb (implicit_partition=gather: x1 line solve across
blocks/ranks BITWISE; global Picard test; bmom_half default; implicit He column restart bitwise;
FOUND: box_convection not decomposition-invariant along x1 (3e-4 even explicit) and its M1 hst
columns are summed over MeshBlocks). 3c committed (next hash): implicit_flux=blend. USER'S tau-only
weight FAILS (atmosphere 3.7e-1, He column x4.4, no Picard convergence: top of an atmosphere is
thin AND diffuse); implicit_blend=tau_f (thin AND beamed, f smoothstep 0.6-0.9, fmode=max) is
bitwise central on all diffuse gates and berthon in free streaming (pulse 0.92/0.72/0.79 with
plm_dc), repairs plm_dc; Marshak not improved -> DEFAULT STAYS central. NEXT: transverse transport
(3b phases 2-4) + the two box fixes.

09-22 night (user asleep; nothing pushed; fork still cde8940e): committed d224ba5c (3c blend),
5ccb298a (box M1 history = GLOBAL planes; x1-decomposition non-invariance is cell-centre round-off,
no bug), 183ff397 (transport = implicit: 3-D face-eliminated solve, line Jacobi, 7-point M-matrix,
halo via MeshBoundaryValuesCC, true-residual test; 2-D/3-D pulse rates 0.9955-1.0000, isotropy
1.7e-9, 1 vs 2x2 blocks bitwise, restart byte-identical; pass counts scale with the transverse
diffusion CFL c dt/(3 tau_cell dx): 131 passes at CFL 1e4). IN FLIGHT (main checkout, uncommitted):
BiCGStab with x1-line preconditioner + the 2-D He slab smoke test (static, then production entropy
seed) -> tests_m1/runs_3b4/. LESSON for briefs: agents given multi-phase milestones burn the
budget on preliminaries; give ONE narrowly scoped deliverable, "code first, gates second".
Still to do after that: results page + handover update; GPU run; transverse Marshak BCs.

09-22 ~morning (user asleep): committed 21971877 (BiCGStab w/ x1-line preconditioner: 131-pass
cases -> 5-7 outer x 5-6 inner, 8-12x faster; STATIC 2-D He slab = 1-D column to 8 digits; line
Jacobi cannot do the slab) and a5c23db9 (implicit_offdiag=operator; implicit_closure_lag=step).
KEY: the seeded slab's Picard divergence is the CLOSURE iteration (chi,n), not the off-diagonal
terms; freezing the closure over the step converges in 4 passes and is 12x faster (costs 7th digit
on the static slab, 6.7e-2 on an oblique thin pulse). Recommended 2-D/3-D: bicgstab + operator +
closure_lag=step. BUT the seeded slab is still UNPHYSICAL (KE_2 grows 75x the convective rate,
F1top/Fin 0.39, positivity fallback nearly every step): the transverse gas-radiation coupling has
never been unit-tested. IN FLIGHT: radiation-modified acoustic wave test along x1 / x2 / 45 deg,
implicit_x1 vs implicit vs explicit (tests_m1/runs_3b6), then bug hunt, then the seeded slab again.
I stopped a lingering finished agent with TaskStop (agents that leave `sleep` waits keep waking).

09-21 afternoon (session crashed 13:54: per-user 108 GB cgroup cap on the viper login node was hit,
claude + a serial athena got SIGBUS; nothing lost). SEEDED SLAB DIAGNOSED (tests_m1/runs_3b7/
RESULTS.txt): the 0.16-0.19 /s is a ONE-OFF RAMP saturating by ~5 s, not a mode; 99 % of v^2 sits at
tau < 0.1 (thin top); cutting the box at tau = 1 (A5) removes every symptom (F1top/Fin 0.394 ->
0.9997, positivity fallback 1388/1393 -> 0, max|v| 740 -> 14 v_MLT). H1 (frozen a_rad_ref makes the
gas-only stratification unstable) is REAL but slow (max 0.0163 /s, hydro-only arm reaches 0.21 v_MLT
in 200 s) and NOT the driver; H2 (WB Phi_eff), H4 (seed consistency) refuted. New default-off pgen
switch problem/m1_seed_consistent. NEXT: cure the thin top (transverse flux form / positivity), then
a long hydro-only A1 (5-10 turnovers) for H1. Caveat for the transverse AP/HLL build: implicit_blend
= tau_f is bitwise central on DIFFUSE fields, and the thin top is diffuse, so a transverse blend must
be keyed on tau (harmless there: a plane-parallel state has zero transverse gradients).

09-21 evening (uncommitted in the main checkout; tests_m1/runs_3b8/RESULTS.txt + addendum): CONFIRMED
the seeded slab is quiet with transport=implicit_x1 (KE 2e25, F 1.000) -> transverse transport in
the thin top is the culprit; dumps show |F2|/(cE) = 1.000 and E fragmenting 1e2-1e5 horizontally in
the top 12 cells within 2 s. BUILT <rad_m1>/implicit_trans_limit = lp (+ implicit_trans_fmax): lagged
limiter opacity |div P|/(phi E) on transverse faces; off = bitwise, static slab bitwise, T10 passes;
PARTIAL cure (fmax 0.5: KE_2 1.7e29 -> 6.8e26, F1top/Fin 0.95-1.0, dt no longer erodes). Negative E /
positivity fallback = the OFF-DIAGONAL P_12 terms (offdiag=operator): with lp 0.5 + offdiag=none
0 fallbacks, Picard 6.3. My F^n-memory hypothesis REFUTED (dbg_trans_memory=0 changes nothing).
OPEN: KE_1 ~ 2e27 (100x x1-only), mostly vertical motion at tau < 10, cause unknown.
CLOSURES (docs/dev/rad_m1_closure_survey.md): exact Hopf surface f=0.577, chi=0.410; Levermore 0.512,
Minerbo 0.485, Kershaw 0.556, Eddington 0.333; uniform radiating sphere chi=(4f^2-2f+1)/3 as a
BENCHMARK only. USER: do NOT invent new closures; add published Minerbo + Kershaw as options.

09-21 night (committed 4e56dc97 = limiter + notes; addendum 2 in runs_3b8/RESULTS.txt uncommitted):
THE THIN-TOP INSTABILITY IS THE LAGGED M1 CLOSURE. Seeded 2-D slab with closure = eddington, no
limiter, offdiag = operator: completely clean (dE/E < 1e-6 at the top, max f2 1e-5, KE_2 5.8e21,
F1top/Fin 1.000, 0 fallbacks, Picard 4.1); with Levermore + lp limiter + offdiag none the top still
fragments (dE/E 0.4-0.5, 5 decades above the seed; the residual KE_1 is its consequence). The limiter
only caps it. Eddington's KE_1 2.7e26 is a coherent 1-D adjustment (IC built with Levermore): rebuild
the IC with chi = 1/3 (bench/m1_stage2/ic/build_ic.py). Note chi = 1/3 is also CLOSER to the exact
Hopf surface value (0.410) than Levermore (0.512). Open: whether Minerbo/Kershaw are stable here; a
stable treatment of a nonlinear closure at CFL 7e3 (implicit/Newton closure, or under-relaxed in time).
09-21 night, later: added closure = minerbo | kershaw (implicit only, uncommitted at writing). BOTH are
as unstable as Levermore on the seeded slab (F1top/Fin 0.36/0.43, fallbacks 200/206, top dE/E ~1);
only eddington is clean -> the instability belongs to the per-step LAGGED nonlinear closure at
c dt/dx ~ 7e3, not to a particular chi(f). COST serial CPU 84x32: hydro 2.1e5 zone-cycles/s, implicit
M1 (x1 or 2-D eddington) 2.0e4 = ~10x hydro; broken Levermore arms 4-5e3 (Picard 12-14).
COST PROFILE 09-21 (tests_m1/runs_3f_prof/PROFILE.txt; perf works on the login node, the Kokkos
kernel timer wrote nothing with the serial binary): implicit M1 2-D eddington = ~10x hydro; of the run
~32 % EOS table lookups in the per-cell T solve (every Picard pass), ~16 % the off-diagonal Krylov
operator (zero work for chi=1/3 -> skip), ~18 % halos (14-component transverse halo every pass),
only ~16 % the linear solver itself. So multigrid is NOT the first lever; first: cheap T solve
(rho frozen -> 1-D e(T) fit per cell per step / Newton with c_v), skip offdiag for eddington, halo
only what changes, Newton instead of Picard (4 -> ~2 passes). MG (src/multigrid, semi-coarsening in
x2/x3 + x1-line smoother) matters later for 3-D production scaling. In flight at writing: agents
runs_3d_edd (Eddington IC V3edd + slab runs) and runs_3e_newton (Anderson-accelerated closure).
OPTION 1 DONE 09-21 night (tests_m1/runs_3d_edd/RESULTS.txt, uncommitted): build_ic.py --closure
eddington -> V3edd files (bench/m1_stage2/ic/ic_m1_V3edd_pgen.txt, arad_V3edd.txt;
tests_m1/runs_3d_edd/m1_rad_ic_V3edd.txt; input he_slab_m1_2d_V3edd.athinput); top f = 1/2 (the
implicit top BC is F = c E/2), V3 vs V2: column |drho/rho| <= 2.8 %, FeCZ tau 4.36-64.27. Verified by
me from the hst: static slab KE_1(30 s) 2.7e26 (V2) -> 2.2e25 (V3edd) = the implicit_x1 level; SEEDED
200 s: F1top/Fin 1.0000 throughout, KE_2 6e22 -> 1.8e21 (DECAYS), max|v1| ~1 v_MLT (coherent 1-D
drift, shared with the 1-D arm), Picard 4.16/5 max, 0 fallbacks. 2-D implicit Eddington is CLEAN.
OPEN: the seed decays over 200 s (0.42 turnover) -- agent says grad_rad never exceeds the GAS-ONLY
grad_ad; conflicts with runs_3b7 H1 table (gas-only unstable at tau 14-45 under g_eff) -> check with a
longer run before believing either. USER approved: after the agents, optimisation items 1 (cheap T
solve), 3 (halo only what changes), 4 (Newton with local T elimination, iterated to convergence).
S3_long 09-21 night: seeded V3edd slab to 2400 s = 5 turnovers: robust (14894 solves, 0 failures,
F1top/Fin 1.0000, dt constant, 33 cpu-min serial). KE_2 decays to 200 s then grows exponentially,
d ln KE_2/dt = 1.2e-3 /s (velocity rate 6e-4 /s = 0.3 v_MLT/H_p, ~30x below the ideal-buoyancy
estimate: radiatively damped) -> the column IS unstable; onset will take tens of turnovers from a
1e-3 seed, like box_w8 (0.9 v_MLT at 13 turnovers). KE_1 1.1e26 flat = the ~1 v_MLT 1-D drift (open).
OPTION 3 FAILED 09-21 night (committed with option 1; tests_m1/runs_3e_newton/RESULTS.txt):
implicit_accel = anderson (default off, inert) cuts the static-slab per-pass closure iteration 80 -> 23
passes but does NOT converge Levermore in the thin top (every step maxit 200, residual in the top 12
cells, negative E every pass; m=5, m=10/beta 0.5, +lp limiter all fail). Remaining route to chi(f) in
thin cells = true Newton with the closure Jacobian on (E,F) or a positivity-preserving upwind transverse
flux: NOT started, discuss with the user first. PRODUCTION PATH NOW = closure = eddington + V3edd IC.
M1 CLOSED OUT 09-21 night (HEAD 86214be9, unpushed; handover docs/handover/HANDOVER-2026-09-22b.md):
speed-ups committed: implicit_gas_newton + implicit_eos_cache (8d09252b; slab 1.37x, 1-D column 5.3x;
pass count unchanged ~4) and shell-only halo copies (de7a52bb; bitwise on 8 arms; slab 3.96e4
zone-cycles/s = ~5x hydro). VET SCAFFOLDING (runs_3i_tensor): tensor frozen from the IC, tilted
(prescribed off-diagonals) and rebuilt EVERY STEP from the column optical depth: all clean (Picard
4.0-4.7, 0 fallbacks, F 1.0000, dt constant) -> the thin-top instability is specific to (chi,n) from
the cell's own F; a short-characteristics VET tensor should be safe. USER: the M1 module is a
stepping stone to VET with short characteristics (tensor refreshed every hydro step); do not invent
closures. M1 work is now PAUSED: see [[m1-paused-back-to-dhj]].
