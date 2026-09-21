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
