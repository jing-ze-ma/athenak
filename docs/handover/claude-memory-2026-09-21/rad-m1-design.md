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
OPAL ships Rosseland only (unresolved; ask the user, who may know what they actually did).

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
