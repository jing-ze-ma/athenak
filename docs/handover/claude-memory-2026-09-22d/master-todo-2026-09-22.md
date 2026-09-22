---
name: master-todo-2026-09-22
description: MASTER PLAN agreed 2026-09-22 02:15 - order = wait for the running agents, CLEAN UP the code, then the three continuation lists (dhj RT, cs vertex/seam, cs dhj productions); M1 stays paused
metadata:
  type: project
---

User 09-22 ~02:15: "wait for all agents to finish, then clean up first, then proceed" (several lists).

## 0. In flight at writing (do NOT start new work before these land)
- CPU agent: IMPLICIT ck column solve (tests_ck_implicit/) + the repair of the ck_spherical=false path to
  the production plane-parallel form (regressed hybrid at HEAD, see [[dhj-improvement-list-0921]]).
- GPU jobs (all with ck_spherical = true unless noted): ck_sph_ab sphbeam 11931204 / prodbin 11931220
  (original binary, plane-parallel); recon plm_ctl/ppmx/wenoz 11931395-7; noang 11931398 (control =
  sph 11931203); wb-cache every1/every10 11931399/11931400 (apu, 10 h).
  READ: bench/*/README.md comparison plans. Results to record here when read.

## 1. CLEAN-UP (first, before any list) -> [[code-cleanup-pause-plan]]
Inventory every runtime switch (two-stream/ck, conduction, WB, rad_m1, dhj/box/RG pgens) with default,
production values, gate coverage; user decides per switch: ABANDON / DEFAULT ON / KEEP. Small commits,
regenerate bitwise references when a default changes. Also: repair or remove the HEAD hybrid; stale
bench NOTES.md; restart-format versions; docs/handover update; push to the fork after user OK.

## 2. LIST A - dhj radiative transfer (after clean-up)  -> [[dhj-improvement-list-0921]]
a. Read the GPU A/B: production (prodbin) vs sph vs sphbeam: dt (sph 12.5 vs prod 20.8 s: WHY), floors,
   T(p) day/night at 1e-6..1e-3 bar, cost. Then decide defaults (ck_spherical, ck_beam_sph on?).
b. Implicit ck column solve: gates, cost; fold conduction (removes the 10 bar split); then default?
c. rad_angular = false for dhj if the noang arm is inert (committed switch 213095e4).
d. Beam power reaching the ck cut is dropped (0.003 %): fine, documented.
e. Long-term: short characteristics per band/g-point (heating, beam, VET tensor). M1 stays PAUSED
   ([[m1-paused-back-to-dhj]]).

## 3. LIST B - cs vertex / seam accuracy  (tests_cs_regions/ is the gate)
a. Angular ppmx or wenoz instead of plm (vertex error x0.31 L1, x0.15 Linf, cost 1.28x); ppm4 is bad.
   GPU recon arms decide for the real atmosphere.
b. Make the RADIAL (x1) reconstruction selectable on cs (always GridPLM today).
c. REAL cube-vertex fill (sampling the third panel; prototype beat extrapolation 3.4-5.5x) - the one
   real build; plus the corner EMF at the 3-valent vertex (field-loop vertex L1 unaffected by any
   reconstruction -> CT/EMF suspect).
d. cs_wellbalanced_src: does NOT help the vertex (8 % worse), unreadable from <hydro>: fix or abandon.
e. Receiver-side seam resample (3rd-order block-end residual): low priority.
f. Closed: seams (2nd order), fc halo (no defect), cc corner slots (fixed ec0a9724), angular momentum
   (cs not worse than sp, tests_cs_angmom/).

## 4. LIST C - cs dhj production practice
a. wb_cache_every = 1 if WB is used (GPU arms test the rot-12.3 death); cs_mhd_prod3_wb ran WITHOUT
   cs_wellbalanced_src (its NOTES.md is wrong).
b. dt_min 1e-2 -> 1e-3.
c. <mhd>/fofc was never on; dhj.log fofc column is meaningless; use eos_dfloor/efloor/tfloor/eos_fail.
d. Restarts: MHD + cs caches now bitwise (0eb7e2c4, 544cb36b); pgen state (dump clocks, WB cache
   flags) not audited; cells flooring in the restart cycle.
e. Floors at the night-side vacuum top (dfloor vs WB background, keep-temperature).
f. Next production: new input with the chosen defaults; chain from scratch or from rot 283?
   (user decision). bench/cs_mhd_prod3/rst = 46 GB (keep a subset? user).
g. Physics questions parked: shallow jet / deep westward flow (sp agrees), day-night contrast at low p.

## 5. Also parked
- He4 global star (he4-session-state-2026-09-21), B-star/He boxes: see their notes.
- Upstream merge (72 commits behind origin/main; 6 overlapping files).
- origin remote URL embeds a token (user action).

## Update 09-22 02:40: implicit ck agent DONE (committed after 76359975). The HEAD hybrid is REPAIRED
(ck_spherical=false = production plane-parallel divergence; sweep-to-gas gap -8/-4 %, not -36 %).
ck_implicit: gap 1e-9..3e-8 in all combos, but 5.1-5.4x RT cost and DIVERGES at 10-100x dt (thin-cell
linearisation; needs the grey mode-3 two-level split rt_impl_tau_min) -> NOT production-ready; list A item b.
ALL CPU AGENTS DONE. Remaining in flight: the GPU jobs only. Next = CLEAN-UP (section 1).

## Update 09-22 ~05:00: clean-up committed (47eb060e G1/2/4/5, 912ef43c G3/6, ac42334d inputs; -4779 src lines);
default flips in progress (Opus). Isolation arm offfix: HEAD plane-parallel = production dt 19.9 s ->
the dt cost is ck_spherical itself. USER ORDER after the flips: start LIST B items c (real cube-vertex
fill + corner EMF at the 3-valent vertex) and b (selectable radial reconstruction on cs); gate =
tests_cs_regions/run_regions.sh. Then C-f (prod4 launch, after the cache arms reach rot 12).

## Update 09-22 06:45: ck_implicit phase 2 committed 4a8daf88 (2.46x cost only WITH ck_impl_once, which is the
box-rejected once-per-cycle split -> must NOT be used without a multi-crossing f-mode growth gate; keep OFF).
USER-APPROVED next for ck_implicit (Opus, after the two cs agents): (1) FREEZE THE EXCHANGE OPERATOR: store
the sweep's layer coefficients once per step and re-apply M_g to the new B(T) in each Newton pass (no
re-sweep, no table lookups; identical fixed point); (3) warm start from the previous step; (6) GPU: with
stored coefficients the re-application parallelises over g-points. Target <= 2x without ck_impl_once.
Also noted: (2) skip saturated g-points in deep cells, (4) reduced-g Newton + full-g correction.

## Update 09-22 09:30: LIST B closed out. B-b committed (reconstruct_x1 = Mignone-2014 curvilinear ppm4/wenoz;
the existing PLM already IS Mignone's; radial choice does nothing at the vertex; ppm4 10x cost). B-c
committed 212857c1 (cc vertex ghosts sampled: 11x better statically, nothing dynamically). VERDICT: the
cs vertex error is NOT ghosts, NOT the radial scheme, NOT WB; angular ppmx/wenoz help hydro tests only
and are null+1.44x on the real dhj atmosphere -> STOP pursuing the vertex for dhj; remaining leads
(fc x2f vertex ghost not converging; vertex EMF across SMR levels; div F vs source cancellation) parked.
HLLD tolerance fix a544a761. Live: top-slab overcooling fix (Opus), GPU profile jobs 11932897-9.

## 09-22 10:30 cs VERTEX CLOSED (tests_cs_vertex_cell): geometry exact closed-form; discrete Gauss identity
holds (vertex 30 % BETTER than interior); pure hydro balances to round-off at the vertex; 100 % of the
1:260 residual is the MAXWELL STRESS -(B.grad)B, converging at the same order as the interior but 3-5x
larger (|d(tangent basis)/d angle| maximal at the 120-deg junction). The field-loop "vertex order 1.57"
is RETRACTED: at fixed tlim it sampled a transient at different phases; at tlim 1.0 the vertex order is
1.83-1.96 with a saturating ratio ~4.9. => the vertex is a constant-factor hot spot at 2nd order, NOT a
lower-order region; no fix warranted. Open (stability, not accuracy): slow error growth in the loop
after t ~ 0.2 in all regions. Vertex ghosts cc/fc both sampled & 3rd order (cc default off).

## 09-22 ~12:30 USER WENT TO BED, left me autonomous. Standing constraints while unattended: nothing pushed
to origin (fork pushes OK), no production launch (prod4 stays DRAFT), no apu-partition jobs beyond what
is running (apudev/apu1 short jobs OK, the approved counts), no chaining of the wb-cache arms, no M1
work, commit only gated results, keep the handover current. IN FLIGHT: probe-free sweep (sd/tm forms),
field-loop growth, MPI gate, Ohmic-cap-in-sphbeam analysis; then: bisect the 1.30x regression
a544a761..HEAD+phase4 before committing phase 4; commit ck_sweep_cache + phase 4 once the header is free;
GPU confirmation arm for pfloor 1e-5 (rot-283 restart, ck defaults); update prod4 draft with pfloor;
handover 22d + memory copy at the end.
12:45 field-loop growth CLOSED (verdict i: O(h^2) truncation accumulating on an exact static equilibrium,
converges 1.9-2.1 everywhere, steady spurious v ~ h, vertex quieter; bs_emf = true is an instability,
never use). LIST B fully closed. MPI gate interim: CPU 1/2/4 ranks BITWISE on cs hydro+MHD blasts with
corner-slot and vertex_fill_cc on and off; GPU 4-rank/2-node arm queued on apu (mgate_gr4).
14:00 MPI GATE PASSED (tests_mpi_gate_0922): CPU 1/2/4 ranks bitwise incl. the new seam/vertex code; restarts
bitwise across rank counts; GPU 1 / 2 / 4-ranks-on-2-nodes write BYTE-IDENTICAL restarts and dumps (the
fixed-order global sums still hold on GPU). Scaling per GPU 1.00 / 0.80 / 0.51 (24/12/6 blocks per GPU;
off-node costs 1.56x) -> the production's 2-GPU/1-node layout is the efficient one; more GPUs need
more blocks. DEFECT: per-rank (rank_%08d) restart files abort on read, build_tree.cpp:429 missing
single_file_per_rank in GetPosition -- not used by the productions; fix when convenient.
18:50 USER: prod4 runs FROM SCRATCH (not from rot 283); approved (7) a rocprof arm on the tm binary and
(8) the semi-implicit vs implicit multi-rotation growth test; wait for the pfloor agent's final verdict
before anything else on the floors. tm = the sweep only; the gas coupling stays semi-implicit.

## 09-22 20:15 SESSION CLOSE-OUT. ALL prod4 PREREQUISITES DONE: MPI gate passed; pfloor 1e-5 confirmed;
tm default (0172f3bc); 1.30x regression fixed by the FOP tag (cfe92ac7; 15.26 cycles/s > old prod 12.84);
cold-start split-lag gate PASSED (growth ratio 1.06 between CFL 0.3/0.15 -> physical); wb-cache arms both
clean to rot 14.6-14.9 (stale cache did NOT kill the old WB arm; keep every 0). prod4 draft final:
inputs/production/deep_hot_jupiter_cs_prod4.athinput, FROM SCRATCH (user). NEEDS ONLY THE USER'S GO.
Launch recipe = bench/cs_mhd_prod3/submit.sh + chain.sh, HEAD binary from a clean snapshot (strings check
ck_sweep_form, hlld_bx_zero_tol), 2 GPUs / 1 node (the efficient layout), apu, 4 chained links.
NEXT COST LEVERS (after launch): (1) sweep rank imbalance 20 % (rank 1 vs 0, new with tm); (2) template
tags IMP (ckjacp_/ckskip_) and LEG (rt_layer_legacy): probe build halves the kernel; (3) occupancy.
OPEN SMALL: eos_efloor nonzero (8.8e6 events / 0.17 rot) from cold start with pfloor 1e-5, harmless,
look at the first production dumps; per-rank restart read defect; grey rt_split discrepancy.
