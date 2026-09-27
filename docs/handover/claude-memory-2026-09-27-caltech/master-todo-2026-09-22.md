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

## 09-22 ~21:30 (Opus 5.5 session) STATE BEFORE THE PLANNED RESTART
- Committed 886bd677: per-rank restart read fix (build_tree.cpp GetPosition flag), gated bitwise
  (bench/rst_perrank_0922). Local only, not pushed.
- IN FLIGHT: (1) rank-imbalance agent -> tests_rank_imbalance/README.md, bench/rank_imb_0922; user's
  remedy to test first = LEADING vs TRAILING hemisphere per GPU (cut on the substellar meridian), not
  day/night. (2) IMP/LEG template-tag agent -> tests_ck_kernel_tags/README_IMP_LEG.md,
  bench/tags_imp_leg_0922/imp_leg.patch (not committed; verify then commit). (3) sweep-form same-binary
  A/B job 11942118, bench/swform_0922 (prod4 binary, tm vs ck_sweep_form=0, spherical+beam on; first try
  11942072 failed: key missing from input, fixed).
- IMPLICIT ck (user: pursue, incl. spherical + beam): design note tests_ck_implicit/DESIGN_tm.md (tm
  T0 brief saved: tests_ck_implicit/T0_BRIEF.md -- launch it with the worker type right after the restart.
  refuses ck_implicit; tm is affine in B at frozen opacity; plan T0-T5). NEXT = launch T0 (measurement
  only: CFL 0.3 vs 0.15 dT above 0.01 bar from tm_prof_growth dumps + one apudev gap run) with the
  `worker` agent type. GO only if gap > ~1 % or dT > ~5 K.
- `worker` agent type = .claude/agents/worker.md (Opus 5.5, effort medium); needs a session restart to
  load. USER: use it for ALL new agents. See [[model-setup-opus55-medium]].
- 09-22 ~22:15 IMP/LEG tags DONE (tests_ck_kernel_tags/README_IMP_LEG.md, bench/tags_imp_leg_0922/imp_leg.patch,
  UNCOMMITTED): production kernel -26 % instr but NO GPU gain (20.20 vs 20.16 s cpu, 3 runs each; verified from
  logs); plane-parallel -15 %. CPU bitwise in 7 configs; GPU dumps differ at float32 round-off (<= 4.8e-6 rel;
  FOP commit's GPU "bitwise" was likewise CPU-only). NOT a prod4 swap candidate (no speedup). Commit = user
  decision; it touches the Jacobian path implicit-ck T1 will build on.
- 09-22 sweep-form same-binary A/B (bench/swform_0922, job 11942118, prod4 binary, sph+beam, rot-283, 300 cyc): tm 22.31/22.41 s cpu vs four-pass 28.48/28.02 -> tm is 21 % FASTER; the 'tm penalty of a few %' estimate in tests_tm_prof_growth is wrong (closed).
- 09-22 REAL-RAY check of the beam (tests_ck_sph/REAL_RAYS.md, prod3 dump rot 283): new ck_beam_sph rule closer to 3-D rays than the old one everywhere except a tie at mu0 0.40-0.45; terminator band (25 % of beam power) new +1.7 % vs old -16 %; planet total new +0.4 %, old -4.8 %; TWILIGHT (mu0 < 0, 2 % of power) new rule gives only 0.5-0.7 of the real rays and misses a high-altitude heating bump near 1.95 a_p -> horizontal inhomogeneity dominates there; fix = rays through neighbouring columns (short characteristics, list A e).
- 09-23 ~00:10 T1 DONE, NOT committed (bench/impl_t1_0922/t1.patch, tests_ck_implicit/README_T1.md): JAC tag, production kernel unchanged (bitwise CPU, identical counts GPU, 0 spills); tm+implicit allowed. CONVERGENCE FAILS: 0/600 calls converge at maxit 8, residual stalls ~4e-6, dstep ~0.11, gap median 5e-6 (semi -9.5 %); four-pass and plain Newton stall too -> not the tm Jacobian; GPU cost implicit tm 5.47x semi. Stall diagnosis agent launched (bench/impl_stall_0922 -> README_T1_stall.md): cells, dropped negative off-diagonals / farther entries (mode-3 exjac lesson), norm floor, non-differentiable switches, sph/beam isolation.
- 09-23 ~01:30 USER: hydro-only prod4 twin (inputs/production/deep_hot_jupiter_cs_hyd4.athinput, smoke test agent -> tests_hyd4/README.md, bench/cs_hyd4_smoke). APPROVED: if the smoke test works, LAUNCH a hydro production (prod4 recipe: clean binary = prod4 binary, 2 GPUs/1 node, apu, chained links like prod4; new run dir bench/cs_hyd4_prod). Also asked: resolution study (radial stretched + horizontal), optimal nodes / MeshBlocks, walltime -> planning agent.
- 09-23 ~01:45 IMPLICIT ck STALL CAUSE (tests_ck_implicit/README_T1_stall.md, verified): ck_impl_arat thin/thick split flips day-side TOP cells every pass (period-2 limit cycle): d ln E/d ln e ~ 53 there (H2 dissociation in e, Wien-side bands E~T^9.8) vs 4 assumed by the thin solve / ~10 by the Newton row (c_v = e/T wrong by ~5.4). Fix = ck_impl_arat 1e30 (input only): non-converged 400/400 -> 122/400, gap max 1.7e-4 -> 5.4e-6 (median 6.6e-7), cost 5.41x -> 4.66x semi (GPU 11943285). Remaining: ONE night-side corner column (gid 13, k=ks j=je, 2-17e-6 bar) with an odd-even T pattern ALREADY in the prod3 restart and cold cells with A = S + E < 0. Also: c_v = e/T is wrong in the H2-dissociation regime -- may matter for the semi-implicit linearisation too (unchecked).
- 09-23 ~02:15 HYDRO SMOKE PASSED (tests_hyd4/README.md, job 11943303: 7000 cycles to rot 0.57, dt identical to MHD to 4 digits, eos_fail 0, hydro 15.4-16.0 vs MHD 14.2 cycles/s) BUT cs_wellbalanced_src is read only from <mhd> (coordinates.cpp:42) -> hydro twin would differ from prod4. Fix+gate agent running (bench/cs_hyd4_prod: patch, clean binary, MHD-path bitwise gate vs prod4 binary, hydro smoke with the switch on, submit/chain prepared NOT submitted). LAUNCH the hydro production after that agent passes (user approval stands). Restarts 78 MB / 0.5 rot = 44 GB per run.
- 09-23 USER: make cs_wellbalanced_src and rt_newton DEFAULT ON (hydro + MHD) IF measurements show they help; plus any other default that should flip. In flight: hydro-fix agent evaluates cs_wellbalanced_src for hydro (tests_hyd4/README_wbsrc.md); rt_newton test vs the converged implicit ck reference (tests_rt_newton/README.md, bench/rtnewton_0923); switch audit with an explicit DEFAULT-FLIP list (docs/dev/cs_dhj_switch_audit.md). Flip defaults in one commit with regenerated references after the verdicts. Committed tonight: 055ca153 fused BiCGStab (implicit_bcg_sync), d837c355 closure = tau (multi-block).
- 09-23 SWITCH AUDIT (docs/dev/cs_dhj_switch_audit.md, verified item 1): dhj pgen NEVER reads rt_semi_lin / rt_newton / rt_rescue_eq -> prod4 runs the LEGACY linearised semi-implicit step (header default rt_semi_lin = true); box_convection defaults to the new step (false/true/true), red_giant per mode. rt_newton agent redirected: expose the keys in dhj (bitwise default) and test legacy / semi_lin false / + newton / + rescue_eq vs the implicit reference. Other audit items: cs_wellbalanced_src neutral or 7-14 % worse with plm in tests_cs_regions (audit says keep default false; hydro agent measuring); hyd4 hydro floors untested (tfloor 20-60x); silent traps: <mhd> efloor_from_ekin / efloor_as_tfloor ignored, eos_floor_consistent unused, prod4 rt_split=false overridden, rad_cap_ang inert. Flip candidates: rt_semi_lin->false (+ newton, rescue_eq), ck_impl_arat 2->1e30, grav_point_mass/rot_potential on for curved grids (user).
- 09-23 rt_newton VERDICT (tests_rt_newton/README.md, verified table): DO NOT flip. vs the converged implicit reference after 400 cycles (rot-283 restart): legacy (prod4) day rms 9.9/3.1/7.0 K at 1e-5/1e-4/1e-3 bar; rt_semi_lin=false alone within ~2 K of legacy; + rt_newton 26.4/21.7/36.4 K (3-7x WORSE day side) but better night side (1e-3 bar 44 -> 22 K). rescue_eq never fires. Cause: Newton uses EOS c_v but emission still ~T^4 (slope 22 in e vs 53 true). Fix idea (untested): band-summed emission slope from the sweep kernel. prod4/hyd4 legacy step is the best available semi-implicit choice on the day side. cs_wellbalanced_src default stays off (hurts pure hydro). pgen patch exposing the 3 keys (bitwise default): bench/rtnewton_0923/pgen.patch, not committed.
- 09-23 ~01:00 USER: continue implicit ck T3 (with T2); do NOT test on prod3 rot-283 restart (odd-even corner column) -> scan earlier prod3 states for the odd-even pattern and pick a developed clean restart. Worker in bench/impl_t3_0923 -> tests_ck_implicit/README_T3.md, t3.patch (T1+T2+T3). Also: prod4 input smoke with cs_wellbalanced_src off = job 11943851 (bench/prod4_wboff_smoke), check before 02:36.
- 09-23 01:05 prod4 wboff smoke PASS (job 11943851: cycle 5500 t 5.622e4 dt 6.91 vs switch-on 5.649e4 / 7.37; no FATAL/NaN). Multi-block VET sweep committed dccdf506 (1 GPU vet_sc = Eddington cost; 2 GPUs +5-6 %).
- 09-23 ~02:00 SMR EQUATOR SCOPE (docs/dev/cs_dhj_smr_equator_scope.md, uncommitted): NOT worth it -- octree refines radially too (3 levels = 1024 radial cells, global dt /8), band holds the narrowest cells, 16-cell blocks make bands illegal (polar seam); best band 137x vs uniform 256^2 164x per rotation (17 % saving for 4-7 weeks). Latent: two-stream column + rad_implicit_x1 silently wrong on radially split blocks under SMR (uniform guard passes on refined blocks) -> add a FATAL guard. First step if 1024 cells/equator is wanted: time uniform hyd4 at 128^2 and 256^2 per panel on 1-16 nodes (off-node penalty unmeasured). Productions started 01:31 (prod4 dt 4.9 s at rot 0.26; hyd4 dt 26 s at rot 0.78).
- 09-23 ~02:45 IMPLICIT ck T1+T2+T3 COMMITTED c02ddd34 (pushed), default off: ck_impl_lin slim passes; needs ck_impl_arat = 1e30; from clean prod3 rst 00135 (rot 67.5): 0/260 non-converged after +20 cycles, 6.2 passes, gap max 7.9e-7; GPU cost vs semi: T1 4.58x, T2 3.20x, T3 2.73x (design hoped 1.75x). Remaining: the two full passes per call (store + Jacobian, 46.8 ms) and ck_impl_tri 1.95 ms/pass. No prod3 state after rot 0 is fully free of odd-even columns (180-340 mildly flagged always); rst 567 has 9 strong ones.
- 09-23 ~03:00 USER: continue accelerating implicit ck -> T4 worker (bench/impl_t4_0923 -> tests_ck_implicit/README_T4.md): merge store+Jacobian passes, pass count (quadratic? dropped off-diagonals, warm start, colskip), PCR for ck_impl_tri, fewer syncs, >32 kB functors; test state prod3 rst 00135.
- 09-23 ~04:40 IMPLICIT ck T4 COMMITTED 95595afc (default off): ck_impl_fuse, ck_impl_jac_lin, ck_impl_cvsec on top of ck_impl_lin; GPU 2 GPUs: semi 10.76 s, T3 26.42 (2.46x), T4 18.33 (1.70x = design goal). OPEN: T4 vs T3 differ up to 56 % in T in 5 cells of MeshBlock 13 after 150 cycles (the odd-even corner column's block) -- check whether it is that column. Inode quota was down to ~4000 files: deleted 49 build dirs of finished tonight tasks -> >= 40000 headroom.
- 09-23 ~05:00 T4 vs T3 REDONE on hyd4 rst 23 (rot 11.5, primary) + prod4 rst 3 (rot 1.5) (tests_ck_implicit/README_T4_prod4.md): cost T4 1.77x / T3 2.48x semi (hyd4), 1.56x / 2.05x (prod4); 0 non-converged; passes 4.5 vs 6.2; gap max T4 1.1e-6 (hyd4) / 7.2e-6 (prod4) vs semi-implicit -10.4 % (hyd4) / -1.4 % (prod4, early). T4-T3 after 150 cycles: rms 3e-4 / 7e-4, max in kinked fronts; hyd4 worst cell: T4 = tight-tolerance reference, T3 7.4 % off. SCAN: new restarts are NOT cleaner: >0.1-flagged columns hyd4 1337, prod4 1359 vs prod3 rst 135 220 (ragged fronts, <=6-cell runs); prod4 MeshBlock 2 top (1e-9..1e-8 bar) T swings 200-1000 K.
