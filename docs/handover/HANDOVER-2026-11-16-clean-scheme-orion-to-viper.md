# HANDOVER 2026-11-16: clean radiation scheme, Orion -> viper (phased)

USER (2026-11-16): "slowly move the clean scheme writing and testing to viper, since the basic scheme is there and
gpu testing is important now". The clean cell-centred implicit grey M1 scheme (`scheme = clean`, `clean_layout =
cell`) moves from the Orion session back to viper in three phases. Times CET. Orion's full log:
/orion/ptmp/jinma/handover/HANDOVER_2026-11-15_clean.md (later lines supersede earlier ones).

## 1. Phases and protocol

| phase | when | viper | Orion |
|---|---|---|---|
| 1 | now | GPU tests (stage-1 G1-G5 on A100/MI300A, decisions arms later on Raven) and the sp/cs GPU divergence debugging (TASK cc_gpu_debug) | owns clean-cc code; audits the sp/cs GPU paths; finishes in-flight items (section 6) |
| 2 | as items land | reads the code, reruns the regression on Raven CPU once to own the harness | finishes items 6.1-6.6, merges them with a green regression, pushes; hands over harnesses (section 8) |
| 3 | after Orion's last push + a NOTE "phase 3" | owns clean-cc development and testing | supports on CPU on request (Orion farm, Python prototypes, references) |

Protocol (unchanged):
- Messages: NOTE/TASK/QUESTION-<date>-<from>-<to>-<topic>.md in docs/handover on `bsg-files-1009` (fork).
- ONE published branch: `origin/clean-cc` (fork), advanced only after a green regression (section 7.1), never main.
  Until phase 3 only Orion pushes it; in phase 3 only viper. `origin/clean-cc-int` is frozen at 4680ce500 (same
  code as push 9), do not use it.
- Restyle: viper merges restyle batches into `origin/clean-1114`; the clean-cc owner merges clean-1114 into
  clean-cc after each batch, gates, pushes with the next batch (agreed fd5106e4b). Restyle never touches the clean
  files (rad_m1_clean*, clean paths of rad_m1_vet*/vetgd*/time2, pgens rad_m1_rhse/clump/tests2). The production
  rad_m1 implicit/IDORT code is removed by the clean owner in one module when clean replaces production.
- Build gotchas: fresh `git worktree add` needs `git submodule update --init kokkos`; CPU builds `-ffp-contract=off`
  (A100 device `--fmad=false`); `M1_VGD_NMAX` 432 default (4096 = 32 KB/thread on GPU, avoid).

## 2. Code

| ref | hash | content |
|---|---|---|
| origin/clean-cc (push 12) | 0eeac5fdb | current published head (push 11 094488c1f + WB limiter fix c16a69021) = Orion's clean-cc-int; contains origin/clean-1114 (restyle) |
| origin/clean-cc push 10 | 32af6349e | GPU stage-1 binaries + ref/ (4e4345d67) were built from this (= 4680ce500 bitwise outputs) |
| origin/clean-cc-int | 4680ce500 | frozen (push 9 code) |
| clean-cc-clump (local, NOT merged) | 2db207987 | pgen rad_m1_clump (porous photosphere), see 5.8 |
| clean-cc-wbfix (merged, push 12) | c16a69021 | FV bracket limiter only on cells that went E<0 (5.2) |
| clean-cc-gpufix (local, in flight) | wt_gpufix on 15499232d | sp/cs GPU divergence audit (5.1) |
| clean-cc-torder (local, in flight) | wt_torder on 15499232d + dirty | time-order diagnosis (5.3) |
| clean-cc-milne 948db6fa6, -rdev-v1/v4/v5, -stage1 6fee08f5f, -causal fc6dfc47c, -lopos-c 136f09f19, -rwthin 933f79685, -stage f67aebb1f (debug) | local | superseded or rejected; do NOT merge (kept for reference) |

Integration worktree: /orion/ptmp/jinma/clean_1115/wt_int (= 0eeac5fdb). Files (7020 lines clean-only):
- `src/rad_m1/rad_m1_clean.{hpp,cpp}` driver of the clean step, preflight, steady start, restart state.
- `rad_m1_clean_cell.{hpp,cpp}` cell-layout rows: face fluxes, bracket, kinetic share, limiter, gas coupling.
- `rad_m1_clean_closure.cpp` closure from rays, extrap1 history, face sums, FV/LLD sweeps (CleanFvSweep).
- `rad_m1_clean_solver.cpp` BiCGStab, multigrid V-cycle, line solve (hybrid PCR-Thomas), fused dots.
- Clean paths in `rad_m1.{hpp,cpp}` (key parsing ~l.77-200), `rad_m1_time2.cpp` (hesdirk2 stages),
  `rad_m1_vet.cpp` (SC sweep, PLM foot, face-crossing F_ray), `rad_m1_vetgd*.cpp/.inc` (sp/cs rays, radial split,
  rotation union sweep).
- Test pgens: `src/pgen/tests/rad_m1_rhse.cpp` (radiative-hydrostatic envelopes, wb), `src/pgen/rad_m1_tests2.cpp`
  (rw, beams, thin/dif/lat/pulse, cs tests), rad_m1_clump (branch only).

## 3. The scheme as it is (0eeac5fdb)

### 3.1 Equations and unknowns
- Unknowns per cell: E and y = F/c (cell layout only; `clean_layout = cell` default; staggered frozen).
- Radiation sources: Ma 2025 mixed frame Eqs. 56-57 exactly: comoving absorption E', momentum source with the
  beta term, energy to beta^2 incl. (1 + beta^2/2) (vc fd147237c + a6e46bed4). v = 0 bitwise to the old form; a
  moving uniform medium is in equilibrium to ~beta^4.
- Gas: momentum from the same face fluxes/sources; gas TOTAL energy changes by -S_E; e_int = that minus the exact
  kinetic-energy change of the momentum coupling, the work taken at the time-centred velocity v_old + dm/(2 rho)
  (exact KE gain); gas velocity implicit in the flux rows via Newton (CleanCellDrag). No e_int deficit transfer.
- Opacities: kappa(T) Newton in the cell solve (prod2). OPEN: kappa(T) derivatives of the new v/c terms (5.9).

### 3.2 Closure (only from the rays)
- Rays give J, K, H+-, K+, C1 (and face sums); D = K/J, h+- = H+-/J. Cells the rays do not reach: isotropic fill
  (`clean_fill_passes`, 8; 0 with consistency path F0).
- `clean_closure = extrap1` (cell default): ray moments extrapolated in time from the stage-1 moments; one formal
  solution per step (~1.05/step), time order 2. Others: `iterate` (reference), `extrap`, `eddington` (fixed D=1/3,
  order gate), `sphere` (needs `clean_sphere_radius`).

### 3.3 Face flux (cell layout)
- Moment flux (central) + dissipation alpha (a/2)[U] with alpha = A(tau)^2, A = (1-e^-tau)/tau of the face optical
  depth, + Rhie-Chow with weight rc = (1-alpha)/tau (CLF_TH). `clean_alpha = tau` (default). Jumps from van Leer
  PLM, lagged slopes (`clean_recon = plm`). Temporary key `clean_dissip_power` = 2 (default) | 3 (A^3; decision 5.11).
- Consistency bracket (`clean_consistency`; default on for Cartesian and for sp with vet_gd, opt-in on cs where
  it is implemented with seam face sums (cs3); `rad_m1.cpp:114`): face flux driven to the rays' face flux,
  C = F_ray,face - F_LO; streaming weight w(C1) (`clean_consistency_weight`, default on, constants `clean_w_c0` 0.1,
  `clean_w_c1` 0.5; weight MUST accompany consistency) times a causal share min(1, E/J) (ramped; fades the CLN_G
  curvature share too).
- Kinetic deviation transport: in the streaming share F = (1-w) f(U) + w [F_ray + F_kin(U - Ur)], F_kin = rays'
  upwind half-range flux; y and E move with the same upwind operator (R form v2), momentum round-off.
- Low-order faces (positivity marking): E flux = rays' kinetic flux blended with compact diffusion via face alpha
  (M-matrix E row); stress = the cell's own ray stress (no D = I). Positivity limiter marks cells with E < 0 only
  (`clean_pos_limit`, true; up to `clean_maxpos` 3 re-solves, then BE with marks, then all-low BE).
- FV bracket limiter (FCT on the bracket E outflow): active with lld and fv1; SC has none (hurts beams 1.5-3.5x).
  Being restricted to cells that went E < 0 (5.2).
- Milne sub-cell face model: `clean_face_model = milne` (opt-in, needs consistency; user-approved exception to the
  chi*dx rule).

### 3.4 Rays
- Cartesian: SC (`vet_sc`), PLM foot interpolation (clean only), face-crossing x1/x2/x3 F_ray sums; opt-in
  `vet_seg_cross` (false), `vet_path = cells`.
- sp/cs (vet_gd): `clean_rays = fv1` (default with vet_gd): conservative FV transport, wavefront-level sweep
  (42129b32c; width x28 sp, x350 Cartesian), bitwise vs sequential. Opt-in `fv` (sweep to tolerance), `lld`
  (lumped linear-discontinuous, asymptotic-preserving within 0.8 % at tau 10-1000; 3x worse than fv1 on thin cs),
  `sc`. FV/LLD also on Cartesian (opt-in).
- Rotation: union sweep of the K (+xfade) sets, device per-set tables, restart bitwise (`vet_gd_rotate_dt` must be
  named in production decks: the auto clock breaks time order; `vet_gd_rotate_xfade`, `vet_gd_rotate_every`).
  Geodesic sets `vet_gd_geo = n` (geodesic 4 best on cs).

### 3.5 Geometry and parallelism
- Cartesian 1/2/3-D; spherical polar wedge, full sphere (polar face zero-area reflect), 2-D axisymmetric vet_gd;
  cubed sphere (y in hydro local frame, y_ph angular-momentum form). One code path with geometry hooks
  (CleanCellRows/CleanCellBlock unified).
- Multi-MeshBlock/MPI: AthenaK CC exchange (chunks <= 32 vars), CleanReduce for all decisions, block-Jacobi + MG
  (coarse levels volume-integrated on sp/cs), exact pipelined banded SC sweep, vet_gd radial split (sp and cs).
  Decompositions agree to solver tolerance (USER accepted; bitwise across decompositions not wanted).
- Conservation to round-off (gas + radiation): total E, momentum (Cartesian), radial momentum (sp, cs),
  L_z (sp, cs; L_z only on cs, full vector L not wanted).

### 3.6 Time integration, starts, restarts
- hesdirk2 stages, BE fallback; dt_retry. No face FCT / k1_redist / ApplyClosureLimits under clean.
- `clean_init = given | steady` (read only when named): steady pre-solve (gas frozen, `clean_init_n` 50,
  `clean_init_dtfac` 1e6, `clean_init_tol` 1e-8); out-of-equilibrium starts must use steady.
- Restart BITWISE (all clean state incl. extrap1 history and wb reference force in the restart file).
- Solver: `clean_tol` 1e-10, `clean_tol_nl` 1e-11, `clean_maxit` 4000, `clean_maxnl` 50, `clean_maxpic` 12,
  `clean_pc_sweeps` 1; round-off-aware Newton exit with hard cap; BiCGStab true-rho breakdown restart. Clean
  ignores implicit_tol / implicit_lin_tol. Temporary key `clean_line_solver = hybrid (default) | pcr | thomas`.

### 3.7 Well-balanced and production features
- `force_reference = clean` + `<hydro> wellbalance_static = true`: Cartesian/sp/cs envelopes v/cs ~1e-8..2e-7.
- Production BCs: inner flux BC (implicit_flux_x1min), outflow, esrc (mlt_flux_frozen) with energy ledger;
  efix still FATAL. `CleanPreflight` (rad_m1_clean.cpp:40-80): production keys mapped to clean keys, ignored ones
  listed, physics-changing ones FATAL. ic_rad off under clean; wb_arad reference work paid.

## 4. Binding rules (USER)
- Code style: original AthenaK upstream style (src/hydro on main): 90 columns, concise comments on each function and
  main step, NO history/patch/tag comments, no dead code, cpplint clean.
- Minimal lines: add as few lines as possible, reuse existing pieces/hooks, scripts outside src; report line counts
  per file.
- Obstacles: state the physical cause and derive OUR OWN physical fix first; only then cross-check how the old
  production scheme did it (engineering reference, not optimal, not an accuracy reference); no numerical patches.
- Closure and regime ONLY from the rays; the optical depth chi*dx ONLY in the numerical dissipation. Exception:
  opt-in `clean_face_model = milne` (default only after Cartesian data).
- Conservation to round-off on every test: total energy and momentum (gas + radiation); on sp and cs also radial
  momentum and z angular momentum. cs: L_z only (no 4th spin unknown).
- Restart bitwise: hard gate, all clean state in the restart file.
- Every run: zero "CLEAN WARNING" lines (rc = 0 is not enough); also grep DT_RETRY.
- Cell layout only (staggered frozen: no runs, no fixes).
- Tests at PRODUCTION tolerance (code defaults; no per-test tight tolerances). rw design: +-amp pair, score the
  odd part (u(+a) - u(-a))/2; amp = min(1e-3, 0.03 E1/E0); default tol except tau 0.01 (per-case clean_tol in the
  driver so the solver resolves the signal); `vet_gd_rotate_dt = 0` for multi-D rw vs vetqs_d; margin: measured
  difference >= ~100x effective tolerance; clean_tol floor ~1e-14. Thin-rw failures are test design: the solver
  stopping rule stays (933f79685 rejected).
- Fast tests: all independent arms concurrently in one allocation, max ranks per run, smallest decks (x0.5, x0.25,
  32^3), skip x2/x4 unless targeted; pre-push regression <= 30 min wall; state the critical-path ETA. Pass keys via
  env/file, never "a;b" through `bash -c`. Farm jobs: `export OMP_PROC_BIND=false
  OMPI_MCA_hwloc_base_binding_policy=none` (+ OMP_NUM_THREADS) or singletons bind to core 0.
- One code for all backends (CPU/CUDA/HIP; no backend branches); reuse the same code between Cartesian, sp and cs
  (geometry hooks, no per-geometry copies of rows/closure/sweeps).
- Push to the fork only; never main.

## 5. Open problems, ranked
1. **GPU sp/cs divergence (top).** NOTE 050a79319: Cartesian G1 exact on A100 and MI300A; sp/cs fail with
   different symptoms (MI300A: spthin m1_f1 maxrel 0.27, cs_thin1 0.99, rhse_sp2 FATAL cycle 220; A100: spthin =
   CPU, rhse_sp2 FATAL cycle 4 via `cln_neint > 0` dt_retry at rad_m1_clean.cpp:779 with no dtr_why). Not the
   wavefront sweep (4680ce500 identical). Suspects: uninitialised View read, race, missing deep_copy in sp/cs-only
   paths (fv1 rays, metric/area factors, panel/pole ghosts). TASK cc_gpu_debug (8e010e76d): D1 compute-sanitizer
   initcheck/racecheck/memcheck, D2 bounds check, D3 cycle-0/1 dumps, D4 HSA_XNACK=0. Orion audit worker:
   clean-cc-gpufix (wt_gpufix, gpufix_tests/). Also G1 flag: cb score 0.4707 (GPU) vs 0.4600 (CPU) with identical
   dumps -> scorer artifact: ana/score.py reads the beam angle from run.log ('implicit_beam ... in-plane ang=') into the exact solution; send the GPU run.log implicit_beam lines.
2. **WB level regression** (sp WB 7e-8 -> 7e-6, fv 4e-3): FV bracket limiter firing in steady-start passes
   (h = 1.5e4). Fix: limit only cells that went E < 0 (c16a69021, merged in push 12: WB back to 6e-8..8e-8, lld beams stay positive).
3. **Time-order reduction in stiff radiation-heavy gas** (tau >= 100, P >= 0.03): order 0.1-1 below 1e-4 relative,
   not tolerance. Worker clean-cc-torder (wt_torder, torder_tests/) diagnosing SDIRK stage order / hydro-rad
   coupling / closure lag.
4. **Vacuum pulse / steady rays** (pulse0 L1 2.28 vs 1.41 earlier, F = 0; |F| > cE ahead of light fronts; weight 1
   unstable on dif/lat): steady rays are acausal. Proposal: time-dependent rays (chi + 1/(c dt)). USER DECISION
   PENDING (Orion first recommended defer; after pulse0 and the weight-1 failure it recommends implementing it for FV/lld rays, without rotation at first).
5. Periodic bubbles: stripes on periodic decks (IDORT page bub/bubu sections).
6. hcyl 1.6 % E below J_sc: SC face sums inconsistent with cell-centred opacity (corner segments, w-ramp leak);
   fixes break nltek4/sp dif/rw; candidate = FV rays on Cartesian (opt-in now).
7. Production gaps (/orion/ptmp/jinma/clean_1115/prodready/GAPS_2026-11-16.md; closed since: inner flux BC,
   esrc, ray rotation on sp/cs, preflight, cell default, kappa(T) Newton, bitwise restart): memory ~2.5x prod per
   block (A100: 1 block per 40 GB); he_star_m1 given-IC ledger miss outside radiation (prod misses similarly:
   wb pressure-correction work, gravity, floorcap?); cs vet_gd SC radial blocks; lbath prod-only; merge_lat_rad,
   sp_x2_periodic_image; cgs absolute thresholds; diagnostics (m1_vet on sp/cs, m1_face x2/x3, CLEAN printf per
   solve -> verbosity knob); efix FATAL; no path to continue a production rst under clean.
8. Clumpy photosphere (clean-cc-clump 2db207987, seed 20261116): SC best (Ftop +2.5..+8.8 %), Milne worse, FV
   needs rerun on head (results were pre-posE/lophys); FV limited by 108 directions.
9. kappa(T) derivatives of the new v/c terms (Eq. 56-57) in the Newton Jacobian.
10. Cartesian reflecting walls not mirrored in FV/LLD rays; nltek-type 1-D FV cost; sp phi-periodic 3-D wedge +
    consistency v/cs ~3e-3; theta-periodic vet_gd + consistency unstable (use physical theta walls).
11. Decisions for Raven GPU tests (TASK cc_gpu_decisions, later): A^2 vs A^3 (`clean_dissip_power`), direction
    set / rotation (prod defaults proposed: geodesic 4 + K4 + xfade on cs, explicit `vet_gd_rotate_dt`), line
    solver (`clean_line_solver`); remove the temporary keys once decided.

## 6. Phase 2: Orion in-flight items

| item | owner (Orion worker) | branch / worktree | status 2026-11-16 |
|---|---|---|---|
| 6.1 GPU sp/cs divergence audit (5.1) | gpufix | clean-cc-gpufix / wt_gpufix | code audit + CPU-side bounds/uninit checks; waits viper D1-D4 |
| 6.2 WB limiter fix (5.2) | wb | c16a69021 | DONE, push 12 0eeac5fdb |
| 6.3 time-order diagnosis (5.3) | torder | clean-cc-torder / wt_torder (dirty) | diagnosis runs (jobs 230434-6) |
| 6.4 pulse / time-dependent rays (5.4) | main | - | awaiting USER decision |
| 6.5 clump test merge (5.8) | clump | clean-cc-clump 2db207987 | rerun FV arms on head, then merge |
| 6.6 push 12 | wb (regression) | origin/clean-cc 0eeac5fdb | DONE (reg13 green) |

Orion ends phase 2 with a NOTE naming the last origin/clean-cc hash, the unmerged branches (pushed to the fork as
`clean-cc-<topic>` if worth keeping) and the state of each open problem.

## 7. Testing

### 7.1 Pre-push regression (Orion CPU, 1 node, <= 30 min)
/orion/ptmp/jinma/clean_1115/wb_tests/reg8: `job9.sbatch` (xargs farm, NPAR 30, CMDS=cmdsNN.txt), latest list
`cmds13.txt` (29 runs), analysis `ana12.sh` (its, CLEAN WARNING count, cpu, scores via gtg.py/ct.py/fb.py/
energy.py). Sections: cc/ (IDORT ladder via run_ladder.sh: hcyl x0.5, cbg/cb x0.5, shd x0.25, nltek1), rw/
(+-amp), sp thin_p/_r/_f, dif, lat, pulse0, s88/s81 steady-start, rhse/ (sp WB), cs/ (uni, closed, thin1, feed1),
fb/ flz q/Q (feedback). Arms X (new) vs W (previous push binary in bin/).
Must be: bitwise for every change claimed bitwise (v = 0, keys off, Cartesian when sp/cs change); 0 CLEAN WARNING,
0 DT_RETRY except listed pre-existing; conservation round-off (closed cs E, L_z; rw dE/dP ~1e-15); scores
unchanged or explained.

### 7.2 Restart harness
/orion/ptmp/jinma/clean_1115/wb_tests/rst: rt.sh / rtmpi.sh (run N, restart at N/2, cmpbin.py bitwise),
decks cart, cb, hcyl, feed, sp3, thin1 (cs); cmds_rst*.txt.

### 7.3 Linear waves (order)
rw driver run_radwave.py + rwfarm.py (handover_1115_clean/rw; Orion copies in rwfix/, rw_prod/); +-amp design
(section 4); gates: space self-convergence (N vs 2N, amp 1e-3) 2nd order (1.83-2.55), time finest pair 1.71-2.29.
Page grid: rw/runs2 (288 runs). Known: 5.3.

### 7.4 IDORT accuracy suite and page
- Decks + scorer: task_1115_cc_idort_head + handover_1115_clean/ana (score.py, nlte_exact.py; new NLTE reference
  nlte_ref.py (Feautrier) in the page tree; hcyl/shd exact npz).
- Page pipeline: /orion/ptmp/jinma/clean_1115/idort_page (suite/run_all3.sh, mklinks_clean.py, extract.py,
  plot.py, mkcontent_clean.py, build_clean.py; inputs from idort_page_inputs/). Published artifact (all 10
  sections, push 11): https://claude.ai/artifact/WTbarCvoRbP4oNhmTGCA4f (the older 5A7UvXPiAAdcnEoxBVzgjW belongs
  to the previous account). Republishing needs the account that owns it; otherwise publish a new one.
- Last logged scores (x0.5 unless noted; rerun on the head before quoting): cb 0.460 (GPU G1 flag: 0.4707),
  cbg 1.308, ws66 1.167, hcyl ~0.012-0.013 (x1 outside 0.0156 after the PLM foot), shd x0.25 0.021-0.023,
  shdg 0.032, nltek1 x1 9.1e-5 (nmu 8 quadrature limit).

### 7.5 GPU packages (bsg-files-1009 docs/handover)
- task_1116_cc_gpu1 (TASK cc_gpu1 + -a100 amendment): decks G1 Cartesian, G2 multi-block, G3 sp, G4 cs, G5 timing;
  ref/ on 4680ce500 (4e4345d67), compare.py, build_cuda_raven.sh. Status: G1 pass both GPUs, G3/G4 fail (5.1).
- task_1116_cc_gpu_decisions (TASK cc_gpu_decisions): arms A^2/A^3, geodesic/rotation, lld, line solver; Raven
  LATER (USER), not part of stage 1. CPU thin-sphere L1(E) (cpu_values.txt, 4680ce500):

| arm | thin1 | thin1l |
|---|---|---|
| base (fv1, A^2) | 4.46e-2 | 1.16e-1 |
| a3 | 4.46e-2 | 1.16e-1 |
| geo4 | 2.50e-2 | 6.31e-2 |
| geo4rot | 1.84e-2 | (not finished) |
| rot | 2.96e-2 | (not finished) |
| lld | 1.50e-1 | 3.09e-1 |

- task_1116 debug TASK cc_gpu_debug (8e010e76d) on 094488c1f: D1-D4 (5.1).
- viper binaries: Raven /raven/ptmp/jinma/ccgpu_1116/bin, viper /viper/ptmp2/jinma/ccgpu_1116/bin; runs
  runs_a100, runs_mi300a, runs_ab4680.

## 8. Where everything is (Orion)
- Work root /orion/ptmp/jinma/clean_1115: wt_int (integration), wt_<topic> worktrees, <topic>_tests dirs.
- Running log /orion/ptmp/jinma/handover/HANDOVER_2026-11-15_clean.md; queue HANDOVER_2026-11-16_queue.md; gaps
  clean_1115/prodready/GAPS_2026-11-16.md; cs design clean_1115/cs_design.md; non-physical-device audit in the log
  (A1-C2, WPs 1-8).
- Test harnesses to hand over in phase 2 (Orion will copy them into a package on bsg-files-1009 with Raven paths
  parametrised): wb_tests/reg8 (regression), wb_tests/rst (restart), cs_tests/ct.py, wb_tests/gtg.py, fb/fb.py,
  rw driver (rwfix/), idort_page tree, clump_tests, torder_tests.
- Prototypes (Python): beamreal/, sp_proto/ (overlap remap), corner/, handover_1115_clean/proto2.
