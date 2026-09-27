- **[CALTECH SESSION START HERE: docs/handover/SESSION-2026-09-27-caltech.md](caltech-cluster-facts.md)** — merges, decisions, next ck/M1 levers, how to test here (rt-integration faf14314)
- **[NEVER RUN ON hpc-sm-01-09 (bad H200 node, user 09-26)](caltech-cluster-facts.md)** — ~6.6x slow + OOM; every GPU sbatch gets `--exclude=hpc-sm-01-09`; pending jobs: `scontrol update JobId=<id> ExcNodeList=hpc-sm-01-09`; brief every agent
- [ck-next resume note (RESOLVED 09-27: merged 6a9ffb0d, pushed)](ck-next-resume-0927.md) — historical
- [Delegate to agents](delegate-to-agents.md) — Opus 5.5 medium agents; AVOID haiku/sonnet (user 09-25); never high effort
- [Caltech cluster facts](caltech-cluster-facts.md) — Slurm, gpu 4/node (18x H200); SCRATCH PURGED after 14 d; H200 ~$1.9/h billed; login 8 GB/proc; srun --mpi=pmix; account carnegie_poc
- [caltech-port MERGED + pushed 09-25](caltech-port-hold-merge.md) — rt-integration 83876228; viper: submodule update --init once + HIP bitwise gate
- [Carnegie HPC usage policy](carnegie-hpc-usage-policy.md) — per quarter: GPU 7,500 h (Tier 1) / 15,000 h (Tier 2); storage 50 / 100 TB; required acknowledgement text
- [Caltech = TESTING ONLY (user 09-25)](caltech-testing-only.md) — no production here unless really needed; short tests on scratch
- [H200 validation 09-25](caltech-h200-validation-0925.md) — port validated (restart, 2-GPU bitwise); H200 ~1.9x viper; 3 porting bug classes; compare GPU vs CPU on new setups
- [M1 real sp wedge -> HANDED TO VIPER 09-25](m1-real-wedge-0925.md) — task doc TASK-2026-09-26-m1-real-wedge.md pushed c1dc8226; no Caltech code
- [Default flips 09-25](default-flips-0925.md) — kkt_row + xstep 8 merged 504d8e30; mg_gf DROPPED 09-26 (+33 % on 2 H200)
- [Kokkos pin bumped to 4.6.02 (09-26)](kokkos-pin-4400-not-4602.md) — 2e028a38; all gates pass; hst 3-mom column is a cancelling sum, cmp of GPU hst shows DIFF even when state is bitwise
- [OVERNIGHT PLAN 09-26 Caltech (user asleep)](overnight-plan-0926-caltech.md) — merge m1-port if gated; ck report-only; Kokkos 4.6.02 bump + re-gate; mg_gf held; one morning summary
- [ck H200 experiments 09-26](ck-h200-experiments-0926.md) — no 25 % lever; beam_par -6.2 % (11 % of ck, within noise); host gaps 1.5 %; dif_dtau +106 %; jlin rewrite proposal
- [M1 port merged + H200 M1 timing 09-26](m1-port-h200-0926.md) — 48321155; mg_gf +33 % on 2 H200 (only wins 1 GPU); one_pass off -3.5 %; launch-latency bound
- [Drop unhelpful default-on switches (user 09-26)](drop-unhelpful-defaults-0926.md) — one_pass kept + auto-disable MERGED dff69bc7 (accuracy 5-22x round-off accepted by user)
- [ck-jlin MERGED 09-26](ck-jlin-merged-0926.md) — 587825ef; ~21-22 % of ck; next lever = store kernel registers; final timing CONFIRMED -13.6 % (2 H200) / -16.1 % (1 H200)
- [ck-next MERGED 6a9ffb0d 09-27 (-9.4 %/-11 %; cumulative ck -27 %/-32 %)](ck-store-0926.md); [ck-store MERGED 55762ada 09-26](ck-store-0926.md) — store pass without chain kernel (ck_store_split); total ~35-38 % of ck with jlin; noise 1.21x at edge
- [Agent collision lesson (09-26)](agent-collision-lesson.md) — agents with background work can resume; TaskStop the old one before a replacement on the same branch
- [Kokkos 5.2.2: NOT adopted (09-26)](kokkos5-decision-0926.md) — 1-3 % dhj, 0 % M1 on good node; big rename cost; branch kokkos5 kept
## RULE ZERO

- **[VIPER = 2 GPUs PER NODE, RECURRING MISTAKE (user 09-26)](viper-2-gpus-per-node.md)** — apu AND apudev nodes are gpu:2; apudev max 1 node; >2 GPUs = apu N x 2; NEVER gpu:4; check every sbatch header incl. agents

- **[NEXT SESSION START HERE: docs/handover/HANDOVER-2026-09-26.md (viper -> Caltech handover)](merge-finished-branches.md)** — data tarball /viper/ptmp2/jinma/caltech_handover_0926; production keys, open decisions and viper-only facts inside

- **[cs_mhd_prod4 STOPPED 09-25 (was running since 09-22)](cs-mhd-prod4-run.md)** — jobs 11941995 (link 1) + 11944873-5 (links 2-4, HSA_NO_SCRATCH_RECLAIM=1), bench/cs_mhd_prod4; what to check at rot 0.5; do not relaunch/chain without the user
- **[cs_hyd4_prod STOPPED 09-25 (hydro twin of prod4)](cs-hyd4-prod-run.md)** — jobs 11943719 (link 1) + 11944876-8 (links 2-4, HSA_NO_SCRATCH_RECLAIM=1), bench/cs_hyd4_prod; do not relaunch/chain without the user

- **[MASTER PLAN 09-22: wait for agents -> CLEAN UP -> lists A (dhj RT) / B (cs vertex) / C (production practice)](master-todo-2026-09-22.md)** — START HERE; in-flight jobs and agents listed inside

- **[M1 REOPENED 09-22 ~22:10 (was paused 09-21)](m1-paused-back-to-dhj.md)** — user: continue implicit M1 for massive stars; VET by short characteristics in progress (bench/m1_vet_0922)

- **[CLEANUP PAUSE after the agents finish (user 09-22)](code-cleanup-pause-plan.md)** — inventory of switches first; user decides abandon / default-on / keep

- **[dhj RANKED TO-DO 09-21 night](dhj-improvement-list-0921.md)** — seam fc halo, r^2 in the ck kernel (confirmed missing), MHD restart caches: 3 agents running, uncommitted

- **[START HERE 09-22: M1 radiation module on rt-integration (unpushed beyond cde8940e); read docs/handover/HANDOVER-2026-09-22.md; seeded 2-D He slab DIAGNOSED (runs_3b7/RESULTS.txt): the damage is made in the optically thin top (tau < 0.1) by the multi-D implicit transport, NOT by the frozen a_rad_ref; next = cure the thin top; login-node RAM cap 108 GB crashed the session 09-21](rad-m1-design.md)**

- **[MODEL SETUP 09-22: Opus 5.5 medium for session AND agents](model-setup-opus55-medium.md)** — trial; ALL agents use subagent_type "worker" (.claude/agents/worker.md, Opus 5.5 medium); supersedes the Opus 5 / Sonnet delegation lines below

- **[LOGIN NODE RAM CAP 108 GB KILLS THE SESSION (09-25, 2nd time)](login-node-ram-cap-pool-size.md)** — claude-work.slice 70 GB hook now guards all Bash; still cap analysis RSS <= 40 GB

Standing user rules (full text in [index-standing-rules](index-standing-rules.md)):

- **[NO BLOCKING WAITS](no-blocking-waits.md)** — no until/sleep loops on jobs in the foreground; submit, one quick look, return control
- **[NO INTERIM AGENT ACKS (user 09-23)](no-interim-agent-acks.md)** — ignore "agent still waiting" notifications; report only results, decisions, failures
- **[NO HEAVY MONITORING](no-heavy-monitoring.md)** — verify a launch once, then stop; no rolling watches
- **[HOUSEKEEPING = MINIMAL TOKENS (user 09-23)](housekeeping-minimal-tokens.md)** — copies/cleanup: launch once, act on the completion notice only
- **[JOB WATCHERS POLL EVERY 10 MIN (user 09-24)](job-watch-10min.md)** — one background squeue loop, sleep 600, notify on completion
- **[MERGE FINISHED BRANCHES INTO rt-integration (user 09-23)](merge-finished-branches.md)** — as soon as gates pass; then push to fork
- **[SBATCH SNAPSHOTS THE SCRIPT](sbatch-snapshots-script.md)** — cancel + resubmit after editing submit.sh; verify with scontrol write batch_script
- **[SAVE TOKENS](save-tokens.md)** — no unneeded work; narrow prompts, cheap models for lookups
- **[SAVE TOKENS everywhere](save-tokens-everywhere.md)** — especially on Fable; no surveys, agents one deliverable each, no timers left armed
- **[TEST STATE = prod4 RESTARTS (user 09-23)](test-state-use-prod4-restarts.md)** — dhj tests start from the newest bench/cs_mhd_prod4/rst restart (read-only copy), not prod3 (odd-even columns)
- **[GPU ENV IN EVERY JOB (user 09-23)](gpu-env-settings-all-jobs.md)** — every GPU sbatch exports the validated set: HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1 (+ whatever the env sweep validates)
- **[MEASURE COST ON GPU (user 09-22, tightened 23:40)](measure-cost-on-gpu.md)** — ALL timing comparisons on GPU (apudev, same binary, interleaved); CPU = correctness gates only; productions are all GPU
- **[TARGET PLANET = WASP-121b (user 09-25)](target-planet-wasp121b.md)** — all parameters from literature; R_p at the computed transit pressure; setup in wasp121_0925
- **[KELT-20b setup ready, not run (09-25)](kelt20b-setup-0925.md)** — mass 2.0 M_J (unmeasured); UV < 0.26 um 14 % folded into band 10; albedo 0.32 disputed
- **[CLOUDS desk study 09-25](clouds-desk-study-0925.md)** — ADAM = SPARC renamed; we lack scattering; cloud-free baseline OK; minimal nightside-deck test
- **[OVERNIGHT PLAN 09-25 (user asleep)](overnight-plan-0925.md)** — smoke -> submit sponge arms if clean; merge gated branches; morning summary only
- **[NO ACCURACY SACRIFICE for speed (user 09-25)](no-accuracy-sacrifice.md)** — only accuracy-neutral speed-ups; multi-rate etc. opt-in, never recommended
- **[THREAD TRIAGE 09-25: production path only (user)](thread-triage-0925.md)** — no new side threads without asking; mg phase 2 RESUMED as m1-coarse2 (user); ck-fast2 running
- **[DELEGATE to SONNET when adequate (user 09-22)](delegate-to-sonnet-when-adequate.md)** — ONLY really simple, single-step, binary-checked tasks -> sonnet (tightened 09-22 after the 3-attempt profiling job); everything else -> opus
- **[apudev FOR SUB-15-MIN JOBS](apudev-for-short-jobs.md)** — user rule 09-15; chained/production jobs go to apu
- **[VERIFY SPECTRAL VERDICTS](verify-spectral-verdicts.md)** — look at the maps; argmax of a broad spectrum is not "box-bounded"
- **[COMPARE WITH v_MLT AND F_conv/F](compare-with-mlt-not-old-runs.md)** — not with box_w4, which was 19x v_MLT
- **[USE rt-integration FOR EVERYTHING](use-rt-integration-branch.md)** — incl. He4; he4-presn-global merged as e688efd4 and frozen
- **[viper INODE QUOTA](viper-inode-quota-builddirs.md)** — worktree add fails while big writes work; delete build dirs in frozen worktrees

- **[ck x1 ONE-BLOCK CONSTRAINT: DEFERRED (user 09-24)](ck-x1-one-block-deferred.md)** — design A (column transpose) / C (affine interface) recorded; weak-scaling test at nx1 256 in /viper/ptmp2/jinma/cs_weak_0924
- **[NEXT dhj PRODUCTION: T4 + c2 + ck_impl_every=4 + kkt_row + xstep 8 (user 09-24/09-26)](next-prod-ck-c2.md)** — GOAL 300 rotations (~1060 Earth days); keys inside; all blockers cleared 09-24; setup needs user choices (grid, MHD/hydro, nx1, start, nodes)
- **[MHD RELAX: hydro spin-up + 10-20 rot C256 MHD is enough (09-25)](mhd-relax-hydro-spinup-ok.md)** — bbot=3 MHD = hydro within noise after rot 6; grid256 (8-coef, nx1 256/320 GPU) merged 587accf0
- **[GRID RESOLUTION COUNTS ONLY BELOW 1e-6 bar (user 09-25)](resolution-only-below-1e-6-bar.md)** — cells/H targets for p > 1e-6 bar only; top is free
- **[MASSIVE STARS: TOP ATMOSPHERE DOES NOT MATTER (user 09-26)](accuracy-region-interior-only.md)** — judge accuracy below the photosphere only (tau >~ 1, no top cells/sponge); crude top OK but it must not limit dt, crash, or cost much; = dhj p < 1e-6 bar
- **[SPONGE RUNS on WASP-121b 1x: KEEP BOTH SPONGES, recheck bottom at rot 20 (user 09-26)](sparc-sponge-campaign-0925.md)** — none stopped (dt collapse); top sponge buys dt 1.8x; floor cells at p>1e-6 bar same in all arms
- **[ck_nquad = 2 FOR dhj PRODUCTION (user 09-25)](ck-nquad2-production.md)** — exact diffusion limit; nq1 is 10 % low; validation branch ck-nq2
- **[ck CADENCE RETEST 09-26](ck-cadence-retest-0926.md)** — xstep<=4 no-op under every 4; xstep 8 -11 % OK; e8x16 -28 % borderline (needs 2nd noise member)
- **[M1 precond = mg merged 091d1422 (09-25)](m1-precond-mg.md)** — box -7 %/cycle, wedge ~0; boxes only, levels 3
- **[dhj DEEP IS CONVECTIVE below ~3 bar; old IC transient ~60 rot (09-25)](dhj-deep-convective-verdict.md)** — fix = exact-adiabat IC on the relaxed adiabat; no MLT needed

- **[NEXT M1 RUN: hesdirk2 cfl 0.9 (user DECIDED 09-26)](m1-hesdirk2-cfl-recommendation.md)** — be@0.3 accuracy, 1.57x cheaper; sp+vet_col supported (m1-sph2) but 0.9 validated only on the Cartesian box
- **[vet_col IS ENOUGH; no SC on sp (He 09-24, dhj 09-25)](vet-col-enough-no-sc-on-sp.md)** — dhj night photosphere f_rr 8 %, flux 1.5 %; SC ~60x cost, only f_rr usable; NO-GO
- **[M1 THIN-CELL CURE: implicit_closure_thin_relax = 1.5 (merged 09-24, default off)](m1-thin-cell-cure.md)** — same mechanism as the 09-22 He-slab thin-top damage; lp 0.5 + offdiag none does NOT cure it
- **[dhj KINK: horizontal rad diffusion RULED OUT (09-24, real profile)](dhj-kink-horizontal-rad-candidate.md)** — rad_angular adds zero flux where the kinks are; irrelevant at 128 and 1024 grids

## Active work (M1 radiation, He4)

Full original lines in [index-active-m1-he4](index-active-m1-he4.md).

- **[M1 RT PLAN 09-21](rad-m1-design.md)** — design note docs/dev/rad_m1_design.md; comoving-frame sources, AP-HLL, RSLA useless at depth, EOS gas-only
- **[M1 STAGE 3 SURVEY](rad-m1-implicit-gpu-survey.md)** — no published implicit two-moment RT on GPUs; matrix-free BiCGStab on Schur-reduced E + line preconditioner
- **[He4 session state 09-21](he4-session-state-2026-09-21.md)** — porous envelope drains in every variant; next = rt_kappa_hsmooth test in tests_r28
- **[He4 radiation redesign](he4-rad-eos-force-split-design.md)** — EOS weight density-only with a low window + optical-depth-gated FORCE in the two-stream
- **[He4 porosity](he4-porosity-luminosity-excess.md)** — 3-D wedge cools after ~2.9 turnovers; L_out/L 1.15; resolved; steady vs runaway open
- **[He4 Strang vs unsplit](he4-strang-vs-unsplit-mlt-seed.md)** — the MLT closure is seeded before the two-stream (unsplit only); Strang is correct; not fixed
- **[He4 r11 floor arms verdict](he4-r11-floor-arms-verdict.md)** — expansion thermally driven; sp wedge + ADI ported; rt_rad_force ghost bug fixed; 6 retractions
- **[He4 PRESN GLOBAL MODEL PLAN](he4-presn-global-plan.md)** — base red_giant.cpp; 2 blockers (cs ADI, no r^2 in the deposit); grid/cost/build order
- **[PRESUPERNOVA He STARS sized](he-presn-sizing.md)** — none fits a box (FeCZ 0.4-0.7 R, super-Eddington); needs global spherical RHD

## Sub-indexes (read the one for your topic)

- [Standing rules](index-standing-rules.md) — 12 entries: the full text of the user rules above.
- [Active M1 / He4](index-active-m1-he4.md) — 9 entries: full lines for the active radiation and He4 work.
- [Convection boxes: He, B star, solar](index-convection-boxes-he-bstar.md) — 41 entries: FeCZ boxes, two-stream/column RT development, f-modes, solar convection.
- [Cubed sphere](index-cubed-sphere.md) — 45 entries: seams, MHD, resistivity, validation, GPU/MPI history.
- [Hot Jupiter](index-hot-jupiter.md) — 92 entries: dhj campaigns, sp poles, correlated-k, grid/EOS/atmosphere.
- [Red giant](index-red-giant.md) — 56 entries: envelope runs, dt collapses, floors, opacity and RT bugs, FOFC.
- [Bugs and performance](index-bugs-performance.md) — 31 entries: code bugs, kernels, GPU profiles, build traps.
- [Working practice](index-working-practice.md) — 22 entries: machines, builds, job submission, output locations.
- [General EOS](index-general-eos.md) — 12 entries: the non-ideal/tabulated EOS project.
- [Old session states](index-old-session-states.md) — 16 entries: orion session notes and superseded state paragraphs.

## Session state

- [Session state 2026-09-18 (viper)](session-state-2026-09-18.md) — the 09-18 CURRENT STATE and NEXT STEPS paragraphs, verbatim.
