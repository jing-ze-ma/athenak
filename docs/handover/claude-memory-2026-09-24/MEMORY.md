## RULE ZERO

- **[NEXT SESSION START HERE: docs/handover/HANDOVER-2026-09-24.md (top block)](merge-finished-branches.md)** — 09-23 session: M1 defaults/fixes merged (rt-integration c65233b4), unmerged m1-launch / ck-fast / m1-thin, kink cause found, open decisions

- **[cs_mhd_prod4 RUNNING since 09-22 20:45](cs-mhd-prod4-run.md)** — jobs 11941995 (link 1) + 11944873-5 (links 2-4, HSA_NO_SCRATCH_RECLAIM=1), bench/cs_mhd_prod4; what to check at rot 0.5; do not relaunch/chain without the user
- **[cs_hyd4_prod LAUNCHED 09-23 (hydro twin of prod4)](cs-hyd4-prod-run.md)** — jobs 11943719 (link 1) + 11944876-8 (links 2-4, HSA_NO_SCRATCH_RECLAIM=1), bench/cs_hyd4_prod; do not relaunch/chain without the user

- **[MASTER PLAN 09-22: wait for agents -> CLEAN UP -> lists A (dhj RT) / B (cs vertex) / C (production practice)](master-todo-2026-09-22.md)** — START HERE; in-flight jobs and agents listed inside

- **[M1 REOPENED 09-22 ~22:10 (was paused 09-21)](m1-paused-back-to-dhj.md)** — user: continue implicit M1 for massive stars; VET by short characteristics in progress (bench/m1_vet_0922)

- **[CLEANUP PAUSE after the agents finish (user 09-22)](code-cleanup-pause-plan.md)** — inventory of switches first; user decides abandon / default-on / keep

- **[dhj RANKED TO-DO 09-21 night](dhj-improvement-list-0921.md)** — seam fc halo, r^2 in the ck kernel (confirmed missing), MHD restart caches: 3 agents running, uncommitted

- **[START HERE 09-22: M1 radiation module on rt-integration (unpushed beyond cde8940e); read docs/handover/HANDOVER-2026-09-22.md; seeded 2-D He slab DIAGNOSED (runs_3b7/RESULTS.txt): the damage is made in the optically thin top (tau < 0.1) by the multi-D implicit transport, NOT by the frozen a_rad_ref; next = cure the thin top; login-node RAM cap 108 GB crashed the session 09-21](rad-m1-design.md)**

- **[MODEL SETUP 09-22: Opus 5.5 medium for session AND agents](model-setup-opus55-medium.md)** — trial; ALL agents use subagent_type "worker" (.claude/agents/worker.md, Opus 5.5 medium); supersedes the Opus 5 / Sonnet delegation lines below

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
- **[DELEGATE to SONNET when adequate (user 09-22)](delegate-to-sonnet-when-adequate.md)** — ONLY really simple, single-step, binary-checked tasks -> sonnet (tightened 09-22 after the 3-attempt profiling job); everything else -> opus
- **[DELEGATE simple tasks to Opus 5](delegate-simple-tasks-to-opus.md)** — from the first tool call; Fable decides, Opus executes
- **[DELEGATE heavy work to Opus 5](delegate-heavy-work-to-opus.md)** — agents read/edit/build/launch/monitor; I diagnose, brief, and VERIFY every diff and number
- **[apudev FOR SUB-15-MIN JOBS](apudev-for-short-jobs.md)** — user rule 09-15; chained/production jobs go to apu
- **[VERIFY SPECTRAL VERDICTS](verify-spectral-verdicts.md)** — look at the maps; argmax of a broad spectrum is not "box-bounded"
- **[COMPARE WITH v_MLT AND F_conv/F](compare-with-mlt-not-old-runs.md)** — not with box_w4, which was 19x v_MLT
- **[USE rt-integration FOR EVERYTHING](use-rt-integration-branch.md)** — incl. He4; he4-presn-global merged as e688efd4 and frozen
- **[viper INODE QUOTA](viper-inode-quota-builddirs.md)** — worktree add fails while big writes work; delete build dirs in frozen worktrees

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
