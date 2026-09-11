# Memory Index
- **[VIPER HIP CONVENTIONS (standing, 09-11): DualView modify_device/sync_host only; kernel splits change GPU round-off; CPU bitwise != GPU bitwise; build on viper before calling GPU code done](viper-hip-code-conventions.md)**
- [VIPER floors_legacy GATE 09-11 (caad9247): orion merge changed GPU answers 1 ulp via the general-EOS c2p; legacy kernels in separate .cpp under floors_legacy (all floor switches off); MERGED locally 7cee2eb7 on 09-12](viper-floors-legacy-gate.md)
- **[STANDING RULE: delegate to Opus 5 agents](delegate-heavy-work-to-opus.md) — EVERY session from the first tool call: Opus agents do reading/editing/building/launching/MONITORING/analysis; I do diagnosis, briefs, VERIFY every diff and number; ONE deliverable per agent, name the hunks to read, "STOP, no timers" (token blow-up 09-11)**
- [SESSION STATE 2026-09-12 01:00](session-state-2026-09-12-0100.md) — START HERE: FOFC committed locally + viper merged 7cee2eb7 (NOT pushed); gates pass except dhj ck 1-ulp from the RT guard e1db81d8; RG_fofc running, no events at 1.5e5
- [RG_fofc LAUNCHED 09-11 15:30: red giant from scratch nghost=3, fofc on, NO vceil, floor fixes on, chain 196662-4; deaths to beat 5.7e5/5.87e5/6.13e5](red-giant-rg-fofc-run.md)
- [SESSION STATE 2026-09-11 12:50](session-state-2026-09-11-0800.md) — START HERE: a10e367d pushed; FOFC 0+1 committed c0163568 (not pushed); I8 done 9e5; CHECK C1 first; uncommitted RT guard (keep) + vacuous dfloor_keep_temperature (drop); decision: top-cell negative density -> FOFC step 3
- [TOP-CELL NEGATIVE DENSITY 09-11 12:40 = root of the RT -nan (rescue was an accidental energy sink) and likely of the 1e16 K vacuum cells; RT guard in tree UNCOMMITTED; with it the implicit conduction fails at 6.1265e5 -> needs FOFC step 3 (cs) or a vacuum/top-cell treatment; user decides](red-giant-top-cell-negative-density.md)
- [DFLOOR-ENERGY DIAGNOSIS WRONG 09-11 11:10: the 1e16 K dt-dip cells are VACUUM cells (rho_old<=0, v=-0) whose energy GROWS at v=0; dfloor_keep_temperature (e*=fv) vacuous and UNCOMMITTED; proposal: vacuum reset + instrument the energy source; user decides](red-giant-dfloor-energy-diagnosis-wrong.md)
- **[SAVE TOKENS (standing rule 09-11): especially on Fable; no unnecessary tool calls, no surveys, agents one deliverable each, no timers left armed](save-tokens-everywhere.md)**
- [FOFC COMPATIBILITY COMPLETE 09-11 incl. MHD on cs + static/dynamic WB tests (new pgen wb_atm; dynamic WB: fofc raises the Mach ceiling 39 -> >80) + MPI/SMR tests; COMMITTED b3345f55 + merged 7cee2eb7 on 09-12, not pushed; 49 fofc tests pass, bit-identical off; left: red-giant fofc run, GPU check, wiki, commit split from red-giant hunks](fofc-compatibility-plan.md)
- [VCEIL WORKS 09-11 07:00: hydro/vceil=5e7 + RT NaN guard; R4 passed the 5.87e5 death; the "hang" was runaway_scan's rank-local line budget skipping an Allreduce (FIXED, pin11, ghost-cell counter fix too); I8_vceil 196221 clean at 6.26e5 -> 9e5; next physics = denser corona](red-giant-vceil-r4-passes-then-hangs.md)
- [R3 RESULT 09-11 07:20: dfloor_keep_velocity correct but NOT the loop (R3_bface died like R2 at 5.8667e5); cells at 1e5 x escape speed, KE-cancellation energy creation + NaN T in RT from negative ei -> NEXT: Newtonian hydro/vceil + RT NaN guard, R4](red-giant-r3-dfloor-keep-velocity-not-the-loop.md)
- [VIPER MERGE LANDED 09-11 08:05: a10e367d on top of 9 topical commits; gates bit-identical both sides; not pushed](viper-merge-trial-2026-09-11.md)
- [Session state 2026-09-10 22:00](session-state-2026-09-10-2200.md) — superseded by the 09-11 03:40 note: implicit radial diffusion + angular cap gated; I2 alive past both deaths; I3 (clamp off) and I4 (400 K active medium, from scratch) launched with watcher agents; tree uncommitted (~29 files)
- [STALE-PRIMITIVE AUDIT (09-11 02:40): 6 more sites of the RT bug class; LIVE: arm rt_use_cons, RedGiantBC ghosts seeded from stale w0 (at the dying faces), rg_wallflux unguarded ei; LATENT: rg_relax, open-inner BC (incl. orthogonal KE on cs), MLT caps, perturb-reconstruction ordering](red-giant-stale-primitive-audit.md)
- [RUNAWAY SOURCE (09-11 03:30): RT read STALE w0 -> fixed (rt_use_cons, floor energy creation gone, R2); loop reopens via the DENSITY FLOOR keeping momentum (v 4e11 in a floored cell) -> dfloor_keep_velocity in progress (R3/I7); I6 DIED the same way at 5.877e5 (v 5e16 in a dfloor cell)(red-giant-runaway-source-rt-stale-w0.md)
- [RT SWITCHES DONE 09-11 06:30: rt_bface + 4 RT defaults now OLD by default, rt_top_clamp new; ck dhj bit-identical to HEAD; red giant inputs need SIX opt-ins from pin10 on (repo input has them, ptmp run inputs do not)](red-giant-rt-bface-switch-and-head-regression.md)
- [I5 DENSITY GATE DONE (09-11 01:40): ran to 9e5, cutoff artefact CONFIRMED (L_out 1.2-3.3 vs 2.6-5.3 with the cutoff, L_cut 0.06 vs 1.36); 0.35 L at restart = clamp off, not the gate; corona inert; NEW: photospheric cells reach 1e16-1e17 K under flat dt -> runaway source (hydro or RT) still open(red-giant-i5-density-gate.md)
- [I4 ACTIVE 400 K MEDIUM: passed V10 death, top never floored, then DIED 6.31e5 of a hydro collapse at an ejected-shell front (all open-top runs launch a shell at 5.7e5; the pressure-supported corona lets I2/I3 ride it out); next = DENSITY-GATED opacity (corona inert, no radius cutoff)](red-giant-i4-active-medium.md)
- [3 L EMERGENT FLUX = OPACITY-CUTOFF ARTEFACT (09-11 01:00): the inflating surface crosses rad_kappa_rmax where kappa=0 and radiates 9000 K gas as a blackbody; lidded top is a uniform 9200 K/1e-8 at 1.3-1.5 L; "+34 L dE/dt" is float32 dump noise; I4 (no cutoff) is the test](red-giant-open-top-kappa-cutoff-artefact.md)
- [OPAC CLAMP UNNECESSARY OPEN-TOP TOO (I3, 09-10 23:30): I2 vs I3 indistinguishable to 6.2e5, top 50-90 K warmer without it, fewer cold cells; open: 31.6 K floor cells in both after 6e5, L_out/L ~2.7 after the cap burst](red-giant-clamp-off-open-top-ok.md)
- [COMMON DEATH 5.71e5 (V9f/V10/V9_dense): EXPLICIT radiative conduction runaway in an underdense supersonic downdraft at 1.058 R on a cube EDGE, inside the star; corona/top innocent; next test = sub-cycled/implicit radial conduction from rst 5e5](red-giant-common-dt-collapse-571e5.md) — START HERE 09-10 18:00
- [V8f 6e5 K CORONA DROPPED (09-10 18:20): it is a Parker wind (T 3x above the critical T at the join, H = 1.9 R); the 14-scale-height atmosphere drains, over-pressured corona shocks down, dfloor cells with no velocity ceiling spike -> hydro dt 0.95 s at 3.2115e5; enable eos_vceil in any open-top run](red-giant-hot-corona-is-a-parker-wind.md)
- [two-stream CANNOT own the interior (09-10 19:00): deposition throttled to src*t_rad deep; core L is a conduction wall term; killer cell is tau_cell~1 in a tau>100 column -> diffusion operator invalid there; cure = implicit radial diffusion (one block radially, no MPI) or local-tau blend](red-giant-two-stream-cannot-own-interior.md)
- [IMPLICIT RADIAL DIFFUSION + ANGULAR CAP DONE 09-10 23:00 (hydro/rad_implicit_x1, rad_cap_ang=0.5): I2 ran 5e5 -> 9e5 with dt flat, tot-E flat, capping only 5.69-5.98e5 then silent; +21% cost; binary build_impl/athena_pin2](red-giant-implicit-radial-diffusion.md)
- [Session state 2026-09-10 09:10](session-state-2026-09-10-0910.md) — superseded by the 22:00 note: WB-walk overflow + RT neighbour-Planck bugs FIXED and gated (build_guard); lidded run needs no clamp (B13); V9f/V10 open-top tests running to ~12:40; commit plan; user decides prod11 revival / open-top config
- [Session state 2026-09-10 05:30](session-state-2026-09-10-0300.md) — superseded by the 09:10 note; keeps the FL/V4/V5/V6/B-series detail
- [Session state 2026-09-10 early](session-state-2026-09-10-early.md) — superseded by the 03:00 note; keeps the FL/V4/catch job list and the agent transcript paths
- [WB POLYTROPIC WALK OVERFLOW = ROOT CAUSE of the prod11 NaN AND the "WB kick at a cold cell": FIXED in the working tree 09-10 (4 hunks), bit-identical on healthy data, B11 passes the death with WB ON and the sawtooth is gone; binary build_guard](red-giant-wb-polytropic-walk-overflow.md)
- [RT NEIGHBOUR-PLANCK BUG FIXED+GATED 09-10 09:00: BFace kappa-weighting at 16 sites; smooth star bit-identical; B13 (lidded, NO clamp, NO wb_rmax) holds dt flat past the old collapse with cold cells 0-27 vs 156 -> the cold-collapse family WAS this bug; V8f join cell holds 3300 K](red-giant-rt-neighbour-planck-bug.md)
- [CORONA V8 DIED 2.96e5 (09-10 07:50): the LAST ACTIVE CELL below the kappa cutoff cools to the 31.6 K floor in 2e4 s, the corona accretes onto it and evaporates; the RT "Newton rescue = 99.9% drop" message precedes every open-top/corona failure (V7, V8) -> suspect; V9 (2e4 K) follows R7's atmosphere-cooling path](red-giant-corona-join-cell-cooling.md)
- [OPEN-TOP prod12 RECIPE FAILS FROM SCRATCH (V7, 09-10 06:50): cold collapse of the thin atmosphere at 5.8-6.0e5, EARLIER than V4; opac_tmin=3200 thins the open atmosphere 19x; corona V8 is the remaining open-top option; lidded B12 is healthy](red-giant-open-top-prod12-recipe-fails.md)
- [FLOOR DEATH 1.43114e6 CLOSED 09-10 02:40: efloor_from_ekin alone cures it (FL2) AND wb_rmax=3.3e12 alone cures it (FL4, floor energy still created but harmless); ship BOTH; ported to main tree + random seed added; smoke tests V5/V6](red-giant-floor-energy-creation-fix.md)
- [CHART-FREE SEED: vertex asymmetry GONE (V4 copies bit-identical) but V4 still DIED at 9.45e5 by an INTERIOR cold-collapse -> conduction-dt (no floor event); seed buys 2x time only; random component added (V6)](red-giant-vpert-seed-chart-imprint.md)
- [prod11: ROOT CAUSE FIXED and the LIDDED REVIVAL GATE PASSED (B12, 09-10 06:55: fixed binary + opac_tmin=3200 from rg.00124 to 3.45e7, dt flat, 0 vertex cold cells; cold population grows slowly, recycles, no dt response); revive from B12_fix_tmin3200/rst](red-giant-prod11-died-3e7.md)
- [MOLECULAR-OPACITY KNEE RUNAWAY (09-10 00:45): open-top cells cool radiatively through the kappa_R knee at 3200-2500 K (emission grows as T falls), park at 200-500 K, become the Riemann precursor; lid keeps the atmosphere 1000 K hotter; TEST opac_tmin~3200](red-giant-molecular-opacity-knee-runaway.md) — START HERE
- [HYDRO TRIGGER CAUGHT (09-09 23:30) + CENSUS 00:15: a COLD-COLLAPSED cell in the thin layer above the 1-cell front (any angle; vertex only first) + vacuum-Riemann energy creation; radiation only amplifies; open-top atmosphere is grid-scale noisy (x3 adjacent-cell rho) and grid-locked](red-giant-vertex-floor-hydro-trigger.md) — START HERE
- [Session state 2026-09-09 night](session-state-2026-09-09-night.md) — superseded by the 09-10 EARLY note: prod12 gate FAILED (see the explosion note), nothing running except prod11 chain 194625; the user must choose: bigger domain, stay lidded, or implicit diffusion
- [VERTEX CHIMNEY dissected 09-09 22:30: sponge+open top lethal, sponge-off survives, wall caps; T15 = prod12 gate sponge-off to 3e6](red-giant-vertex-chimney.md) — START HERE for prod12
- [PHOTOSPHERIC RUNAWAY: photosphere IS inside (1.077 R, 20 cells under the top; the 'outside' claim was a misread); 1-cell ionization front; every blend/two-stream-semi variant dies ~1.84e6, EXPLICIT two-stream sails through and dies at the top cell instead; D1_catch budget pending](red-giant-explicit-conduction-explosion.md) — START HERE 09-10 evening
- [Session state 2026-09-09 evening](session-state-2026-09-09-evening.md) — superseded by the NIGHT note; keeps the T11/T10 details and the earlier agent transcripts
- [cs orthogonal-KE AUDIT: red_giant ghost fills + IC + inner passes FIXED in the working tree 09-09 (T14: correct, not the killer); latent sites listed in docs/handover/NOTE-2026-09-09-cs-covariant-basis.md (pushed 40aef0ad)](cs-orthogonal-ke-audit-2026-09-09.md)
- [SPONGE KE BUG on the cubed sphere: real, fixed (887c241e), but NOT the killer — T11 died the same way at 8.72e5](red-giant-sponge-cs-kinetic-energy-bug.md) — START HERE 09-09 17:50: the lid-free deaths at ~1e6 were the sponge draining the cube-vertex cells; not RT, not the ghost, not dust (that was the EARLIER deaths)
- [T7 (grains off + WB fix) STILL died at 1.04e6 with no precursor: the 100-cycle NaN check hides the origin (dt ignores NaN); T8 194637 reruns with the check every cycle — OPEN](red-giant-nan-check-hides-origin.md) — START HERE 09-09 17:00
- [RESTART x730 radial-KE kick: FOUND+FIXED (WB cache zero for ncycle%N cycles after restart; affects ALL well-balanced restarts, dhj too), COMMITTED 5c0b98e4, pushed; viper must pull](red-giant-restart-radial-ke-injection.md) — START HERE 09-09 16:20: prod11 chain 194624/5 with wb_cache_every=1; T7_nodust_fix 194627 is the prod12 gate
- [DUST OPACITY kills lid-free runs: PROVEN, opac_tmin=2500 survives; explicit RT is a negative control](red-giant-dust-opacity-kills-lidfree.md) — START HERE 09-09 14:40: T6 extended to 3e6 (194594); prod12 candidate = lid-free+open+grains off; prod11 lidded control at 4.1e6
- [PRODUCTION prod11 LAUNCHED 06:40: 1.1 R, all fixes, chain 194515-7](red-giant-prod11-launched.md) — START HERE next session: check dt, L_out, deep quiet, surface onset; old chain cancelled; option A (hot corona) proposed, not approved
- [GREY OPAQUE-LID BUG: grav=0 on the grey path -> dtau_top = inf in every grey run incl. prod11](red-giant-grey-opaque-lid-bug.md) — fix doubles L_out at t=0 but the thin atmosphere cools 3364->2200 K; being isolated on the plain star (R8); corona itself was inert and fine
- [Session state 2026-09-07 (orion)](session-state-2026-09-07-orion.md) — START HERE: pulled viper's handover (5e9b1273); read docs/handover/HANDOVER-2026-09-07.md; steps 1-2, the determinism test and the EOS table all DONE; left: GPU build with d3d74f2b (needs viper)
- [Orion build traps: OpenMP-off binary + module load](orion-build-openmp-and-module-traps.md) — a Kokkos_ENABLE_OPENMP=OFF build is SILENT and 3.7x slow at 96x7; `module load` does not survive between tool calls, so make falls back to gcc 7.5
- [Red giant: the 0.34 L deficit is a SPIN-UP transient, not a bug](red-giant-flux-deficit-is-spinup.md) — CHECKED 09-08 23:20 at t=5.5e6: the convection front IS advancing outward ~25 cells/1e6 s, deep runaway saturated; chain 194259+194262-5 healthy, just needs wall time
- [Red giant: the ANALYSIS used the wrong opacity](red-giant-analysis-opacity-trap.md) — eoslib.get_kapr is solar_convection's analytic formula, 30-44x too high in the interior; the RUN is clean; retracts the "wall cannot shed the flux" story
- [RT source in a transparent cell: RESOLVED](rt-transparent-cell-cancellation.md) — direct source in all 3 kernels + the REAL cause: the opacity lookup CLAMPED rho to the table edge (1e-14) and returned molecular opacity; table now extended to 1e-24 at the floor; background holds 411 K
- [Red giant DEEP motion = RCB pile-up transient](red-giant-deep-onset-is-rcb-pileup.md) — L switched on with no convection piles heat above the RCB; sup grows linearly from t=0, ignites at 7e5 at 200x MLT driving, overshoots 10x; surface convects first and correctly; remedy under test: mlt_alpha=3
- [Red giant: the ambient medium is killed by the WELL-BALANCED scheme](red-giant-wb-kills-ambient-medium.md) — WB on: NaN in 2 cycles, ~400 g spurious force; WB off: background free-falls at exactly g and survives. RT/BC/opacity fixes were all red herrings for THIS failure. Fix: <hydro>/wb_rmax
- [two_stream_rt has THREE sweep kernels (grey/ck/generic)](two-stream-rt-three-kernels-trap.md) — a sweep edit must hit all three; red_giant uses rt_chain_grey; verify with rt_apply_debug. Also: <hydro>/wb_rmax added+verified
- [MLT closure lessons: per-column MLT flickers deep and cannot hand over at the top](red-giant-mlt-closure-lessons.md) — R1/R1b/R1c all failed at the photosphere; fix = shell-averaged 1D MLT with f = min(F_MLT, L - F_rad,diff) (problem/mlt_mean)
- [PROVEN: the vpert seed across the RCB causes the deep transient](red-giant-seed-across-rcb-proven.md) — vpert_rmin=1.6e12 flattens the dipole with/without MLT; nquad=2 killed both no-MLT runs at 1.44e6 (bisecting); MLT gap = chi limiter
- [1.5 R grid STARVED the deep: 0.92 cells/H_p at the wall -> dies at 3.2e5](red-giant-15R-grid-starved-deep.md) — ambient medium was fine; refit with >=5 cells/H_p in the star and more nx1; measure on the actual grid, not the fitter
- [Quadrature 0.904x, RCB dipole = seed?, pin binaries](red-giant-quadrature-seed-binary.md) — ck_nquad=1 is 10% below diffusion (use 2); the RCB dipole is closure-blind (R0==MLT to 4 digits) -> testing vpert_rmin; never rebuild under a running job
- [Red giant: open outer BC sealed the star, FIXED bd2b974d](red-giant-open-outer-seals-the-star.md) — `problem/rt_top_re` makes the column above the domain RADIATE and the 6650 K seal is gone; its "the defect is the handover" conclusion is SUPERSEDED by the spin-up note above
- [Red giant: the GREY two-stream, rebuilt](red-giant-grey-two-stream.md) — bd156d77 `problem/rt_grey`: one band, the conduction module's opacity; removes the thermal runaway and the runs live past every ck death point. The drain/outer-BC follow-ups are CLOSED by rt_top_re
- [Red giant: the TOP cools away and kills the run](red-giant-top-cooling-runaway.md) — historical: all three jobs died of a HYDRO dt collapse; the optically thin top cell falls 3364 -> 1100 K and the cooling GROWS as it cools; the emergent-flux decay is this, not the interior
- [Red giant: dt collapse SOLVED](red-giant-dt-collapse-solved.md) — the t=1.9e5 collapse was the EXPLICIT radiative diffusion at a cooling photospheric cell; a DEEPER tau blend (10/100) hands it to the semi-implicit two-stream and dt stays 30.7 s
- [Red giant session 2026-09-08](red-giant-session-2026-09-08.md) — the 10 commits, the working production config, and the orion SLURM core-vs-hyperthread trap; its jobs are all finished/cancelled
- [Correlated-k shared module](correlated-k-shared-module.md) — DONE 2026-09-07 in TWO stages: correlated_k.hpp (opacity+Rosseland/conduction coupling) and two_stream_rt.hpp + atm_column.hpp (the solver); pgen 6082 -> 3648 lines; gated BITWISE identical
- [Red giant envelope project](red-giant-envelope-project.md) — 2026-09-07: pgen builds and RUNS on the cubed sphere with the CORRELATED-K two-stream (problem/rt_ck); stellar opacities merged (OPLIB+AESOPUS); grey 1-D validation passed; COMMITTED+PUSHED (35ea57ba, 1d99e22e, 361be207)
- [sp polar bisection on orion](sp-pole-bisect-orion.md) — jobs 193408-11; the production grid RECONSTRUCTED (nx2=64 f_stretch_theta=3 verified against the 8.1 deg polar cell); dt still 7x off, so the CONTROLS decide
- [cs determinism test on orion](cs-determinism-test-orion.md) — DONE 193377: CPU restarts are BIT-IDENTICAL, even 16x7 vs 28x4; the viper 5e-6 divergence is GPU-specific or a binary difference

## Orion-side notes (this machine's history)
- [Session state 2026-09-04](session-state-2026-09-04.md) — START HERE: artifact PUBLISHED with the time player + ray-marched convection plumes; nothing pending
- [solar_convection restart gravity bug](solar-convection-restart-gravity-bug.md) — FIXED/pushed ef9561e2: restarts ran with zero gravity, dt collapsed to 1e-135
- [solar_convection scaling wall](solar-convection-scaling-wall.md) — MEASURED: saturates at ~36 ranks; run it on ONE node, multinode buys 6%
- [solar_convection binary provenance](solar-convection-binary-provenance.md) — the pinned ideal192/table192 binaries cannot be rebuilt from git

- [Session state 2026-08-25](session-state-2026-08-25.md) — superseded: job 190939 was left running that day; the viper doc's blow-up does NOT reproduce on orion CPU; input pushed as cd2e8815

- [Session state 2026-08-24](session-state-2026-08-24.md) — superseded: the ckrepro reproducers ran clean; hst is unusable here (cadence + Cartesian volume)

- [Session state 2026-08-18](session-state-2026-08-18.md) — superseded entry point: no code changed today; all cost questions answered; the ONE open thread is the viper dt discrepancy
- [dhj viper dt discrepancy](dhj-viper-dt-discrepancy.md) — LARGELY DISSOLVED: viper bottoms at ~2.6 not ">3", and V_old (viper commit 2af153a3) reproduces it; only the t~1e6 dip is still unchecked
- [Session state 2026-08-17](session-state-2026-08-17.md) — previous entry point (origin/general-eos = d1c289dd)
- [dhj ideal vs general cost](dhj-ideal-vs-general-cost.md) — general is 3.09x ideal per simulated second (10 G, done); max_eta=1e13 beats 1e14 at both fields. Its 3 G runs were CANCELLED and superseded by the viper-base set
- [dhj dt limited by Alfven floor](dhj-dt-limited-by-alfven-floor.md) — MEASURED: dt is the radial Alfven CFL at r/Rp~1.28 in a cell sitting ON dfloor; raising dfloor x10 buys 2.7x
- [Freya build procedure](freya-build-procedure.md) — how to build AthenaK on Freya (MPI+OpenMP, SPR); clean only build/ contents
- [Solar convection test](solar-convection-test.md) — tuning solar_convection pgen (two-stream RT) for solar-like convection + stable atmosphere
- [Freya/Orion job submission](freya-job-submission.md) — p.shared vs p.exclusive, keep --cpus-per-task with OMP, and an orion node has 112 PHYSICAL cores (the 224 SLURM reports are hyperthreads)
- [Work pace preference](work-pace-preference.md) — move faster / less deliberation on routine tasks
- [run/ is untouchable](run-directory-untouchable.md) — ~40GB of simulation data; never add, clean, or modify it
- [General EOS project](general-eos-project.md) — non-ideal EOS for Newtonian hydro/MHD; branch `general-eos`; ALL stages including the analytic EOS are done and verified. START HERE.
- [General EOS: Stage 3 tabulated EOS](general-eos-stage3-table.md) — the analytic EOS itself: H2 + Saha + radiation via a (log rho, log T) table; done and verified, commit 03febafa
- [General EOS: Stage 3 perf cleanup](general-eos-stage3-perf-cleanup.md) — one T solve per cell instead of ~11-30; fully verified and committed as 7d1fed3f on 2026-08-12
- [AMR null-tree bug](amr-broken-on-branch.md) — fixed; the cubed-sphere refactor left `Mesh::ptree` unconstructed, and the same half-migration may bite elsewhere
- [FOFC 1D segfault](fofc-1d-segfault.md) — pre-existing crash on the branch: FOFC + 1D hydro + PLM, not caused by the EOS work
- [Scratchpad invisible to compute nodes](scratchpad-not-visible-to-compute-nodes.md) — /tmp is node-local; stage SLURM jobs on /orion
- [GPFS quota wall](gpfs-quota-wall.md) — df lies; a per-user quota SILENTLY truncates AthenaK dumps without any error
- [Fork git setup](athenak-fork-git-setup.md) — origin = jing-ze-ma/athenak fork over SSH (all git ops); upstream = IAS-Astrophysics; the user's terminal cannot copy OUT
- [General EOS: MHD+polar difference](general-eos-mhd-polar-bug.md) — resolved, not a bug: the polar boundary forces HLLD→HLLE at the pole, and HLLE is the solver whose Roe average has no general-EOS analogue
- [dhj restart bugs](dhj-restart-loses-gravity-potential.md) — both FIXED (30d21859, 4cfdc329): restarts lost the gravitational potential and the cell-centered field bcc0
- [eostest input inventory](eostest-input-inventory.md) — what the ideal/general comparison inputs cover and how to run them; .bin output is single precision
- [General EOS: table linear-wave tests](general-eos-table-linear-wave.md) — how general_eos=table is regression tested; the tabulated EOS cannot be called from host code
- [General EOS: Stage 3 loose ends](general-eos-stage3-loose-ends.md) — the last three TODO(stage3)s closed; why the WB background hands back a TEMPERATURE and not a pressure
- [Pre-existing branch breakage](branch-preexisting-breakage.md) — non-MPI build fails in bfield_bcs.cpp; hydro/mhd_linwave test scripts pass a `vflow` param that no longer exists
- [solar_convection with the general EOS](solar-convection-general-eos.md) — OPEN: the atmospheric runaway survives the corrected sponge; sub-photospheric comparison is done and solid. RESUME HERE for this thread.
- [EOS picked by run name](eos-selection-by-run-name-bug.md) — analysis scripts guessed the EOS from the directory name and silently used the wrong one; fixed, plus the audit that found the scope
- [Sponge inert at 192^3](sponge-inert-at-192.md) — RESOLVED: the sponge read pcoord->x1v, a 1x1 placeholder on Cartesian meshes; the 64^2 sponge results are invalid
- [pcoord arrays unallocated in Cartesian](pgen-cartesian-coordinate-arrays.md) — which geometry arrays are placeholders on a Cartesian mesh, and the solar/cooling_convection audit for unguarded uses
- [EOS x_e inert in the resistivity](eos-xe-resistivity-capped.md) — RESOLVED: not the EOS at all; two out-of-bounds reads of eta_b, one of which made ohmic_resistivity=constant a silent no-op
- [EOS electron regression test](eos-electron-regression-test.md) — DONE: mhd_eos_electrons.py, plus how to dry-run a test without letting run_tests.py delete build/
- [Perna resistivity for UHJs](resistivity-perna-uhj.md) — good to ~3x below 5000 K; the min_xe→max_eta rename that aborts the dhj runs; a tabulated x_e would be 6x FASTER
- [dhj high-B crash is the outer BC](dhj-highB-outer-bc.md) — FIXED bc6b5774: RKG super-time-stepping already fixes the diffusive dt; the high-B crash was the outer-x1 Maxwell term, not the Alfven speed
- [General EOS table cost](general-eos-table-cost.md) — RE-MEASURED 2026-08-17: 3.07x per cycle / 2.32x per simulated second, NOT 4.4x; profile (49% libm) still stands. Restrict any such comparison to the hydro-limited window
- [General EOS Stage 4: (rho,e) table](general-eos-stage4-rho-e-table.md) — SCOPED, NOT STARTED; re-baselined 2026-08-17 to ~1.5x (3.07x -> ~2.0x), proxy ceiling measured twice; task 1 is the (rho,e) domain shape
- [x_e table for ideal-gas runs](resistivity-xe-table-for-ideal.md) — DONE 5d1c3436: ohmic_resistivity=eos under eos=ideal; 7% FASTER than perna, plus two latent bugs fixed
- [dhj ideal + EOS x_e: floors, max_eta, STS](dhj-ideal-xe-floor-relaxation.md) — CLOSED, shipped as 6003c1ce: max_eta=1e13 + STS off is 1.54x faster at the same dt; STS is NOT a standing win at 1e14
- [Event log MPI deadlock](eventlog-mpi-deadlock.md) — FIXED 27da6380: file_type=log hangs ANY clean MPI run (early return skipped the schedule advance -> Allreduce every cycle); likely upstream too
- [Exo-FMS correlated-k tables](exo-fms-ck-tables.md) — installed and md5-verified in the repo at data/exo_fms_ck/; ROMIO cannot write AthenaK output on /tmp
- [Test output location](test-output-location.md) — all run output/analysis goes in /orion/ptmp/jinma/Athenak/ (separate GPFS, no quota); never /orion/u or /tmp
- [dhj conservation check](dhj-conservation-check.md) — mass +0.21%, internal E +0.56% over t=4.32e5; the hst cadence and weighting both make hst useless here
- [hst uses Cartesian cell volume](hst-cartesian-volume-on-spherical.md) — .hst mass/tot-E are NOT physical integrals on the spherical dhj mesh
- [dhj general EOS + correlated-k blow-up](dhj-general-eos-ck-blowup.md) — OPEN: read docs/HANDOFF_dhj_ck_eos.md IN THE REPO first (9 refuted hypotheses); TWO blow-ups, only the general-EOS one is live
- [RT clip warning is rank-local](rt-srclim-warn-rank-local.md) — BUG, unfixed: the clip count is rank 0's share, and the warning is silent when rank 0 never clips

## Viper-side notes (imported 2026-09-07 from docs/handover/claude-memory-2026-09-07/)
## Cubed sphere — start here

- **[cs RADIAL = sp RADIAL, DONE 979edada: centroid everywhere, x1 Grid-PLM always, full audit](cs-radial-unification.md) — every cs baseline before it is slightly stale; rcmfix resubmitted as 11442863**
- **[cs STRETCHED-grid SOURCE TERM off by 0.57-2.15x: FIXED 13a97399 + the full sp-vs-cs source AUDIT](cs-stretched-source-term-bug.md) — index-space dr and r in the angular-momentum curvature source; cs_test was unstretched too. Affects every stretched cs dhj run, hydro included**
- **[cs STRETCHED-grid RESISTIVITY was ANTI-DIFFUSIVE: FIXED c5c85e3b](cs-stretched-resistive-rcm-bug.md) — r_cm never got the stretch, sign flip near the top; the production dhj config. Every pre-fix cs resistive MHD result is suspect. START HERE for the deep sheet**
- **[cs seam CO-LOCATION: the binding term, FOUND and FIXED](cs-seam-colocation-fixed.md) — the clamped-window cubic; 5-11x on the evolved field. START HERE for the seam halo**
- **[cs seam order: NOTHING limits it, 2nd order CONFIRMED](cs-seam-order-limiter.md) — global L1 2.01; never read an order off the fixed-cell region bins. The seam is DONE**
- **[cs seam: de-staggering fix REFUTED](cs-seam-destag-refuted.md) — superseded by [[cs-seam-colocation-fixed]]; kept for the trail of what was refuted**
- **[resistive seam 1st order: MECHANISM CLOSED](cs-resistive-seam-order.md) — the operator is innocent (1.98), the halo INPUTS were guilty; six committed instruments**
- **[cs cross-level SEAM halo: CLOSED by two fixes](cs-crosslevel-seam-halo-first-order.md) — a01ace75 + 180a9b3e; a level boundary no longer degrades the seam. START HERE for cs+SMR halos**
- **[cs CUBE-VERTEX corner halo x RADIAL GHOST: FIXED](cs-cube-vertex-corner-radial-ghost.md) — 397b4ad3; the whole "resistivity across a radial block interface" defect. START HERE for cs halo bugs**
- **[Cubed sphere: seam conservation CLOSED](cubed-sphere-seam-conservation.md) — 985faa22 +; mass and energy exact to round-off, and history.cpp used the CARTESIAN volume on every grid**
- **[Cubed sphere: seam EMF CLOSED](cubed-sphere-seam-emf.md) — 42323a66 +; exact at 1/4/16 blocks per panel. div B is the WRONG gate. START HERE for cs+MHD**
- **[Cubed sphere: RESISTIVITY](cubed-sphere-resistivity.md) — 4bfacdd8; gnomonic two-pass curl, cs_test iprob=11, 2nd order; refinement supported since c4735f18. START HERE**
- **[Cubed sphere SMR: refined MHD CONVERGES and is MPI-clean](cubed-sphere-smr.md) — 74cbc8df; the edge FLUX buffers had no seam transform. "Rank dependence" RETRACTED. START HERE for cs+SMR**
- **[Cubed sphere MHD convergence](cubed-sphere-mhd-convergence.md) — e63b571a; mhd_corner_e had no cs form. Residual is the radial-BC phase lag, one law in dt**
- **[SHOCKS through cs seams](cs-shocks-through-seams.md) — structure ON a seam NaN'd everything SILENTLY; the unlimited along-seam resample, now clamped. FOFC DECIDED not needed. START HERE**
- **[cs blast vs a CARTESIAN grid](cs-blast-vs-cartesian.md) — the seam costs NOTHING (4.1% vs Cartesian's own 4.5%); only the cube VERTEX is ~2x worse. START HERE for "is cs good enough for shocks"**

## Cubed sphere — validation and limits

- **[cs_test FF decay + rot_axis (4990eb41): cs matrix vs sp, no 1.5-order region on cs](cs-test-ffdecay-rotaxis.md) — START HERE for cs-vs-sp simple tests**
- **[cs PURE HYDRO validation](cs-hydro-validation.md) — mass and energy to MACHINE PRECISION even under a shock; space 2.3-2.4 L1. Records the FALSE first-order reading. START HERE**
- **[cs MHD + RESISTIVE validation, BY REGION](cs-mhd-validation.md) — RE-MEASURED: ideal MHD 2nd order EVERYWHERE (the old 1.84/1.76 was a radial floor); the resistive energy drift is PHYSICAL Ohmic heating, NO BUG**
- **[cs ANGULAR MOMENTUM: why sp is exact and cs cannot be; RE-MEASURED 09-06](cs-angular-momentum.md) — rigid rotation 6.9e-4/rot at nx2=32 (~2.7th order); in the dhj runs the rotating-frame Lz tracks sp to 1e-4 of the frame L over 100 rot**
- **[cs: radiation + srcterms are now REFUSED](cs-unsupported-physics-guards.md) — b283ad3a; they used to run and return a wrong answer. Two traps that make a guard which can never fire**
- **[cs NARROW-BLOCK resample degeneracy](cs-narrow-block-resample-degeneracy.md) — the along-seam stencil INVERTS below 3 cells, i.e. any MeshBlock under 6. Test at the MINIMUM legal block size**
- **[cs GENERAL/TABULATED EOS stale cache: FIXED](cs-general-eos-stale-cache.md) — p/Gamma_1/T cached BEFORE the gnomonic correction, 6-15% wrong. The gate is VACUOUS without a flow**
- **[cs MHD blast from a VECTOR POTENTIAL](cs-mhd-blast.md) — b283ad3a; div B zero by construction. Gate the construction with TWO numbers — div B alone is satisfied by a ZERO field**
- **[cube-vertex corner fill: PROMOTED, default ON](cs-wire-fill-wip.md) — branch cs-wire-wip, `<mesh>/cs_vertex_fill`; ownership is the exact chart DIAGONAL. Merged; widening coarse/fine MEASURED not worth it**
- **[cs cube-vertex REAL fill: prototyped, NOT built](cs-cube-vertex-real-exchange.md) — sampling beats extrapolation 3.4-5.5x and needs NO new exchange, just a widened destination range**

## Cubed sphere — history

- **[Cubed sphere ON GPU](gpu-this-capture-device-lambda.md) — af539941; hydro and MHD match CPU. The constant-memory CLOSURE LIMIT; incremental builds hide warnings**
- **[cs GPU multi-block fault: FIXED](cs-gpu-multiblock-fault.md) — 4302a008; one kernel captured both send and recv buffers. Second over-large-functor bug that session**
- **[Cubed sphere + MPI: FIXED](cubed-sphere-mpi-hang.md) — 1c2e29d6; unmatched receive, cube-vertex cause. Run the RECIPROCITY AUDIT first for any boundary bug**
- **[Cubed sphere: MHD seam FIXED](cubed-sphere-mhd-seam.md) — 4f19a244; three plumbing bugs, the transform was always exact**
- **[Cubed sphere: MHD](cubed-sphere-mhd.md) — e0c74357; dxedge was zero, mhd_fluxes had no gnomonic rotation, bcc was a non-orthogonal triple**
- **[Cubed sphere: COMMITTED as 9492a946](cubed-sphere-committed.md) — five sessions of hydro work is in git**
- **[Cubed sphere: along-seam resample, 2nd ORDER](cubed-sphere-seam-interp.md) — the "flat interior residual" is RETRACTED**
- **[Cubed sphere: seam basis transform, FIXED](cubed-sphere-seam-basis.md) — a tangent-BASIS transform, not a signed permutation**
- [Cubed sphere panel frames](cubed-sphere-panel-frames.md) — the six frames are forced by panel_neighbors; a duplicate copy had 3/4 swapped
- [Cubed-sphere hydro fix](cubed-sphere-hydro-fix.md) — halo axis bug FIXED; its numbers predate the panel 3/4 swap
- [Cubed sphere: x1 as radial](cubed-sphere-x1-radial.md) — DONE, bit-identical; axis roles cyclically permuted
- [Cubed sphere: hydro state](cubed-sphere-hydro-state.md) — 4 early bugs FIXED; its rigid-rotation section is RETRACTED
- [Cubed sphere for the hot Jupiter](cubed-sphere-for-hot-jupiter.md) — the dt argument is DEAD (radial binds)

## Hot Jupiter: physics, runs, campaigns

- **[cs_mhd_prod DIED at rot 41.6: NaN everywhere in one interval, no precursor, ran on NaN to rot 54](cs-mhd-prod-nan-rot41.md) — NaN-hunt restart 11529359 (0.05-rot dumps from rot 40). The code does not stop on NaN. START HERE for the old cs production**
- **[cs_hyd_rs: pure-hydro cs twin x {hllc, lhllc, ausmpup, ppmx} for the jet comparison](cs-hyd-rs-run.md) — launched 09-07 08:25, jobs 11529253/56/59 + chains. START HERE for the solver comparison**
- **[cs_mhd_prod2: cs MHD production FROM SCRATCH on the 09-07 defaults (WB+rot_potential+deep RT)](cs-mhd-prod2-run.md) — lead 11526708 (05:19, 6.5 h), chain 11526681/2; A/B against cs_mhd_prod. START HERE for the current production**

- **[Deep RT: radiative conduction + the tau_R BLEND fad2f5db](radiative-conduction-deep-interior.md) — `rad_tau_lo/hi`; sp column: dt+cost equal, RT chain SHORTER, 1-cell cut spread; the sp dt collapse (angular widths) FIXED; `rad_kappa_src = table` (kappa_R tabulated from the ck table at start-up) makes the handover CONSISTENT to 10 %; Freedman blend was off 1.7x. START HERE for deep RT**
- **[cs MHD dhj VALIDATED at rot 20: T structure = sp, jet WEAK (~1 km/s) on every grid, B does NOT slow it (beta>=100)](cs-dhj-diagnostics-rot20.md) — report URL + the EOS-inversion / cs lat-lon / wind-projection pipeline and its traps. START HERE for any cs dhj science plot**
- **[sp_mhd_prod (18f5dd21) DIES at rot 10.3: polar-row RADIAL field at the BOTTOM; BISECTED: 2a64e2c7 CLEAN, ideal 18f5dd21 DIRTY, resistivity innocent](sp-pole-bottom-radial-blowup.md) — round 2 arms wall_r (4990eb41) + x3shift_r (863e8337), jobs 11526902/3, gate ME1 at rot 4.7. START HERE for sp MHD**
- **[cs MHD dhj BOTTOM-BOUNDARY DRIFT: FOUND + FIXED 248f1b77](cs-mhd-bottom-inflow.md) — the radial BC built cs ghost bcc as plain face averages, the ghost energy drifted every call (8 % by rot 16), bottom cells drained, NaN at rot 32. EVERY cs MHD dhj run before it is suspect; redo from scratch. START HERE**
- **[sp MHD ENERGY EXCESS vs hydro (+1 %/2 rot, 3-4x horizontal KE), ABSENT on cs](sp-mhd-energy-excess.md) — floors on the upper NIGHTSIDE where sp carries 7x the field of cs; deep field identical; every sp MHD variant has it; cause OPEN. READ before comparing sp MHD with anything**
- **[Correlated-k design + build](correlated-k-design.md) — THE MAIN THREAD. COMPLETE: physics, docs, tests, MPI. START HERE**
- **[dhj RT on the CUBED SPHERE: WORKING](dhj-cubed-sphere-port.md) — the full production physics tracks spherical polar; four bugs found, three PRE-EXISTING**
- **[cs NON-ORTHOGONAL audit: BOTH FIXED](cs-nonorthogonal-audit.md) — e415b91a; history.cpp KE/ME were ~1% wrong on cs, newdt understated the fast speed by 1/sin_cell. Pre-fix cs history files carry that 1%**
- **[A real cs MHD bug, FIXED: the C2P floor corrupts u.e](cs-mhd-c2p-floor-corrupts-ue.md) — ConsToPrim's magnetic energy uses the NON-ORTHOGONAL triple; a floor then OVERWRITES u.e. Vertex order -0.04 -> +1.97. But it does NOT save the dhj run (0.2023 -> 0.2129) — that blow-up stays OPEN**
- **[cs MHD low beta: per-step operator 2nd ORDER everywhere incl. the vertex](cs-mhd-low-beta-divergent.md) — the instability is NONLINEAR on a consistent scheme; HLLE fallback is a legitimate stabiliser. START HERE**
- **[cs MHD low-beta: fixed ON THE TEST PROBLEMS ONLY](cs-mhd-lowbeta-fix.md) — cb1afd28, `<mhd>/cs_lowbeta_llf`, a per-face HLLD->HLLE fallback below beta 0.5, default ON for cs. Real but INSUFFICIENT: the production dhj run still dies at 0.21 rot**
- **[cs SEAM GHOST e_int: RETRACTED as a mechanism](cs-seam-ghost-eint.md) — the seam degrades EVERY quantity by the same ~15-19x; eint is not special, the 1/beta signature was in the null row too. Real E-vs-B inconsistency is ~14 %. Fix b9efbbfa single-rank only. READ THE RETRACTION FIRST**
- **[GENERAL-EOS AUDIT](general-eos-audit.md) — the (rho,e)->T inversion is BRACKETED, monotone, warm-start-independent; no new live bug. Off-table states are silently log-linear-extrapolated**
- **[3D RESISTIVITY BUG since 2026-06-20: curl B missing its theta terms on sp + Cartesian, FIXED (uncommitted)](resistivity-3d-curl-missing-terms.md) — every sp resistive dhj run since June is wrong; cs unaffected. START HERE for any sp resistive claim**
- **[sp POLE-EDGE CURRENT: far-side loop segment had the WRONG SIGN, J_r grew ~1/dtheta, FIXED](sp-pole-edge-current-sign.md) — the m=2 polar pressure pattern; new edge-by-edge EMF check in sp_test iprob=11**
- **[sp POLE: resistive J_r divided by a ZERO dual area, FIXED](sp-pole-edge-area-zero.md) — NaN in cycle 1 on any uniform-theta grid; production hid it with the theta stretch**
- **[sp POLE FIXED: third-difference polar_emf_diss + x3 face-state theta shift](sp-pole-fixes.md) — de667f32/863e8337; polar-row L1(B) rate 0.76 -> 2.1, v_r 2.3; tangential stays ~1.5 (cart REPAIRED but no better, 4th-order faces REJECTED). dhj checkerboard gate PASSED to rot 1.2. START HERE for the pole**
- **[sp_test: RIGID ROTATION + RESISTIVE DECAY tests + rates](sp-test-rigidrot-resist.md) — hydro 2.6, resistive B 2.1; the POLAR ROWS do not converge for a tangential field (O(1))**
- **[RESISTIVITY AUDIT](resistivity-audit.md) — eta_b ghosts, resistive dt, and the pole-face curl are CLEAN (with file:line); the one live gap is the polar EMF average excluding the resistive EMF, fixed behind `use_polar_average_eresist`**
- **[sp POLAR FIELD BLOW-UP: FOUND + FIXED 1cabe85c `<mhd>/polar_emf_diss`, DEFAULT ON 3147de4f](sp-polar-field-blowup.md) — GS05 corner-EMF upwind terms grow a (-1)^(i+k) face checkerboard in the polar row. READ BEFORE ANY cs-vs-sp MHD CLAIM**
- **[cs ROTATION SOURCE BUG: cs ran the CARTESIAN beta-plane, sp ran full Coriolis](cs-rotation-source-bug.md) — f75ad783; the switch was REMOVED in ae831665, full rotation is now unconditional. "Same physics, cs dies and sp does not" was FALSE. READ THIS BEFORE TRUSTING ANY cs-vs-sp COMPARISON**
- **[cs 1-ULP amplification: nx=64 was DYNAMICALLY UNSTABLE on PRE-FIX binaries; RETEST cs_n64_fixed 11517116 ALIVE past 0.5 rot, dt flat 2.3 s -- old claim STALE](cs-ulp-amplification.md) — round-off grows to O(1) and KILLS nx=64 at 0.309 rot while it saturates ~1% at nx=32; the perturbation script must skip the TREE metadata. START HERE**
- **[WHAT RECONSTRUCTION ACTUALLY RUNS: sp ignores `reconstruct` in ALL directions (stretched PLM hard-wired); cs honours it but x1 ignores the radial stretch](reconstruction-what-actually-runs.md) — READ before any reconstruction claim**
- **[cs VERTEX = LIMITER CLIPPING: ppmx (the existing extremum-preserving PPM) and WENO-Z cut the worst-cell vertex force 6.5x; the Cartesian update (d166238b) does nothing there](cs-vertex-limiter-clipping.md) — same root cause as the sp pole; ppmx is the production candidate (nghost 4), gates in memory**
- **[cs VERTEX ORACLE: exact ghosts change the vertex residual <1 %; it is the CELL BALANCE, no ghost stencil can help](cs-vertex-oracle-halo-innocent.md) — the through-vertex seam-continuation idea DECIDED against; Cartesian-momentum update is the only lever left**
- **[sp POLAR-ROW RESIDUAL FIXED 7.7x: `mesh/polar_quadratic_recon` 6410be99, default OFF](sp-polar-row-reconstruction.md) — limiters see the pole as an extremum; `reconstruct` is IGNORED on sp. READ before any sp accuracy claim**
- **[sp GEOMETRIC-SOURCE RESIDUAL: polar rows 14x interior; face-sum WB source (default OFF) does NOT touch it; Cartesian polar-row update 34882bd7 (default OFF) fixes the transverse-field part only](sp-geometric-source-residual.md) — my "WB will fix it" RETRACTED**
- **[cs WELL-BALANCED SOURCE CACHED: 46 % -> 4 % of a run, bitwise identical, restart-clean](cs-wb-source-cached.md) — 88064f67; records the tlim-clipping restart-test TRAP and that cs_test cannot restart**
- **[GS07 corner EMF is UNCONDITIONAL on cs since b71ef392 (switch REMOVED)](cs-gs07-emf-gate.md) — gate: same order as the plain average, 1.8x its smooth-field L1(B); `bs_emf` is the only way back, for diagnosis**
- **[cs DEEP TOROIDAL SHEET: CLOSED, it was the stretched resistive bug c5c85e3b](cs-deep-toroidal-sheet.md) — rcmfix on 979edada has NO sheet; the old cs_prod_mhd_rot (sheet, eruption, death) is archived. START HERE**
- **[cs dhj production RETRY: ARM MATRIX + production reruns](cs-dhj-production-retry.md) — nx=32 clean to 5.5 rot, nx=64 dies at 0.3 rot (1-ULP seed amplifies). Pre-fix binaries. START HERE**
- **[cs MHD instability CHARACTERIZED](cs-mhd-instability-characterized.md) — a6401406, df2edc8e, c3992145. NOT gnomonic-specific (spherical polar does it too). The WELL-BALANCED source is BUILT and does NOT fix it, refuting the quadrature-mismatch hypothesis. TWO RETRACTIONS: the vertex order deficit and "the rate rises with resolution". START HERE**
- **[cs MHD MINIMAL REPRODUCER: it blows up from NOISE](cs-mhd-minimal-reproducer.md) — f3d35a96, cs_test iprob=13. An atmosphere AT REST + a force-free field dies in 2.5 min on ONE CORE, no RT/rotation/EOS/resistivity/shear. Growth rate RISES with resolution. Halo and stratification EXONERATED by A/B. START HERE**
- **[cs MHD dhj blow-up: OPEN, at a CUBE VERTEX](cs-mhd-dhj-blowup.md) — 10 mechanisms ELIMINATED incl. timestep, CT/monopoles, seam transform+resample (exact), block decomposition, general EOS, radial stretch. Refining makes it WORSE. Minimal reproducer: ideal gas + ideal MHD + uniform grid, cs dies 0.15 rot, sp clean. START HERE**
- **[cs dhj NaN: CAUSE FOUND](cs-raisevel-missing-floor.md) — GnomonicEquiangleRaiseVel never re-applied the floors after the metric correction. START HERE for cs+dhj**
- **[cs_dhj_long went NaN in ONE rotation](cs-dhj-long-run.md) — RESOLVED; kept for the 2x2 that pointed at the grid. A job that COMPLETES with exit 0 is not a job that ran**
- **[ck_grav_prod: THE PRODUCTION CAMPAIGN](ck-grav-prod-run.md) — point-mass gravity, nx1 234 stretched, 2 GPUs, ~70-110 h. HELD, NOT LAUNCHED. START HERE**
- **[ck_mhd_b3: MHD + EOS resistivity, bbot 3 G](ck-mhd-b3-run.md) — ck_grav_prod's input, 6 lines changed. SMOKE PASSED, NOT LAUNCHED**
- [ck_limb: the CONSTANT-g run](ck-limb-run.md) — superseded as science by ck_grav_prod; still the calibration reference (1.76 H buffer)
- [ck_grav_size: the sizing run](ck-grav-size-run.md) — ANSWERED: r99 = 1.794e10, H = 22 cells
- [ck_hydro_long: RETIRED](ck-hydro-long-run.md) — sim "days" are Earth days; its mass drift was a transient
- [Ideal-gas + EOS x_e resistive runs](xe-resistivity-long-runs.md) — COMPLETED but carry the stellar-heating bug; do not restart on a fixed binary
- **[dhj blow-up: CLOSED, it was the race](dhj-ck-eos-blowup.md) — no thermodynamic blow-up exists; every pre-fix arm comparison is retracted**
- [Exo-FMS cross-validation](exofms-cross-validation.md) — LW agrees 1.7%; found and fixed a ~25% stellar-heating loss (b4e0953c changes production answers)

## Hot Jupiter: grid, EOS, atmosphere

- **[HYDROSTATIC well-balanced scheme SUPPORTS cs + point-mass gravity; NEW wb_option=polytropic 140dbf9e](wb-hydrostatic-scheme-cs.md) — cs == sp to 4 digits; 40-100x residual gain, floor = the IC's own discrete inconsistency; avg-anchor correction REFUTED; table cost 1.03x of plain (polytropic + wb_cache_every=10); VERIFIED hydro/MHD/resistive on cs, sp, Cartesian x1/x3, RT on, rotation on (WB = the CONVERGED answer), GPU == CPU to 1e-11. Recipe: polytropic + wb_x1 + wb_cache_every=10 + problem/rot_potential=true (barotropic rotating IC). Not yet in production**
- **[Radial grid stretch](radial-grid-stretch.md) — 1e19d4e7, `mesh/use_grid_stretch_r_poly`. Post-processing MUST apply the map to x1v**
- **[Radial stretch REFIT](radial-stretch-refit.md) — the production coefficients leave a 3.7x spread in cells/H; refit vs H_rho gives 1.34x and 1.46x dt on cs**
- **[dt binding direction](dt-binding-direction.md) — dt is a straight MIN over x1/x2/x3; the radial stretch flips binding from azimuthal to radial, costing 1.76x**
- **[Point-mass gravity flag](grav-point-mass-flag.md) — `problem/grav_point_mass` (c37ebe75); constant g understates H by 2.7x at 1e-6 bar**
- **[Stellar tide flag](stellar-tide-flag.md) — c4aa730f, `problem/stellar_tide`; validated, deliberately OFF (mirror-symmetric between limbs)**
- [Stellar tide at the domain top](stellar-tide-at-domain-top.md) — the missing terms are 0.5% at the limb, +6.8 cells substellar; r_L1 = 4.28 R_p
- **[Isobar, not shell](dhj-isobar-vs-shell.md) — the 1e-6 bar level must be measured on the ISOBAR. Supersedes the floor conclusions**
- [Sizing x1max from the isobar](isobar-buffer-calibration.md) — buffer in SCALE HEIGHTS against the accepted run (1.76 H), plus two mu validations
- [Floors for the 1e-6 bar level](dhj-floors-for-1e-6-bar.md) — lower pfloor to 1e-3 (free); do NOT extend x1max or raise dfloor; tfloor_kelvin=50 crashes
- [Grid resolution design](dhj-grid-resolution-design.md) — measure cells/H AT the feature; the cost is the timestep, not the cells
- **[Upper-atmosphere mottling](upper-atm-mottling.md) — the speckled T maps are the H2 dissociation front plus sinking plumes, NOT noise or floors**
- **[H2 chemistry quenches above 1e-3.6 bar](h2-chemistry-quench.md) — the EOS assumes instant dissociation; recombination is 8 decades too slow, so that front is not physical**
- **[sp HYDRO vs MHD, controlled at 100 rot](sp-hydro-vs-mhd-comparison.md) — the field takes 31% of the zonal KE, 89% of it below 20 bar, and shows at the photosphere as ~2% in radius. START HERE for "what does B change"**
- **[COMPOSITION + FIELD maps in the atlas](dhj-composition-maps.md) — dayside atomic + iron ionized, nightside molecular + neutral; the atlas now carries BOTH runs, with log|B|, B_r and plasma beta (min 0.96 at tau=2/3)**
- **[tau=2/3 PHOTOSPHERE diagnostic](dhj-photosphere-diagnostic.md) — a4d30d43, `problem/photosphere_dump`; the emitting surface spans 7040 km between near- and far-IR. Records the x2v CENTROID trap**
- **[EOS table dump](eos-table-dump.md) — f882159f, `<hydro>/eos_table_dump`; recovers p and T from binary dumps. Table is NaN below 71 K. REGENERATED on orion 2026-09-07 at /orion/ptmp/jinma/Athenak/eos_table/, validated to 7e-4; the column dump j,k are GHOST-INCLUSIVE**
- **[EOS inversion NaN trap](eos-inversion-nan-trap.md) — a naive root find on the dumped table is SILENTLY 0.35 dex low in T. Validate against ck_dump_file**
- [General EOS optimization](general-eos-optimization.md) — 1.41x on CPU but only 1.062x on GPU; further EOS work capped at ~8%
- [FastChem vs the general EOS](fastchem-vs-general-eos.md) — MEASURED: buys nothing, but validates the composition. DECIDED against
- [UHJ band structure in the literature](uhj-band-structure-literature.md) — Parmentier+2018 and Tan+2024 both use 11 Kataria bins x 8 g; our grid choice is right

## Solar / stellar convection

- **[run/sun viz + WAVES NOT SHOCKS](sun-convection-viz.md) — the "Buoyancy Box" artifact; tau=2/3 granulation from the run's OWN opacity; four measures say gravity waves (max compression 0.048, 87% solenoidal). The 512^3 convection data is GONE from viper (it is on ORION). START HERE**

## Bugs and performance

- **[All-Mach solvers lhllc/ausmpup: algebra matches the papers; rk2 == rk3, NO CFL~M restriction on the stratified rest test](lowmach-solver-audit.md) — two design choices (normal-velocity chi; AUSM Mref=1 hard-wired). wb_column tgrad=-0.3 is SUPERADIABATIC (trap)**

- **[REFLECTING WALLS leaked mass under a blast on sp (always) + cs (since 979edada): FIXED, mirror the wall state](reflect-wall-mass-leak.md) — invisible in every v_r = 0 test; closed blasts back to round-off on both grids**
- **[GPU nondeterminism, FIXED 6e600f12](dhj-run-to-run-nondeterminism.md) — a missing team_barrier before the x1 interface-pressure floor; new _gpu regression test**
- [Polar MPI host-mirror bug, FIXED](polar-mpi-host-mirror-bug.md) — 3882e37f; every multi-rank polar run died on cycle 1
- [Restart bug with tfloor_kelvin, FIXED](athenak-restart-tfloor-bug.md) — ba2f0943; general-EOS runs could never be restarted
- [dhj ideal input out-of-bounds, FIXED](dhj-ideal-input-oob.md) — 0d1f6f6a; use a Debug build (Kokkos bounds checking) for this class
- [MHD flux stale x1 limits](mhd-fluxes-stale-x1-limits.md) — 3429f59f; LATENT, changes no answer. Also records the full race audit
- [RT chain-parallel split kernel](rt-chain-parallel-split.md) — 38311a8a, the live architecture; its cost numbers predate the coalescing fix
- [RT source is semi-implicit](rt-source-semi-implicit.md) — FIXED 048dff30: exponential relaxation at rate 4E/e; cured the red giant dt collapse; the wall control is what ruled out the boundary
- **[red giant dt COLLAPSE is RADIATIVE CONDUCTION, still OPEN](red-giant-conduction-dt-collapse.md) — RESUME HERE: the RT fix moved it 2.9e4 -> 1.86e5 but did not remove it; MEASURED: the limiting cell is IN THE BLEND WINDOW at the photosphere (tau 1-3), both times; tau blend 10/100 NaNs; STS for conduction is the promising cure**
- **[MAGNETIC ENERGY is not in the dumps](mhd-energy-not-from-dumps.md) — bcc is the mean of the FIELDS, the code integrates the mean of the ENERGIES; 8.8x off in B_r. Quote the history**
- [Overall GPU profile of a dhj run](dhj-overall-gpu-profile.md) — RT 35.7%, fluxes 27.8%, ConsToPrim 12.2%, boundaries 9.6%
- [Meshblock decomposition on one GPU](meshblock-decomposition-gpu.md) — 32 blocks beats 2 by 1.35-1.45x; the second APU buys ~1%
- [nx1 ceiling from MHD LDS](nx1-ceiling-lds.md) — radial resolution capped at nx1=264; x1 cannot be split because RT is a column solve
- [par_for_outer team size](par-for-outer-team-size.md) — ~3% of wall. DECIDED not worth it, do not re-propose
- RT kernel history (all HISTORICAL, superseded by [[rt-chain-parallel-split]]): [interleaving](inflight-rt-kernel-optimization.md), [band scaling](inflight-rt-band-scaling.md), [occupancy](rt-kernel-occupancy-limit.md)
- [IN FLIGHT: nx2=64 + GPU timings](inflight-nx64-timings.md) — paused 2026-08-16; sections 5/6 still to re-measure on GPU

## Working practice

- **[VALIDATE THE INSTRUMENT](validate-the-instrument.md) — a gate reporting "nothing" usually means "I did not look". READ THIS BEFORE TRUSTING ANY CLEAN RESULT**
- **[Measure impact before claiming it](measure-impact-before-claiming.md) — report the defect, not its consequences, until an A/B has shown them; two retractions**
- [Localise by dilution](localise-by-dilution.md) — a FLAT error profile does not exonerate a boundary; vary the DOMAIN SIZE at fixed dx and look for 1/L
- [Don't propose confounded tests](confounded-tests-rejected.md) — the user vetoes experiments whose result could not discriminate
- [Use fork, not origin](use-fork-not-origin.md) — all git fetch/push goes to jing-ze-ma/athenak
- **[Restart output dt trap](restart-output-dt-trap.md) — last_time comes from the restart (a dump silently never fires), and output dt is SIMULATED time (38 GB in two minutes)**
- [Never write in run/](never-write-in-run-dir.md) — 165 GB of output; read-only, and never `git add -A`
- [Viper HIP build recipe](viper-hip-build-recipe.md) — module loads and cmake flags for the MI300A APU nodes
- [Reusable CPU baseline binary](bench-baseline-worktree.md) — a worktree at bench/base_wt reproduces every baseline digit-for-digit; the other bench binaries are GPU builds

---

**CURRENT STATE, 2026-09-07 (~09:45), HANDOVER TO ORION.** HEAD 7c652768 on `polar-average-perf`, PUSHED
to the fork (jing-ze-ma/athenak). Viper is in MAINTENANCE 09-07 12:00 -> 09-12 12:00; the next
session runs on ORION, which cannot see viper's bench/ or scratch. Everything needed is in git:
docs/handover/HANDOVER-2026-09-07.md (READ FIRST), docs/handover/scripts/ (the analysis scripts),
docs/handover/claude-memory-2026-09-07/ (a copy of this memory directory as of the handover).

### Running on VIPER through/after the maintenance (check when back on viper, NOT from orion)
1. cs_mhd_prod2 (11526708 -> chain 11526681/2): cs MHD production on the new defaults, rot ~11 at 09:00.
2. cs_hyd_rs/{hllc,lhllc,ausmpup,ppmx} (11529253/56/59/11529330 + chains): pure-hydro solver and
   reconstruction comparison; the jet at rot 15-20 is the question. LHLLC keeps 70x the deep radial KE.
3. sp_pole_bisect wall_r (11526902) and x3shift_r (11526903): gate ME1 at rot 4.7 / polerow at rot 4.
4. cs_mhd_prod_nan (11529359): passed rot 41.6 without NaN -> the original NaN is not reproducible.
5. Login-node serial: lowmach/isoL_N64_* (2000 t_BV Edelmann slow-mode test), ~50 % at 09:00.

### NEXT STEPS on ORION (in order)
1. Clone/pull the fork branch; `git submodule update --init`; build CPU (PROBLEM=deep_hot_jupiter_rt) and GPU
   per viper-hip-build-recipe's orion equivalent; the ck tables are at /orion/u/jinma/ATHENAK/athenak/data/exo_fms_ck.
2. Regression: `cd tst && python run_test_suite.py --test test_suite/rad/test_rad_dhj_ck_cpu.py` and
   test_rad_cs_raddiff_cpu.py (both PASS on viper at ebd57244).
3. GPU binary needs the NaN guard (d3d74f2b) -- rebuild before any GPU run.
4. Determinism test (cs MHD restart twice from one rst, 200 cycles, bitwise) -- the rot-41.6 NaN hunt showed
   restarts diverge at 5e-6 in 0.1 rot ([[cs-mhd-prod-nan-rot41]]).
5. Open science: the jet comparison (results on viper after 09-12), the deep interior with the radiative
   outlet, and the sp polar blow-up bisection result.
