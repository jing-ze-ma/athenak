---
name: inflight-2026-09-09-viper
description: IN FLIGHT, STOPPED 2026-09-11 ~03:30 - sp_mhd_prod3 clean at rot 23.8 (chain on apu), cs collapse CLOSED, instrument + rad_tmax_kappa + handover bundle COMMITTED and PUSHED (c26b01ac), explicit ensemble arms CANCELLED; next = dt_min 1e-3, density floor vs WB, longer si ensemble. READ THE BOTTOM SECTION FIRST
metadata:
  type: project
---

State when the session ended (2026-09-09 ~21:45, viper):

1. **DONE**: `time/dt_min` committed (ca616d1a); orion's 33 commits merged (387c5794, clean);
   HIP build fix cb8c6c2f (orion's `dt_diag.template sync<HostMemSpace>()` in conduction.cpp
   failed Kokkos' static assert on HIP -> `modify_device()/sync_host()`); build_dhj_cpu and
   build_dhj_gpu rebuilt on the merged tree (GPU binary 21:29); regression tests
   test_rad_dhj_ck_cpu and test_rad_cs_raddiff_cpu both PASS (note: `--test` is not
   cumulative, run them one at a time; use `python3`). Push to fork was delegated; verify
   with `git status -sb`.
2. **MERGE WAS NOT NEUTRAL** [[rt-semi-implicit-changes-dhj-answer]]: flag `problem/rt_semi_implicit`
   committed ae300379 (default true; false = pre-merge physics, byte-identical to ca616d1a),
   GPU binary rebuilt 22:31, pushed. The ablation was RE-STAGED with the flag false:
   - from scratch, bench/cs_ablate/{ctl2,nocond,nowb,norot,cache1}: 11539561/62/63/64/65 (5 h,
     apu1). Old ctl 11539431 keeps running as the SEMI-IMPLICIT control (flag true); *_si dirs
     are cancelled partial runs.
   - RESULT: all 10 valid restart arms (f_/fh_) SURVIVED past the original death times; see the
     update in [[cs-vertex-dt-collapse-0907-defaults]]. Read ONLY the from-scratch arms next.
   - RESTART BUG [[restart-rot-potential-bug]] fixed 42c0a6a2 (GPU binary 23:03); EVERY earlier
     restart arm is INVALID (they all 'survived' on a doubled centrifugal force). Valid restart
     arms on apudev (14 min each), bench/cs_ablate_r: f_{ctl,nocond,nowb,norot,cache1}
     11539616-20 (old ids 11539567-71 cancelled) restart lhllc rot 4.0 (dhj.00002.rst, original collapse at t 1.345e6, tlim 1.45e6);
     fh_* 11539621-25 restart hllc rot 15.11 (with output*/last_time=4.608e6 to defeat the
     restart-output trap) (dhj.00008.rst, collapse at 4.6330e6, tlim 4.70e6,
     dumps every 5e3 s to catch the onset). OLD-binary controls old_cache1 11539512 and
     h_old_cache1 11539513 (f6f0d6f4, wb_cache_every=1). First-round restart arms (ctl, nocond,
     nowb, norot, cache1, h_*) ran with the flag TRUE: ctl and nocond did NOT collapse -> not
     interpretable. cs_ablate_r/NOTES.md has the map. rst_vmax2.py (hydro-capable) is in the
     session scratchpad precursor/ -- copy it to bench/ if still there.
   - Precursor analysis: NO precursor in the last 2-rot dump; the vertex column goes from
     unremarkable (62-96th pct) to 1-2e7 cm/s within <0.02 rot. All four cs_hyd_rs arms AND
     cs_mhd_prod2 collapsed (hllc at rot 15.19, 0.08 rot after its last dump).
   Original plan for reference: bench/cs_ablate/{ctl,nocond,nowb,norot,cache1}, jobs
   11539431/32/33/34/35 (5 h walltime, 2 GPUs each, apu1). Base = cs_hyd_rs/lhllc input
   (dies at 4.4 rot), tlim 6 rot, dt_min 1e-2, rst+bin every rot. Gate and post-mortem recipe
   in bench/cs_ablate/NOTES.md. All five passed a 5-cycle CPU smoke on the merged binary.
   READ THIS FIRST next session: which arm survives 6 rot with dt flat in dhj.hydro.hst?
   nowb keeps rot_potential=true (the explicit radial centrifugal term comes back).
3. **sp bisection round 3** [[sp-pole-bisect-round2-null]]: 11539009/10/11 running (rot 0.9
   at 21:20), 11539012 (eb8e3cb1) pending; bench/sp_pole_bisect/r3_*; gate 1-ME at rot 4.7.
4. Orion's basis note (docs/handover/NOTE-2026-09-09-cs-covariant-basis.md) says the dhj
   runs are NOT affected; latent sites for viper (ismcooling drag, solid_body_rot,
   derived_variables) + shared Gnomonic kinetic/lower helpers + two dhj sponge physics
   flags (thermalised KE, explicit damping). Not started; avoid red_giant.cpp (orion has
   uncommitted edits).
5. [[wb-restart-cache-bug]] fix 5c0b98e4 is now in the viper binaries; any chained WB run
   restarted before this build got the zero-background kick.

### STATE AT STOP, 2026-09-10 ~02:00 (read this first)
- Decision (user): discriminate the 09-07 defaults with an ENSEMBLE (`problem/seed`,
  `seed_amp=1e-10`, 10 seeds x 5 arms x 8 rot, ~1.6 h each on 2 GPUs, bench/cs_ens/ planned).
  NOT LAUNCHED. The user required the perturbation to be a function of the GLOBAL cell only
  (include the PANEL index in the hash; identical for any meshblock/rank layout, incl. seam and
  vertex cells) and to verify that CELL-WISE (bin dumps at t=0 for 16x16 vs 8x8 meshblocks,
  numpy array_equal; no duplicated noise across panels) BEFORE running.
- Seed code COMMITTED and PUSHED: 3f8d2415 (hash of panel + global i/j/k, modify_host/
  sync_device, prints seed/amplitude). build_dhj_gpu was rebuilt and `strings` shows the
  rt_semi_implicit read IS compiled now (count 1). NOT DONE (subagent died on an API 500):
  (a) the cell-wise decomposition check (16x16 vs 8x8 meshblocks, t=0 bin dumps,
  array_equal; no duplicated noise across panels) -- the user REQUIRED this before running;
  (b) the startup print of rt_semi_implicit; (c) the seed's effect on the very first cycle
  (seed 1 vs 2 differ, seed 0 == no seed) on CPU. Do (a)-(c) first next time, via Opus.
- Every GPU arm of 09-09/10 ran semi-implicit ON. Before the ensemble, rerun ONE arm (ctl2
  input, seed 0) 0.2 rot on apudev on the fixed binary and confirm its first hst row equals
  det_cyc_old's (1-KE 1.80756e+28 at cycle 1 => explicit source active).
- Results so far: all 6-rot arms survived (ctl, ctl2, nocond, nowb, norot, cache1); all 20
  restart arms survived; GPU deterministic; collapse = chaos-selected rare event.
- Old ablation dirs: bench/cs_ablate (6-rot arms, det_*, det_cyc_*), bench/cs_ablate_r
  (restart arms, all INVALID or non-discriminating). Worktrees in the session scratchpad
  (ab/wt_pre, wt_ext, wt_si) are gone with the scratchpad; det_cyc_* keep their binaries.
- sp bisection round 3 (11539009-12): gate 1-ME at rot 4.7, NOT yet read.


## UPDATE 2026-09-10 ~03:30
- Seed gates (a)-(c) PASSED (bench/cs_seedchk): 16x16 vs 8x8 MeshBlocks bitwise equal at
  amp 1e-3 (the bin dumps are FLOAT32, so a 1e-10 kick is invisible in a t=0 dump -- gate
  at 1e-3, run at 1e-10); seed 1 vs 2 differ in 99.996% of cells; seed 0 amp 1e-10 differs
  from amp 0; panel-panel correlation of the noise < 4e-3. The seed line prints ONLY when
  seed_amp > 0.
- 95f2b07c: startup print "RT two-stream source is SEMI-IMPLICIT/EXPLICIT". GPU arm
  ctl2_si_chk (11543953, apudev) reproduced det_cyc_old bit for bit (1-KE 1.80756e+28 at
  cycle 1): the flag IS honoured now; the build race is closed.
- Ensemble: bench/cs_ens/<arm>/s01..s10, arms ctl2/nocond/nowb/norot/cache1, all with
  rt_semi_implicit=false, tlim 2.44e6 (8 rot), bin dumps 1/rot, self-chaining submit.sh,
  jobs in cs_ens/jobs.txt. Read: which seeds/arms hit dt_min (log) vs reach tlim; a rate
  per arm.
- sp: see [[sp-pole-bisect-culprit-863e8337]]; worktree bench/wt_x3shift, branch
  polar-x3-shift-flag (NOT merged), arms bench/sp_pole_bisect/hd_{shift,rotate,scalar}.
  Gate 1-ME col 11 at rot 4.7. If rotate is clean on HEAD: merge the branch (default
  rotate) and restart sp production.


## UPDATE 2026-09-10 ~09:05 -- STOPPED at the user's request. READ THIS FIRST.
### Ensemble (bench/cs_ens, explicit RT unless noted; collapse = "FATAL: timestep collapsed" in log.out)
| arm | dead | survived 8 rot | note |
|---|---|---|---|
| ctl2 (explicit) | 10/10, rot 4.3-5.7 | 0 | a narrow window: NOT a rare chaotic event |
| nocond | 3/5 (4.3, 5.5, 5.6) | 2 (s04, s05) | conduction is a secondary amplifier |
| si (rt_semi_implicit=true, 6th arm added) | 1/8 (s05 at 5.3) | s01, s02 done; s03-s08 reached tlim 8 rot (FINAL: 7/8 survived) | THE EXPLICIT TWO-STREAM SOURCE IS THE DOMINANT TRIGGER; s02 ended with dt 3.1 s (precursor) |
| nowb, norot, cache1 s01-s05; si s09, s10 | not started, HELD (JobHeldUser) | | s06-s10 of nowb/norot/cache1/nocond CANCELLED to free apu1 |
Read next: `for d in cs_ens/*/s*; do ...grep 'timestep collapsed' log.out.*; tail -1 dhj.hydro.hst` (rot = t/3.05e5).
Holds: a watcher agent releases them once spb_hd_rotate_apu1 (11552424) runs; if any ens_* still shows
JobHeldUser next session, `scontrol release` them by hand. Account cap: MaxJobs=8 RUNNING on apu.
Next physics question: WHY the explicit source collapses dt at a cube vertex (stiff heating at the
vertex? check the collapse cell's T/tau history in the last dumps of ctl2/s01, bin every 1 rot) and
whether semi-implicit should become the production default (it changes the answer 3% KE / 20 cycles,
see [[rt-semi-implicit-changes-dhj-answer]]).

### sp A/B of <mesh>/polar_x3_shift (branch polar-x3-shift-flag c5740ec4 in bench/wt_x3shift, NOT merged)
Gate: dhj.mhd.hst column 11 (1-ME) at rot 4.7: clean ~8e31, dirty ~1.6e34. All arms clean at rot <2 (expected).
- apudev chains (15-min links, --dependency=afterany, ~0.19 rot/link, STOP file stops a chain):
  bench/sp_pole_bisect/hd_shift (rot 1.74), hd_rotate (rot 1.60); hd_scalar's chain was cancelled.
- apu1 single 6-h jobs: hd_shift_apu1 11552423 (from rst rot 0.5, at 0.8), hd_scalar 11552169 (at 0.6),
  hd_rotate_apu1 11552424 PENDING (Priority). Whichever route reaches 4.7 first is read.
- Traps: `mesh/polar_x3_shift=` on the command line needs the line present in the athinput (it is);
  sbatch from inside a job is DENIED on viper (hence pre-submitted chains); scontrol top/Nice denied.
If rotate clean and shift dirty on HEAD: merge polar-x3-shift-flag into polar-average-perf (default
rotate), push to the fork, restart sp MHD production. scalar tells whether the vector or scalar shift
is the destabiliser (mechanism only).

## UPDATE 2026-09-10 ~13:40 (resumed)
- sp A/B FINISHED: shift DIRTY 6.4e34, rotate CLEAN 7.2e31, scalar CLEAN 5.4e31 at rot 4.7 -> merging the flag branch (default rotate), relaunching sp production on apu. See [[sp-pole-bisect-culprit-863e8337]].
- Ensemble si arm COMPLETE: 7/8 reached 8 rot (s05 dead at 5.3); ctl2 10/10 dead 4.3-5.7; nocond 3/5 dead. The 17 held arms (nowb/norot/cache1 s01-05, si s09-10) are being released by hand (watcher never fired).
- apudev is now ONE-TIME TESTS ONLY (user rule): chains and production go to apu. [[apudev-one-time-tests-only]]
- Next: WHY the explicit source collapses dt at the vertex (collapse-cell T/tau/heating history in ctl2/s01 last dumps vs si/s01).

## STATE AT STOP, 2026-09-10 ~19:30 (read this first next session)
- sp: polar-x3-shift-flag MERGED (9a9396f7, default rotate), pushed to fork; GPU binary build_dhj_gpu at 9a9396f7 (strings-verified).
  sp_mhd_prod3 (bench/sp_mhd_prod3, apu, 4-link chain 11569324-27, STOP-file guard) PENDING on Priority at stop; the first link is also
  the smoke test of the 09-07 defaults on sp. CHECK FIRST: is it running, did link 1 start cleanly (log, "RT two-stream source is
  SEMI-IMPLICIT" line, first hst rows), gate 1-ME col 11 at rot 4.7 should stay ~1e32.
- cs ensemble: ens_si_s09/s10 RUNNING (finish ~5 h); nowb/norot/cache1 s01-05 (15 jobs) HELD BY USER ORDER so production goes first.
  Release them only after sp_mhd_prod3 is running (scontrol release; the user wanted production first). No watcher exists.
- cs collapse: MECHANISM + ONSET recorded in [[cs-vertex-dt-collapse-0907-defaults]] (bottom). The NEXT INSTRUMENT is not built:
  (1) add a per-cycle diagnostic dump of the RT source term (de/dt), kappa_rad and the conduction dt owner for one meshblock (m=6)
  in the dhj pgen or two_stream_rt.hpp; (2) run ctl2/s01 from rst rot 5 to write an rst at t=1.7385e6 (rot 5.700, output4/dt override),
  ~10 min on apudev (one-time test OK); (3) restart from it twice, explicit and rt_semi_implicit=true, ~400 cycles, per-cycle dumps;
  compare the vertex column i=40-47 on m=6. Question: does the explicit source or conduction inject the energy in the drained cells.
- Working rule: apudev ONE-TIME tests only (MaxTime 15 min); chains/production on apu. [[apudev-one-time-tests-only]]

## UPDATE 2026-09-10 late evening (resumed ~22:00)
- sp_mhd_prod3 (11569324, apu) RUNNING clean: rot 1.19 at 20:16, dt ~2.07 s, 1-ME flat ~7e31, startup echoes
  `polar_x3_shift = rotate` and `RT ... SEMI-IMPLICIT`, BUILD_COMMIT 9a9396f7. Gate 1-ME col 11 at rot 4.7 not yet reached.
- The 15 held ensemble arms (nowb/norot/cache1 s01-05) RELEASED (scontrol release), PENDING Priority on apu behind prod3.
  ens_si_s09/s10 completed. Read the ensemble rates when they finish: bench/cs_ens/<arm>/s*/log.out.*
- The rot-5.700 restart EXISTS: bench/cs_ens/ctl2_s01_rst570/rst/dhj.00006.rst (cycle 138353, t 1.7385e6, explicit RT).
  Lesson: Driver::Finalize writes EVERY output at tlim regardless of dt/last_time (a bin dump cannot be fully suppressed).
- THE INSTRUMENT IS BUILT (uncommitted in the working tree, CPU-verified: style-neutral, regression PASS, bitwise
  identical with the flag off): `problem/diag_gid` (+ diag_cycle_min/max) dumps ONE meshblock every cycle to
  cyclediag/<base>.cyclediag.<cycle>.dat, 23 double arrays: w, u, rt_src (rate, pre-limiter), rt_de (applied, LAST
  RK stage only), rt_clip, cond_f1/f2/f3 (conduction-only face fluxes), cond_kappa/keff/dtcell, rt_T, rt_icut,
  rad_w, rad_tauf. New generic hook pgen->user_cycle_func called in Driver::Execute after NewTimeStep. Reader:
  docs/handover/scripts/read_cyclediag.py. Only the split RT path (rt_ck=true) is instrumented.
- IN FLIGHT: GPU rebuild of build_dhj_gpu with the diag, then twins bench/cs_ens/ctl2_s01_twin/{exp,si} restarting
  from rst570 with diag_gid=6 (single apudev one-time jobs each), first analysis in ctl2_s01_twin/analysis/.
  Question: at the abort column (gid 6, k=17, j=2, i=42-47) which term supplies the runaway energy: rt_de or the
  conduction divergence, and does rt_clip fire. If exp reproduces the abort and si does not, the answer is per-cell.
- User rule REPEATED: no chained/long jobs on apudev ([[apudev-one-time-tests-only]]).

## STATE 2026-09-11 ~01:00: cs collapse mechanism CLOSED (see [[cs-vertex-dt-collapse-0907-defaults]] bottom). Instrument
(problem/diag_gid) UNCOMMITTED in the working tree -- commit it (6 files + docs/handover/scripts/read_cyclediag.py). sp_mhd_prod3
at rot 2.9 clean at 00:30; gate rot 4.7 still ahead. 15 ens arms queued on apu (Priority). Analyses: bench/cs_ens/ctl2_s01_twin2/analysis,
bench/cs_ens/analysis/coldmap. Open follow-ups: (1) implicit/sub-cycled RT source for x>1 cells; (2) read the nowb/norot/cache1 ens
rates when done; (3) origin remote embeds a plaintext GitHub token -- rotate it, switch to SSH.
- 2026-09-11 ~02:00: COMMITTED + PUSHED to fork: 1d988879 (cyclediag instrument), 5c58426e (rad_tmax_kappa cap, default off,
  a no-op on semi-implicit, insufficient on explicit). The 15 explicit ensemble arms (nowb/norot/cache1) CANCELLED by the user's
  decision: their question is moot (explicit is not a production candidate). Running: si/s05 post-mortem (bench/cs_ens/analysis/si_s05,
  bench/cs_ens/si_s05_pm) -- the only semi-implicit death (rot 5.3). User asked whether semi-implicit == explicit: answered from the
  record (semi-implicit is 3 days old, 048dff30; the 09-09 switch to explicit was for ablation reproduction; the difference is the
  explicit scheme's overcooling error).
- 2026-09-11 ~03:00: HANDOVER BUNDLE pushed as c26b01ac (docs/handover/HANDOVER-2026-09-11.md, scripts/cs_collapse, claude-memory-2026-09-11). Fork = c26b01ac.

## STATE AT STOP, 2026-09-11 ~03:30 (read this first next session)
- Code: HEAD = fork/polar-average-perf = c26b01ac (docs/handover bundle) on top of 5c58426e (rad_tmax_kappa) and 1d988879
  (cyclediag instrument). Working tree clean. docs/handover/HANDOVER-2026-09-11.md is the handover; memory snapshot in
  docs/handover/claude-memory-2026-09-11/.
- Running on viper: sp_mhd_prod3 only (11569324 RUNNING, chain 11569325-27 pending on apu). At 03:00: rot 23.8, dt 13.6 s,
  1-ME 8.4e30 (peak 8.5e31 at rot 4.3, decaying), zero NaN/FATAL/collapse. It passed the old death point rot 10.3.
  CHECK FIRST: is link 2 running, is 1-ME still decaying, any dt COLLAPSE lines.
- Nothing else running. The 15 explicit ensemble arms are CANCELLED (user decision: moot).
- Science: cs vertex collapse CLOSED (chain in [[cs-vertex-dt-collapse-0907-defaults]] bottom); semi-implicit RT is the production
  default; rescue 1 (rad_tmax_kappa) is a no-op on semi-implicit, insufficient on explicit; si/s05 death = marginal self-healing
  excursion under dt_min 1e-2.
- NEXT (user to choose): (1) production dt_min 1e-2 -> 1e-3; (2) option 2 density floor vs the WB background (not built);
  (3) longer si ensemble 4 seeds x 20 rot; (4) rotate the GitHub token in the origin URL (user only); (5) implicit/sub-cycled RT
  source for x>1 cells (design only). Also open: the sp polar 1-ME long-term behaviour and the jet comparison on cs_hyd_rs.

## UPDATE 2026-09-11 ~05:30 (viper, HELD)
- User decision: DO NOT start any of the next steps (no dt_min change, no density floor, no ensemble). WAITING for
  bug-fixes arriving from ORION. When they land: fetch fork, merge into polar-average-perf, then re-evaluate the next-step list.
- sp_mhd_prod3 checked at 05:10: link 1 (11569324) RUNNING, rot 35.8, dt 12-14 s, 1-ME ~1e31 flat, zero NaN/FATAL. Links 2-4 queued.
- Reminder: bench/ is /viper/u2/jinma/ATHENAK/bench (sibling of athenak/), logs are log.out.<jobid>.

## UPDATE 2026-09-11 ~06:30: orion fixes pulled, HEAD a10e367d (ff from c26b01ac, 10 commits)
- HEAD does NOT build on HIP: src/hydro/hydro_newdt.cpp:223 `dt_diag.template sync<HostMemSpace>()` is illegal when the
  DualView device is HIPSpace. One-line fix (modify_device()/sync_host(), diagnostic path only) sits UNCOMMITTED in the tree.
- A/B old (c26b01ac, bench/wt_c26b01ac/build_gpu) vs new (athenak/build_hip) on apudev, bench/ab_orion/{cs,sp}/{old,new}, 400 cycles:
  NOT bitwise. cs: cycles 0-11 identical, cycle 12 one cell vely 1 float32-ULP; 1e-12 rel in dens/eint to ~cycle 100; 1e-3..1e-2
  by cycle 400 (chaotic amplification). sp MHD: 1e-12 at cycle 104, 1e-10 dens at 400, bcc3 1e-3. old-vs-old repeat bitwise identical,
  so the seed is a real code difference at round-off level (likely FP reassociation in restructured hunks), not a scheme change.
  Not bisected. ab_orion holds 3.7 GB (per-cycle dumps).

## UPDATE 2026-09-11 ~07:30: HEAD = fork = cdd7d2a5 (HIP fix committed+pushed). Bisect DONE, see [[orion-merge-roundoff-8da093f5]].
ab_orion trimmed 3.7 -> 1.5 GB. Worktrees bench/wt_{k3,k6,nEOS,nRT,c26b01ac} still present. Nothing running except sp_mhd_prod3.
## UPDATE 2026-09-11 ~09:00: HEAD = fork = caad9247 (floors_legacy gate). cs+sp gates bit-identical to c26b01ac. Idle.

## UPDATE 2026-09-11 ~11:00: orion-switch arms on the cs si setup (bench/cs_ens/orion, apu)
- HEAD = fork = 18b9c563 (dhj pgen now reads problem/rt_use_cons; it was INERT before, only red_giant read it).
- Arms = si/s01 input (8 rot, same seed) on the new binary: ix1 11591438 (rad_implicit_x1 + rad_cap_ang 0.5),
  cons 11591809 (rt_use_cons), both 11591810, ctl 11591811 (unchanged input, new binary). Notes + recipe in orion/NOTES.md.
- Conduction NEVER bound dt in s01 (zero conduction-dt reports in 207k cycles), so ix1 can only matter in a draining column.
- Gate: "FATAL: timestep collapsed", min dt (hst col 2), drains at rot 5 (analysis/extract.py), KE vs si/s01. ctl must == s01.

## UPDATE 2026-09-11 ~13:00: orion-switch arms READ (all four COMPLETED, apu)
- ix1/cons/both/ctl all reached tlim 8 rot, min dt 5-8.7 s (s01 8.8), zero FATAL/NaN. Logs: bench/cs_ens/orion/<arm>/log.out.<jobid>.
- BUT s01 itself survived 8 rot, so survival on s01 CANNOT discriminate; the switches are shown harmless, NOT shown to rescue.
  The discriminating seed is s05 (the only si death, rot 5.3). Next: rerun s05 input with cons / ix1 / both.
- ctl vs s01 is not bitwise (orion merge): dt diff 19% at rot 8, KE -1.3%. cons +14% KE, both +11%, ix1 -4% vs ctl:
  larger than the chaos bound; rt_use_cons may change the answer -> needs the seed ensemble before any claim.
- Gate item 3 (drained columns at rot 5, analysis/extract.py on <arm>/bin/) NOT run yet.
- sp_mhd_prod3: link 2 (11569325) RUNNING, rot 107.9, dt 12-14 s, clean; 1-ME still decaying (3.6e30). Old "rot 35.8" is stale.
- 2026-09-11 ~14:00: s05 arms SUBMITTED on apu (user granted job-submission permission): s05_ctl 11610533, s05_cons 11610534,
  s05_ix1 11610535, s05_both 11610536 in bench/cs_ens/orion. All on the 18b9c563 binary. Read: does s05_ctl still die ~rot 5.3;
  does any switch outlive it. cs_mhd_prod3 (cs twin of sp_mhd_prod3, no WB, ix1 + cons + ck RT) being staged in bench/cs_mhd_prod3.
- 2026-09-11 ~15:00 cs_mhd_prod3 STAGED (bench/cs_mhd_prod3, NOT submitted). Two blockers found while staging:
  (1) rad_implicit_x1 is HYDRO-ONLY (conduction.cpp fatal ~l.173; no MHD::ImplicitConduction task).
  (2) BUG: rt_use_cons in MHD is WRONG -- EintFromCons (utils/eint_from_cons.hpp) subtracts KE and rho*Phi but NOT the
      magnetic energy, so the eiN accessor in two_stream_rt.hpp feeds ME into T (beta~1 at the photosphere -> O(1) error).
      No MHD run has used rt_use_cons yet (sp_mhd_prod3 binary predates the flag), so no result is contaminated.
  Opus agent implementing: MagEnergyCC + ME subtraction in both callers, MHD::ImplicitConduction after id.ct, HIP rebuild,
  apudev smoke (bench/cs_mhd_prod3/smoke). Submit cs_mhd_prod3 only after the smoke passes and the diff is reviewed.
- 2026-09-11 ~15:45: BOTH BLOCKERS FIXED in bc972cc8 (committed, NOT pushed): MagEnergyCC + MHD EintFromCons overload
  (bcc0 is orthonormal-frame, plain sum of squares is right on cs); MHD::ImplicitConduction after MHDSrcTerms, before the
  u0 send (after CT the ghosts would miss it; O(dt) ME mismatch accepted); MHD::Fluxes now sets stage_beta_dt, without
  which rad_cap_ang was INERT in MHD while NewTimeStep dropped dt2/dt3. Smoke apudev 11611107 (bench/cs_mhd_prod3/smoke):
  cons+implicit vs explicit agree to O(dt) over 20 cycles. cs_mhd_prod3 SUBMITTED, 4 links on apu: 11611223/29/30/31.
  Gate: compare with sp_mhd_prod3 (rot ~108) and cs_mhd_prod2; first check ~4 rot/h -> rot 5 by ~17:00, rot 42 (old
  cs_mhd_prod death) by tomorrow. The smoke did NOT exercise a low-beta cell; watch the rt_eiclamp warning and 1-ME.
- 2026-09-11 ~19:00 cs SPEED: PROFILED (bench/prof, rocprofv3 --kokkos-trace, apudev 11614468; stats in
  bench/prof/{cs,sp}/prof_kk/rank_*_kernel_stats.csv). cs 62.0 ms/cycle vs sp 38.5 on 0.75x the cells = 2.15x per zone;
  GPU busy 78%/72% so NOT host-bound, NOT launch-bound. Two kernels carry the gap: SendBuff (cc+fc pack share the label)
  14.9x per call on cs = 25.6% -- the seam TransformMomentum/TransformField geometry + resample weights are recomputed
  for EVERY radial cell i though they depend only on (kk,jj); radimpx1 21.6% -- one thread per column (12k/rank), serial
  EOS T inversion + Cv + RadFaceKappa per cell. SeamFlux* 6.5%, same per-i redundancy. Opus implementing: hoist the
  geometry per (kk,jj) (bitwise-identical arithmetic) in bvals_cc/bvals_fc/flux_seam_cc; split radimpx1 into a
  per-cell coefficient par_for + per-column sweeps. Validation = cmp of hst on CPU and GPU, timing in bench/prof/cs_new.
- 2026-09-11 ~22:30 cs SPEED-UP COMMITTED + PUSHED: 1bd31298 (seam geometry hoist in cc SendBuff + SeamFluxSend*,
  BITWISE on CPU and GPU) and 77de5618 (radimpx1 split into coef/face/sweep kernels, 8.6x; CPU bitwise, GPU differs
  1 ulp in 1 cell after 3 cycles from hipcc contraction, deterministic, NO race -- revert alone if needed). Aggregate on
  the cs prod restart: 62.0 -> 43.8 ms/cycle (1.42x) with switches, 40.0 without. rt_chain_ck is now the top kernel (30%).
  fc SendBuff (1432 ms/300 cyc) UNTOUCHED: the register hoist was 6x slower (spills + serial i); a precomputed geometry
  TABLE is being built by Opus in worktree /viper/u2/jinma/ATHENAK/wt_fcseam (branch fc-seam-table, on top of Task A).
  cs_mhd_prod3: binary SWAPPED to 77de5618 at ~22:45 on user instruction (links 2+; old kept as athena.bc972cc8).
- 2026-09-11 ~23:55 fc SEAM TABLE COMMITTED + PUSHED 27ca5b13 (bitwise CPU+GPU). Final cs prod restart: 40.5 ms/cycle with
  switches (was 62.0, 1.53x), 37.2 without vs sp 38.5 -> the cs GRID gap is CLOSED. Remaining: rt_chain_ck 34% (grid-neutral),
  radimpx1 365 ms. cs_mhd_prod3 binary swapped to 27ca5b13 for links 2-4. Worktrees wt_fcseam / wt_fcref can be removed.
## STOPPED 2026-09-12 ~00:30 (viper). Session summary
- Pushed fb69c554. Production: cs_mhd_prod3 link 1 on bc972cc8 image, links 2-4 on 27ca5b13 (bench/cs_mhd_prod3/NOTES.md,
  BUILD_COMMIT.txt). sp_mhd_prod3 link 2 running. s05 arms all COMPLETED clean at 8 rot (bench/cs_ens/orion/s05_*).
- Speed work: profile method = bench/prof (rocprofv3 --kernel-trace --kokkos-trace --stats via wrap.sh, 300-cycle restart of
  the production run; job ~2 min on apudev). Next lever rt_chain_ck (2.8 s / 8.2 s GPU on cs; also 40% on sp). radimpx1 365 ms.
- OPEN science: cs_mhd_prod3 vs sp_mhd_prod3 (no WB, ix1+cons+ck): watch rt_eiclamp count (was 3 cells at startup),
  1-ME, dt band, drained columns (cs_ens/analysis/extract.py); the ME subtraction has not yet been exercised at beta~1.
- The switches (ix1, cons) are harmless on 2 seeds x 4 arms; rescue is unproven because no control dies on the new binary.
  A rescue test needs a seed that dies on the CURRENT binary: run a 4-seed si ensemble on 27ca5b13 first, then arm the dead seed.
## UPDATE 2026-09-12 ~05:00: cs_mhd_prod3 READ at rot 8.5 (link 1, 11611223, still on the bc972cc8 image, 60.7 ms/cyc)
- CLEAN: no FATAL/NaN, dt 12.6-14.4 s (sp 13.0-13.7), rt_use_cons clamp fired ONCE (3 cells, cycle ~1300), rt_de_max once.
- 1-ME flat 4-6e31 through rot 8 = sp band; no radial-field runaway. ZERO drained columns at every dump 0-8 rot
  (coldmap threshold rho(i=52)<1e-9; min rho 4e-9, stopped falling rot 6->8) vs 35-65 on the surviving cs_ens si arms.
- OPEN: cs KE 7x below sp at rot 8 (1.7e34 vs 1.2e35); cs tot-E -0.12%/8 rot, sp +4.5% (sp excess is [[sp-mhd-energy-excess]]).
- sp_mhd_prod3 link 2 at rot 131, 31.6 ms/cyc, clean; 1-ME decayed to 1.1e30. Runs live in /viper/u2/jinma/ATHENAK/bench/.
- Drain tool = bench/cs_ens/analysis/coldmap/coldmap.py (density-only regrid); extract.py is NOT it. P = 3.05e5 s.
- NEXT: reread cs_mhd_prod3 at rot ~20 (~3 h) and rot 42 (old cs_mhd_prod death); confirm link 2 picks up 27ca5b13 speed.
- 2026-09-12 ~05:40: cs_mhd_prod3_wb SUBMITTED (bench/cs_mhd_prod3_wb, links 11618999/11619000/01/02, apu), = cs_mhd_prod3
  + sp's <mhd> WB block (wellbalance_dynamic, wb_x1, wb_option=polytropic, wb_cache_every=10), binary 27ca5b13 hard link.
  Smoke apudev 11618946 passed. Purpose: isolate WB as the cause of cs KE 7x < sp and the opposite tot-E drifts.
  Riemann solver/recon/integrator identical across cs and sp (hlld/plm/rk2); cs x1 recon is stretch-aware since 4a16f07b.
  Gate: KE(rot), tot-E(rot), drains (coldmap.py), 1-ME vs cs_mhd_prod3 and sp_mhd_prod3 at rot 2, 4, 8.
- 2026-09-12 ~09:30: pulled orion's FOFC handover (e3410692); memory snapshot merged (122 new files). Handover item 2 DONE:
  see [[fofc-gpu-verified]]; fix b18f3205 pushed. Open from orion: (1) RT-guard gating for dhj bitwise (user's call),
  (3) wiki entries, (4) RG_fofc outcome on orion. cs_mhd_prod3_wb links 11618999-11619002 queued/running on apu.
- 2026-09-12 ~12:30 RG_fofc ON GPU SUBMITTED: bench/RG_fofc links 11625391/11625392 (apu, 1 node, 2 MI300A, 24 blocks
  480x16x16 = 109 ms/cycle vs 121 at 96 blocks; ~30k cycles to tlim 9e5 => ~1 h). Fixes on the way: da9734fd (red_giant
  sponge host-global capture, HIP compile), 955c39df (hydro/scratch_level=1, LDS cap at nx1=480). Gate = rg.log fofc
  column first nonzero, deaths to beat t=5.7e5/5.87e5/6.13e5, i8.py sections 1-3/4a. GPU = different chaotic realization.
- 2026-09-12_05:18 cs_mhd_prod3_fofc SUBMITTED (bench/cs_mhd_prod3_fofc, 4 links, see NOTES.md for ids): cs_mhd_prod3 +
  fofc=true/llf, cs_lowbeta_fallback=0.0, nghost=3 (from scratch). Binary build_hip_dhj at 955c39df (first dhj GPU binary
  with FOFC + orion RT guard). Gate vs cs_mhd_prod3: dt band, 1-ME, KE, drains, dhj.log fofc column (zero => vacuous).
  ALSO RUNNING: RG_fofc GPU 11625391 (rg.log dfloor from t=1.6e5 = cube-vertex corner ghost extrapolation going negative
  in the inflating layer above the join; active cells clean) and CPU twin RG_fofc_cpu 11626142 (apu host cores).
- 2026-09-12 ~08:00: RG_fofc DONE (survived to 9e5 on GPU and on orion CPU: FOFC replaces vceil, [[rg-fofc-gpu-survived]]);
  corner limiter 2e3ae96e pushed; RG_fofc_clamp long run 11633164/5 to 3e6. cs_mhd_prod3_fofc (11627070) alive, dt recovered
  to 10.5 s by rot 2.8 (2-rot transient vs 1). cs_mhd_prod3 rot 36 clean. sp_mhd_prod3 rot 167. cs_mhd_prod3_wb DEAD rot 12.3.
- 2026-09-12 ~14:00: sp excess SOLVED ([[sp-polar-hlle-swap-is-the-excess]], 7d87f3c5 pushed); sp_mhd_nopole 11639436/7
  running (gate rot 10-12 polar stability). Pushed today: b18f3205 c6cbab62 111a0eed 955c39df da9734fd 2e3ae96e 7d87f3c5.
  cs deaths: fofc arm died rot 11.4, WB arm rot 12.3, BUT restarts (bench/cs_death_ab) survive: fofc_new rot 11->12 clean,
  prod3_new rot 10->13 clean on HEAD binary => rare chaotic event, not a regression; prod3_old control pending.
  RG_fofc_clamp long run 11633164 to 3e6 in progress. rt sweep split measured 1.13x, shelved.
- 2026-09-12 ~15:10: RG_fofc_clamp long run DIED t=1.119e6 (corona T 1e17 K cell, r 3.9e12; see [[rg-fofc-gpu-survived]]);
  killed. Control RG_fofc_ctl 11639634 (unclamped, 9e5->1.3e6) + diag 11639635 (clamped HEAD, 1.1e6->1.125e6 with
  hydro_fofc dumps) running. sp_mhd_nopole 11639436/7 running. cs_death_ab prod3_old control pending.
