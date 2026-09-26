---
name: index-working-practice
description: Memory links for working practice, machine/build/job conventions on viper, orion and freya.
metadata:
  type: reference
---

## Working practice

- **[VALIDATE THE INSTRUMENT](validate-the-instrument.md) — "nothing" usually means "I did not look"**
- **[Measure impact before claiming it](measure-impact-before-claiming.md) — report the defect, not consequences**
- [Localise by dilution](localise-by-dilution.md) — vary domain size at fixed dx
- [Don't propose confounded tests](confounded-tests-rejected.md) — user vetoes non-discriminating tests
- **[apudev: ONE-TIME TESTS ONLY](apudev-one-time-tests-only.md) — chained/production jobs go to `apu`**
- [HIP DualView sync idiom + hipcc contraction](hip-dualview-sync-idiom.md) — modify_host/sync_device only; kernel splits lose GPU bitwise
- [Job submission permitted](job-submission-permitted.md) — from the main session
- [Use fork, not origin](use-fork-not-origin.md) — push to jing-ze-ma/athenak
- **[Restart output dt trap](restart-output-dt-trap.md) — last_time from the restart; dt is sim time**
- [Never write in run/](never-write-in-run-dir.md) — read-only
- [Viper HIP build recipe](viper-hip-build-recipe.md) — cmake flags for MI300A
- [Reusable CPU baseline binary](bench-baseline-worktree.md) — bench/base_wt
- [Fork git setup](athenak-fork-git-setup.md) — origin = jing-ze-ma/athenak fork over SSH (all git ops); upstream = IAS-Astrophysics; the user's terminal cannot copy OUT
- [Freya build procedure](freya-build-procedure.md) — how to build AthenaK on Freya (MPI+OpenMP, SPR); clean only build/ contents
- [Freya/Orion job submission](freya-job-submission.md) — p.shared vs p.exclusive, keep --cpus-per-task with OMP, and an orion node has 112 PHYSICAL cores (the 224 SLURM reports are hyperthreads)
- [GPFS quota wall](gpfs-quota-wall.md) — df lies; a per-user quota SILENTLY truncates AthenaK dumps without any error
- [Orion build traps: OpenMP-off binary + module load](orion-build-openmp-and-module-traps.md) — a Kokkos_ENABLE_OPENMP=OFF build is SILENT and 3.7x slow at 96x7; `module load` does not survive between tool calls, so make falls back to gcc 7.5
- [run/ is untouchable](run-directory-untouchable.md) — ~40GB of simulation data; never add, clean, or modify it
- [Scratchpad invisible to compute nodes](scratchpad-not-visible-to-compute-nodes.md) — /tmp is node-local; stage SLURM jobs on /orion
- [Test output location](test-output-location.md) — all run output/analysis goes in /orion/ptmp/jinma/Athenak/ (separate GPFS, no quota); never /orion/u or /tmp
- **[VIPER HIP CONVENTIONS (standing, 09-11): DualView modify_device/sync_host only; kernel splits change GPU round-off; CPU bitwise != GPU bitwise; build on viper before calling GPU code done](viper-hip-code-conventions.md)**
- [Work pace preference](work-pace-preference.md) — move faster / less deliberation on routine tasks
