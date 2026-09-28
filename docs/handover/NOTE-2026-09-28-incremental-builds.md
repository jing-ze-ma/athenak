# NOTE for the Caltech and DeltaAI sessions: build incrementally (2026-09-28)

The user asked all machines to stop doing full rebuilds (fresh `git archive` + new build dir each time, Kokkos
included, 15-25 min) and to rebuild only what changed.

**Pattern** (viper: `/viper/ptmp2/jinma/builds/build_inc_viper.sh`, being set up 09-28):
- One persistent git worktree + one persistent CMake build dir per TARGET (problem x backend), e.g. dhj_gpu,
  box_gpu, none_gpu, dhj_cpu. Configure once.
- To build a commit: `git -C <worktree> checkout --detach <commit>`, then `reset --hard` and check that
  `git status --porcelain` is empty (never build a dirty tree), then `make -j`. A checkout rewrites only the files
  that changed, so make recompiles only the affected objects; Kokkos and untouched files are reused.
- PITFALL: do NOT `git archive` a new snapshot on top of an existing build tree: archive stamps every file with
  the commit time, which can be older than the object files, and make then misses real changes.
- Take a lock per target (flock) so two agents never build the same target at once.
- Copy the binary out with its commit and md5 (e.g. `bin/athena_<target>_<sha>`), and log it; reuse an existing
  binary of the same commit instead of rebuilding (check the log/md5).
- A header included everywhere (e.g. src/utils/two_stream_rt.hpp) still recompiles everything that includes it,
  but not Kokkos.

**Caltech:** you already have `docs/handover/caltech-2026-09-26/scripts/build_inc.sh`. Make sure it follows the
rules above (in particular the git-archive mtime pitfall: if it unpacks an archive over a persistent tree, switch
to a worktree checkout or an rsync that only touches changed files), and use it for every build, including the
cross-cluster benchmark (`TASK-2026-09-28-cross-cluster-timing.md`, commit 11c9a5be).

**DeltaAI:** set up the same pattern (targets: dhj_gpu for the benchmark first; box_gpu/none_gpu as needed), and
record the script and its first timings in `DELTAAI_FACTS.md`.
