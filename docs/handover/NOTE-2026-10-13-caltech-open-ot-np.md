# NOTE 2026-10-13 from Caltech: open-ot NP (re TASK-2026-10-13-caltech-open-ot-np.md)

## QUEUED (10-12 ~23:50 PDT)
- Binary: 281dd686 + patch/0001-0003 (git apply onto an archive of 281dd686) -> athena_cpu_omp_281dd686_np3, md5 0d486e6a7421d3600d4febc0b27dbe04.
  **Patch 0002 also changes `kokkos` from the submodule to a symlink** (an artefact of the tree it was made in); applied with --exclude=kokkos
  (our kokkos = the repo's pinned submodule). All source hunks of 0001-0003 applied cleanly.
- Debug job 4432031 (after the CV rounds): the 4 runs in parallel at **8 OpenMP threads each** (list.txt says 1; level 64 is ~25-30 min serial,
  which does not fit the 28-min debug cutoff). Results + rc follow here.
