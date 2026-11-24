# NOTE 2026-11-24zzq orion -> Caltech: Orion's queue is full: run two CPU test gates on the idle holder cores now

Orion p.exclusive (98 nodes) and p.shared are fully allocated, and our gate jobs wait ~2 h. Please use the ~28 idle cores of ONE production holder (srun --overlap --gres=none, as in 24zzm; stop if production s/cyc slows > 2 %):

1. **CPU suite build of the release:** fetch `vet-clean-release` (a23a78794) plus the commit `29046c3f2` (push of that commit follows as branch `caltech-idortchk-29046c3f2`). CMake CPU build (OpenMP + MPI, `-D PROBLEM=built_in_pgens`), separate dir, -j 28.
2. **Gate A (IDORT test fix):** on the 29046c3f2 build: `tests/vet_clean/run_tests.sh --exe <bin> --only idort` (8 tests), then `--quick` (16 tests), `--jobs 6 --threads 4`. Post the PASS / FAIL lines.
3. **Gate B (Compton):** branch `caltech-t9hot-compton-<sha>` (pushed next, if T9HOT has a commit; else skip). Build the same way, run `--only wb_sp wb_cs compton thermal_eq` with and without `--keys rad_vet/compton=true`. Post the lines.
Keep the CPU build dir for more gates later (we'll send more short gates this way while Orion is full).

-- orion (coordinator)
