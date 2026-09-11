---
name: orion-build-openmp-and-module-traps
description: "two orion build traps that each cost a wasted job: a build dir with Kokkos_ENABLE_OPENMP=OFF runs 3.7x slow on a 96x7 job with NO error, and `module load` does not survive between tool calls so make silently falls back to system gcc 7.5"
metadata:
  type: feedback
---

Both hit on 2026-09-08 while setting up the red giant production chain.

## 1. An OpenMP-less build is SILENT, and 3.7x slower

`build_rg_topre` was configured for the 6-rank testbed with
**`Kokkos_ENABLE_OPENMP=OFF`**. Submitted at the production shape (96 ranks x
`--cpus-per-task=7`, `OMP_NUM_THREADS=7`) it ran perfectly correctly and 3.7x slower --
0.165 s/cycle against `build_rg_diag`'s 0.0452 -- because each rank was serial, so the
job used 96 of 672 cores. **Nothing warns.** dt, mass and the physics were all fine.

**Why:** a 40 h campaign estimate becomes 150 h, and the only symptom is the wall clock.

**How to apply:** before submitting a multi-node job, check the binary's build:
`grep -E "^(Kokkos_ENABLE_OPENMP|CMAKE_BUILD_TYPE):" <build>/CMakeCache.txt`.
Testbed builds and production builds are not interchangeable on orion.
`build_rg_prod` is the OpenMP-on red_giant build; `build_rg_diag` also has it on.

## 2. `module load` does not persist between tool calls

The Bash tool starts a fresh shell each call, so `module load gcc/13 openmpi/4.1` in one
call and `make` in the next means make runs against **system gcc 7.5**, and Kokkos stops
with `#error "Compiling with GCC version earlier than 8.2.0 is not supported."`.
Configuring with `-D CMAKE_CXX_COMPILER=$(which mpicxx)` is not enough -- mpicxx is a
wrapper and finds whatever gcc is on PATH at compile time.

**How to apply:** put the module load and the build in ONE command:
`module purge && module load gcc/13 openmpi/4.1 && make -C <build> -j 32`.

Related: [[red-giant-flux-deficit-is-spinup]], [[freya-job-submission]].
