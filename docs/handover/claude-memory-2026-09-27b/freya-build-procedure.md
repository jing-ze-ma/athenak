---
name: freya-build-procedure
description: How to build AthenaK on the Freya (MPCDF) cluster with MPI+OpenMP
metadata: 
  node_type: memory
  type: project
  originSessionId: 0df72474-f27b-49d2-9277-124d645c5193
  modified: 2026-08-09T13:08:51.675Z
---

AthenaK build procedure on Freya (MPCDF cluster), Sapphire Rapids nodes, MPI+OpenMP:

```bash
export SYSTYPE="Freya-OpenMPI"
module purge
module load gcc/13 openmpi/4.1 cmake/4.0
cd build && find . -mindepth 1 -delete      # clean ONLY build/ contents
cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_SERIAL=ON \
  -DAthena_ENABLE_MPI=ON -DKokkos_ARCH_SPR=ON \
  -D PROBLEM=<name>
make -j 32
```

- **Use this same recipe on `orion` too, not just Freya** — the user confirmed this explicitly. `/mpcdf` is mounted on orion, so the Freya toolchain is what compiles there as well.
- **`export SYSTYPE="Freya-OpenMPI"` is mandatory and must come BEFORE `module load`.** Without it the module loads silently no-op, `g++` stays at the system 7.5.0, and the build dies with ``bad value ('sapphirerapids') for '-march='``. That error means the modules did not load — it is not a wrong-architecture problem.
- Every Bash call is a fresh shell, so re-export `SYSTYPE` and re-run `module load` in the *same* command as `make`.
- Rebuilding an already-configured `build/` only needs the env + `make -j 32`; a full clean+`cmake` is not required just to recompile.
- `<name>` is the `src/pgen/` file basename (current task uses `solar_convection`).
- Executable lands at `build/src/athena`.
- **When cleaning, delete only the contents of `build/`** — never `build/` itself or anything in the `athenak` root. The user explicitly corrected me on this; use a guarded `find build -mindepth 1 -delete`, not `rm -rf build/*` with a bare `cd`.
- Also documented in [[CLAUDE.md]] under "Building on Freya".

**This recipe applies to the regression-test harness too — do NOT improvise a plain `g++`
build for `tst/run_tests.py`.** (Corrected on this point 2026-08-09 after trying exactly
that and getting `Compiler: GNU 7.5.0 ... Configuring incomplete`.) Three extra things the
harness needs:

1. `tst/scripts/utils/athena.py` invokes **`cmake3`**, which does not exist on
   orion/Freya. Put a shim early on `PATH` that execs the *module's* cmake — resolve it
   via `PATH` (`exec cmake "$@"`), NOT a hardcoded `/usr/bin/cmake`, which is an ancient
   version and fails the same way.
2. Pass the whole flag set through repeated `--cmake=` arguments, e.g.
   `--cmake=-DCMAKE_CXX_COMPILER=mpicxx --cmake=-DCMAKE_C_COMPILER=mpicc
   --cmake=-DCMAKE_BUILD_TYPE=Release --cmake=-DKokkos_ENABLE_OPENMP=ON ...`.
3. `run_tests.py` builds into **`tst/build`** and `rm -rf`s that directory, not the
   top-level `build/` — so running tests does not destroy the main working build.
