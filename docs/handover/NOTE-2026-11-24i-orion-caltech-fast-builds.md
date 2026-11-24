# NOTE 2026-11-24i orion -> Caltech: faster AthenaK builds on H200 (USER: make builds fast everywhere)

Most of our build time was the file system, not the compiler. On Orion, the same rc1 build took 238 s with the build dir on GPFS and 51 s with it in /dev/shm. On GPFS, a one-file rebuild spent 65 s of its 101 s build step in the final link, because the linker reads 260 freshly written objects back from GPFS. The recipe below keeps the toolchain you use now: CUDA 12.9 + g++ 13, never 12.2. It changes only where the build happens and what gets cached.

## Recipe (one H200 node, once per build)
1. **Build dir on node-local storage.** Use a fixed path that is the same on every node, e.g. `L=/dev/shm/$USER/fb` (or the node's local scratch).
   - Copy the frozen Kokkos once: `cp -a <athenak>/kokkos $L/kokkos`. It is the release submodule 4.6.02 (d8e9af03).
   - Sources go in `$L/w/<slot>/src`, the build in `$L/w/<slot>/b`.
   - One slot per arch + flags, e.g. `h200_vet`. Reuse it for every commit: copy new sources in with `rsync -rlc --delete` so unchanged files keep their mtime, and ninja recompiles only the changed TUs.
2. **ccache on the shared file system**, the same for every node, with these settings:
   - `CCACHE_DIR=<shared>/ccache_h200`, `max_size = 20G`;
   - `CCACHE_BASEDIR=$L/w`, `CCACHE_NOHASHDIR=1`;
   - `CCACHE_SLOPPINESS=include_file_mtime,include_file_ctime,time_macros,locale`;
   - `CCACHE_COMPILERCHECK="string:caltech-cuda12.9-gcc13"`, changed whenever a module changes;
   - `CCACHE_TEMPDIR=/tmp/ccache-$USER`.

   A new slot or node then compiles an already-built commit as 100 % hits: on Orion, 254/254 hits in 9 s, with a byte-identical binary.
3. **CMake** (Ninja, Release, only the athena target):
   ```
   export NVCC_WRAPPER_DEFAULT_COMPILER=mpicxx      # the g++ 13 MPI wrapper
   cmake -S $L/w/h200_vet/src -B $L/w/h200_vet/b -G Ninja -D CMAKE_BUILD_TYPE=Release \
     -D CMAKE_CXX_COMPILER=$L/kokkos/bin/nvcc_wrapper -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
     -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_ENABLE_MPI=ON \
     -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
     -D Athena_VET_VGD_NMAX=512 -D PROBLEM=he_star_vet
   ninja -C $L/w/h200_vet/b -j $(nproc) athena
   ```
   - Re-run cmake only when these arguments change.
   - Check that `--fmad=false` appears in build.ninja (`-ffp-contract=off` adds it).
   - Point CMAKE_CXX_COMPILER at the one fixed nvcc_wrapper path. A per-tree path is baked into config.hpp and makes every TU miss the cache.
   - Build one arch only: Hopper 90.
4. **Copy the binary out** of `/dev/shm` to the shared file system, and record next to it:
   - md5, commit, `git diff HEAD` (if dirty), the cmake args, `nvcc --version`, `g++ --version` and the Kokkos version.

   Our .prov format is in /orion/ptmp/jinma/clean_1115/fastbuild/README.md.
5. **Gate every new binary** with the tiny BSG on 1 H200: our bsgs_ia121 deck at mesh/nx2=nx3=32, meshblock 32x32, `time/nlim=0`, plus `rad_vet/vet_gd_theta_seam_sym=true` in the deck.
   - `clean_init steady pass 1` must be **3.686e-11**.
   - Orion CPU, the rc1 build: 3.686e-11.
   - A CUDA 12.2 build fails this kind of gate. The old wngl smoke gave 3.499e-08 instead of 4.020e-06.

## What to expect (Orion CPU, rc1 da92efe70, he_star_vet, same node)
| case | build dir on GPFS (before) | /dev/shm + ccache (after) |
|---|---|---|
| full build, nothing cached | 238 s | 68 s (51 s without ccache) |
| new slot, commit already cached | ~30 s | 9 s |
| one .cpp changed (rad_vet_clean.cpp) | 130 s (link 65 s) | 31 s |
| rad_vet.hpp changed (52 TUs) | 295 s (MERGE123's real 46-TU rebuild) | 46 s |

The GPU numbers will differ. nvcc's critical path is ptxas on rad_vet_vetgd_sweep.cpp (~7 min on A100) and rad_vet_implicit (~4 min), so a cold H200 build stays several minutes. With a warm slot and the ccache, a typical one-commit rebuild should be 1-3 min.

Please post one cold and one incremental wall time with this recipe, and the tiny-gate pass-1 value.

Apply when convenient (between arms); keep your existing rc builds until a new fb-style build passes the tiny gate.

-- orion (coordinator)
