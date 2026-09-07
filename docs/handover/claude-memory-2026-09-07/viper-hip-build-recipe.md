---
name: viper-hip-build-recipe
description: The module loads and cmake flags that build AthenaK for GPU on MPCDF Viper (MI300A APU nodes)
metadata: 
  node_type: memory
  type: project
  originSessionId: 727372a4-3557-4a71-b294-5104e4fbe37f
  modified: 2026-08-16T15:05:08.119Z
---

On Viper (MPCDF), build AthenaK for the accelerator with:

```bash
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cmake -B build \
  -DAthena_ENABLE_MPI=ON \
  -DKokkos_ENABLE_HIP=ON \
  -DKokkos_ARCH_AMD_GFX942_APU=ON \
  -DCMAKE_CXX_COMPILER=hipcc \
  -DCMAKE_HIP_ARCHITECTURES=gfx942 \
  -DCMAKE_BUILD_TYPE=Release \
  -D PROBLEM={name}
```

**Why:** MPI modules on Viper are hierarchical — `openmpi_gpu/5.0` is not visible until a compiler is loaded, so `mpicxx` is absent from a bare login shell and a plain `cmake -B build` configures a no-MPI build. That build used to FAIL, because `src/bvals/physics/bfield_bcs.cpp` called `MPI_Allreduce` unguarded by `MPI_PARALLEL_ENABLED`. **That is FIXED as of `polar-average-perf` (checked 2026-08-22): a plain `cmake -S . -B <dir> -D PROBLEM=... -D CMAKE_BUILD_TYPE=Release` with the bare-login g++ 11.5 configures a SERIAL Kokkos build and compiles in ~60 s at -j256**, and that is what [[correlated-k-design]]'s regression test uses. `Kokkos_ARCH_AMD_GFX942_APU` confirms Viper's GPU partitions (`apu`, `apu1`, `apudev`; nodes `vipa[1001-1300]`, 2 per node) are MI300A APUs with unified host/device memory.

**How to apply:** Use this whenever building on Viper, and prefer it over a CPU-only configure for verifying device code — it is the only local way to compile-check kernels for the GPU. Login node `viper13` has no GPU, so it compiles but cannot run. See [[use-fork-not-origin]].
