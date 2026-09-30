---
name: rocm72-adopted-0929
description: user 09-29 adopted the ROCm 7.2 stack on viper (gcc/16 rocm/7.2 openmpi_gpu/5.0 + Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF): ~10 % faster; build targets *_gpu72; jobs must load the same modules
metadata:
  type: project
---
gputune_0928: ROCm 7.2 + HIP malloc async OFF = -9.3 % (2 GPUs) / -10.3 % (4 GPUs) per cycle on WASP-121b; the gain is
the compiler (6.3 + async off = bitwise, no gain). Plain 7.2 segfaults (Kokkos 4.6.02 writes View headers into
hipMallocAsync memory). Equal-time noise gate: pass at 1e-3..100 bar, 1.21-1.28x above 100 bar (no bias). MPCDF lists
gcc/14 rocm/6.3 as valid only until 08/2026. Other candidates (ZEN4, multi-kernel instantiation, gpu-bind, UCX) = noise.
**How to apply:** build with /viper/ptmp2/jinma/builds/build_inc_viper.sh dhj_gpu72|box_gpu72|none_gpu72|rg_gpu72;
job scripts `module load gcc/16 rocm/7.2 openmpi_gpu/5.0` (a 7.2 binary cannot run under 6.3 modules); keep
HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1. The viper WASP-121b production (w121prod_0929) uses it: jobs 12018387/89.
