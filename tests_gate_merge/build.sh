#!/bin/bash -l
# GPU (MI300A APU) red_giant build: build.sh <srcdir> <builddir>
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cmake -S "$1" -B "$2" -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_COMPILER=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4/bin/hipcc \
  -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
  -D Athena_ENABLE_MPI=ON -D PROBLEM=red_giant || { echo CMAKE_FAILED; exit 1; }
make -C "$2" -j 12
echo "MAKE_EXIT=$?"
[ -x "$2/src/athena" ] && echo BUILD_DONE_OK || echo BUILD_DONE_FAIL
