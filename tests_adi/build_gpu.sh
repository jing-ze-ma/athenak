#!/bin/bash -l
# GPU (MI300A APU) build of one PROBLEM into wt_he4_adi/build_gpu_<name>.
# usage: build_gpu.sh <pgen-name-or-"default"> [srcdir]
set -u
name=$1
src=${2:-/viper/u2/jinma/ATHENAK/bench/wt_he4_adi}
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
b=$src/build_gpu_$name
args=(-D CMAKE_BUILD_TYPE=Release
      -D CMAKE_CXX_COMPILER=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4/bin/hipcc
      -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On
      -D Athena_ENABLE_MPI=ON)
if [ "$name" != "default" ]; then args+=(-D PROBLEM=$name); fi
cmake -S "$src" -B "$b" "${args[@]}" || { echo "CMAKE_FAILED"; exit 1; }
make -C "$b" -j 16
echo "MAKE_EXIT=$?"
[ -x "$b/src/athena" ] && echo "BUILD_DONE_OK $b" || echo "BUILD_DONE_FAIL $b"
