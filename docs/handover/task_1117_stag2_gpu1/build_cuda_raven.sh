#!/bin/bash
# build_cuda_raven.sh <athenak worktree at clean-cc-stag2-gpu1> [nmax ...]: Raven A100 binaries
# Recipe = your Raven A100 production builds (gcc/13 cuda/12.6 openmpi_gpu/5.0, nvcc_wrapper,
# Kokkos AMPERE80, CUDA malloc async off, MPI, Release) with -ffp-contract=off: the fmad-1031
# CMake logic then appends --fmad=false for the device code (check flags.make).  Built-in pgens.
WT=$(readlink -f $1); shift; NM=${*:-432}
source /etc/profile.d/modules.sh; module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake
CFG="-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF
     -D CMAKE_CXX_COMPILER=$WT/kokkos/bin/nvcc_wrapper -D Athena_ENABLE_MPI=ON
     -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off"
for N in $NM; do
  B=$WT/build_a100_n$N
  X=""; [ $N != 432 ] && X="-D Athena_M1_VGD_NMAX=$N"
  cmake -S $WT -B $B $CFG $X > $B.cmake.log 2>&1 || { echo "CMAKEFAIL $N"; continue; }
  grep -m1 'fmad-1031' $B.cmake.log
  nice -n 10 make -C $B -k -j16 > $B.make.log 2>&1; echo "make n$N rc=$?"
  grep -l -- '--fmad=false' $B/src/CMakeFiles/athena.dir/flags.make > /dev/null && echo "flags.make: --fmad=false present"
  grep -E 'error|warning' $B.make.log | grep -E 'rad_m1_clean|rad_m1_rhse|rad_m1_vet|rad_m1_vgd' | sort | uniq -c > $B.clean_diag.txt
  [ -x $B/src/athena ] && md5sum $B/src/athena
done
