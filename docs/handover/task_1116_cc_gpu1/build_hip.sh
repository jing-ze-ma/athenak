#!/bin/bash
# build_hip.sh <athenak worktree at 5ddf0e886>: the two MI300A binaries (NMAX 432 default, 4096)
# Recipe = your MI300A he_star_m1 builds (gcc/16 rocm/7.2 openmpi_gpu/5.0, hipcc, GFX942_APU,
# HIP malloc async off) + -ffp-contract=off (clang applies it to host AND device code; the
# fmad-1031 CMake logic is CUDA-only and does nothing here).  Built-in pgens (no PROBLEM).
WT=$(readlink -f $1)
source /etc/profile.d/modules.sh; module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
HIPCC=$(which hipcc)
CFG="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=$HIPCC
     -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release
     -D CMAKE_CXX_FLAGS=-ffp-contract=off"
for N in 432 4096; do
  B=$WT/build_hip_n$N
  X=""; [ $N = 4096 ] && X="-D Athena_M1_VGD_NMAX=4096"
  cmake -S $WT -B $B $CFG $X > $B.cmake.log 2>&1 || { echo "CMAKEFAIL $N"; continue; }
  nice -n 10 make -C $B -k -j16 VERBOSE=1 > $B.make.log 2>&1; echo "make n$N rc=$?"
  grep -c -- '-ffp-contract=off' $B.make.log | sed "s/^/compile lines with -ffp-contract=off: /"
  grep -E 'error|warning' $B.make.log | grep -E 'rad_m1_clean|rad_m1_rhse|rad_m1_vet|rad_m1_vgd' | sort | uniq -c > $B.clean_diag.txt
  [ -x $B/src/athena ] && md5sum $B/src/athena
done
