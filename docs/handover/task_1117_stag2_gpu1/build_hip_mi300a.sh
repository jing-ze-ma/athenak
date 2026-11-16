#!/bin/bash
# build_hip_mi300a.sh <athenak worktree at clean-cc-stag2-gpu1> [nmax ...]: MI300A binaries
# Recipe of the stage-1 MI300A TASK (TASK-2026-11-16-viper-cc_gpu1.md): gcc/16 rocm/7.2
# openmpi_gpu/5.0, hipcc, Kokkos GFX942_APU, HIP malloc async off, MPI, Release,
# -ffp-contract=off (clang applies it to host and device code).  Built-in pgens.
WT=$(readlink -f $1); shift; NM=${*:-432}
source /etc/profile.d/modules.sh; module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
CFG="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF
     -D CMAKE_CXX_COMPILER=$(which hipcc) -D Athena_ENABLE_MPI=ON
     -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off"
for N in $NM; do
  B=$WT/build_hip_n$N
  X=""; [ $N != 432 ] && X="-D Athena_M1_VGD_NMAX=$N"
  cmake -S $WT -B $B $CFG $X > $B.cmake.log 2>&1 || { echo "CMAKEFAIL $N"; continue; }
  nice -n 10 make -C $B -k -j16 > $B.make.log 2>&1; echo "make n$N rc=$?"
  grep -E 'error|warning' $B.make.log | grep -E 'rad_m1_clean|rad_m1_rhse|rad_m1_vet|rad_m1_vgd' | sort | uniq -c > $B.clean_diag.txt
  [ -x $B/src/athena ] && md5sum $B/src/athena
done
