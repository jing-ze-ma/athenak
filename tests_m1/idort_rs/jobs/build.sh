#!/bin/bash -l
#SBATCH -J rsbld
#SBATCH -p interactive
#SBATCH -N 1
#SBATCH -c 36
#SBATCH --mem=100000
#SBATCH --time 00:50:00
#SBATCH -o /raven/ptmp/jinma/idortrs_1017/logs/build.%j.out
# fin5-1014: fin4_1014 / mem_1014 build recipe (same cmake flags as athena_he_a100_5c504e49) on git-archive snapshots src/<H>
# (kokkos -> wt-prodfinal/kokkos). WHAT items: gpu:<H> (he a100), nccl:<H> (he a100 -DM1_VGD_NCCL), cpu:<H> (plain cpu), cpuh:<H> (he cpu)
# binaries bin/athena_he_a100_<H>, athena_he_a100_nccl_<H>, athena_cpu_<H>, athena_he_cpu_<H>
W=/raven/ptmp/jinma/idortrs_1017; NJ=${NJ:-16}
NH=/mpcdf/soft/SLE_15/packages/skylake/nccl/cuda_12.6-12.6.1/2.23.4
for it in ${WHAT:?}; do t=${it%%:*}; H=${it#*:}; S=$W/src/$H; (
  case $t in
  gpu|nccl) module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake/3.30; export NVCC_WRAPPER_DEFAULT_COMPILER=g++
    XF=""; O=$W/bin/athena_he_a100_$H
    [ $t = nccl ] && { XF="-D CMAKE_CXX_FLAGS=-DM1_VGD_NCCL\ -I$NH/include -D CMAKE_EXE_LINKER_FLAGS=-L$NH/lib\ -lnccl\ -Wl,-rpath,$NH/lib"; O=$W/bin/athena_he_a100_nccl_$H; }
    B=$W/build/${t}_$H; [ -f $B/Makefile ] || eval cmake -S $S -B $B -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On -D Kokkos_ENABLE_CUDA_LAMBDA=ON -D CMAKE_CXX_COMPILER=$S/kokkos/bin/nvcc_wrapper -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=he_star_m1 $XF > $W/logs/cmake_${t}_$H.log 2>&1;;
  cpu|cpuh) module purge; module load gcc/13 openmpi/5.0 cmake/3.30; XP=""; O=$W/bin/athena_cpu_$H
    [ $t = cpuh ] && { XP="-D PROBLEM=he_star_m1"; O=$W/bin/athena_he_cpu_$H; }
    B=$W/build/${t}_$H; [ -f $B/Makefile ] || cmake -S $S -B $B -D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx $XP > $W/logs/cmake_${t}_$H.log 2>&1;;
  esac
  make -C $B -j $NJ > $W/logs/make_${t}_$H.log 2>&1 && cp -p $B/src/athena $O && md5sum $O || { echo FAIL $it; grep -E "error" $W/logs/make_${t}_$H.log | head; } ) &
done
wait
echo BUILD_DONE $(date)
