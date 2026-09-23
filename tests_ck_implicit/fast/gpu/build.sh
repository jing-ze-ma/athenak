#!/bin/bash -l
# usage: build.sh <commit> <tag> <cpu|gpu> <jobs> [patch]
set -e
P=/viper/ptmp2/jinma/ckfast_0923
C=$1; T=$2; K=$3; J=$4; PAT=$5
S=$P/src_$T
module purge
if [ $K = gpu ]; then module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0; else module load gcc/14 openmpi/5.0 cmake/4.0; fi
if [ ! -d $S ]; then
  mkdir -p $S; git -C /viper/ptmp2/jinma/wt_ckfast archive $C | tar -x -C $S
  rmdir $S/kokkos 2>/dev/null || true; ln -sfn /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
  [ -n "$PAT" ] && (cd $S && patch -p1 < $PAT > /dev/null)
fi
cd $S
if [ $K = gpu ]; then
  cmake -B build_gpu -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
    -D CMAKE_CXX_COMPILER=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4/bin/hipcc \
    -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=deep_hot_jupiter_rt > $P/cmake.$T.$K.log 2>&1
  cd build_gpu
else
  cmake -B build_cpu -D CMAKE_CXX_COMPILER=mpicxx -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
    -D PROBLEM=deep_hot_jupiter_rt > $P/cmake.$T.$K.log 2>&1
  cd build_cpu
fi
nice -n 10 make -j$J > $P/make.$T.$K.log 2>&1
cp src/athena $P/athena.$K.$T
echo BUILD_OK $T $K
