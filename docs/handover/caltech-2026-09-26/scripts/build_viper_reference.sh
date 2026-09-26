#!/bin/bash -l
# build.sh <tag> <cpu|gpu> <commit>   git archive snapshot of <commit> (from the ck-ncloc worktree)
set -e
P=/viper/ptmp2/jinma/ckcad_0926
A=$1; D=$2; C=$3
S=$P/src_${A}_$D
rm -rf $S; mkdir -p $S; git -C /viper/ptmp2/jinma/wt_ckncloc archive $C | tar -x -C $S
rmdir $S/kokkos 2>/dev/null || true; ln -sfn /viper/u2/jinma/ATHENAK/athenak/kokkos $S/kokkos
git -C /viper/ptmp2/jinma/wt_ckncloc rev-parse $C > $P/COMMIT_$A
module purge
cd $S
if [ $D = gpu ]; then
  module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  cmake -B build -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
    -D CMAKE_CXX_COMPILER=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/6.3.4/bin/hipcc \
    -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=deep_hot_jupiter_rt > $P/log/cmake_${A}_$D.log 2>&1
else
  module load gcc/14 openmpi/5.0 cmake/4.0
  cmake -B build -D CMAKE_CXX_COMPILER=mpicxx -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
    -D PROBLEM=deep_hot_jupiter_rt > $P/log/cmake_${A}_$D.log 2>&1
fi
cd build
nice -n 10 make -j12 > $P/log/make_${A}_$D.log 2>&1
cp src/athena $P/athena_${A}_$D
rm -rf $S
echo BUILD_OK $A $D
