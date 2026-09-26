#!/bin/bash -l
# build_caltech.sh <tag> <cpu|gpu> [commit=HEAD] [problem=deep_hot_jupiter_rt]
# Caltech (Resnick HPC) build from a git archive snapshot, so a dirty tree never leaks into a binary.
# GPU target: H100/H200 (HOPPER90). Output: $P/athena_<tag>_<cpu|gpu>, commit in $P/COMMIT_<tag>.
set -e
R=/resnick/home/jingze/ATHENAK/athenak
P=${BUILD_ROOT:-/resnick/home/jingze/ATHENAK/builds}
A=$1; D=$2; C=${3:-HEAD}; PROB=${4:-deep_hot_jupiter_rt}
S=$P/src_${A}_$D
mkdir -p $P/log
rm -rf $S; mkdir -p $S; git -C $R archive $C | tar -x -C $S
rm -rf $S/kokkos; ln -sfn $R/kokkos $S/kokkos
git -C $R rev-parse $C > $P/COMMIT_$A
module purge
cd $S
if [ $D = gpu ]; then
  # HPC-X OpenMPI 4.1.7 is CUDA-aware; the spack openmpi modules are not (ompi_info, 2026-09-25)
  module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
  cmake -B build -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On \
    -D CMAKE_CXX_COMPILER=$S/kokkos/bin/nvcc_wrapper \
    -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=$PROB > $P/log/cmake_${A}_$D.log 2>&1
else
  module load gcc/13.2.0-gcc-13.2.0-w55nxkl openmpi/5.0.1-gcc-13.2.0-xpoh5uw
  cmake -B build -D CMAKE_CXX_COMPILER=mpicxx -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
    -D PROBLEM=$PROB > $P/log/cmake_${A}_$D.log 2>&1
fi
cd build
nice -n 10 make -j ${NJ:-16} > $P/log/make_${A}_$D.log 2>&1
cp src/athena $P/athena_${A}_$D
rm -rf $S
echo BUILD_OK $A $D $(cat $P/COMMIT_$A)
