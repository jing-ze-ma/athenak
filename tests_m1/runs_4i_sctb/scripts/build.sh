#!/bin/bash -l
# usage: build.sh <srcdir> <builddir> gpu|cpu
set -e
S=$1; B=$2
module purge
if [ $3 = gpu ]; then
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cmake -S $S -B $B -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
  -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=box_convection \
  -D CMAKE_CXX_COMPILER=hipcc > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 12 > $B.make.log 2>&1
else
module load gcc/14 openmpi/5.0 cmake/4.0
cmake -S $S -B $B -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=box_convection > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 12 > $B.make.log 2>&1
fi
md5sum $B/src/athena
