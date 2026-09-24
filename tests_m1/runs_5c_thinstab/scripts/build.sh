#!/bin/bash -l
# usage: build.sh <srcdir> <tag> <cpu|gpu>   (problem = none, built-in pgens)
set -e
W=/viper/ptmp2/jinma/thinstab_0924; S=$1; B=$W/b_$2_$3
module purge
if [ "$3" = gpu ]; then
  module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=hipcc"
else
  module load gcc/14 openmpi/5.0 cmake/4.0; X=""
fi
cmake -S $S -B $B $X -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release > $B.cmake.log 2>&1
cmake --build $B -j 32 > $B.make.log 2>&1
cp $B/src/athena $W/bin/athena_$2_$3
md5sum $W/bin/athena_$2_$3
