#!/bin/bash -l
# usage: build.sh <tag> <git-rev> <cpu|gpu> [repo]   snapshot + build, binary -> bin/athena_<tag>_<cpu|gpu>
set -e
W=/viper/ptmp2/jinma/h2div_0924; T=$1; R=$2; D=$3; G=${4:-/viper/ptmp2/jinma/wt_h2div}
S=$W/src_$T; P=${5:-none}; B=$W/b_${T}_${P}_$D
if [ ! -d $S ]; then mkdir -p $S; git -C $G archive $R | tar -x -C $S; rmdir $S/kokkos 2>/dev/null || true; ln -sfn /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos; fi
PP=""; [ "$P" = box ] && PP="-D PROBLEM=box_convection"
module purge
if [ "$D" = gpu ]; then
  module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=hipcc -D CMAKE_HIP_ARCHITECTURES=gfx942"
else
  module load gcc/14 openmpi/5.0 cmake/4.0; X=""
fi
cmake -S $S -B $B $X -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release $PP > $B.cmake.log 2>&1
nice cmake --build $B -j 24 > $B.make.log 2>&1
cp $B/src/athena $W/bin/athena_${T}_${P}_$D
md5sum $W/bin/athena_${T}_${P}_$D
