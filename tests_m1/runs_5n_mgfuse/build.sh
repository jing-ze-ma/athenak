#!/bin/bash -l
# usage: [REV=<commit>] build.sh <tag> <problem|none> <cpu|gpu>: git archive REV (default
# m1-mgfuse HEAD) into src_<tag> (kokkos symlinked), build b_<tag>_..., copy to bin/
set -e
W=/viper/ptmp2/jinma/mgf_0925; S=$W/src_$1; B=$W/b_$1_$2_$3; R=${REV:-HEAD}
if [ ! -d $S ]; then
  mkdir -p $S
  git -C /viper/ptmp2/jinma/wt_mgfuse archive $R | tar -x -C $S
  rm -rf $S/kokkos; ln -s /viper/ptmp2/jinma/wt_mgfuse/kokkos $S/kokkos
fi
module purge
if [ "$3" = gpu ]; then
  module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=hipcc"
else
  module load gcc/14 openmpi/5.0 cmake/4.0; X=""
fi
P=""; [ "$2" != none ] && P="-D PROBLEM=$2"
cmake -S $S -B $B $X -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 32 > $B.make.log 2>&1
cp $B/src/athena $W/bin/athena_$1_$2_$3
md5sum $W/bin/athena_$1_$2_$3
