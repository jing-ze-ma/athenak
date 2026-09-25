#!/bin/bash -l
# usage: [REV=<commit>|WT] build.sh <tag> <problem|none> <cpu|gpu>: git archive REV (default
# m1-fast5-box HEAD, or the working tree if REV=WT) into src_<tag>, build, copy to bin/,
# delete the build dir
set -e
W=/viper/ptmp2/jinma/fast5box_0925; S=$W/src_$1; B=$W/b_$1_$2_$3; R=${REV:-HEAD}; WT=/viper/ptmp2/jinma/wt_fast5box
if [ ! -d $S ]; then
  mkdir -p $S
  if [ "$R" = WT ]; then
    git -C $WT archive HEAD | tar -x -C $S; rm -rf $S/kokkos
    (cd $WT && git diff --name-only HEAD | xargs -r cp --parents -t $S/)
  else
    git -C $WT archive $R | tar -x -C $S; rm -rf $S/kokkos
  fi
  ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
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
rm -rf $B
md5sum $W/bin/athena_$1_$2_$3
