#!/bin/bash -l
# usage: build.sh <cpu|gpu> <tag> <git-rev> [problem]
# snapshot of <git-rev> of the m1-wedge worktree into src/<tag>, build, copy to bin/
set -e
W=/viper/ptmp2/jinma/sprhd_0926; T=$2; REV=$3; P=${4:-none}; S=$W/src/$T; B=$W/b_${T}_${P}_$1
if [ ! -d $S ]; then
  mkdir -p $S; git -C $W/wt archive $REV | tar -x -C $S
  rm -rf $S/kokkos; ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
  git -C $W/wt rev-parse $REV > $S/REV
fi
module purge
if [ "$1" = gpu ]; then
  module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=hipcc"
else
  module load gcc/14 openmpi/5.0 cmake/4.0; X=""
fi
PP=""; [ "$P" != none ] && PP="-D PROBLEM=$P"
cmake -S $S -B $B $X -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release $PP > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 24 > $B.make.log 2>&1
cp $B/src/athena $W/bin/athena_${T}_${P}_$1; md5sum $W/bin/athena_${T}_${P}_$1
