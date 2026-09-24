#!/bin/bash -l
# usage: build.sh <tag> <commit> <kinds...>   kinds: box dhj mpi
W=/viper/ptmp2/jinma/h2fast_0924; tag=$1; c=$2; shift 2
module purge; module load gcc/14 openmpi/5.0 cmake/4.0
S=$W/src_$tag; rm -rf $S; mkdir -p $S
git -C /viper/ptmp2/jinma/wt_h2fast archive $c | tar -x -C $S
rm -rf $S/kokkos; ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
for k in "$@"; do
  case $k in box) P="-D PROBLEM=box_convection";; dhj) P="-D PROBLEM=deep_hot_jupiter_rt";; mpi) P="";; esac
  B=$S/b_$k
  cmake -S $S -B $B -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
  nice -n 10 cmake --build $B -j 12 > $B.make.log 2>&1 && cp $B/src/athena $W/bin/athena_${tag}_$k
  echo "$tag $k exit $? $(md5sum $W/bin/athena_${tag}_$k 2>/dev/null)"
done
