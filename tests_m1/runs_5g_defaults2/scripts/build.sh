#!/bin/bash -l
# usage: build.sh <base|new> <commit> <box|none> <cpu|gpu>
W=/viper/ptmp2/jinma/m1def2_0924; tag=$1; c=$2; p=$3; d=$4
S=$W/src_$tag
if [ ! -d $S ]; then mkdir -p $S; git -C /viper/ptmp2/jinma/wt_m1def2 archive $c | tar -x -C $S
  rm -rf $S/kokkos; ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos; fi
module purge
if [ $d = gpu ]; then module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942_APU=ON -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_HIP_ARCHITECTURES=gfx942"
else module load gcc/14 openmpi/5.0 cmake/4.0; X=""; fi
P=""; [ $p = box ] && P="-DPROBLEM=box_convection"
B=$W/b_${tag}_${p}_$d
cmake -S $S -B $B $X -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 12 > $B.make.log 2>&1 && cp $B/src/athena $W/bin/athena_${tag}_${p}_$d
echo "$tag $p $d exit $? $(md5sum $W/bin/athena_${tag}_${p}_$d 2>/dev/null)"
