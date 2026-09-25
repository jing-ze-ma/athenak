#!/bin/bash -l
# usage: build.sh <ref|new> <box|none> <cpu|gpu>
W=/viper/ptmp2/jinma/sporder2b_0925; a=$1; p=$2; d=$3
S=$W/src_$a
module purge
if [ $d = gpu ]; then module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942_APU=ON -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_HIP_ARCHITECTURES=gfx942"
else module load gcc/14 openmpi/5.0 cmake/4.0; X=""; fi
P=""; [ $p = box ] && P="-DPROBLEM=box_convection"
B=$W/b_${a}_${p}_$d
[ -f $B/CMakeCache.txt ] || cmake -S $S -B $B $X -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 16 > $B.make.log 2>&1 && cp $B/src/athena $W/bin/athena_${a}_${p}_$d
echo "$a $p $d exit $? $(md5sum $W/bin/athena_${a}_${p}_$d 2>/dev/null)"
