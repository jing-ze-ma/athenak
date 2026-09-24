#!/bin/bash -l
# usage: build.sh <box|none> <cpu|gpu>   (new = snapshot of the wt_sph2 working tree)
W=/viper/ptmp2/jinma/sph2_0924; p=$1; d=$2
S=$W/src_new
module purge
if [ $d = gpu ]; then module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942_APU=ON -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_HIP_ARCHITECTURES=gfx942"
else module load gcc/14 openmpi/5.0 cmake/4.0; X=""; fi
P=""; [ $p = box ] && P="-DPROBLEM=box_convection"
B=$W/b_new_${p}_$d
[ -f $B/CMakeCache.txt ] || cmake -S $S -B $B $X -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j 16 > $B.make.log 2>&1 && cp $B/src/athena $W/bin/athena_new_${p}_$d
echo "new $p $d exit $? $(md5sum $W/bin/athena_new_${p}_$d 2>/dev/null)"
