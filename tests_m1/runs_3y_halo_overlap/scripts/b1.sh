#!/bin/bash -l
# usage: b1.sh <ref|new> <boxcpu|boxgpu|dhjcpu|dhjgpu|rwcpu> <jobs>
W=/viper/ptmp2/jinma/m1int_0924; s=$1; t=$2; j=$3
S=$W/$s; B=$W/$s/b_$t
module purge
case $t in *gpu) module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
  X="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=hipcc";;
  *) module load gcc/14 openmpi/5.0 cmake/4.0; X="";; esac
case $t in box*) P="-D PROBLEM=box_convection"; M=ON;; dhj*) P="-D PROBLEM=deep_hot_jupiter_rt"; M=ON;;
  rw*) P=""; M=OFF;; esac
cmake -S $S -B $B $X -D Athena_ENABLE_MPI=$M -D CMAKE_BUILD_TYPE=Release $P > $B.cmake.log 2>&1
nice -n 10 cmake --build $B -j $j > $B.make.log 2>&1
rc=$?; echo "$s $t rc=$rc"
[ $rc = 0 ] && cp $B/src/athena $W/bin/athena_${s}_$t && md5sum $W/bin/athena_${s}_$t
