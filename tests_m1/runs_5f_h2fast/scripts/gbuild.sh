#!/bin/bash -l
# usage: gbuild.sh <tag> <commit>   GPU MPI box_convection build (HIP gfx942)
W=/viper/ptmp2/jinma/h2fast_0924; tag=$1; c=$2
module purge; module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
S=$W/gsrc_$tag; rm -rf $S; mkdir -p $S
git -C /viper/ptmp2/jinma/wt_h2fast archive $c | tar -x -C $S
rm -rf $S/kokkos; ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
cmake -S $S -B $S/b_gpu -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942_APU=ON \
  -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_HIP_ARCHITECTURES=gfx942 -DCMAKE_BUILD_TYPE=Release \
  -DPROBLEM=box_convection > $S/cmake.log 2>&1
nice -n 10 make -C $S/b_gpu -j12 > $S/make.log 2>&1 && cp $S/b_gpu/src/athena $W/bin/athena_${tag}_gpu
echo "$tag gpu exit $? $(md5sum $W/bin/athena_${tag}_gpu 2>/dev/null)"
