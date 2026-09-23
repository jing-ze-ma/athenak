#!/bin/bash
# usage: bgpu.sh <srcdir>   GPU MPI box_convection build (-j48), log in <srcdir>/make_build_gpu.log
source /etc/profile.d/modules.sh
module purge; module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cd $1
[ -f build_gpu/CMakeCache.txt ] || cmake -B build_gpu -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_AMD_GFX942_APU=ON \
    -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_HIP_ARCHITECTURES=gfx942 -DCMAKE_BUILD_TYPE=Release -DPROBLEM=box_convection > cmake_gpu.log 2>&1
nice -n 10 make -C build_gpu -j48 > make_build_gpu.log 2>&1; echo "exit $?" >> make_build_gpu.log
