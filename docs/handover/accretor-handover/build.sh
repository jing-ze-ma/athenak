#!/bin/bash
# usage: build.sh cpu|gpu72
SRC=/viper/ptmp2/jinma/wt_accretor
B=/viper/ptmp2/jinma/accretor_1006/build_$1
source /etc/profile.d/modules.sh
module purge >/dev/null 2>&1
if [ "$1" = gpu72 ]; then
  module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
  CFG="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/7.2.4/bin/hipcc -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF"
else
  module load gcc/14 openmpi/5.0 cmake/4.0
  CFG=""
fi
cmake -S $SRC -B $B $CFG -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=ry_per_accretor > /viper/ptmp2/jinma/accretor_1006/logs/cmake_$1.log 2>&1 || { echo CMAKE_FAIL; exit 1; }
nice -n 10 make -C $B -j8 > /viper/ptmp2/jinma/accretor_1006/logs/make_$1.log 2>&1 || { echo MAKE_FAIL; tail -30 /viper/ptmp2/jinma/accretor_1006/logs/make_$1.log; exit 1; }
md5sum $B/src/athena
echo BUILD_OK $1
