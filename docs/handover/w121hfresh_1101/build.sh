#!/bin/bash
cd /viper/ptmp2/jinma/w121hfresh_1101 || exit 1
source /etc/profile.d/modules.sh; module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
HIPCC72=/mpcdf/soft/RHEL_9/packages/x86_64/rocm/7.2.4/bin/hipcc
CFG="-D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On -D CMAKE_CXX_COMPILER=$HIPCC72 -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=deep_hot_jupiter_rt"
T0=$(date +%s)
cmake -S wt -B build $CFG > build_cmake.log 2>&1 || { echo CMAKEFAIL; exit 1; }
nice -n 10 make -C build -j16 VERBOSE=1 > build_make.log 2>&1 || { echo MAKEFAIL; exit 1; }
mkdir -p bin; cp build/src/athena bin/athena_dhj_gpu72_18247755_hfresh
echo "DONE wall=$(( $(date +%s)-T0 ))s md5=$(md5sum bin/athena_dhj_gpu72_18247755_hfresh | cut -c1-32)"
