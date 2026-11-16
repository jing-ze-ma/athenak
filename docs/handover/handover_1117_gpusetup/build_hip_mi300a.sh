#!/bin/bash
# MI300A build of clean-cc-push17 (worktree wt) in the persistent tree build_hip: configured
# once with the push-16 recipe, then incremental make; binary -> bin/athena_p17_hip
P=/viper/ptmp2/jinma/ccpush17; B=$P/build_hip; L=$P/build_hip.log
source /etc/profile.d/modules.sh; module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
T0=$(date +%s)
if [ ! -f $B/Makefile ]; then
  cmake -S $P/wt -B $B -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON \
    -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
    -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF -D CMAKE_CXX_COMPILER=$(which hipcc) \
    -D CMAKE_CXX_FLAGS=-ffp-contract=off > $P/build_hip.cmake.log 2>&1 || { echo CMAKE_FAIL; exit 1; }
fi
make -C $B -j${NJ:-96} athena > $L 2>&1 && { cp -p $B/src/athena $P/bin/athena_p17_hip; echo "BUILD_OK $(( $(date +%s)-T0 )) s $(git -C $P/wt rev-parse --short HEAD) $(md5sum $P/bin/athena_p17_hip)"; } || { echo BUILD_FAIL; grep -m8 -E 'error' $L; exit 1; }
