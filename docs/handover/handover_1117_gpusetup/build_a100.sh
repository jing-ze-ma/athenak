#!/bin/bash -l
# A100 build of clean-cc-push17 (src synced from viper ccpush17/wt) in the persistent tree
# src/build_a100 (push-16 recipe): configure once, then make; run in a -c 54 overlap step
S=/raven/ptmp/jinma/ccpush17/src; B=$S/build_a100; P=/raven/ptmp/jinma/ccpush17
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake
if [ ! -f $B/Makefile ]; then
  cmake -S $S -B $B -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On \
    -D Kokkos_ARCH_AMPERE80=On -D CMAKE_CXX_COMPILER=$S/kokkos/bin/nvcc_wrapper \
    -D CMAKE_CXX_FLAGS=-ffp-contract=off > $P/cmake.log 2>&1 || { echo CMAKE_FAIL; exit 2; }
fi
make -C $B -j54 athena > $P/build.log 2>&1 || { echo BUILD_FAIL; grep -m5 error $P/build.log; exit 3; }
cp -p $B/src/athena $P/bin/athena_p17_a100; md5sum $P/bin/athena_p17_a100; echo BUILD_DONE
