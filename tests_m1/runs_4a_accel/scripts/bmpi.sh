#!/bin/bash
# usage: bmpi.sh <srcdir>   CPU MPI build without PROBLEM (radwave), -j32
source /etc/profile.d/modules.sh
module purge; module load gcc/14 openmpi/5.0
cd $1 && { [ -f build_mpi/CMakeCache.txt ] || cmake -B build_mpi -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release > cmake_mpi.log 2>&1; } && nice -n 10 make -C build_mpi -j32 > make_mpi.log 2>&1; echo "exit $?" >> make_mpi.log
