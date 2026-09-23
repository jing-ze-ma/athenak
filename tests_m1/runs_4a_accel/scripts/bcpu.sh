#!/bin/bash
# usage: bcpu.sh <srcdir>   CPU MPI box_convection build (-j32)
source /etc/profile.d/modules.sh
module purge; module load gcc/14 openmpi/5.0
cd $1 && { [ -f build_boxcpu/CMakeCache.txt ] || cmake -B build_boxcpu -D Athena_ENABLE_MPI=ON -D PROBLEM=box_convection -D CMAKE_BUILD_TYPE=Release > cmake_boxcpu.log 2>&1; } && nice -n 10 make -C build_boxcpu -j32 > make_boxcpu.log 2>&1; echo "exit $?" >> make_boxcpu.log
