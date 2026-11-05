#!/bin/bash -l
W=/raven/ptmp/jinma/athpp_audit_1015/h2h; H=9a68a5e3; S=$W/src/$H
module purge; module load gcc/13 openmpi/5.0 cmake/3.30
B=$W/build/cpuh_$H; [ -f $B/Makefile ] || cmake -S $S -B $B -D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D PROBLEM=he_star_m1 > $W/logs/cmake_$H.log 2>&1
make -C $B -j 36 > $W/logs/make_$H.log 2>&1 && cp -p $B/src/athena $W/bin/athena_he_cpu_$H && md5sum $W/bin/athena_he_cpu_$H || { echo FAIL; grep error $W/logs/make_$H.log | head; }
echo BUILD_DONE
