#!/bin/bash -l
module purge
module load gcc/14 cmake/4.0
cmake -S "$1" -B "$2" -D CMAKE_BUILD_TYPE=Release -D PROBLEM=box_convection || { echo CMAKE_FAILED; exit 1; }
make -C "$2" -j 12
echo "MAKE_EXIT=$?"
[ -x "$2/src/athena" ] && echo BUILD_DONE_OK || echo BUILD_DONE_FAIL
