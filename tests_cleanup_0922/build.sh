#!/bin/bash -l
# build.sh <pgen> : configure+build /viper/u2/jinma/ATHENAK/athenak into build_clean_<pgen>
module purge >/dev/null 2>&1; module load gcc/14 cmake/4.0 >/dev/null 2>&1
A=/viper/u2/jinma/ATHENAK/athenak
P=$1; B=$A/build_clean_$P
if [ "$P" = "base" ]; then POPT=""; else POPT="-D PROBLEM=$P"; fi
cmake -S $A -B $B -D CMAKE_BUILD_TYPE=Release $POPT > $B.cmake.log 2>&1 || { echo "CMAKE_FAILED $P"; exit 1; }
make -C $B -j 8 > $B.make.log 2>&1
[ -x $B/src/athena ] && echo "BUILD_OK $P" || { echo "BUILD_FAIL $P"; tail -30 $B.make.log; }
