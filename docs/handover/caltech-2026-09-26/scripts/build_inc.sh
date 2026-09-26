#!/bin/bash -l
# build_inc.sh <tag> <cpu|gpu> [commit=HEAD] [problem=deep_hot_jupiter_rt]
# INCREMENTAL build from a git archive snapshot (no dirty tree ever reaches a binary).
# A persistent tree $P/inc_<cpu|gpu>_<problem>/src is kept per device+problem; each call
# exports <commit> to a temp dir and rsync's it in by CONTENT (-c, no -t): unchanged files keep
# their mtime, changed files get a fresh one, so make recompiles only what changed.
# Output: $P/athena_<tag>_<cpu|gpu>, commit in $P/COMMIT_<tag>.
# GPU builds belong on a compute node:
#   sbatch -A carnegie_poc -p expansion -c 32 --mem=64G -t 01:00:00 -o <log> \
#     --wrap "NJ=32 bash build_inc.sh <tag> gpu <commit>"
set -e
R=/resnick/home/jingze/ATHENAK/athenak
P=${BUILD_ROOT:-/resnick/home/jingze/ATHENAK/builds}
A=$1; D=$2; C=${3:-HEAD}; PROB=${4:-deep_hot_jupiter_rt}
I=$P/inc_${D}_$PROB
mkdir -p $P/log $I/src
exec 9>$I/.lock; flock -w 3600 9            # one build per tree at a time
T=$(mktemp -d $P/snap_XXXXXX)
git -C $R archive $C | tar -x -C $T
rm -rf $T/kokkos
rsync -rlc --delete --exclude=/kokkos --exclude=/build $T/ $I/src/
rm -rf $T
ln -sfn $R/kokkos $I/src/kokkos
# a kokkos version change (e.g. the 4.4.00 -> 4.6.02 bump) needs a fresh CMake configure
KC=$(git -C $R/kokkos rev-parse HEAD)
[ "$(cat $I/KOKKOS_COMMIT 2>/dev/null)" = "$KC" ] || rm -rf $I/src/build
echo $KC > $I/KOKKOS_COMMIT
git -C $R rev-parse $C > $P/COMMIT_$A
module purge
cd $I/src
if [ $D = gpu ]; then
  # HPC-X OpenMPI 4.1.7 is CUDA-aware; the spack openmpi modules are not (ompi_info, 2026-09-25)
  module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
  [ -f build/CMakeCache.txt ] || cmake -B build -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On \
    -D CMAKE_CXX_COMPILER=$I/src/kokkos/bin/nvcc_wrapper \
    -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D PROBLEM=$PROB > $P/log/cmake_${A}_$D.log 2>&1
else
  module load gcc/13.2.0-gcc-13.2.0-w55nxkl openmpi/5.0.1-gcc-13.2.0-xpoh5uw
  [ -f build/CMakeCache.txt ] || cmake -B build -D CMAKE_CXX_COMPILER=mpicxx -D Athena_ENABLE_MPI=ON \
    -D CMAKE_BUILD_TYPE=Release -D PROBLEM=$PROB > $P/log/cmake_${A}_$D.log 2>&1
fi
cd build
nice -n 10 make -j ${NJ:-16} > $P/log/make_${A}_$D.log 2>&1
cp src/athena $P/athena_${A}_$D
echo BUILD_OK $A $D $(cat $P/COMMIT_$A) "recompiled: $(grep -c 'Building CXX' $P/log/make_${A}_$D.log)"
