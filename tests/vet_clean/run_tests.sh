#!/bin/bash
# run_tests.sh [--quick] [--exe PATH] [--build] [other run_tests.py options]
#
# The clean VET test suite on a CPU machine.  Needs python3 with numpy, cmake, a C++17
# compiler and MPI (mpicxx).  See tests/vet_clean/README.md.
#   --exe PATH   test this athena binary
#   --build      configure and build a CPU binary in <repo>/build_vettest first (the
#                recipe of docs/vet_clean/BUILD.md); without --exe or --build the script
#                uses <repo>/build_vettest/src/athena if it exists, else builds it
#   --quick      the short subset (25 tests, ~1 min on 48 cores)
#   --star-exe PATH   a -D PROBLEM=star_vet binary, for star_smoke and wb_star (else skipped)
# Exit status 0 if every test passes.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
EXE=""; BUILD=0; ARGS=()
while [ $# -gt 0 ]; do
  case $1 in
    --exe) EXE=$2; shift 2 ;;
    --build) BUILD=1; shift ;;
    *) ARGS+=("$1"); shift ;;
  esac
done
B=$REPO/build_vettest
if [ -z "$EXE" ]; then
  EXE=$B/src/athena
  if [ $BUILD = 1 ] || [ ! -x "$EXE" ]; then
    [ -f "$REPO/kokkos/CMakeLists.txt" ] || git -C "$REPO" submodule update --init kokkos
    cmake -S "$REPO" -B "$B" -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx \
      -D Athena_ENABLE_MPI=ON -D Athena_ENABLE_OPENMP=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off \
      -D Athena_VET_VGD_NMAX=4096 > "$B.cmake.log" 2>&1 || { echo "cmake failed: $B.cmake.log"; exit 2; }
    make -C "$B" -j "$(nproc)" athena > "$B.make.log" 2>&1 || { echo "build failed: $B.make.log"; exit 2; }
  fi
fi
[ -x "$EXE" ] || { echo "no binary $EXE"; exit 2; }
exec python3 "$HERE/run_tests.py" --exe "$EXE" "${ARGS[@]}"
