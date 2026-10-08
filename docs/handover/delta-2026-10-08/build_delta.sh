#!/bin/bash
# Incremental AthenaK build on NCSA Delta (A100, CUDA) -- see docs/handover/NOTE-2026-09-28-incremental-builds.md
# usage: bash build_delta.sh <target> <commit>       targets: he_gpu (PROBLEM=he_star_m1, CUDA+MPI)
# One persistent worktree + build dir per target; never builds a dirty tree; binary -> bin/athena_<target>_<sha>.
set -eo pipefail
B=/work/nvme/bivj/jma20/delta_1008
REPO=/u/jma20/ATHENAK/athenak
T=$1; C=$2
[ -n "$T" ] && [ -n "$C" ] || { echo "usage: $0 <target> <commit>"; exit 1; }

# --- environment (Delta Cray PE, 2026-10): default PrgEnv-gnu (gcc-native/14), cray-mpich/9.1.0,
# craype-accel-nvidia80 (GPU-aware MPI via GTL). CUDA 12.9 instead of the default cudatoolkit 13.2:
# Kokkos 4.6.2 (in-tree) does not compile with CUDA 13 (cudaMemAdvise(int device) removed).
module unload cudatoolkit
module load cuda/12.9
export NVCC_WRAPPER_DEFAULT_COMPILER=CC     # Cray wrapper -> brings in cray-mpich + GTL
module -t list 2>&1 | tr '\n' ' '; echo

WT=$B/wt_$T; BD=$B/build_$T
exec 9>$B/.lock_$T; flock -n 9 || { echo "target $T is being built by someone else"; exit 1; }
if [ ! -d $WT ]; then
  git -C $REPO worktree add --detach $WT $C
fi
git -C $WT fetch -q origin 2>/dev/null || true
git -C $WT checkout -q --detach $C
git -C $WT reset -q --hard
git -C $WT submodule update --init --recursive --depth 1
[ -z "$(git -C $WT status --porcelain)" ] || { echo "dirty tree in $WT"; exit 1; }
SHA=$(git -C $WT rev-parse --short=8 HEAD)
OUT=$B/bin/athena_${T}_${SHA}
mkdir -p $B/bin
[ -x $OUT ] && { echo "exists: $OUT $(md5sum < $OUT | cut -d' ' -f1)"; exit 0; }

case $T in
  he_gpu) OPTS="-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On
                -D Kokkos_ARCH_AMPERE80=On -D Kokkos_ARCH_ZEN3=On
                -D CMAKE_CXX_COMPILER=$WT/kokkos/bin/nvcc_wrapper" ;;
  *) echo "unknown target $T"; exit 1 ;;
esac
[ -f $BD/CMakeCache.txt ] || cmake -S $WT -B $BD $OPTS
t0=$(date +%s)
nice make -C $BD -j ${NJ:-32}
cp $BD/src/athena $OUT
M=$(md5sum < $OUT | cut -d' ' -f1)
echo "$(date '+%F %T') $T $SHA $M $(( $(date +%s)-t0 ))s" | tee -a $B/bin/BUILD_LOG
