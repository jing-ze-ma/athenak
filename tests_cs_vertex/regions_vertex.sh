#!/bin/bash
# GATE (c): the per-region error gate of tests_cs_regions, rerun with the cube-vertex
# fills switched.  Same configuration as the 76359975 baseline: plm,
# cs_wellbalanced_src off, nghost 3, one MeshBlock per panel, conv_nband scaled with n
# (2/4/8 at n = 16/32/64) so the band has a fixed PHYSICAL width, and the per-test
# physical time (loop 0.25, strat 0.04) rather than a fixed cycle count.
#
#   ./regions_vertex.sh <binary> <tag> [extra athinput overrides...]
set -u
B=${1:-./athena_new}
TAG=${2:-new}
shift 2 || true
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1
mkdir -p logs out
for t in loop strat; do
  case $t in
    loop)  tlim=0.25; cerr=0 ;;
    strat) tlim=0.04; cerr=1 ;;
  esac
  for n in 16 32 64; do
    nb=$(( 4 * n / 32 )); [ $nb -lt 1 ] && nb=1
    log="logs/${t}_${TAG}_n${n}.log"
    if [ -s "$log" ] && grep -q REGION "$log"; then echo "skip $log"; continue; fi
    echo "run $log"
    ( cd out && /usr/bin/time -f "WALLCLOCK %e s" "../$B" -i "../in_${t}.athinput" -d . \
        job/basename="${t}_${TAG}_n${n}" \
        mesh/nghost=3 mesh/nx1=$n mesh/nx2=$n mesh/nx3=$n \
        meshblock/nx1=$n meshblock/nx2=$n meshblock/nx3=$n \
        time/tlim=$tlim time/nlim=-1 \
        mhd/reconstruct=plm mhd/cs_wellbalanced_src=false \
        problem/conv_errors=$cerr problem/conv_nband=$nb "$@" ) > "$log" 2>&1
  done
done
