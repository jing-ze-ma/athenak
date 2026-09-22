#!/bin/bash
set -u
ROOT=/viper/u2/jinma/ATHENAK/athenak
HERE=$ROOT/tests_cs_vertex_cell
ATH=${ATH:-$HERE/athena_ref}
TAG=${TAG:-ref}
EXTRA=${EXTRA:-}
NS=${NS:-"16 32 64"}
TS=${TS:-loop}
for t in $TS; do
  case $t in
    loop) inp=$ROOT/inputs/tests/cs_regions_loop.athinput; ce=0; tl=${TL:-0.25} ;;
    strat) inp=$ROOT/inputs/tests/cs_regions_strat.athinput; ce=1; tl=${TL:-0.04} ;;
  esac
  for n in $NS; do
    nb=$(( 4 * n / 32 )); [ $nb -lt 1 ] && nb=1
    log=$HERE/logs/${TAG}_${t}_n${n}_t${tl}.log
    [ -s "$log" ] && grep -q REGION "$log" && { echo "skip $log"; continue; }
    ( cd $HERE/out && /usr/bin/time -f "WALLCLOCK %e s" $ATH -i $inp -d . \
        job/basename=${TAG}_${t}_${n} \
        mesh/nghost=3 mesh/nx1=$n mesh/nx2=$n mesh/nx3=$n \
        meshblock/nx1=$n meshblock/nx2=$n meshblock/nx3=$n \
        time/tlim=$tl time/nlim=-1 mhd/reconstruct=plm \
        mhd/cs_wellbalanced_src=false \
        problem/conv_errors=$ce problem/conv_nband=$nb $EXTRA ) > $log 2>&1 \
      || echo "FAILED $log"
  done
done
