#!/bin/bash
# ONE-STEP (local truncation error) per-region scan.  Both exact solutions are STATIC,
# so the error after nlim steps is nlim*dt*(spatial truncation error) and, since dt ~ h,
# the measured order is (spatial order + 1).  What matters is VERTEX vs INTERIOR.
set -u
ROOT=/viper/u2/jinma/ATHENAK/athenak
HERE=$ROOT/tests_cs_vertex_cell
ATH=${ATH:-$HERE/athena_ref}
TAG=${TAG:-ref}
NLIM=${NLIM:-1}
EXTRA=${EXTRA:-}
for t in loop strat; do
  case $t in
    loop) inp=$ROOT/inputs/tests/cs_regions_loop.athinput; ce=0 ;;
    strat) inp=$ROOT/inputs/tests/cs_regions_strat.athinput; ce=1 ;;
  esac
  for n in 16 32 64; do
    nb=$(( 4 * n / 32 )); [ $nb -lt 1 ] && nb=1
    log=$HERE/logs/${TAG}_${t}_n${n}_s${NLIM}.log
    ( cd $HERE/out && $ATH -i $inp -d . job/basename=${TAG}_${t}_${n} \
        mesh/nghost=3 mesh/nx1=$n mesh/nx2=$n mesh/nx3=$n \
        meshblock/nx1=$n meshblock/nx2=$n meshblock/nx3=$n \
        time/tlim=1e30 time/nlim=$NLIM mhd/reconstruct=plm \
        mhd/cs_wellbalanced_src=false \
        problem/conv_errors=$ce problem/conv_nband=$nb $EXTRA ) > $log 2>&1 \
      || echo "FAILED $log"
  done
done

# The CONTROLS of section 3 (strat with the gravity and/or the field switched off):
#   A: problem/grav=0 problem/b0c=0        uniform gas, no field
#   C: problem/grav=1 problem/b0c=0        hydrostatic atmosphere, no field
#   B: problem/grav=0 problem/b0c=0.3162   uniform gas + uniform Cartesian field
# pass them through EXTRA, e.g.
#   EXTRA="problem/grav=0 problem/b0c=0" TAG=ctlA ./run1.sh
