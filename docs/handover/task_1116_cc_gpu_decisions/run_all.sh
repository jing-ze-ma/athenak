#!/bin/bash
# run_all.sh <athena> [group ...]   run the runs of runs.txt (default G1-G4, or the named groups)
#   one run per line: runs/<tag>/ with in.athinput = decks/<tag>.athinput, nranks = GPUs.
#   Fresh starts only: an existing runs/<tag>/ is moved to runs/<tag>.old.<date>.
#   LAUNCH (env) is the launcher, %N = ranks; default for 1 rank per GPU under Slurm:
#     LAUNCH='srun -n %N --gpus-per-task=1'   (add your usual MI300A binding options)
#   PAR (env, default 1) = runs at the same time (keep 1 for the timing group G5).
#   Stops at the first run with rc != 0, FATAL or NaN in its log unless KEEPGOING=1.
#   TIMEOUT (env, e.g. 30m) caps each run; NLIM (env) overrides <time>/nlim (G5 memory runs);
#   KOKKOS_TOOLS_LIBS is passed through (kernel timer / space-time-stack write into the run dir).
#   G5, G5L, G6 and BIT run only when named.
H=$(cd "$(dirname "$0")" && pwd); B=$(readlink -f $1); shift
GROUPS_SEL="$*"
L=${LAUNCH:-'srun -n %N --gpus-per-task=1'}
PAR=${PAR:-1}
R=${RUNS:-$H/runs}   # run root (env RUNS)
mkdir -p $R
md5sum $B > $R/binary_md5_$(basename $B).txt
one() {
  local tag=$1 np=$2 d=$R/$1
  [ -d $d ] && mv $d $d.old.$(date +%m%d%H%M%S)
  mkdir -p $d && cd $d || return 1
  cp $H/decks/$tag.athinput in.athinput
  [ -n "$NLIM" ] && sed -i "s/^nlim *=.*/nlim = $NLIM/" in.athinput
  local cmd=${L//%N/$np}
  [ -n "$TIMEOUT" ] && cmd="timeout $TIMEOUT $cmd"
  echo "$cmd $B -i in.athinput -d ." > cmd.txt
  local t0=$(date +%s.%N)
  $cmd $B -i in.athinput -d . < /dev/null > run.log 2>&1
  local rc=$?
  echo "wall $(awk -v a=$t0 -v b=$(date +%s.%N) 'BEGIN{printf "%.2f", b-a}')" >> run.log
  echo "rc=$rc" >> run.log
  local bad=$(grep -c -E '### FATAL|(^|[^A-Za-z_])-?(nan|NaN)([^A-Za-z_]|$)|DIVERGED' run.log)
  echo "$tag rc=$rc fatal/nan=$bad warn=$(grep -c 'CLEAN WARNING' run.log)"
  [ $rc -eq 0 ] && [ $bad -eq 0 ]
}
grep -v '^#' $H/runs.txt | { while read -r g tag np; do
  [ -z "$tag" ] && continue
  [ -z "$GROUPS_SEL" ] && case $g in BIT|G5|G5L|G6) continue;; esac
  [ -n "$GROUPS_SEL" ] && ! echo " $GROUPS_SEL " | grep -q " $g " && continue
  if [ $PAR -le 1 ]; then
    one $tag $np || { [ "$KEEPGOING" = 1 ] || { echo "STOP at $tag"; exit 1; }; }
  else
    while [ $(jobs -rp | wc -l) -ge $PAR ]; do sleep 2; done
    one $tag $np &
  fi
done
wait; }
