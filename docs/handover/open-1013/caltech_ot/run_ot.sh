#!/bin/bash
# open-1013 OT: order diagnostics.  bash run_ot.sh <athena CPU binary (plain, built-in pgens, 281dd686)> [parallel=56]
# list.txt: "<run dir> <OpenMP threads>", costliest first.  Each dir holds a self-contained in.athinput (no extra args).
HERE=$(cd $(dirname $0) && pwd); B=$(readlink -f ${1:?}); P=${2:-56}
one() { d=$HERE/runs/$1; cd $d || return; [ -f rc.txt ] && [ "$(cat rc.txt)" = 0 ] && return
  rm -rf tab bin; OMP_NUM_THREADS=$2 OMP_PROC_BIND=false /usr/bin/time -f "%e s" -o time.txt $B -i in.athinput -d $d > run.log 2>&1
  echo $? > rc.txt; echo "$1 rc=$(cat rc.txt) $(cat time.txt)"; }
export -f one; export HERE B
md5sum $B; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
