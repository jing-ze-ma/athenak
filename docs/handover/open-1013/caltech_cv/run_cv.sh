#!/bin/bash
# open-ot-1013 CV: FIN3 gas-feedback IDORT set with implicit_cs2_vimp = false.
# bash run_cv.sh <athena CPU binary built from 281dd686 + patch/0001-*.patch> [parallel=16]
# list.txt: "<run dir> <OpenMP threads>"; each dir: in.athinput + args.txt (command-line keys); the new key is IN the deck (a command-line key must already exist in the deck)
HERE=$(cd $(dirname $0) && pwd); B=$(readlink -f ${1:?}); P=${2:-16}
one() { d=$HERE/runs/$1; cd $d || return; [ -f rc.txt ] && [ "$(cat rc.txt)" = 0 ] && return
  rm -rf tab bin; OMP_NUM_THREADS=$2 OMP_PROC_BIND=false /usr/bin/time -f "%e s" -o time.txt $B -i in.athinput -d $d $(cat args.txt) > run.log 2>&1
  echo $? > rc.txt; echo "$1 rc=$(cat rc.txt) $(cat time.txt)"; }
export -f one; export HERE B
md5sum $B; xargs -P $P -L 1 bash -c 'one "$@"' _ < $HERE/list.txt; echo "done $(date)"
