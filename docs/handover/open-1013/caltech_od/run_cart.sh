#!/bin/bash
# open-od-1013: the full Cartesian IDORT set with vet_path = cells (+ 5 vet_path = none controls).  bash run_cart.sh <athena CPU binary (plain, built-in pgens)> [threads per run=8] [parallel=8]
HERE=$(cd $(dirname $0) && pwd); B=$(readlink -f ${1:?}); NT=${2:-8}; P=${3:-8}
one() { cd $HERE/runs/$1; rm -rf bin tab rc.txt; OMP_NUM_THREADS=$NT OMP_PROC_BIND=false /usr/bin/time -f "%e s" -o time.txt $B -i in.athinput -d $PWD $(cat args.txt) > run.log 2>&1; echo $? > rc.txt; echo "$1 rc=$(cat rc.txt) $(cat time.txt)"; }
export -f one; export HERE B NT
md5sum $B; xargs -P $P -I{} bash -c 'one {}' < $HERE/list.txt; echo "done $(date)"
