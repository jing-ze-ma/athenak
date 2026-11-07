#!/bin/bash
# parred.sh rsgE|hggl "idx idx ..." : reduce the listed dumps in parallel (one reduce.py per dump), for parpage fan-out (USER 10-22)
J=/resnick/groups/carnegie_poc/jingze; N=$1
case $N in rsgE) PG=$J/rsgE_page;; hggl) PG=$J/hgglob_page;; esac
eval "$(grep -E '^(export|HERE=|R=)' $PG/update.sh)"; export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2; cd $PG/out
echo $2 | tr ' ' '\n' | grep . | xargs -P ${W:-8} -I{} env ONLY={} $PY $SRC/reduce.py $RUN data
