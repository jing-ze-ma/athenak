#!/bin/bash
# run_list.sh <taskfile> <nparallel>
export OMP_NUM_THREADS=1
cat $1 | xargs -P $2 -L 1 /resnick/groups/carnegie_poc/jingze/lowbeta_1102/run_one.sh
echo ALLDONE $1 >> /resnick/groups/carnegie_poc/jingze/lowbeta_1102/done.txt
