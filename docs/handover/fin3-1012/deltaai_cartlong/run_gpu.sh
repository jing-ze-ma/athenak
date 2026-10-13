#!/bin/bash
# fin3-1012 final: 4 long Cartesian IDORT runs, one per GH200 GPU, on ONE node (4 GPUs).
# bash run_gpu.sh <athena plain CUDA binary>   (inside the job allocation; export MPICH_GPU_SUPPORT_ENABLED=1)
HERE=$(cd $(dirname $0) && pwd); B=$(readlink -f ${1:?}); md5sum $B
i=0
for d in $(cat $HERE/list.txt); do
  ( cd $HERE/runs/$d; rm -rf bin tab rc.txt
    CUDA_VISIBLE_DEVICES=$i OMP_NUM_THREADS=1 /usr/bin/time -f "%e s" -o time.txt $B -i in.athinput -d $PWD $(cat args.txt) > run.log 2>&1
    echo $? > rc.txt; echo "$d rc=$(cat rc.txt) $(cat time.txt)" ) &
  i=$((i+1))
done; wait; echo "done $(date)"
