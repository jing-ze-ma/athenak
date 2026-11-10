#!/bin/bash
# waits for a free test GPU (mkdir lock, owner fehgpu, <100 MiB used), runs queue_gpu.sh there, releases the lock
FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj
declare -A J=( [hpc-sm-01-11]=5274196 [hpc-sm-02-03]=5544871 )
while true; do
  for c in hpc-sm-02-03:3 hpc-sm-02-03:1 hpc-sm-01-11:3 hpc-sm-01-11:1; do
    n=${c%:*}; g=${c#*:}; d=$L/lock_${n}_$g
    [ -e $d ] && continue
    mkdir $d 2>/dev/null || continue
    echo "fehgpu queue $$ $(date -Iseconds)" > $d/owner
    mem=$(timeout 60 srun --overlap --jobid=${J[$n]} -N1 -n1 -c1 nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $g 2>/dev/null | tr -d ' ')
    if [ -n "$mem" ] && [ "$mem" -lt 100 ]; then
      echo "$(date) got $n GPU $g (mem $mem MiB)"
      srun --overlap --jobid=${J[$n]} -N1 -n1 -c8 bash -c "export CUDA_VISIBLE_DEVICES=$g; $FF/queue_gpu.sh"
      rm -rf $d; echo "$(date) released $n GPU $g"; exit 0
    fi
    rm -rf $d
  done
  sleep 60
done
