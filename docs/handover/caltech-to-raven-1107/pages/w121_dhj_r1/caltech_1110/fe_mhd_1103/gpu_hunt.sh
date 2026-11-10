#!/bin/bash
# gpu_hunt.sh: every 120 s try a free test GPU on hpc-sm-01-11 (alloc 5274196; mkdir lock + < 100 MiB), run a GPU worker on it,
# release the lock; stop when every queue line is claimed.
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj
while true; do
  n=$(grep -c . $O/queue.txt); c=$(ls $O/claim | wc -l); [ $c -ge $n ] && break
  for g in 3 1 0 2; do
    d=$L/lock_hpc-sm-01-11_$g; mkdir $d 2>/dev/null || continue
    echo "fe-mhd-1103 agent (fast fe_hires on MHD dumps, $O) alloc 5274196 $(date -Iseconds)" > $d/owner
    mem=$(timeout 60 srun --overlap --jobid=5274196 -N1 -n1 -c1 --gres=none nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $g 2>/dev/null | tr -d ' ')
    if [ -n "$mem" ] && [ "$mem" -lt 100 ]; then
      echo "$(date +%T) GPU $g" >> $O/logs/hunt.log
      srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c8 -J fe_mhd bash $O/worker2.sh $g >> $O/logs/worker_g$g.log 2>&1
      rm -rf $d; continue 2
    fi
    rm -rf $d
  done
  sleep 120
done; echo "hunt done $(date +%T)" >> $O/logs/hunt.log
