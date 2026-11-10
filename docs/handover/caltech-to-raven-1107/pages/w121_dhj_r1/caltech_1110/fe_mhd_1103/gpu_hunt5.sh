#!/bin/bash
# gpu_hunt2.sh: every 120 s look for an UNLOCKED, idle (< 5000 MiB, user 11-03 via ssh nvidia-smi) GPU on hpc-sm-01-11 (alloc 5274196).
# Takes its lock, stops the 2nd CPU worker after its current dump, runs a GPU worker (-c4), releases the lock.
# Never takes a locked GPU.  Stops when every queue line is claimed.
O=/resnick/groups/carnegie_poc/jingze/fe_mhd_1103; L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj
while true; do
  bash $O/addm2.sh >> $O/logs/addm2.log
  n=$(grep -c . $O/queue3.txt); c=$(ls $O/claim3 | wc -l); [ $c -ge $n ] && break
  for g in 0 2 1 3; do
    d=$L/lock_hpc-sm-01-11_$g; [ -e $d ] && continue
    mem=$(timeout 30 ssh -n hpc-sm-01-11 nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $g 2>/dev/null | tr -d ' ')
    [[ "$mem" =~ ^[0-9]+$ ]] && [ "$mem" -lt 5000 ] || continue
    mkdir $d 2>/dev/null || continue
    echo "fe-mhd-1103 agent (fast fe_hires on MHD dumps, $O) alloc 5274196 $(date -Iseconds)" > $d/owner
    echo "$(date +%T) GPU $g" >> $O/logs/hunt.log; touch $O/stop_cpu2
    srun --overlap --jobid=5274196 -w hpc-sm-01-11 -N1 -n1 -c4 -J fe_mhd bash -c "bash $O/worker3.sh $g" >> $O/logs/worker_g$g.log 2>&1
    rm -rf $d; continue 2
  done
  echo "$(date +%T) no free GPU" >> $O/logs/hunt_poll.log
  sleep 120
done
echo "hunt done $(date +%T)" >> $O/logs/hunt.log
