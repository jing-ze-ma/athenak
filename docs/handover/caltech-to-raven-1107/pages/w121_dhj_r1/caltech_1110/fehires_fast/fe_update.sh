#!/bin/bash
# fe_update.sh [--stride R] [--max N] [--img DIR] [--wait SEC]
#   incremental fe_hires blueshift series for the dhj page; prints the image paths at the end.
#   1) list the dumps not done yet (newest 5 per planet + one every R rotations, default 10);
#   2) run them on a free TEST GPU (mkdir lock, owner fehgpu, < 100 MiB in use; one dump ~1.7 min,
#      both cuts in one trace); if no GPU is free (after --wait seconds, default 0): CPU one-pass
#      fallback for at most 1 dump (--max applies to the GPU only; CPU ~30-60 min per dump);
#   3) plots: <img>/fe_vs_rot.png, <img>/fehires_vs_seidel.png, <img>/fe_takeaway.txt.
#   Never production GPUs (01-11 GPU 0/2, 02-03 GPU 0/2), never 01-05 / 01-09.
FF=/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast; PY=/resnick/groups/carnegie_poc/jingze/fewind_1102/venv_fehgpu/bin/python
L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj
STRIDE=10; MAX=8; IMG=$FF/series/img; WAIT=0
while [ $# -gt 0 ]; do case $1 in --stride) STRIDE=$2; shift;; --max) MAX=$2; shift;; --img) IMG=$2; shift;; --wait) WAIT=$2; shift;; esac; shift; done
declare -A JOB=( [hpc-sm-01-11]=5274196 [hpc-sm-02-03]=5544871 )
CAND="hpc-sm-01-11:3 hpc-sm-01-11:1 hpc-sm-02-03:3 hpc-sm-02-03:1"
mkdir -p $FF/series/logs; LOG=$FF/series/logs/fe_update_$(date +%m%d_%H%M%S).log
TODO=$FF/series/todo_$$.txt
$PY $FF/fe_series.py todo --stride $STRIDE --max $MAX > $TODO
N=$(grep -c . $TODO)
if [ "$N" -gt 0 ]; then
  GOT=""; T0=$(date +%s)
  while [ -z "$GOT" ]; do
    for c in $CAND; do
      n=${c%:*}; g=${c#*:}; d=$L/lock_${n}_$g
      mkdir $d 2>/dev/null || continue
      echo "fehgpu fe_update $$ $(date -Iseconds)" > $d/owner
      mem=$(timeout 60 srun --overlap --jobid=${JOB[$n]} -N1 -n1 -c1 nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i $g 2>/dev/null | tr -d ' ')
      if [ -n "$mem" ] && [ "$mem" -lt 100 ]; then GOT=$c; break; fi
      rm -rf $d
    done
    [ -n "$GOT" ] || [ $(( $(date +%s) - T0 )) -ge $WAIT ] && break
    sleep 30
  done
  if [ -n "$GOT" ]; then
    n=${GOT%:*}; g=${GOT#*:}; d=$L/lock_${n}_$g
    trap "rm -rf $d" EXIT
    echo "$(date) GPU $n:$g, $N dumps" >> $LOG
    srun --overlap --jobid=${JOB[$n]} -N1 -n1 -c8 bash -c "export CUDA_VISIBLE_DEVICES=$g NPROC=8 NUMBA_NUM_THREADS=8; $PY $FF/fe_series.py run $TODO --gpu" >> $LOG 2>&1
    rm -rf $d; trap - EXIT
  else
    head -n 1 $TODO > $TODO.cpu
    echo "$(date) no free GPU: CPU one-pass for $(cat $TODO.cpu)" >> $LOG
    srun --overlap --jobid=5274196 -N1 -n1 -c16 bash -c "export NPROC=16 NUMBA_NUM_THREADS=16; $PY $FF/fe_series.py run $TODO.cpu" >> $LOG 2>&1
    rm -f $TODO.cpu
  fi
fi
rm -f $TODO
$PY $FF/fe_series.py plot --img $IMG 2>> $LOG
