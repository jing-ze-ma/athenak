#!/bin/bash
# profc.sh TAG BIN RST NC KEYSFILE [EXTRA key=val ...]   (orion H200PROF 2026-11-21): H200 profile of one production step.
#   Same calling convention/env as ../armc.sh (HOLDER NODE J OUT required; NR default 4 = 1 node x 4 H200; TL srun minutes).
#   Phase N: all NR ranks under nsys (cuda,nvtx,mpi), N0 -> N0+NC, no final restart; per-rank sqlite -> kname.py (kernel table),
#            gaps.py (GPU idle attribution: MPI / memcpy / launch / host), nsys stats (rank 0).
#   Phase C: NCU=1 (default): picktop.py chooses the top NTOP (default 10) kernels of rank 0 (cycles 2..), and for each one run
#            N0 -> N0+3 with rank 0 under "ncu --set full" on ONE launch of that kernel in cycle 2 (other ranks plain).
#            Stops after the first ncu failure (e.g. ERR_NVGPUCTRPERM) and records it. NCU binary: env NCU, else ncu in PATH, else
#            next to nsys, else /usr/local/cuda*/bin/ncu.
#   Output: $OUT/$TAG/{info.txt, n/ (run.log, kname_r*.txt, gaps_r*.txt, stats_r0.txt), c/ (k<i>.txt details, ncusum.txt), res.txt}
#   Only text files need to come back (nsys-rep/ncu-rep/sqlite may stay at Caltech; please keep them until 11-23).
TAG=$1; BIN=$2; RST=$3; NC=$4; KF=$5; shift 5; EX="$*"
: ${HOLDER:?} ${NODE:?} ${J:?} ${OUT:?}; NR=${NR:-4}; NTOP=${NTOP:-10}
H=$(dirname "$(readlink -f "$0")"); HP=$(dirname $H)
case $(basename $RST) in rsg.01295.rst) N0=589788;; hegiant.00139.rst) N0=624961;; *) N0=${N0:?set N0};; esac
D=$OUT/$TAG; rm -rf $D; mkdir -p $D/n $D/c; cd $D
K=$(grep -v '^#' $KF | sed "s#@J@#$J#g" | tr '\n' ' ')
NSYS=${NSYS:-$(command -v nsys)}; NCU=${NCU:-$(command -v ncu || ls $(dirname "$NSYS")/ncu /usr/local/cuda*/bin/ncu 2>/dev/null | head -1)}
echo "TAG=$TAG BIN=$BIN md5=$(md5sum < $BIN | cut -c1-32) RST=$RST N0=$N0 NC=$NC HOLDER=$HOLDER NODE=$NODE NR=$NR EX=$EX NSYS=$NSYS ($($NSYS --version 2>&1 | head -1)) NCU=$NCU ($($NCU --version 2>&1 | tail -1))" > info.txt
nvidia-smi --query-gpu=name,memory.total,clocks.max.sm,clocks.max.mem --format=csv,noheader >> info.txt 2>&1
run() { local d=$1 n=$2 w=$3; shift 3; cd $D/$d; t0=$(date +%s)
  OMPI_MCA_io=romio321 srun --overlap --jobid=$HOLDER -w $NODE -N 1 -n $NR --mpi=pmix -t ${TL:-12} -J p_$TAG \
    env KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1 W_NSYS=$NSYS W_NCU=$NCU W_OUT=$D/$d NCU_RE="$NCU_RE" NCU_SKIP="$NCU_SKIP" NCU_OUT="$NCU_OUT" $w $BIN "$@" $K time/nlim=$n \
    job/final_restart=false $EX -d $D/$d < /dev/null > $D/$d/run$RUNID.log 2>&1
  echo "run$RUNID rc=$? wall=$(( $(date +%s)-t0 ))" >> $D/$d/rc.txt; cd $D; }
# ---- phase N
RUNID= run n $((N0+NC)) $H/wrap_nsys.sh -r $RST
for f in $D/n/nsys_r*.nsys-rep; do $NSYS export -t sqlite -f true -o ${f%.nsys-rep}.sqlite $f > /dev/null 2>&1; done
for s in $D/n/nsys_r*.sqlite; do r=$(basename $s .sqlite | sed 's/nsys_//')
  python3 $HP/batch4/kname.py $s 2 99 60 2 > $D/n/kname_$r.txt 2>&1; python3 $H/gaps.py $s 2 99 2 > $D/n/gaps_$r.txt 2>&1; done
for rp in cuda_gpu_kern_sum cuda_gpu_mem_time_sum cuda_gpu_mem_size_sum cuda_api_sum mpi_event_sum osrt_sum; do
  echo "=== $rp"; $NSYS stats -q -r $rp -f csv $D/n/nsys_r0.sqlite 2>&1 | head -45; done > $D/n/stats_r0.txt
# ---- phase C
if [ "${NCU_ON:-1}" = 1 ] && [ -x "$NCU" ]; then
  python3 $H/picktop.py $D/n/nsys_r0.sqlite 2 $NTOP 2 > $D/c/top.txt 2>&1
  i=0; while IFS=$'\t' read -r key skip re; do i=$((i+1)); [ -n "$re" ] || continue
    echo "k$i $key skip=$skip" >> $D/c/rc.txt
    RUNID=_k$i NCU_RE="$re" NCU_SKIP=$skip NCU_OUT=$D/c/k$i run c $((N0+3)) $H/wrap_ncu.sh -r $RST
    if [ -f $D/c/k$i.ncu-rep ]; then
      $NCU --import $D/c/k$i.ncu-rep --page details --print-units base > $D/c/k$i.txt 2>&1
      $NCU --import $D/c/k$i.ncu-rep --page raw --csv --print-units base > $D/c/k$i.csv 2>&1
    else grep -h -i -E "==ERROR==|==WARNING==|ERR_|No kernels" $D/c/run_k$i.log | head -5 >> $D/c/rc.txt
      grep -q -E "ERR_NVGPUCTRPERM|==ERROR==" $D/c/run_k$i.log && { echo "ncu failed, phase C stopped" >> $D/c/rc.txt; break; }
    fi
  done < <(grep -v '^#' $D/c/top.txt)
  python3 $H/ncusum.py $D/c > $D/c/ncusum.txt 2>&1
else echo "phase C skipped (NCU=$NCU NCU_ON=${NCU_ON:-1})" >> $D/c/rc.txt; fi
# ---- summary
{ cat info.txt; echo "CW $(grep -c 'CLEAN WARNING' $D/n/run.log) FATAL $(grep -c FATAL $D/n/run.log)"
  grep -E "main loop time" $D/n/run.log; grep -oE 'elapsed=\S+ cycle=[0-9]+' $D/n/run.log | tail -6
  head -3 $D/n/kname_r0.txt; head -12 $D/n/gaps_r0.txt; cat $D/n/rc.txt $D/c/rc.txt; head -30 $D/c/ncusum.txt; } > $D/res.txt 2>&1
find $D -name "*.rst" -delete; echo DONE >> $D/res.txt
