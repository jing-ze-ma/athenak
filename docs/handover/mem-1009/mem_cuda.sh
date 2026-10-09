#!/bin/bash
# BSG GPU-memory inventory + s/cycle (mem-1009), CUDA port of viper /viper/ptmp2/jinma/mem_1009/inv.sh.
# usage (inside a 1-node allocation with 4 GPUs):
#   bash mem_cuda.sh <BO = 98835d99 binary> <BN = mem-1009 binary> <files dir from SETUP.sh> <out dir> [KT .so]
# KT = kp_memory_events.so (Kokkos Tools memory-events); without it only nvidia-smi peaks are taken.
# arms (in this order; each is independent, an OOM in one does not stop the others):
#   mem_c8  BN, FULL 16-block BSG mesh, 2 ranks x 8 blocks (THE CALTECH LAYOUT), 6 cycles  -> peak/GPU, s/cycle
#   mem_r4  BN, reduced mesh (4 blocks 256x64x64), 4 ranks x 1 block, 5 cycles (dumps at 0 and 5)
#   base_r4 BO, same as mem_r4 -> bitwise cmp mem_r4 vs base_r4, peak per block before/after
#   mem_c4  BN, FULL mesh, 4 ranks x 4 blocks, 6 cycles
#   base_c4 BO, FULL mesh, 4 ranks x 4 blocks, 6 cycles (OOM is a valid result: report it)
# env: ARMS (subset/order), SRUN2/SRUN4 (launchers for 2 / 4 ranks, 1 GPU per rank)
BO=$1; BN=$2; F=$3; S=$4; KT=${5:-}
SRUN2=${SRUN2:-srun -n 2 --gpus-per-task=1 --gpu-bind=closest --cpu-bind=cores}
SRUN4=${SRUN4:-srun -n 4 --gpus-per-task=1 --gpu-bind=closest --cpu-bind=cores}
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
IH=$F/bsg_hr_dc5.athinput
GR="mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472 mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64"
for x in "$BO" "$BN"; do [ -x "$x" ] || { echo "no binary $x"; exit 1; }; done
[ -f "$IH" ] || { echo "no input $IH"; exit 1; }
go(){ # name bin launcher mesh(red|full) nlim
  D=$S/$1; mkdir -p $D && cd $D || exit 1
  G=""; [ "$4" = red ] && G="$GR"
  echo "== $1 $(md5sum $2 | cut -c1-32) $4 nlim=$5 start $(date +%T)"; t0=$(date +%s)
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -l 1 > smi.log 2>&1 &
  SP=$!
  KOKKOS_TOOLS_LIBS=$KT $3 $2 -i $IH -d $D $G time/nlim=$5 time/ndiag=1 output5/dt=1e30 \
    < /dev/null 2>&1 | awk '{print systime(), $0; fflush()}' > out.log
  rc=${PIPESTATUS[0]}; kill $SP
  pk=$(awk -F, '{if($2>m[$1])m[$1]=$2} END{for(i in m) printf "gpu%s=%dMiB ", i, m[i]}' smi.log)
  c1=$(grep "cycle=1 " out.log | head -1 | awk '{print $1}'); cl=$(grep "cycle=" out.log | tail -1 | awk '{print $1}')
  nc=$(grep "cycle=" out.log | tail -1 | sed 's/.*cycle=\([0-9]*\).*/\1/')
  spc=$(awk -v a="$c1" -v b="$cl" -v n="$nc" 'BEGIN{if(n>1&&a!="")printf "%.2f", (b-a)/(n-1); else print "na"}')
  echo "rc=$rc fatal=$(grep -c FATAL out.log) oom=$(grep -ci 'out of memory\|failed to allocate' out.log)" \
       "wall=$(( $(date +%s)-t0 ))s s/cycle=$spc peak: $pk $(date +%T)"; }
for a in ${ARMS:-mem_c8 mem_r4 base_r4 mem_c4 base_c4}; do
  case $a in
    mem_c8)  go $a $BN "$SRUN2" full 6;;
    mem_r4)  go $a $BN "$SRUN4" red 5;;
    base_r4) go $a $BO "$SRUN4" red 5;;
    mem_c4)  go $a $BN "$SRUN4" full 6;;
    base_c4) go $a $BO "$SRUN4" full 6;;
  esac
done
if [ -d $S/mem_r4 ] && [ -d $S/base_r4 ]; then
  n=0; d=0
  for f in $S/base_r4/bin/*.bin $S/base_r4/*.hst; do
    g=$S/mem_r4/${f#$S/base_r4/}; n=$((n+1)); cmp -s $f $g || { d=$((d+1)); echo "  DIFF $(basename $f)"; }
  done
  echo "bitwise base_r4 vs mem_r4: $n files, $d differ"
fi
