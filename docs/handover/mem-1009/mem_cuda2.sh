#!/bin/bash
# BSG GPU-memory + s/cycle + bitwise, round 2 (mem-1009 vet_gd_twin_lowmem).
# usage (inside a 1-node allocation with 4 GPUs):
#   bash mem_cuda2.sh <BO = 98835d99 binary> <BN = mem-1009 6c5d8fb2 binary> <files dir> <out dir> [KT .so]
# arms "name:B:nranks:mesh:nlim:lm" (B = O|N, mesh = red|full, lm = 1: input copy with vet_gd_twin_lowmem = true)
# bitwise: a lowmem bin differs from base only in its header (the extra input line): compare the
# data after <par_end> (bitcmp.py) and the hst byte-wise
# default order (the first two answer the Caltech question):
#   lm_c8   N, FULL mesh, 2 ranks x 8 blocks, lowmem      -> fits? peak/GPU, s/cycle (THE CALTECH LAYOUT)
#   lm_c4   N, FULL mesh, 4 ranks x 4 blocks, lowmem      -> bitwise vs base_c4 (dumps at 0 and 5), cost/cycle
#   base_c4 O, FULL mesh, 4 ranks x 4 blocks
#   new_c4  N, FULL mesh, 4 ranks x 4 blocks, keys as input (fused twin) -> bitwise vs base_c4
#   lm_r4   N, reduced mesh, 4 ranks x 1 block, lowmem   -> bitwise vs base_r4, peak per block
#   base_r4 O, reduced mesh, 4 ranks x 1 block
# env: ARMS, SRUN2/SRUN4 (launchers for 2 / 4 ranks, 1 GPU per rank)
BO=$1; BN=$2; F=$3; S=$4; KT=${5:-}
SRUN2=${SRUN2:-srun -n 2 --gpus-per-task=1 --gpu-bind=closest --cpu-bind=cores}
SRUN4=${SRUN4:-srun -n 4 --gpus-per-task=1 --gpu-bind=closest --cpu-bind=cores}
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
IH=$F/bsg_hr_dc5.athinput
GR="mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472 mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64"
for x in "$BO" "$BN"; do [ -x "$x" ] || { echo "no binary $x"; exit 1; }; done
[ -f "$IH" ] || { echo "no input $IH"; exit 1; }
go(){ # name bin nranks mesh nlim lm
  D=$S/$1; mkdir -p $D && cd $D || exit 1
  G=""; [ "$4" = red ] && G="$GR"
  L="$SRUN4"; [ "$3" = 2 ] && L="$SRUN2"
  # lowmem: the key must be IN the input file (AthenaK rejects a command-line key that the
  # input does not have): an input copy with the line after vet_gd_twin_fuse
  K=""; I=$IH
  if [ "$6" = 1 ]; then
    I=$S/bsg_hr_dc5_lm.athinput; K="(in file $I)"
    [ -f $I ] || sed '/^vet_gd_twin_fuse/a vet_gd_twin_lowmem = true' $IH > $I
    grep -q '^vet_gd_twin_lowmem = true' $I || { echo "no lowmem line in $I"; return; }
  fi
  echo "== $1 $(md5sum $2 | cut -c1-32) n=$3 $4 nlim=$5 keys=[$K] start $(date +%T)"; t0=$(date +%s)
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits -l 1 > smi.log 2>&1 &
  SP=$!
  KOKKOS_TOOLS_LIBS=$KT $L $2 -i $I -d $D $G time/nlim=$5 time/ndiag=1 output5/dt=1e30 \
    < /dev/null > out.log 2>&1
  rc=$?; kill $SP
  pk=$(awk -F, '{if($2>m[$1])m[$1]=$2} END{for(i in m) printf "gpu%s=%dMiB ", i, m[i]}' smi.log)
  # steady s/cycle: mean wall of cycles 2-4 from the elapsed= log lines (no dump cycles)
  spc=$(grep -o "elapsed=[0-9.e+-]* cycle=[0-9]*" out.log | sed 's/elapsed=//;s/ cycle=/ /' | \
        awk '{t[$2]=$1} END{if(4 in t && 1 in t) printf "%.3f", (t[4]-t[1])/3; else print "na"}')
  echo "rc=$rc fatal=$(grep -c FATAL out.log) oom=$(grep -ci 'out of memory\|failed to allocate' out.log)" \
       "wall=$(( $(date +%s)-t0 ))s s/cycle(2-4)=$spc peak: $pk $(date +%T)"
  grep -h "ragged per-shell band\|vet_gd_twin_lowmem=" out.log | head -2; }
cmpd(){ # a b
  [ -d $S/$1 ] && [ -d $S/$2 ] || return
  n=0; d=0
  for f in $S/$1/bin/*.bin $S/$1/*.hst; do
    g=$S/$2/${f#$S/$1/}; n=$((n+1)); cmp -s $f $g || { d=$((d+1)); echo "  DIFF $(basename $f)"; }
  done
  echo "bitwise $1 vs $2: $n files, $d differ"; }
for a in ${ARMS:-lm_c8:N:2:full:6:1 lm_c4:N:4:full:6:1 base_c4:O:4:full:6:0 new_c4:N:4:full:6:0 lm_r4:N:4:red:5:1 base_r4:O:4:red:5:0}; do
  IFS=: read n b r m l k <<< "$a"
  B=$BN; [ "$b" = O ] && B=$BO
  go $n $B $r $m $l $k
done
cmpd base_c4 lm_c4; cmpd base_c4 new_c4; cmpd base_r4 lm_r4; cmpd base_r4 new_r4; cmpd new_r4 lm_r4
cmpd base_r2 lm_r2; cmpd base_r2 new_r2; cmpd new_r2 lm_r2
