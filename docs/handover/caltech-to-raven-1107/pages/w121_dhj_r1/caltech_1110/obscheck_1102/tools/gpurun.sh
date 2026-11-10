#!/bin/bash
# gpurun.sh NODE GPU cmd... : take the mkdir lock lock_<node>_<g> (owner file), check idle (<100 MiB), run, release.
N=$1; G=$2; shift 2
L=/resnick/groups/carnegie_poc/jingze/lhllc2_1102/dhj/lock_${N}_${G}
mkdir $L 2>/dev/null || { echo "LOCK BUSY $L: $(cat $L/owner 2>/dev/null)"; exit 3; }
echo "obscheck_1102 $$ $(date -Iseconds)" > $L/owner
trap 'rm -rf $L' EXIT
U=$(nvidia-smi -i $G --query-gpu=memory.used --format=csv,noheader,nounits)
[ "$U" -lt 100 ] || { echo "GPU $G busy ($U MiB)"; exit 4; }
export CUDA_VISIBLE_DEVICES=$G OMP_NUM_THREADS=8
"$@"
