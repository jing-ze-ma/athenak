#!/bin/bash
# TASK-2026-10-10-caltech-vgdfuse-scaling: run ONE arm inside an allocation and print its result line.
# usage: bash vgdfuse_run.sh <out dir> <name> <binary> <input> <nranks> <nlim> "<extra keys>"
# - copies <input> to <out dir>/<name>.athinput with the vgdfuse-1010 keys added to <rad_m1> (at the 1010
#   defaults, so that arms can override them on the command line; BASE ignores unknown keys);
# - runs nranks ranks (1 GPU per rank) with dumps off (bin, rst), hst kept;
# - prints "ARM <name> rc= fatal= oom= s/cycle mean(6-N)= median(6-N)= steady(no rotation)= peakMiB= bin= keys=".
# env: SRUN (default "srun --mpi=pmix -c 8 --cpu-bind=cores $GPUWRAP"), GPUWRAP (sets CUDA_VISIBLE_DEVICES =
#      SLURM_LOCALID), KT=1: KOKKOS_TOOLS_LIBS=$KTLIB (kernel timer; KOKKOS_TOOLS_TIMER_BINARY=1 is exported)
S=$1; N=$2; B=$3; I=$4; NR=$5; NL=$6; X=$7
mkdir -p $S/$N || exit 1
IN=$S/$N.athinput
sed 's/^<rad_m1>.*$/&\nvet_gd_halo_exact = 2\nvet_gd_overlap = 0\nvet_gd_overlap_mb = 1024\nvet_gd_fuse_shells = 1/' $I > $IN
grep -q "^vet_gd_fuse_shells = 1" $IN || { echo "ARM $N SETUP_FAILED (no <rad_m1> block?)"; exit 1; }
# vet_gd_twin_lowmem named false = its default (so that the lowmem arm can set it)
grep -q "^ *vet_gd_twin_lowmem" $IN || sed -i 's/^vet_gd_fuse_shells = 1$/&\nvet_gd_twin_lowmem = false/' $IN
O="output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 time/ndiag=1"
# (keys absent from an input make AthenaK fatal: keep only the output blocks that exist)
OK=""; for k in $O; do b=${k%%/*}; [ "$b" = time ] || grep -q "^<$b>" $IN || continue; OK="$OK $k"; done
export KOKKOS_TOOLS_TIMER_BINARY=1 OMP_NUM_THREADS=1 HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
L=""; [ "$KT" = 1 ] && L=$KTLIB
SRUN=${SRUN:-srun --mpi=pmix -c 8 --cpu-bind=cores $GPUWRAP}
( while true; do nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits >> $S/$N/mem.txt 2>/dev/null; sleep 5; done ) &
MP=$!
cd $S/$N || exit 1
KOKKOS_TOOLS_LIBS=$L $SRUN -n $NR $B -i $IN -d $S/$N $OK time/nlim=$NL $X < /dev/null > out.log 2>&1
rc=$?
kill $MP 2>/dev/null
python3 - "$S/$N/out.log" "$S/$N/mem.txt" > res.txt <<'EOF'
import re, statistics, sys
el = {}
for x in open(sys.argv[1], errors='replace'):
    m = re.match(r'elapsed=([0-9.eE+-]+) cycle=(\d+)', x)
    if m:
        el[int(m.group(2))] = float(m.group(1))
d = {c: el[c] - el[c - 1] for c in el if c - 1 in el and c > 6}
v = list(d.values())
st = [t for c, t in d.items() if (c - 1) % 10 != 0]
pk = 0
try:
    pk = max(int(y) for y in open(sys.argv[2]).read().split())
except (OSError, ValueError):
    pass
if v:
    print('%.4f %.4f %.4f %d' % (sum(v)/len(v), statistics.median(v), sum(st)/max(len(st), 1), pk))
else:
    print('nan nan nan %d' % pk)
EOF
read mean med st pk < res.txt
echo "ARM $N rc=$rc fatal=$(grep -c FATAL out.log) oom=$(grep -ci 'out of memory' out.log)" \
     "s/cycle mean(7-$NL)=$mean median=$med steady(no rotation)=$st peakMiB(max GPU on node)=$pk" \
     "nranks=$NR bin=$(md5sum $B | cut -c1-8) keys=$X"
grep -h "vgdfuse-1010 (rank 0)\|vgdspeed-1009 (rank 0)" out.log | cut -c1-260
