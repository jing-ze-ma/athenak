#!/bin/bash
# TASK-2026-10-19 determinism, REDUCED smoke layout (Raven smoke: cs N=48, 24 blocks 591x24x24) on 1 node x 4 H200 (6 blocks/GPU).
# C = 200 cyc fresh (rst/bins at cycle 100, 200; hst every cycle), C2 = identical rerun, R = restart from C cycle-100 rst to 200.
#SBATCH -A carnegie_poc -p gpu --exclude=hpc-sm-01-09
#SBATCH -N 1 --ntasks-per-node=4 --gres=gpu:nvidia_h200:4 --cpus-per-task=8 --mem=480G
#SBATCH -t 01:30:00 -J rsg_det1 -o /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/det/det1.%j.out
P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; S=$P/det/n48; source $P/run.cfg; IN=$P/det/rsg_det.athinput
source /etc/profile >/dev/null 2>&1
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export OMP_NUM_THREADS=1 UCX_TLS=^xpmem OMPI_MCA_io=romio321
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date)"; md5sum $BIN $IN
K="mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 time/nlim=200 output6/dt=1e30"
R="srun --mpi=pmix -N 1 -n 4 --ntasks-per-node=4 -c 8 --cpu-bind=cores $WRAP $BIN"
run() { d=$1; shift; mkdir -p $S/$d; echo "== $d $(date)"
  nvidia-smi --query-gpu=index,memory.used --format=csv,noheader -l 20 > $S/gpumem_$d.txt 2>&1 & NV=$!
  T0=$SECONDS; $R "$@" -d $S/$d > $S/$d/run.log 2>&1; echo "$d rc=$? wall $((SECONDS-T0)) s"; kill $NV
  L=$S/$d/run.log; echo "  peak MiB $(awk -F, '{gsub(/ MiB/,"",$2); if($2>m)m=$2} END{print m}' $S/gpumem_$d.txt) FATAL $(grep -c FATAL $L) NaN $(tail -60 $L | grep -ciw nan) DIVERGED $(grep -c DIVERGED $L) NC>0 $(grep -cE 'NON-CONVERGED=[1-9]' $L)"
  grep 'Picard iterations' $L | tail -1; grep -E 'cpu time used|zone-cycles/cpu_second|cycle=' $L | tail -3; }
run C  -i $IN $K
run C2 -i $IN $K
RS=""; for f in $S/C/rst/*.rst; do [ "$(python3 $INFO $f | head -1 | awk '{print $2}')" = 100 ] && RS=$f; done
echo "R from $RS"; [ -n "$RS" ] && run R -r $RS $K
echo "== compare"
hc() { python3 - "$@" <<'PY'
import sys
def rows(p):
    d={}
    for l in open(p):
        if l.startswith('#'): continue
        d[l.split()[0]]=l
    return d
a,b=rows(sys.argv[1]),rows(sys.argv[2]); ks=sorted(set(a)&set(b),key=float)
bad=[k for k in ks if a[k]!=b[k]]
print(sys.argv[3], 'hst', 'IDENT' if not bad and ks else 'DIFF', len(ks),'rows', ('first diff t='+bad[0]) if bad else '')
PY
}
hc $S/C/rsg.user.hst $S/C2/rsg.user.hst "C-C2" 2>/dev/null || hc $S/C/rsg.hst $S/C2/rsg.hst "C-C2"
for h in $S/C/*.hst; do n=$(basename $h); hc $h $S/C2/$n "C-C2 $n"; [ -f $S/R/$n ] && hc $h $S/R/$n "R-C $n"; done
for f in $S/C/bin/*.bin; do n=$(basename $f); python3 /resnick/groups/carnegie_poc/jingze/hegiant_large_1014/bincmp.py $f $S/C2/bin/$n | sed 's/^/C-C2 /'
  [ -f $S/R/bin/$n ] && python3 /resnick/groups/carnegie_poc/jingze/hegiant_large_1014/bincmp.py $f $S/R/bin/$n | sed 's/^/R-C /'; done
echo "done $(date)"
