#!/bin/bash
# TASK-2026-10-19-caltech-rsg-restart-diff: R2 (identical restart), CV (continuous, vet_rotate_avg=1), RV (restart, vet_rotate_avg=1)
# same layout/deck/binary as det1.sh (job 5008037): cs N=48, 24 blocks, 1 node x 4 H200 (6 blocks/GPU), rst from n48/C cycle 100.
#SBATCH -A carnegie_poc -p gpu --exclude=hpc-sm-01-09
#SBATCH -N 1 --ntasks-per-node=4 --gres=gpu:nvidia_h200:4 --cpus-per-task=8 --mem=480G
#SBATCH -t 00:30:00 -J rsg_rv2 -o /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/det/rdiff_1019/rv2.%j.out
# RV2 = restart from CV own rst.00001 (cycle 100) with vet_rotate_avg=1 -> RV2 vs CV (RV used C rst = invalid pair)
P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; S=$P/det/rdiff_1019; source $P/run.cfg; IN=$P/det/rsg_det.athinput
source /etc/profile >/dev/null 2>&1
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export OMP_NUM_THREADS=1 UCX_TLS=^xpmem OMPI_MCA_io=romio321
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date)"; md5sum $BIN $IN $S/rst_C/rsg.00001.rst
K="mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 time/nlim=200 output6/dt=1e30"
R="srun --mpi=pmix -N 1 -n 4 --ntasks-per-node=4 -c 8 --cpu-bind=cores $WRAP $BIN"
RS=$S/CV/rst/rsg.00001.rst
run() { d=$1; shift; mkdir -p $S/$d; echo "== $d $(date)"; echo "$R $* -d $S/$d" | tee $S/$d/cmdline.txt
  T0=$SECONDS; $R "$@" -d $S/$d > $S/$d/run.log 2>&1; echo "$d rc=$? wall $((SECONDS-T0)) s"
  grep -E '^elapsed' $S/$d/run.log | tail -1; }
run RV2 -r $RS $K rad_m1/vet_rotate_avg=1
python3 - $S <<'PY'
import sys,glob,os
S=sys.argv[1]
def rows(p):
    d={}
    for l in open(p):
        if l.startswith('#'): continue
        x=l.split()
        if x: d[x[0]]=l
    return d
def data(p):
    b=open(p,'rb').read(); i=b.find(b'<par_end>'); i=b.find(b'\n',i)+1; return b[i:]
for h in ('rsg.hydro.hst','rsg.user.hst'):
    a=rows(f'{S}/CV/{h}'); b=rows(f'{S}/RV2/{h}'); ks=sorted(set(a)&set(b),key=float); d=[k for k in ks if a[k]!=b[k]]
    print('RV2-CV',h,'IDENT' if ks and not d else 'DIFF',len(ks),'rows',('first diff t='+d[0]) if d else '')
for g in sorted(glob.glob(f'{S}/RV2/bin/*.bin')):
    f=f'{S}/CV/bin/'+os.path.basename(g)
    if os.path.exists(f): print('RV2-CV',os.path.basename(g),'SAME' if data(f)==data(g) else 'DIFF')
PY
echo "done $(date)"
