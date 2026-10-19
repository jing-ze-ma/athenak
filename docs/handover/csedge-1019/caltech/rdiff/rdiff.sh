#!/bin/bash
# TASK-2026-10-19-caltech-rsg-restart-diff: R2 (identical restart), CV (continuous, vet_rotate_avg=1), RV (restart, vet_rotate_avg=1)
# same layout/deck/binary as det1.sh (job 5008037): cs N=48, 24 blocks, 1 node x 4 H200 (6 blocks/GPU), rst from n48/C cycle 100.
#SBATCH -A carnegie_poc -p gpu --exclude=hpc-sm-01-09
#SBATCH -N 1 --ntasks-per-node=4 --gres=gpu:nvidia_h200:4 --cpus-per-task=8 --mem=480G
#SBATCH -t 01:00:00 -J rsg_rdiff -o /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/det/rdiff_1019/rdiff.%j.out
P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; S=$P/det/rdiff_1019; source $P/run.cfg; IN=$P/det/rsg_det.athinput
source /etc/profile >/dev/null 2>&1
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export OMP_NUM_THREADS=1 UCX_TLS=^xpmem OMPI_MCA_io=romio321
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date)"; md5sum $BIN $IN $S/rst_C/rsg.00001.rst
K="mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 time/nlim=200 output6/dt=1e30"
R="srun --mpi=pmix -N 1 -n 4 --ntasks-per-node=4 -c 8 --cpu-bind=cores $WRAP $BIN"
RS=$S/rst_C/rsg.00001.rst
run() { d=$1; shift; mkdir -p $S/$d; echo "== $d $(date)"; echo "$R $* -d $S/$d" | tee $S/$d/cmdline.txt
  T0=$SECONDS; $R "$@" -d $S/$d > $S/$d/run.log 2>&1; echo "$d rc=$? wall $((SECONDS-T0)) s"
  grep -E '^elapsed' $S/$d/run.log | tail -1; }
run R2 -r $RS $K
run CV -i $IN $K rad_m1/vet_rotate_avg=1
run RV -r $RS $K rad_m1/vet_rotate_avg=1
echo "== compare"; bash $S/compare.sh; echo "done $(date)"
