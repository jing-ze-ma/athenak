#!/bin/bash
# TASK-2026-10-19-caltech-csedge: one arm per job, 1 H200, 24 blocks. usage: sbatch -J cse_<ARM> arm.sh <ARM> [extra keys...]
#SBATCH -A carnegie_poc -p gpu --exclude=hpc-sm-01-09
#SBATCH -N 1 -n 1 --gres=gpu:nvidia_h200:1 --cpus-per-task=8 --mem=200G
#SBATCH -t 16:00:00 -o /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/csedge_1019/%x.%j.out
P=/resnick/groups/carnegie_poc/jingze/rsg_g25_prod; source $P/run.cfg; IN=$P/rsg_g25_prod.athinput
ARM=$1; shift; D=$P/csedge_1019/$ARM; mkdir -p $D
source /etc/profile >/dev/null 2>&1
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export OMP_NUM_THREADS=1 UCX_TLS=^xpmem
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST $(date) arm $ARM"; md5sum $BIN $IN
M="mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24"
[ "$ARM" = B24 ] && M="mesh/nx2=24 mesh/nx3=24 meshblock/nx2=24 meshblock/nx3=24"
K="$M time/tlim=1.04e6 time/ndiag=100 output1/dt=2000 output2/dt=43200 output3/dt=43200 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 $*"
CMD="srun --mpi=pmix -N 1 -n 1 -c 8 --cpu-bind=cores $WRAP $BIN -i $IN $K -d $D"
echo "$CMD" | tee $D/cmdline.txt
nvidia-smi --query-gpu=index,memory.used --format=csv,noheader -l 60 > $D/gpumem.txt 2>&1 & NV=$!
T0=$SECONDS; $CMD > $D/run.log 2>&1; echo "rc=$? wall $((SECONDS-T0)) s"; kill $NV
grep -E '^elapsed' $D/run.log | tail -2; echo "done $(date)"
