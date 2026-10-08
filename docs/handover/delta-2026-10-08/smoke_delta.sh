#!/bin/bash
#SBATCH -J agcsmoke
#SBATCH -A bivj-delta-gpu
#SBATCH -p gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH -o /work/nvme/bivj/jma20/delta_1008/smoke/smoke.%j.out
# AG Car A then B smoke (time/nlim=10), 1 node x 4 A100, 4 ranks = one 480x64x64 MeshBlock per GPU.
# usage: sbatch smoke_delta.sh <binary> [cases, default "A B"]
B=/work/nvme/bivj/jma20/delta_1008
BIN=$1
[ -x "$BIN" ] || { echo "no binary $BIN"; exit 1; }
module unload cudatoolkit
module load cuda/12.9
module -t list 2>&1 | tr '\n' ' '; echo
export MPICH_GPU_SUPPORT_ENABLED=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
export OMP_NUM_THREADS=1
echo "job $SLURM_JOB_ID node $SLURMD_NODENAME submit->start: $(squeue -h -j $SLURM_JOB_ID -o '%V %S')"
md5sum $BIN $B/agcar_files/*.txt $B/agcar_files/*.athinput
nvidia-smi --query-gpu=index,name,memory.total,driver_version --format=csv
for C in ${2:-A B}; do
  D=$B/smoke/smoke$C.$SLURM_JOB_ID; mkdir -p $D
  echo "case $C start $(date +%T)"
  srun -n 4 --cpu-bind=cores $BIN -i $B/agcar_files/agcar_shake${C}_ge.athinput -d $D -t 00:08:00 \
      time/nlim=10 > $D/run.log 2>&1
  echo "case $C rc=$? end $(date +%T)"
done
