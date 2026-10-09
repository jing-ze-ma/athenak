#!/bin/bash
#SBATCH -J radprof
#SBATCH -A bivj-delta-gpu
#SBATCH -p gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH -o /work/nvme/bivj/jma20/delta_1008/prof/j.%j.out
# prof batch1 (TASK-2026-10-08-delta-prof-batch1.md). Partition: gpuA100x4 or gpuA100x4-interactive
# (whichever starts first; interactive bills 2x: still <= 4 GPU-h at the 30 min limit).
module unload cudatoolkit; module load cuda/12.9
W=/work/nvme/bivj/jma20/delta_1008
P=$W/prof
bash $P/bundle/docs/handover/delta-prof-1009/prof_cuda.sh \
  $W/bin/athena_he_gpu_98835d99 \
  $P/agcar_files \
  $W/bsg_hrdet/files \
  $P/runs/$SLURM_JOB_ID \
  $P/kt/profiling/simple-kernel-timer/kp_kernel_timer.so
