#!/bin/bash
#SBATCH -J radprof
#SBATCH -A bivj-dtai-gh
#SBATCH -p ghx4-interactive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:30:00
#SBATCH -o /work/nvme/bivj/jma20/prof_1009/j.%j.out
# DeltaAI prof (TASK-2026-10-09-deltaai-prof.md): ONE job, 1 node x 4 GH200, 8 arms in sequence.
# Module stack: the same one the AG Car smokes / mem rounds ran with (the build stack of
# athena_hes_gpu_98835d99a16f: PrgEnv-gnu 8.6.0, cray-mpich 9.0.1, cudatoolkit 25.5_12.9). If your smoke/mem
# job scripts load modules explicitly, put the same lines here.
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # job-env rule (ROCm-only, no-op on CUDA)
W=/work/nvme/bivj/jma20
P=$W/prof_1009
# binary + md5 exactly as the AG Car smokes used them
BIN=""; BIN_MD5=""
source $W/agcar_1009/smokeB/run.cfg
M=$(md5sum "$BIN" | cut -d' ' -f1)
[ "$M" = c9c6d16705d07dc6e97de6c808d02182 ] || { echo "md5 $M of $BIN != c9c6d167...: exit"; exit 1; }
bash $P/bundle/docs/handover/deltaai-prof-1009/prof_cuda.sh \
  "$BIN" \
  $W/agcar_1009/files \
  $W/bsg_hrdet_1009/files \
  $P/runs/$SLURM_JOB_ID \
  $W/bsg_hrdet_1009/kt/profiling/simple-kernel-timer/kp_kernel_timer.so
