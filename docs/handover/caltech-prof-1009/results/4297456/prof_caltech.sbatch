#!/bin/bash
#SBATCH -J radprof
#SBATCH -A carnegie_poc
#SBATCH -p gpu
#SBATCH --exclude=hpc-sm-01-09
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:nvidia_h200:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=240G
#SBATCH -t 01:00:00
#SBATCH -o /resnick/groups/carnegie_poc/jingze/prof_1009/j.%j.out
# Caltech port of deltaai-prof-1009/prof.sbatch (TASK-2026-10-09-caltech-prof.md): 1 node x 2 H200, 8 arms in sequence.
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1
P=/resnick/groups/carnegie_poc/jingze/prof_1009
# first-result-wins check at job start
G=/home/jingze/ATHENAK/athenak
if git -C $G fetch -q origin 2>/dev/null; then
  if git -C $G cat-file -e origin/bsg-files-1009:docs/handover/NOTE-2026-10-09-deltaai-prof.md 2>/dev/null || \
     git -C $G cat-file -e origin/rad-beam-1008:docs/handover/NOTE-2026-10-09-delta-prof-batch2.md 2>/dev/null; then
    echo "SKIP: DeltaAI/Delta prof NOTE already on the fork"; exit 0; fi
else echo "warning: git fetch failed at job start; proceeding"; fi
BIN=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_98835d99_nofma
M=$(md5sum "$BIN" | cut -d' ' -f1)
[ "$M" = bc6f4e9e11c106a3606246d2ac40d296 ] || { echo "md5 $M != bc6f4e9e: exit"; exit 1; }
nvidia-smi -L; hostname
export SRUN4="srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores /resnick/groups/carnegie_poc/jingze/bsg_1009/gpuwrap.sh"
bash $P/bundle/docs/handover/deltaai-prof-1009/prof_cuda.sh \
  "$BIN" \
  /resnick/groups/carnegie_poc/jingze/agcarb_1009/files \
  /resnick/groups/carnegie_poc/jingze/bsg_1009/files \
  $P/runs/$SLURM_JOB_ID \
  $P/kt/profiling/simple-kernel-timer/kp_kernel_timer.so
