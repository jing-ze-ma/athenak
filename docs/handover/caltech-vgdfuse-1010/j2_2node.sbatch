#!/bin/bash
#SBATCH -J vgdfuse2
#SBATCH -A carnegie_poc
#SBATCH -p gpu
#SBATCH --exclude=hpc-sm-01-09
#SBATCH -N 2
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:nvidia_h200:4
#SBATCH --cpus-per-task=8
#SBATCH --mem=480G
#SBATCH -t 00:45:00
#SBATCH -o /resnick/groups/carnegie_poc/jingze/vgdfuse_1010/j2.%j.out
# TASK-2026-10-10-caltech-vgdfuse-scaling, job 2 (2 nodes x 4 H200): BSG hr 16 blocks on 8 H200 (2 blocks per GPU),
# BASE and NEW, plus NEW repeated (run-to-run).  Same CONFIG block as j1_1node.sbatch.
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
# ---- CONFIG ----
V=/resnick/groups/carnegie_poc/jingze/vgdfuse_1010
BB=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_5b304cf9_nofma     # BASE
BN=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_712ace27_nofma     # NEW
IB=/resnick/groups/carnegie_poc/jingze/bsg_1009/files/bsg_hr_dc5.athinput
export GPUWRAP=/resnick/groups/carnegie_poc/jingze/bsg_1009/gpuwrap.sh
R=$(dirname $(readlink -f $0))/vgdfuse_run.sh
[ -f "$R" ] || R=$V/vgdfuse_run.sh
# ----------------
S=$V/runs/$SLURM_JOB_ID; mkdir -p $S; nvidia-smi -L; hostname; md5sum $BB $BN $IB
BHR="mesh/nx2=256 mesh/nx3=256 meshblock/nx2=64 meshblock/nx3=64"
for a in bsghr_s8_base:$BB bsghr_s8_new:$BN bsghr_s8_new2:$BN; do
  n=${a%%:*}; b=${a#*:}
  bash $R $S $n $b $IB 8 20 "$BHR" | tee -a $S/results.txt
  grep "^ARM $n " $S/results.txt | grep -q "rc=0 fatal=0" || { echo "STOP: $n failed"; exit 1; }
done
echo DONE $(date +%T)
