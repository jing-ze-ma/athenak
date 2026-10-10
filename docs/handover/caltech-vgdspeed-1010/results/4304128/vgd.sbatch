#!/bin/bash
#SBATCH -J vgdspeed
#SBATCH -A carnegie_poc
#SBATCH -p gpu
#SBATCH --exclude=hpc-sm-01-09
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:nvidia_h200:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=240G
#SBATCH -t 01:00:00
#SBATCH -o /resnick/groups/carnegie_poc/jingze/vgdspeed_1010/j.%j.out
# TASK-2026-10-10-caltech-vgdspeed: BASE b2f2e897 vs NEW ec0757c1, interleaved B N B N, 1 node x 2 H200
module purge; module load gcc/13.2.0-gcc-13.2.0-w55nxkl cuda/12.9.0-none-none-pfmzfdv hpcx/2.17.1/hpcx-ompi
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1 KOKKOS_TOOLS_TIMER_BINARY=1
P=/resnick/groups/carnegie_poc/jingze/prof_1009
V=/resnick/groups/carnegie_poc/jingze/vgdspeed_1010
BB=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_b2f2e897_nofma
BN=/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_ec0757c1_nofma
mkdir -p $V/runs/$SLURM_JOB_ID; nvidia-smi -L; hostname; md5sum $BB $BN
export SRUN4="srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores /resnick/groups/carnegie_poc/jingze/bsg_1009/gpuwrap.sh"
A1="agc_plain:agc:plain:60 agc_kt10:agc:kt:10 agc_kt60:agc:kt:60 bsg_plain:bsg:plain:30 bsg_kt5:bsg:kt:5 bsg_kt30:bsg:kt:30"
A2="agc_plain:agc:plain:60 bsg_plain:bsg:plain:30"
for r in "base1:$BB:$A1" "new1:$BN:$A1" "base2:$BB:$A2" "new2:$BN:$A2"; do
  IFS=: read n b a <<< "$r"
  echo "######## $n $b $(date +%T)"
  ARMS="$a" bash $P/bundle/docs/handover/deltaai-prof-1009/prof_cuda.sh "$b" \
    /resnick/groups/carnegie_poc/jingze/agcarb_1009/files /resnick/groups/carnegie_poc/jingze/bsg_1009/files \
    $V/runs/$SLURM_JOB_ID/$n $P/kt/profiling/simple-kernel-timer/kp_kernel_timer.so 2>&1 | tee $V/runs/$SLURM_JOB_ID/$n.log
  [ ${PIPESTATUS[0]} -eq 0 ] || { echo "STOP: $n failed"; exit 1; }
  if grep -q "^rc=[^0]" $V/runs/$SLURM_JOB_ID/$n.log; then echo "STOP: rc != 0 in $n"; exit 1; fi
  if grep -q FATAL $V/runs/$SLURM_JOB_ID/$n/*/out.log; then echo "STOP: FATAL in $n"; exit 1; fi
done
echo DONE $(date +%T)
