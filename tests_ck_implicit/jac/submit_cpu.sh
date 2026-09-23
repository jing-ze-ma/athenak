#!/bin/bash -l
# the ck-jac CPU gate suite (A1, A2, D, E; correctness only, no timing) on the host cores
# of ONE apu node (the login node caps a user at 12 cores).  GPUs unused.
#SBATCH -J ckjac_cpu
#SBATCH -p apu
#SBATCH --constraint="apu"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --gres=gpu:2
#SBATCH -t 04:00:00
#SBATCH -o /viper/ptmp2/jinma/ckjac_0923/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/ckjac_0923/log.err.%j
export OMP_NUM_THREADS=1
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
P=/viper/ptmp2/jinma/ckjac_0923
cd $P
ARMS=${ARMS:-"t4 sub esc aa aas aae"}
for a in $ARMS; do
  bash $P/run_jac.sh A1 $a > A1_$a.log 2>&1 &
  bash $P/run_jac.sh D $a > D_$a.log 2>&1 &
done
( bash $P/run_jac.sh E $ARMS > E.log 2>&1; bash $P/run_jac.sh A2 $ARMS > A2.log 2>&1 ) &
wait
echo ALL DONE
