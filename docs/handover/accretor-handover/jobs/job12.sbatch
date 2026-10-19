#!/bin/bash -l
# ry_per_accretor stage-2 (envelope) GPU tests (apudev, 1 node x 2 MI300A, 2 ranks).
# usage: sbatch job.sbatch ARM [ARM2 ...]; keys for each ARM in $G/$ARM/keys.txt
#SBATCH -p apudev
#SBATCH --constraint=apu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
G=/viper/ptmp2/jinma/accretor_1006/gpu2
B=/viper/ptmp2/jinma/accretor_1006/bin/athena_ryper_gpu72_env7
BIN_MD5=a5ad0a18c1d3a73e4cb41e833f1821d5
INP=/viper/ptmp2/jinma/accretor_1006/bin/plaskett_env12.athinput
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
M=$(md5sum $B | cut -d' ' -f1); echo "$M  $B"
[ "$M" != "$BIN_MD5" ] && { echo "md5 mismatch: exit 1"; exit 1; }
for ARM in "$@"; do
  cd $G/$ARM || exit 1
  KEYS=$(cat $G/$ARM/keys.txt)
  echo "job $SLURM_JOB_ID arm $ARM start $(date +%F_%T)"
  echo "srun -n $SLURM_NTASKS $B -i $INP $KEYS -d $G/$ARM"
  srun -n $SLURM_NTASKS -c 24 $B -i $INP $KEYS -d $G/$ARM > $G/$ARM/run.log 2>&1
  rc=$?
  echo "arm $ARM rc=$rc $(date +%F_%T)"
  grep -q "FATAL" $G/$ARM/run.log && { echo "FATAL in $ARM/run.log"; exit 3; }
  [ $rc -ne 0 ] && exit $rc
done
exit 0
