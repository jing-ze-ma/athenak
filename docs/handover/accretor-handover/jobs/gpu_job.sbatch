#!/bin/bash -l
# ry_per_accretor GPU tests (apudev, 1 node x 2 MI300A, 2 ranks).  usage: sbatch job.sbatch ARM
# keys for ARM are in $G/$ARM/keys.txt (one line of command-line overrides)
#SBATCH -p apudev
#SBATCH --constraint=apu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
G=/viper/ptmp2/jinma/accretor_1006/gpu
B=/viper/ptmp2/jinma/accretor_1006/bin/athena_ryper_gpu72_1e8d0ff8
BIN_MD5=b6ab7220c3dab29fb93bfd03d5bfc6b5
INP=/viper/ptmp2/jinma/accretor_1006/bin/ry_per_accretor_1e8d0ff8.athinput
ARM=$1
cd $G/$ARM || exit 1
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
M=$(md5sum $B | cut -d' ' -f1); echo "$M  $B"
[ "$M" != "$BIN_MD5" ] && { echo "md5 mismatch: exit 1"; exit 1; }
KEYS=$(cat $G/$ARM/keys.txt)
echo "job $SLURM_JOB_ID arm $ARM start $(date +%F_%T)"
echo "srun -n $SLURM_NTASKS $B -i $INP $KEYS -d $G/$ARM -t 00:12:00"
srun -n $SLURM_NTASKS -c 24 $B -i $INP $KEYS -d $G/$ARM -t 00:12:00 > $G/$ARM/run.log 2>&1
rc=$?
echo "rc=$rc $(date +%F_%T)"
grep -q "FATAL" $G/$ARM/run.log && { echo "FATAL in run.log"; exit 3; }
exit $rc
