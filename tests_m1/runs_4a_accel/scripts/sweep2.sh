#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/accel_0923/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/accel_0923/log.err.%j
#SBATCH -J acsweep2
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# timing arms, 3-D He box, 120 cycles (ms/cycle 20-120).  Env: BIN, ARMS ("name|input|overrides;..."), NREP
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/accel_0923
md5sum $BIN
IFS=';' read -ra AL <<< "$ARMS"
NREP=${NREP:-1}
for r in $(seq 1 $NREP); do
  if [ $((r % 2)) -eq 0 ]; then ORD=$(seq $((${#AL[@]}-1)) -1 0); else ORD=$(seq 0 $((${#AL[@]}-1))); fi
  for x in $ORD; do
    IFS='|' read -r name inp ov <<< "${AL[$x]}"
    d=$W/runs/${name}_$r; rm -rf $d; mkdir -p $d; cd $d
    echo "#### ${name}_$r $(date +%T)"
    srun -n 2 $BIN -i $W/inp/$inp.athinput -d . -t 00:01:40 time/nlim=120 $ov > run.log 2> run.err
    echo "rc=$? $(date +%T)"
  done
done
