#!/bin/bash -l
#SBATCH -o log.out.%j
#SBATCH -e log.err.%j
#SBATCH -J h2g5_V30
#SBATCH -p apu
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
# runs_3x_hesdirk2 G5: 3-D seeded He slab, vet_sc full (D^n) + plm + vimp + hesdirk2 (central), cfl 0.3, tlim 1000.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
LOG=run.log.$SLURM_JOB_ID
if [ -f STOP ]; then echo "CHAIN: STOP present"; exit 0; fi
OV="time/cfl_number=0.3 time/tlim=1000.0"
LAST=$(ls -1t rst/*.rst 2>/dev/null | head -1)
if [ -n "$LAST" ]; then
  echo "CHAIN: restarting from $LAST"
  srun -n 1 /viper/ptmp2/jinma/h2_3x/g5/athena_gpu -r "$LAST" -d . -t 01:50:00 $OV > $LOG 2>&1
else
  echo "CHAIN: starting from scratch"
  srun -n 1 /viper/ptmp2/jinma/h2_3x/g5/athena_gpu -i ../he_slab_m1_3d_h2v.athinput -d . -t 01:50:00 $OV > $LOG 2>&1
fi
echo "rc=$?"
if grep -q "Terminating on wall clock limit" "$LOG" && ! grep -q "FATAL" "$LOG"; then
  echo "CHAIN: wall clock limit reached, next link may run"
else
  echo "CHAIN: stopping (tlim reached, or FATAL)"; touch STOP
fi
