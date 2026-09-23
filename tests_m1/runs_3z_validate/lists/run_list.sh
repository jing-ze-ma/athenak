#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/validate_0923/gpu2d/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/validate_0923/gpu2d/log.err.%j
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# usage: sbatch run_list.sh <listfile>; each line: <name> <input> [overrides]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
V=/viper/ptmp2/jinma/validate_0923; X=$V/wt/build_gpu_gen/src/athena
while read -r n inp ov; do
  [ -z "$n" ] && continue
  d=$V/gpu2d/$n; rm -rf $d; mkdir -p $d; cd $d
  srun -n 1 $X -i $inp -d $d $ov > run.log 2>&1 < /dev/null; echo "$n rc=$?"
done < $1
