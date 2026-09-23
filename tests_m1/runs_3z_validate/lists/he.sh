#!/bin/bash -l
#SBATCH -o log.out.%j
#SBATCH -e log.err.%j
#SBATCH -p apu
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=01:00:00
# runs_3z item 4: continue the common restart kehdt_0924/R0/rst/m1slab.00001.rst (t=20,
# written by eddington + be, cfl 0.15) with vet_sc + plm + vimp + <scheme> at <cfl> to t=500.
# usage: sbatch -J he_<s><c> he.sh <h2|be> <cfl>
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
V=/viper/ptmp2/jinma/validate_0923; X=$V/wt/build_gpu/src/athena
RST=/viper/ptmp2/jinma/kehdt_0924/R0/rst/m1slab.00001.rst
D=$V/he3d/$1_c$2; mkdir -p $D; cd $D
srun -n 1 $X -r $RST -i $V/inp/he_vet_$1.athinput -d $D time/cfl_number=$2 time/tlim=500.0 > run.log.$SLURM_JOB_ID 2>&1
echo "rc=$?"
