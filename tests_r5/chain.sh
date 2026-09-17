#!/bin/bash -l
#SBATCH -J r5he4
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM, WALL, TLIM, VPERT, MLT come in through --export
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r5/$ARM
mkdir -p $A && cd $A
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=$MLT problem/vpert=$VPERT $EXTRA"
COMMON="time/tlim=$TLIM output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_he4_${ARM}.txt"
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
if [ -n "$RST" ]; then
  echo "### RESTART from $RST  (job $SLURM_JOB_ID)"
  srun -n 2 $B/tests_r4/athena_v2 -r $RST \
       -i $B/inputs/hydro/he4_presn_cs.athinput -t $WALL $COMMON $OV
else
  echo "### FRESH START (job $SLURM_JOB_ID)"
  srun -n 2 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t $WALL $COMMON $OV
fi
echo "### athena exit $?"
