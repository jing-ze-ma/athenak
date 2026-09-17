#!/bin/bash -l
#SBATCH -J he4_3d
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM and WALL come in through --export
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_3d/arms/$ARM
mkdir -p $A && cd $A

# arm-specific overrides
case $ARM in
  A) OV="problem/inner_bc=open problem/rt_bottom_flux=false hydro/rad_flux_inner=0.0" ;;
  B) OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15" ;;
  C) OV="$ARMC_OV problem/vpert=1.0e-2" ;;
esac
[ -n "$EXTRA" ] && OV="$OV $EXTRA"

# cadences: hst/rt_profile 47 s, rt_surface 94 s, bin 0.5 turnover, no periodic rst
COMMON="time/tlim=1.4115e4 output1/dt=47.0 output2/dt=2352.5 output3/dt=1.0e30
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_he4_${ARM}.txt"

RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
if [ -n "$RST" ]; then
  echo "### RESTART from $RST  (job $SLURM_JOB_ID)"
  srun -n 2 $B/build_gpu_rg/src/athena -r $RST \
       -i $B/inputs/hydro/he4_presn_cs.athinput -t $WALL $COMMON $OV
else
  echo "### FRESH START (job $SLURM_JOB_ID)"
  srun -n 2 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t $WALL $COMMON $OV
fi
rc=$?
echo "### athena exit $rc"
