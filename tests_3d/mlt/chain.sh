#!/bin/bash -l
#SBATCH -J he4_mlt
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM (D|E), WALL and TLIM come in through --export
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_3d/mlt/$ARM
mkdir -p $A && cd $A

# arm B's inner boundary (wall + the base flux on the two-stream), which tests_3d/arms
# settled on, plus the sub-grid MLT flux at the alpha the ic was built with.
BASE="problem/inner_bc=wall problem/rt_bottom_flux=true
 hydro/rad_flux_inner=1.305278e15 problem/mlt_alpha=1.5"
case $ARM in
  # D: held sub-grid flux, one turnover, the arm-B seed.  Does the MLT flux close the
  #    budget?
  D) OV="$BASE problem/vpert=1.0e-3 problem/mlt_ramp_down_time=0.0" ;;
  # E: stronger seed, flux held for one turnover then withdrawn over three, five
  #    turnovers.  Does resolved convection take the load?
  E) OV="$BASE problem/vpert=1.0e-2 problem/mlt_hold_time=4705.0
       problem/mlt_ramp_down_time=1.4115e4" ;;
esac
[ -n "$EXTRA" ] && OV="$OV $EXTRA"

COMMON="time/tlim=${TLIM} output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_he4_${ARM}.txt problem/mlt_dump=mltfaces_${ARM}.txt"

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
