#!/bin/bash -l
#SBATCH -J he4_hov
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM (BH|DH) and WALL come in through --export
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_3d/handover/$ARM
mkdir -p $A && cd $A

# THE H1 HANDOVER: wall inner boundary, the luminosity handed to the CONDUCTION wall
# face (not the two-stream's own lower boundary -- with a live blend the two would
# double count), tau blend 20 -> 300 so Rosseland owns the deep interior.
BASE="problem/inner_bc=wall problem/rt_bottom_flux=false
 hydro/rad_flux_inner=1.305278e15 hydro/rad_tau_lo=20 hydro/rad_tau_hi=300
 problem/vpert=1.0e-3"
case $ARM in
  # BH: pure transport, arm B's seed, NO sub-grid MLT flux.
  BH) OV="$BASE problem/mlt_alpha=0.0" ;;
  # DH: arm D's held sub-grid flux on top of the handover.
  DH) OV="$BASE problem/mlt_alpha=1.5 problem/mlt_ramp_down_time=0.0" ;;
esac
[ -n "$EXTRA" ] && OV="$OV $EXTRA"
TL=${TLIM:-4705.0}

# DISK: no bin (output2), no cbin (output4).  hst + event log + rt_profile at 47 s;
# the restart comes from the -t wall-clock dump.
COMMON="time/tlim=${TL} output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0
 problem/rt_surface_dt=94.0 problem/column_dump=column_he4_${ARM}.txt"

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
