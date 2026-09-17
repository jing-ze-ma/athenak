#!/bin/bash -l
#SBATCH -J r6he4
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM, NLEG, LPJ, VPERT, MLT, EXTRA come in through --export
# LEGS: the run is cut into 0.5-turnover legs so that problem/mlt_dump (a ONE-SHOT per
# process) is written at the start of every leg -> a Gamma_rad time series.  Leg k runs
# from (k-1)*0.5 to k*0.5 turnover and dumps mltfaces_<k-1 in halves>.txt at its start.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r6/$ARM
mkdir -p $A && cd $A
TURN=4705.0
HALF=2352.5
: ${LPJ:=3}
: ${TCAP:=23525.0}
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=$MLT problem/vpert=$VPERT $EXTRA"
COMMON="output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_he4_${ARM}.txt"
# which leg is next: one marker file per completed leg
K=1
while [ -f $A/leg_$(printf %02d $K).done ]; do K=$((K+1)); done
NRUN=0
while [ $K -le $NLEG ] && [ $NRUN -lt $LPJ ]; do
  TAG=$(printf %02d $K)
  TLIM=$(python3 -c "print(min($K*$HALF, $TCAP+47.0))")
  MDUMP="mltfaces_$(printf %02d $((K-1))).txt"
  RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
  echo "### LEG $TAG  tlim=$TLIM  mlt_dump=$MDUMP  job $SLURM_JOB_ID"
  if [ -n "$RST" ]; then
    srun -n 2 $B/tests_r4/athena_v2 -r $RST \
         -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:04:00 \
         $COMMON $OV time/tlim=$TLIM problem/mlt_dump=$MDUMP
  else
    srun -n 2 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput \
         -t 00:04:00 $COMMON $OV time/tlim=$TLIM problem/mlt_dump=$MDUMP
  fi
  EX=$?
  echo "### athena exit $EX (leg $TAG)"
  # keep only the newest restart; the .bin/.cbin outputs are off
  ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
  rm -rf $A/bin $A/cbin
  if [ $EX -ne 0 ]; then echo "### STOP: leg $TAG failed"; exit $EX; fi
  touch $A/leg_$TAG.done
  K=$((K+1)); NRUN=$((NRUN+1))
done
echo "### JOB DONE, next leg $K of $NLEG"
