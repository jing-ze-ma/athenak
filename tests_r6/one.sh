#!/bin/bash -l
#SBATCH -J r6one
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
# ARM, TCAP, VPERT, MLT, RSTDT, WALL, MDUMP, RSTF come in through --export
# ONE athena invocation, NO restart chaining (the 0.5-turnover restarts of chain.sh
# destabilised both arms within 400 s of the first restart -- see tests_r6/README).
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r6/$ARM
mkdir -p $A && cd $A
: ${RSTDT:=1.0e30}
: ${WALL:=00:40:00}
: ${MDUMP:=mltfaces_00.txt}
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=$MLT problem/vpert=$VPERT $EXTRA"
COMMON="output1/dt=47.0 output2/dt=1.0e30 output3/dt=$RSTDT
 output4/dt=1.0e30 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_he4_${ARM}.txt problem/mlt_dump=$MDUMP time/tlim=$TCAP"
if [ -n "$RSTF" ]; then
  echo "### HARVEST restart $RSTF -> $MDUMP"
  srun -n 2 $B/tests_r4/athena_v2 -r $RSTF -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t $WALL $COMMON $OV
else
  echo "### SINGLE RUN tlim=$TCAP  rst dt=$RSTDT  job $SLURM_JOB_ID"
  srun -n 2 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput \
       -t $WALL $COMMON $OV
fi
echo "### athena exit $?"
rm -rf $A/bin $A/cbin_hydro_w_2
