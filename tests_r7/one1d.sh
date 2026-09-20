#!/bin/bash -l
#SBATCH -J r7_1d
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
# ARM, TCAP, EXTRA come in through --export
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r7/$ARM
mkdir -p $A && cd $A
: ${TCAP:=6120.0}
# 1-D column: the angular grid collapsed to one 4x4 block per panel, seed off.
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
# arm Gnr's configuration (the floors and vceil are already in the input file)
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=23.5 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30 output5/dt=47.0
 problem/rt_profile_dt=47.0 problem/rt_surface_dt=1.0e30
 problem/column_dump=column_${ARM}.txt problem/mlt_dump=mltfaces_${ARM}.txt"
echo "### ARM $ARM tlim=$TCAP EXTRA='$EXTRA' job $SLURM_JOB_ID"
srun -n 2 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput \
     -t 00:13:00 $ONED $GNR $OUT time/tlim=$TCAP $EXTRA
echo "### athena exit $?"
rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
