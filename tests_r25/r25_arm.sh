#!/bin/bash -l
#SBATCH -J r25arm
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
# tests_r25: <hydro>/vceil_thermalise.  Body copied from tests_r23/r23_arm.sh; only
# the run directory (tests_r25) and the face-budget cadence (FACEB) differ.
# usage: sbatch --export=ALL,XBIN=...,TLIM=...,FACEB=... -o <log> r25_arm.sh <arm> [ovr]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v24}
IC=${ICF:-$B/tests_r14/ic_he4_tall_neww.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
# 1 turnover = 4705 s.
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=${PROFDT:-235.0} problem/rt_surface_dt=470.5
 problem/face_budget=${FACEB:-0}"
d=$1; shift
A=$B/tests_r25/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : tlim=${TLIM:-300.0} : $@"
srun -n 1 $X -i $IN -t ${WALLI:-00:13:00} $ONED $GNR $OUT time/tlim=${TLIM:-300.0} \
  problem/ic_profile=$IC \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" \
  > $A/full.log 2>&1
grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|face budget' \
  $A/full.log | tail -20
rm -rf $A/bin $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
