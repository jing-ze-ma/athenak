#!/bin/bash -l
#SBATCH -J r17arm
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
# tests_r17: 1-D discriminator for the standing inflow at the IC density cliff
# (physical i ~ 99/100, r/R ~ 1.025) in the tall He4 column.  Same body as
# tests_r14/r14_long.sh, but the run directory is under tests_r17, and tlim /
# rt_profile_dt / partition come from the environment.
# usage: sbatch [-p apudev] [--time=hh:mm:ss] \
#          --export=ALL,XBIN=...,ICF=...,TLIM=...,PROFDT=...,WALLI=hh:mm:ss \
#          -o <log> r15_arm.sh <armname> [overrides ...]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v21}
IC=${ICF:-$B/tests_r14/ic_he4_tall_neww.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
# 1 turnover = 4705 s.
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=${RSTDT:-1.0e30} output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=${PROFDT:-235.0} problem/rt_surface_dt=470.5
 problem/face_budget=0"
d=$1; shift
A=$B/tests_r17/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : IC=$IC : tlim=${TLIM:-7058.0} : $@"
srun -n 1 $X -i $IN -t ${WALLI:-00:28:00} $ONED $GNR $OUT time/tlim=${TLIM:-7058.0} \
  problem/ic_profile=$IC \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" \
  > $A/full.log 2>&1
grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|TAPERED|GATED|e_ledger' \
  $A/full.log | tail -40
rm -rf $A/bin $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
