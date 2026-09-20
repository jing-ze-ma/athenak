#!/bin/bash -l
#SBATCH -J r14long
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
# tests_r14: a LONG (3-turnover) 1-D sp-column arm of the radiation redesign.
# Same body as tests_r13/r13_arm.sh, but on apu (apudev is for <15 min), with the full
# stdout kept in <arm>/full.log (the grep-filtered tail alone loses a 3-turnover history),
# with the event log (output5) on a FINITE dt, and with the restart cadence and the IC
# file overridable through the environment.
# usage: sbatch [--time=hh:mm:ss] --export=ALL,XBIN=...,ICF=...,RSTDT=...,WALLI=hh:mm:ss \
#               -o <log> r14_long.sh <armname> [overrides ...]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v21}
IC=${ICF:-$B/tests_r14/ic_he4_neww.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
# 1 turnover = 4705 s.  hst + event log every 0.01, rt_profile every 0.25, rst per RSTDT.
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=${RSTDT:-1.0e30} output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=1176.25 problem/rt_surface_dt=470.5
 problem/face_budget=0"
d=$1; shift
A=$B/tests_r14/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : IC=$IC : $@"
srun -n 1 $X -i $IN -t ${WALLI:-00:55:00} $ONED $GNR $OUT time/tlim=14115.0 \
  problem/ic_profile=$IC \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" \
  > $A/full.log 2>&1
grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|TAPERED|GATED|e_ledger' \
  $A/full.log | tail -80
rm -rf $A/bin $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
