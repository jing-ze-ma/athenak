#!/bin/bash -l
#SBATCH -J r11arm
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
# tests_r9 task 2: ONE bisection arm against the tests_r8 x_gate 1-D baseline.
# usage: sbatch -o <log> r9_arm.sh <armname> [extra athinput overrides ...]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=$B/tests_r11/athena_v12
IN=$B/inputs/hydro/he4_presn_cs.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=235.25 problem/rt_surface_dt=470.5
 problem/face_budget=0"
d=$1; shift
A=$B/tests_r12/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : $@"
srun -n 1 $X -i $IN -t 01:00:00 $ONED $GNR $OUT time/tlim=14115.0 \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" 2>&1 \
  | grep -E 'cycle=|COLLAPSE|dt is set|^    r=|FATAL|Terminating|^time=|TAPERED|GATED|e_ledger' \
  | tail -60
rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
