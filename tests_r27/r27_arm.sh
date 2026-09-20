#!/bin/bash -l
#SBATCH -J r27arm
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
# tests_r27: 1-D validation of the DEEPER INNER WALL grid (r_in = 0.400 R).
# Body copied from tests_r23/r23_arm.sh; the grid, the inner flux and the bottom
# sponge are the tests_r27 ones, and mlt_alpha = 0 / vceil_thermalise are baked in
# so the arm is directly comparable with tests_r23/V100 (the OLD-grid reference).
# usage: sbatch [--export=ALL,XBIN=..,ICF=..,TLIM=..,VDB=..] -o <log> r27_arm.sh <arm> [ov ...]
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v24}
IC=${ICF:-$B/tests_r23/ic_rad_f100.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
# THE NEW GRID: inner wall at 0.400 R, 216 radial cells in one MeshBlock, stretch
# refitted by tests_r27/fit2.py; rad_flux_inner = L/(4 pi r_in^2).
GRID="mesh/x1min=9.486800e10 mesh/x1max=2.964625e11 mesh/nx1=216 meshblock/nx1=216
 mesh/f_stretch_r_c1=-0.391990 mesh/f_stretch_r_c2=+0.936503
 mesh/f_stretch_r_c3=+4.196941 mesh/f_stretch_r_c4=-6.428076
 hydro/rad_flux_inner=2.039496e15 problem/vdamp_bot_cells=${VDB:-74}"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true problem/mlt_alpha=0.0
 hydro/vceil_thermalise=true time/cfl_number=0.3 problem/rt_top_vacuum=true
 hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0
 hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11
 hydro/rad_gate_rho=0.0 problem/rt_rad_force=false
 problem/rt_force_tau_gate=false problem/e_ledger=100"
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=${RSTDT:-1.0e30} output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=${PROFDT:-235.0} problem/rt_surface_dt=470.5
 problem/face_budget=0"
d=$1; shift
A=$B/tests_r27/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : IC=$IC : tlim=${TLIM:-7058.0} : VDB=${VDB:-74} : $@"
srun -n 1 $X -i $IN -t ${WALLI:-00:13:00} $ONED $GNR $GRID $OUT \
  time/tlim=${TLIM:-7058.0} problem/ic_profile=$IC \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" \
  > $A/full.log 2>&1
grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|TAPERED|GATED|e_ledger' \
  $A/full.log | tail -40
rm -rf $A/bin $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
