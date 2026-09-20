#!/bin/bash -l
#SBATCH -J r16arm
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
# tests_r16: Strang vs unsplit radiation-operator discriminator.  Body copied from
# tests_r15/r15_arm.sh; run directory under tests_r16, mlt_alpha from the environment.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=${XBIN:-$B/tests_r11/athena_v21}
IC=${ICF:-$B/tests_r14/ic_he4_tall_neww.txt}
IN=$B/inputs/hydro/he4_presn_sp.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=${MLTA:-1.5}"
OUT="output1/dt=${HSTDT:-5.0} output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 output5/dt=1.0e30 problem/rt_profile_dt=${PROFDT:-1.0} problem/rt_surface_dt=1.0e30
 problem/face_budget=0 problem/e_ledger=${ELED:-50}"
GRID="problem/rt_top_vacuum=true hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0
 mesh/x1max=2.964625e11 mesh/nx1=144 meshblock/nx1=144 mesh/f_stretch_r_c1=0.992525
 mesh/f_stretch_r_c2=-0.404821 mesh/f_stretch_r_c3=-0.372255
 mesh/f_stretch_r_c4=-1.499759 hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11
 hydro/rad_gate_rho=0.0 problem/rt_rad_force=false problem/rt_force_tau_gate=false
 time/cfl_number=${CFL:-0.15}"
d=$1; shift
A=$B/tests_r16/$d; rm -rf $A; mkdir -p $A; cd $A
echo "=========== ARM $d : X=$X : IC=$IC : tlim=${TLIM:-600.0} : $@"
srun -n 1 $X -i $IN -t ${WALLI:-00:12:00} $ONED $GNR $OUT $GRID time/tlim=${TLIM:-600.0} \
  problem/ic_profile=$IC \
  problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" \
  > $A/full.log 2>&1
grep -E 'cycle=|COLLAPSE|dt is set|FATAL|Terminating|^time=|e_ledger' $A/full.log | tail -30
rm -rf $A/bin $A/cbin_hydro_w_2
echo "=========== ARM $d DONE"
