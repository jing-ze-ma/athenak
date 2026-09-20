#!/bin/bash -l
#SBATCH -J r8_1dL
#SBATCH -p apu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r8/r8_1d_long.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r8/r8_1d_long.err
# GATE 2 of tests_r8: the 1-D x_base configuration of tests_r7 section 3, WITH the
# temperature-gated taper and rt_rad_force = true.  tests_r7's x_base dies at t = 5541 s
# (1.18 turnover) and only rt_rad_force = false survived (> 2.42); this must survive to 5.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=$B/tests_r8/athena_v4
IN=$B/inputs/hydro/he4_presn_cs.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
OUT="output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=235.25 problem/rt_surface_dt=470.5
 problem/face_budget=0"
run () {
  d=$1; shift
  A=$B/tests_r8/$d; rm -rf $A; mkdir -p $A; cd $A
  echo "=========== ARM $d : $@"
  srun -n 1 $X -i $IN -t 02:50:00 $ONED $GNR $OUT time/tlim=23525.0 \
    problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt "$@" 2>&1 \
    | grep -E 'cycle=|COLLAPSE|dt is set|^    r=|FATAL|Terminating|^time=|TAPERED|GATED' \
    | tail -25
  rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
}
# the gate ON (the deliverable), and the same thing with the force off as the control
run x_gate5
