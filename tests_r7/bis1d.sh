#!/bin/bash -l
#SBATCH -J r7_bis
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/bis.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r7/bis.err
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5"
run () {
  d=$1; shift
  A=$B/tests_r7/$d; mkdir -p $A && cd $A
  OUT="output1/dt=23.5 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
   output5/dt=1.0e30 problem/rt_profile_dt=94.0 problem/rt_surface_dt=1.0e30
   problem/column_dump=column_${d}.txt problem/mlt_dump=mltfaces_${d}.txt
   problem/face_budget=0"
  echo "=========== ARM $d : $@"
  srun -n 1 $B/tests_r4/athena_v2 -i $B/inputs/hydro/he4_presn_cs.athinput \
    -t 00:01:30 $ONED $GNR $OUT time/tlim=6120.0 "$@" 2>&1 \
    | grep -E 'ARM|cycle=|COLLAPSE|dt is set|^    r=|FATAL|Terminating' | tail -22
  rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
}
run b_mlt0      problem/mlt_alpha=0.0
run c_noforce   problem/rt_rad_force=false
run d_nostrang  problem/rt_strang=false
run e_noimplx1  hydro/rad_implicit_x1=false
run e2_nouform  hydro/rad_x1_uform=false
run f_noang     hydro/rad_implicit_ang=false
run g_percol    problem/mlt_mean=false
