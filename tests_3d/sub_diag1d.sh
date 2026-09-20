#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_3d/diag1d.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_3d/diag1d.err
#SBATCH -J he4_d1d
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
# every arm: 1-D, vpert 0, inner wall + base flux, 120 cycles, hst every 2 s
COMMON="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0
 problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 output1/dt=2.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
 problem/rt_surface_dt=1.0e30 problem/rt_profile_dt=1.0e30 time/nlim=120"
run () {  # $1 = arm dir, rest = overrides
  local d=$1; shift
  mkdir -p $B/tests_3d/diag/$d && cd $B/tests_3d/diag/$d
  echo "########## ARM $d : $* ##########"
  srun -n 2 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/he4_presn_cs.athinput \
       $COMMON "$@" 2>&1 | tail -40
}
run a_base
run a_m0        problem/rt_implicit_column=0
run a_m3_plain  problem/rt_impl_mixed=0 problem/rt_impl_warm=0
run a_noforce   problem/rt_rad_force=false
run a_noang     hydro/rad_implicit_ang=false
run a_m0_noforce problem/rt_implicit_column=0 problem/rt_rad_force=false
