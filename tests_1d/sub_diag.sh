#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/diag.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/diag.err
#SBATCH -J he4_diag
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
COM="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0 \
 time/nlim=300 output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30"
run () {
  d=$1; shift
  mkdir -p $B/tests_1d/$d && cd $B/tests_1d/$d
  echo "=========== ARM $d : $@"
  srun -n 1 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/he4_presn_cs.athinput $COM "$@" \
     2>&1 | grep -E 'ARM|cycle=|COLLAPSE|FATAL|WARNING|guard|face budget' | tail -25
}
run d_base
run d_noforce   problem/rt_rad_force=false
run d_nowb      hydro/wellbalance_dynamic=false
run d_nostrang  problem/rt_strang=false
run d_wall      problem/outer_bc=wall
run d_nosponge  problem/vdamp_top_tau=0.0 problem/vdamp_bot_cells=0 problem/sponge=false
run d_norad     problem/rt_grey=false
