#!/bin/bash -l
#SBATCH -J r8_3d
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:14:00
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r8/r8_3d_dev.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r8/r8_3d_dev.err
# GATE 3, the short version: the tests_r6 arm Gnr configuration WITH the gated taper, in
# ONE piece, through the 1.04-1.08 turnover window where every tests_r6 3-D arm died.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r8/Gnr8dev; rm -rf $A; mkdir -p $A; cd $A
srun -n 2 $B/tests_r8/athena_v4 -i $B/inputs/hydro/he4_presn_cs.athinput -t 00:12:00 \
  problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15 \
  problem/mlt_alpha=1.5 problem/vpert=1.0e-3 \
  output1/dt=47.0 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30 \
  output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0 \
  problem/face_budget=0 problem/column_dump=column_he4_Gnr8dev.txt \
  problem/mlt_dump=mltfaces_00.txt time/tlim=7057.5 2>&1 \
  | grep -E 'cycle=|COLLAPSE|FATAL|Terminating|^time=|GATED|TAPERED' | tail -30
rm -rf $A/bin $A/rst $A/cbin_hydro_w_2
