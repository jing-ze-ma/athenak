#!/bin/bash -l
#SBATCH -J he4w6d
#SBATCH -p apu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:30:00
# ONSET REPLAY of tests_r11/wedge6: restart from the kept rst at t = 11280 s and dump the
# 3-D primitives every 23.5 s up to t = 11760, through the dt collapse at t = 11714.6.
# Restart + cbin output are OFF.  Nothing in wedge6/ is touched.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ARM=wedge6d; A=$B/tests_r11/$ARM; mkdir -p $A; cd $A
RST=$B/tests_r11/wedge6/rst_keep/he4.00012.rst
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5 problem/vpert=1.0e-3"
COMMON="output1/dt=47.0 output2/dt=23.5 output2/last_time=11280.0 output3/dt=470.0
 output4/dt=1.0e30 output5/dt=23.5 problem/rt_profile_dt=23.5
 problem/rt_surface_dt=1.0e30
 problem/column_dump=column_$ARM.txt problem/mlt_dump=mltfaces_$ARM.txt"
NT=$((SLURM_JOB_NUM_NODES*2))
echo "### nodes $SLURM_JOB_NUM_NODES ranks $NT restart $RST"
srun -n $NT $B/tests_r11/athena_v21 -r $RST -i $B/inputs/hydro/he4_presn_sp.athinput \
     -t 00:20:00 $COMMON $OV time/tlim=12100 \
     time/cfl_number=0.15 problem/rt_top_vacuum=true problem/mlt_split_deposit=true \
     mesh/nx2=192 mesh/nx3=192 meshblock/nx2=24 meshblock/nx3=24 \
     problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true \
     problem/vdamp_all_mean_only=true problem/vdamp_all_until=18820.0 \
     problem/vdamp_all_time=500.0 \
     hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11 \
     hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0 \
     problem/rt_force_tau_gate=true hydro/rad_gate_rho=1.0e-11 \
     problem/ic_profile=$B/tests_r14/ic_he4_tall_neww.txt \
     mesh/x1max=2.964625e11 mesh/nx1=144 meshblock/nx1=144 \
     mesh/f_stretch_r_c1=0.992525 mesh/f_stretch_r_c2=-0.404821 \
     mesh/f_stretch_r_c3=-0.372255 mesh/f_stretch_r_c4=-1.499759
echo "### athena exit $?"
