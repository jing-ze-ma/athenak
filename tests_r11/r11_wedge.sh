#!/bin/bash -l
#SBATCH -J he4wedge
#SBATCH -p apu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
# tests_r11: the 90x90 deg spherical-polar wedge, 96x256x256, recipe of w3d, toward 5
# turnovers; resubmitting the same command chains from the newest restart file.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
ARM=$1; shift; A=$B/tests_r11/$ARM; mkdir -p $A; cd $A
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5 problem/vpert=1.0e-3"
COMMON="output1/dt=47.0 output2/dt=2352.5 output3/dt=940.0 output4/dt=1.0e30
 output5/dt=47.0 problem/rt_profile_dt=47.0 problem/rt_surface_dt=94.0
 problem/column_dump=column_$ARM.txt problem/mlt_dump=mltfaces_$ARM.txt"
RST=$(ls -t $A/rst/*.rst 2>/dev/null | head -1)
srun -n 8 $B/tests_r11/athena_v16 ${RST:+-r $RST} -i $B/inputs/hydro/he4_presn_sp.athinput \
     -t 05:50:00 $COMMON $OV time/tlim=23525 "$@"
echo "### athena exit $?"
ls -t $A/rst/*.rst 2>/dev/null | tail -n +2 | xargs -r rm -f
# the 3-D dumps (bin every 0.5 turnover) are KEPT: they are the convection diagnostics
