#!/bin/bash -l
#SBATCH -J r11smsp
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
# smoke of the 64^2-per-panel layout (96 MeshBlocks of 96x16x16), 40 cycles, no outputs
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
A=$B/tests_r11/smokeseed; rm -rf $A; mkdir -p $A; cd $A
OV="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15
 problem/mlt_alpha=1.5 problem/vpert=1.0e-3"
srun -n 2 $B/tests_r11/athena_v18 -i $B/inputs/hydro/he4_presn_sp.athinput $OV \
  time/nlim=20 problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true output2/dt=30.0 mesh/nx2=192 mesh/nx3=192 meshblock/nx2=24 meshblock/nx3=24 output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 \
  output4/dt=1e30 output5/dt=1e30 problem/rt_profile_dt=1e30 problem/rt_surface_dt=1e30 \
  time/cfl_number=0.15 problem/rt_top_vacuum=true problem/mlt_split_deposit=true \
  problem/cs_max=1.0e8 2>&1 | grep -E 'cycle=|FATAL|periodic random|ADI|transverse|trs|Conduction|zone-cycles|cpu time|MeshBlocks|Error|error' | tail -25
echo "### exit $?"
