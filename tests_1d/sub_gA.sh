#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/gA.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/gA.err
#SBATCH -J he4_gA
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
cd $B/tests_1d/gA
srun -n 2 $B/build_gpu_rg/src/athena -i $B/inputs/hydro/he4_presn_cs.athinput \
  mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 \
  problem/vpert=0.0 \
  output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30
