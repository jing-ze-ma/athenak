#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/g34.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_1d/g34.err
#SBATCH -J he4_g34
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
I=$B/inputs/hydro/he4_presn_cs.athinput
COM="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0 \
 problem/rt_rad_force=false output1/dt=1.0e-30 output2/dt=1.0e30 output4/dt=1.0e30"
# --- gate (iii): mode 3 vs rt_implicit_column = 0, 500 cycles, emergent flux
for M in 3 0; do
  d=$B/tests_1d/g3_m$M; mkdir -p $d; cd $d
  echo "======== gate(iii) rt_implicit_column=$M"
  srun -n 1 $B/build_gpu_rg/src/athena -i $I $COM time/nlim=500 \
    problem/rt_implicit_column=$M output3/dt=1.0e30 2>&1 | grep -E 'face budget|cycle=500|FATAL' | tail -3
done
# --- gate (iv): restart, 200 + 50 vs 250 straight
d=$B/tests_1d/g4_cont; mkdir -p $d; cd $d
echo "======== gate(iv) continuous 250"
srun -n 1 $B/build_gpu_rg/src/athena -i $I $COM time/nlim=250 output3/dt=1.0e30 \
   2>&1 | tail -2
d=$B/tests_1d/g4_rst; mkdir -p $d; cd $d
echo "======== gate(iv) 200 then restart +50"
srun -n 1 $B/build_gpu_rg/src/athena -i $I $COM time/nlim=200 output3/dt=1.0e30 \
   2>&1 | tail -2
ls -la $d/rst
R=$(ls -t $d/rst/*.rst 2>/dev/null | head -1)
echo "restart file: $R"
srun -n 1 $B/build_gpu_rg/src/athena -r $R time/nlim=250 2>&1 | tail -2
