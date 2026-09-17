#!/bin/bash -l
#SBATCH -o /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r10/he4_relax/relax.log
#SBATCH -e /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r10/he4_relax/relax.err
#SBATCH -J he4relax
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
B=/viper/u2/jinma/ATHENAK/bench/wt_he4
X=$B/build_gpu_rg/src/athena
IN=$B/inputs/hydro/he4_presn_cs.athinput
ONED="mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0"
GNR="problem/inner_bc=wall problem/rt_bottom_flux=true hydro/rad_flux_inner=1.305278e15 problem/mlt_alpha=1.5"
OUT="output1/dt=1.0e30 output2/dt=1.0e30 output3/dt=1.0e30 output4/dt=1.0e30 output5/dt=1.0e30 problem/rt_profile_dt=1.0e30 problem/rt_surface_dt=100.0 problem/face_budget=0"
cd /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r10/he4_relax
rm -rf rst
srun -n 1 $X -i $IN $ONED $GNR $OUT time/nlim=200 \
  problem/mlt_dump=mltfaces_t0.txt \
  < /dev/null > run_stage1.log 2>&1
echo "stage1 exit $?"; tail -5 run_stage1.log
RST=$(ls -t rst/*.rst | head -1)
echo "restart from $RST"
srun -n 1 $X -r $RST problem/mlt_dump=mltfaces_relax.txt time/nlim=201 \
  < /dev/null > run_stage2.log 2>&1
echo "stage2 exit $?"; tail -5 run_stage2.log
rm -rf rst *.bin *.cbin*
