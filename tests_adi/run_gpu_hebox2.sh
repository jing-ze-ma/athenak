#!/bin/bash -l
#SBATCH -o gpu_hebox2.out.%j
#SBATCH -e gpu_hebox2.err.%j
#SBATCH -J adi_hebox2
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# Is the He box bitwise REPRODUCIBLE at all on this machine?  Same binary, twice.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
I=$PWD/he_box_w8_smoke.athinput
for r in a b; do
  d=heboxr_$r; rm -rf $d; mkdir -p $d
  (cd $d && srun -n 2 /viper/u2/jinma/ATHENAK/bench/wt_he4_adi/tests_adi/athena_head_box_gpu -i $I time/nlim=10 > log.txt 2>&1)
done
echo GPU_HEBOX2_DONE
