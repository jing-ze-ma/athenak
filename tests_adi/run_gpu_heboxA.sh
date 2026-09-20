#!/bin/bash -l
#SBATCH -o gpu_heboxA.out.%j
#SBATCH -e gpu_heboxA.err.%j
#SBATCH -J adi_heboxA
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# The PRODUCTION Cartesian configuration -- the He-star FeCZ box (box_w8): general
# tabulated EOS, rad_implicit_x1, rad_implicit_ang with rad_ang_solver = adi and
# rad_adi_scheme = lod2, 104/52 = 2 MeshBlocks per ADI line (so the partitioned line
# solve and its reduced interface system run), on 2 ranks.  nlim = 50.
# Nothing this branch adds may move a bit of it.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
I=$PWD/he_box_w8_smoke.athinput
for v in head new; do
  B=$PWD/athena_head_box_gpu
  [ $v = new ] && B=$PWD/athena_bisA_box_gpu
  d=heboxA_$v; rm -rf $d; mkdir -p $d
  (cd $d && srun -n 2 $B -i $I time/nlim=30 > log.txt 2>&1)
done
echo GPU_HEBOXA_DONE
