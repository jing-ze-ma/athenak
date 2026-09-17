#!/bin/bash -l
#SBATCH -o gpu_bw.out.%j
#SBATCH -e gpu_bw.err.%j
#SBATCH -J adi_bw_gpu
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# The CARTESIAN BITWISE regression on the GPU: the 2D Gaussian transverse-diffusion test
# (pgen rad_diff2d) through the unmodified head binary and the new one, for both
# transverse solvers.  Nothing this branch adds may move a Cartesian bit.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
cd "$SLURM_SUBMIT_DIR"
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
HEADB=$1
I=$PWD/rad_transverse_gauss_bin.athinput
for v in head new; do
  B=$HEADB
  [ $v = new ] && B=$W/build_gpu_default/src/athena
  for s in sts adi; do
    # 1 rank, the input's single 64x64 MeshBlock
    d=gbw_${v}_${s}; rm -rf $d; mkdir -p $d
    (cd $d && srun -n 1 $B -i $I time/nlim=30 hydro/rad_ang_solver=$s \
        hydro/rad_adi_scheme=lod2 > log.txt 2>&1)
    # 2 ranks, 2x2 MeshBlocks: the ADI line then crosses a block AND a rank boundary,
    # so the partitioned line solve and its ring gather are in the comparison too
    d=gbwm_${v}_${s}; rm -rf $d; mkdir -p $d
    (cd $d && srun -n 2 $B -i $I time/nlim=30 hydro/rad_ang_solver=$s \
        hydro/rad_adi_scheme=lod2 meshblock/nx2=32 meshblock/nx3=32 > log.txt 2>&1)
  done
done
echo GPU_BW_DONE
