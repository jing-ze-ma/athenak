#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/h2_3x/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/h2_3x/gpu/log.err.%j
#SBATCH -J h2g0
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# G0 on the GPU (runs_3x_hesdirk2): time_scheme default (be) of the new build vs the
# 2c3c4178 build, 3-D He box (plm + vimp) 40 cycles, hst + bin + rst compared with cmp.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/h2_3x
OV="time/nlim=40 output2/dt=2.0 output5/dt=2.0 output3/dt=1.0e9"
for b in base new; do
  md5sum $W/$b/build_gpu/src/athena
  d=$W/gpu/g0_$b; rm -rf $d; mkdir -p $d; cd $d
  srun -n 1 $W/$b/build_gpu/src/athena -i $W/inp/box3d_plm_vimp.athinput -d . $OV > run.log 2>&1
  echo "$b rc=$?"
done
cd $W/gpu
for f in $(cd g0_new && find . -name "*.bin" -o -name "*.hst" | sort); do
  cmp -s g0_new/$f g0_base/$f && echo "SAME $f" || echo "DIFF $f"
done
