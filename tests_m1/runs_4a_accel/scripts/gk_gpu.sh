#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/accel_0923/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/accel_0923/log.err.%j
#SBATCH -J acgk
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
# runs_4a_accel: implicit_fast_kernels on vs off (bitwise expected), 3-D He box, 40 cycles,
# hesdirk2 + the other levers, Eddington and vet_sc.  Env: BIN
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/accel_0923
md5sum $BIN
OV="time/nlim=40 output2/dt=2.0 output5/dt=2.0 output3/dt=1.0e9"
L="rad_m1/implicit_vimp_fold=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2"
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
for c in E V; do
  [ $c == V ] && X="$V" || X=""
  for k in false true; do
    d=$W/gpu_gk/${c}_$k; rm -rf $d; mkdir -p $d; cd $d
    srun -n 1 $BIN -i $W/inp/box3d_plm_vimp_h2_x.athinput -d . $OV $L $X rad_m1/implicit_fast_kernels=$k > run.log 2>&1
    echo "gk $c $k rc=$?"
  done
  cd $W/gpu_gk
  for f in $(cd ${c}_true && find . -name "*.hst" -o -name "rt_profile.bin" | sort); do
    cmp -s ${c}_true/$f ${c}_false/$f && echo "SAME $c $f" || echo "DIFF $c $f"
  done
done
