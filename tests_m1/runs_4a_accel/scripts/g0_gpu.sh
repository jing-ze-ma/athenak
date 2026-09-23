#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/accel_0923/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/accel_0923/log.err.%j
#SBATCH -J acg0
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# runs_4a_accel G0 on the GPU: switches off, bitwise vs d0c59f7c (bin_base_gpu).  3-D He box,
# plm + vimp, 40 cycles, be / hesdirk2 / hesdirk2 + vet_sc.  Env: BIN (new binary)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/accel_0923
md5sum $W/bin_base_gpu $BIN
OV="time/nlim=40 output2/dt=2.0 output5/dt=2.0 output3/dt=1.0e9"
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
for c in be h2 hv; do
  case $c in be) I=box3d_plm_vimp; X="";; h2) I=box3d_plm_vimp_h2; X="";; hv) I=box3d_plm_vimp_h2; X="$V";; esac
  for b in base new; do
    [ $b == base ] && E=$W/bin_base_gpu || E=$BIN
    d=$W/gpu_g0/${c}_$b; rm -rf $d; mkdir -p $d; cd $d
    srun -n 1 $E -i $W/inp/$I.athinput -d . $OV $X > run.log 2>&1
    echo "g0 $c $b rc=$?"
  done
  cd $W/gpu_g0
  for f in $(cd ${c}_new && find . -name "*.bin" -o -name "*.hst" -o -name "rt_profile.bin" | sort); do
    cmp -s ${c}_new/$f ${c}_base/$f && echo "SAME $c $f" || echo "DIFF $c $f"
  done
done
