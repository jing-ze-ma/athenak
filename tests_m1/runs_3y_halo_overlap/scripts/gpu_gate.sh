#!/bin/bash -l
#SBATCH -J m1int_gpug
#SBATCH -p apudev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:15:00
#SBATCH -o /viper/ptmp2/jinma/m1int_0924/gpu_gate.%j.out
# merge gate (b): defaults on GPU, ref = rt-integration d80c84ca vs new = m1-int, bitwise
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
W=/viper/ptmp2/jinma/m1int_0924; B=/viper/u2/jinma/ATHENAK/bench
md5sum $W/bin/athena_*gpu
g() {  # <name> <exe> <np> args...
  local n=$1 x=$2 np=$3; shift 3
  local d=$W/gpu/$n; rm -rf $d; mkdir -p $d; cd $d
  srun -n $np --ntasks-per-node=2 bash -c 'export ROCR_VISIBLE_DEVICES=$SLURM_LOCALID; exec "$0" "$@"' \
    $x -d $d "$@" > log.txt 2>&1
  echo "$n rc=$? $(date +%T)"
}
for w in ref new; do
  g box1_$w $W/bin/athena_${w}_boxgpu 1 -i $W/gpu/he3d_fast.athinput meshblock/nx3=26 time/nlim=30
  g box2_$w $W/bin/athena_${w}_boxgpu 2 -i $W/gpu/he3d_fast.athinput meshblock/nx3=26 time/nlim=30
  g boxe2_$w $W/bin/athena_${w}_boxgpu 2 -i $W/gpu/he3d_fast.athinput meshblock/nx3=26 time/nlim=30 rad_m1/closure=eddington
  g dh_$w $W/bin/athena_${w}_dhjgpu 2 -r $B/cs_hyd4_prod/rst/dhj.00113.rst \
    -i $B/cs_hyd4_prod/deep_hot_jupiter.athinput time/nlim=873600 output1/dt=1.0 \
    output2/dt=1.0 output3/dt=1.0e20 output4/dt=1.0e20
done
for n in box1 box2 boxe2 dh; do
  for f in $(cd $W/gpu/${n}_ref && find . -name "*.hst" -o -name "*.bin" -o -name "*.rst" | sort); do
    cmp -s $W/gpu/${n}_ref/$f $W/gpu/${n}_new/$f && echo "$n $f IDENTICAL" || echo "$n $f DIFFERS"
  done
  grep -h "NON-CONVERGED\|cycle=" $W/gpu/${n}_new/log.txt | tail -2
done
