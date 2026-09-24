#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/thinstab_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/thinstab_0924/gpu/log.err.%j
#SBATCH -J thinstab_gpu
#SBATCH -p apudev
#SBATCH --ntasks=1
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:00
# runs_5c_thinstab GPU gate: (1) switch off = bitwise the base binary (s0 = rt-integration
# 304fb38f + the pgen seed, s1 = m1-thinstab), (2) the cure on the GPU, (3) the cost of
# the switch vs off vs closure_lag = pass + anderson, same binary, interleaved, 2 repeats.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/thinstab_0924; G=$W/gpu; I=$W/inp/uni_ctr_rst.athinput
RW=$W/inp/rw_gpu.athinput
run() {  # run <dir> <s0|s1> <args...>
  local t=$1 c=$2; shift 2
  local d=$G/$t; rm -rf $d; mkdir -p $d; cd $d
  echo "== $t $c $(date +%T)"
  local t0=$(date +%s.%N)
  srun -n 1 $W/bin/athena_${c}_gpu -d $d "$@" > log.txt 2>&1
  echo "rc=$? wall=$(echo "$(date +%s.%N) - $t0" | bc)" >> log.txt
}
IMPL="rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"
for c in s0 s1; do
  run off_uni_$c $c -i $I time/nlim=60
  run off_rwxy_$c $c -i $RW mesh/nx1=64 mesh/nx2=64 meshblock/nx1=64 meshblock/nx2=64 problem/radwave_dir=xy $IMPL time/tlim=0.8374357893586236 output1/dt=0.0261698684
done
run on_uni_k1 s1 -i $I time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5
run on_uni_k025 s1 -i $I time/nlim=500 output3/dcycle=25 rad_m1/implicit_closure_thin_relax=1.5 rad_m1/kappa_s=0.25
T="mesh/nx1=128 meshblock/nx1=128 mesh/nx2=64 meshblock/nx2=64 mesh/x2max=2.0 rad_m1/kappa_s=4.0 time/nlim=200 output3/dcycle=1000"
for r in 1 2; do
  run tim_off_$r s1 -i $I $T
  run tim_ctr_$r s1 -i $I $T rad_m1/implicit_closure_thin_relax=1.5
  run tim_and_$r s1 -i $I $T rad_m1/implicit_closure_lag=pass rad_m1/implicit_accel=anderson
done
for t in off_uni off_rwxy; do
  n=0; m=0
  for f in $(cd $G/${t}_s0 && find . -type f ! -name log.txt); do
    m=$((m+1)); cmp -s $G/${t}_s0/$f $G/${t}_s1/$f || n=$((n+1))
  done
  echo "BITWISE $t: $m files, $n differ"
done
for d in $G/tim_*; do
  echo "TIME $(basename $d) $(tail -n1 $d/log.txt) $(grep -h 'cpu time used' $d/log.txt) $(grep -h 'Picard iterations' $d/log.txt | sed 's/.*Picard/Picard/')"
done
echo GPU GATE DONE
