#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/ckglob_0924/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/ckglob_0924/gpu/log.err.%j
#SBATCH -J ckglob_gpu
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# ck_impl_glob GPU cost (README_T4_prod4.md section 3 layout): semi vs T4 vs T4+ls vs
# T4+ls_sub from the hyd4 production restart dhj.00113.rst (rot 56.50, t = 1.72325e7 s,
# ncycle 873580, dt 18.7 s), READ IN PLACE; ONE binary (athena.gpu.glob), 2 ranks / 2
# GPUs, 150 cycles, arms interleaved s t4 ls sub | sub ls t4 s | s t4 ls sub.  Last run:
# t4 with the BASE binary (HEAD d80c84ca), whose final restart must be bitwise t4_r1's.
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
P=/viper/ptmp2/jinma/ckglob_0924
BIN=$P/athena.gpu.glob
BBASE=$P/athena.gpu.base
RST=/viper/u2/jinma/ATHENAK/bench/cs_hyd4_prod/rst/dhj.00113.rst
INP=$P/gpu/hyd.athinput
N0=873580
NCYC=${NCYC:-150}
N2=$((N0 + NCYC))
IMP="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true"
T3="$IMP problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_debug=-2"
T4="$T3 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true"
args() {
  case $1 in
    s) echo "problem/ck_impl_debug=-2" ;;
    t4) echo "$T4" ;;
    ls) echo "$T4 problem/ck_impl_glob=ls" ;;
    sub) echo "$T4 problem/ck_impl_glob=ls_sub" ;;
  esac
}
arm() {
  local a=$1; local d=$2; local b=${3:-$BIN}
  rm -rf "$P/gpu/hyd/$d"; mkdir -p "$P/gpu/hyd/$d"; cd "$P/gpu/hyd/$d"
  /usr/bin/time -f "WALL_$d %e s" srun -n 2 $b \
      -r $RST -i $INP -t 00:03:00 output3/dt=1e30 \
      time/nlim=$N2 $(args $a) > run.log 2>run.err
  echo "rc_$d=$?"
  grep -E "WALL_|cpu time used" run.err run.log | tail -2
}
echo "RST=$RST md5=$(md5sum $RST | cut -c1-32) BIN md5=$(md5sum $BIN | cut -c1-32)" \
     "N0=$N0 NCYC=$NCYC HSA_XNACK=$HSA_XNACK HSA_NO_SCRATCH_RECLAIM=$HSA_NO_SCRATCH_RECLAIM"
for a in s t4 ls sub; do arm $a ${a}_r1; done
for a in sub ls t4 s; do arm $a ${a}_r2; done
for a in s t4 ls sub; do arm $a ${a}_r3; done
arm t4 t4base_r1 $BBASE
python3 /viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed/cmpdir.py \
    $P/gpu/hyd/t4_r1 $P/gpu/hyd/t4base_r1
echo GLOB_GPU_DONE
