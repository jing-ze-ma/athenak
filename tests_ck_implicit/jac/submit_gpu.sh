#!/bin/bash -l
#SBATCH -o /viper/ptmp2/jinma/ckjac_0923/gpu/log.out.%j
#SBATCH -e /viper/ptmp2/jinma/ckjac_0923/gpu/log.err.%j
#SBATCH -J ckjac_gpu
#SBATCH -p apudev
#SBATCH --ntasks=2
#SBATCH --constraint="apu"
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=24
#SBATCH --time=00:14:50
# ck_impl_esc / ck_impl_aa GPU cost (as glob/submit_gpu.sh): semi vs T4 vs each lever
# from the hyd4 production restart dhj.00127.rst (rot 63.50, t = 1.93675e7 s, ncycle
# 986896, dt 19.47 s), READ IN PLACE; ONE binary (athena.gpu.jac), 2 ranks / 2 GPUs, 150
# cycles, arms interleaved r1 / reversed r2 / r3.  Last run: the first lever arm with
# the BASE binary (HEAD be02c647 + hooks), whose output must be bitwise the jac r1's.
# usage: SET=A (s t4 aa aae, base t4) | SET=B (s t4 sub esc, base sub) sbatch submit_gpu.sh
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
P=/viper/ptmp2/jinma/ckjac_0923
BIN=$P/athena.gpu.jac
BBASE=$P/athena.gpu.base
RST=/viper/u2/jinma/ATHENAK/bench/cs_hyd4_prod/rst/dhj.00127.rst
INP=$P/gpu/hyd.athinput
N0=986896
NCYC=${NCYC:-150}
N2=$((N0 + NCYC))
IMP="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true"
T3="$IMP problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_debug=-2"
T4="$T3 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true"
SUB="problem/ck_impl_glob=ls_sub problem/ck_impl_sub_max=32"
args() {
  case $1 in
    s) echo "problem/ck_impl_debug=-2" ;;
    t4) echo "$T4" ;;
    sub) echo "$T4 $SUB" ;;
    esc) echo "$T4 $SUB problem/ck_impl_esc=1" ;;
    aa) echo "$T4 problem/ck_impl_aa=4" ;;
    aae) echo "$T4 $SUB problem/ck_impl_esc=1 problem/ck_impl_aa=4" ;;
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
case ${SET:-A} in
  A) L1="s t4 aa aae"; L2="aae aa t4 s"; BA=t4 ;;
  B) L1="s t4 sub esc"; L2="esc sub t4 s"; BA=sub ;;
esac
echo "SET=${SET:-A} RST=$RST md5=$(md5sum $RST | cut -c1-32) BIN md5=$(md5sum $BIN |
      cut -c1-32) N0=$N0 NCYC=$NCYC HSA_XNACK=$HSA_XNACK" \
     "HSA_NO_SCRATCH_RECLAIM=$HSA_NO_SCRATCH_RECLAIM"
for a in $L1; do arm $a ${a}_r1; done
for a in $L2; do arm $a ${a}_r2; done
for a in $L1; do arm $a ${a}_r3; done
arm $BA ${BA}base_r1 $BBASE
python3 /viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed/cmpdir.py \
    $P/gpu/hyd/${BA}_r1 $P/gpu/hyd/${BA}base_r1
echo JAC_GPU_DONE
