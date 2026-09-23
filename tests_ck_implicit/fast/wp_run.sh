#!/bin/bash
# ck-fast well-posed gates (CPU, serial, login node nice 10): A1 (hold the A0b steady
# state, 200 calls), A2 (transient from the IC, 100 calls, vs the wellposed t4x dt = 4 s
# reference), D (eos_h2 = true, mu0 = 1, 100 calls, vs t4x of this binary), all at
# dt = 20 s; C (energy budget) is read from the A1/A2 records.
# usage: wp_run.sh <bin> <root> <test A1|A2|D> <arm> [<arm> ...]   (arms run in parallel)
set -u
BIN=$1; ROOT=$2; TEST=$3; shift 3
F=/viper/ptmp2/jinma/wt_ckfast/tests_ck_implicit/fast
IN=$F/wp_fast.athinput
WP=/viper/ptmp2/jinma/wellposed_0923
GEO="problem/rt_ck=true problem/rt_use_cons=true problem/ck_spherical=true \
problem/ck_beam_sph=true problem/ck_sweep_form=1"
BASE="problem/rt_test_freeze=true problem/rt_test_out=rec.txt output1/dt=1e30 \
output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 time/tlim=1e30"
# T4 as profiled in prof_0923 (arm t8): no chord, no seed, no lw
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true \
problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true \
problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true problem/ck_impl_verbose=true \
problem/ck_impl_reuse_jac=0 problem/ck_impl_seed=0"
TX="problem/ck_impl_tol=1e-10 problem/ck_impl_dtol=1e-10 problem/ck_impl_maxit=20"
# levers 3-5 of the combination (the pred_chk diagnostic makes res/ckdesum/ck_src the
# final state's; it is not part of the timed configuration)
CMB="problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 \
problem/ck_impl_pred_chk=true problem/ck_impl_nosync=true"
arm() {
  case $1 in
    semi)  echo "" ;;
    t4)    echo "$T4" ;;
    t4x)   echo "$T4 $TX" ;;
    x2)    echo "$T4 problem/ck_impl_xstep=2" ;;
    x4)    echo "$T4 problem/ck_impl_xstep=4" ;;
    x8)    echo "$T4 problem/ck_impl_xstep=8" ;;
    x8t)   echo "$T4 problem/ck_impl_xstep=8 problem/ck_impl_xstep_thr=0.02" ;;
    j)     echo "$T4 problem/ck_impl_jreuse=0.2" ;;
    w)     echo "$T4 problem/ck_impl_warm=true" ;;
    n)     echo "$T4 problem/ck_impl_nosync=true" ;;
    x4t)   echo "$T4 problem/ck_impl_xstep=4 problem/ck_impl_xstep_thr=0.01" ;;
    x8t1)  echo "$T4 problem/ck_impl_xstep=8 problem/ck_impl_xstep_thr=0.01" ;;
    j5)    echo "$T4 problem/ck_impl_jreuse=0.5" ;;
    p)     echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 \
problem/ck_impl_pred_chk=true" ;;
    p1)    echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_chk=true" ;;
    p5)    echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 \
problem/ck_impl_pred_chk=true" ;;
    ws)    echo "$T4 problem/ck_impl_warm=true problem/ck_impl_warm_step=true" ;;
    c2|c4|c8) echo "$T4 problem/ck_impl_xstep=${1#c} $CMB" ;;
    c4t|c8t) echo "$T4 problem/ck_impl_xstep=${1:1:1} problem/ck_impl_xstep_thr=0.01 $CMB" ;;
    k)     echo "$T4 problem/ck_impl_cvkeep=true" ;;
    ck2|ck4|ck8) echo "$T4 problem/ck_impl_xstep=${1#ck} $CMB problem/ck_impl_cvkeep=true" ;;
    ck4t|ck8t) echo "$T4 problem/ck_impl_xstep=${1:2:1} problem/ck_impl_xstep_thr=0.01 $CMB \
problem/ck_impl_cvkeep=true" ;;
    cw4)   echo "$T4 problem/ck_impl_xstep=4 $CMB problem/ck_impl_warm=true \
problem/ck_impl_warm_step=true" ;;
    X*)    echo "$T4 ${XARGS[$1]:-}" ;;
    *)     echo "UNKNOWN_ARM_$1" ;;
  esac
}
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
mkdir -p $ROOT/$TEST; cd $ROOT/$TEST
for a in "$@"; do
  d=$ROOT/$TEST/$a; rm -rf $d; mkdir -p $d
  case $TEST in
    A1) R=$(ls -1 $WP/A/a0b_t4/rst/*.rst | tail -1)
        CMD="$BIN -r $R -i $IN $GEO $BASE $(arm $a) problem/rt_test_mu0=0.5 \
problem/rt_test_dt=20 time/nlim=7200 problem/rt_test_every=10" ;;
    A2) CMD="$BIN -i $IN $GEO $BASE $(arm $a) problem/rt_test_mu0=0.5 \
problem/rt_test_dt=20 time/nlim=100 problem/rt_test_every=100" ;;
    D)  CMD="$BIN -i $IN $GEO $BASE $(arm $a) problem/rt_test_mu0=1.0 hydro/eos_h2=true \
problem/rt_test_dt=20 time/nlim=100 problem/rt_test_every=10" ;;
  esac
  ( cd $d && nice -n 10 $CMD > run.log 2>&1; echo "rc=$? $TEST $a nonconv=$(grep -c NOT-CONV run.log)" ) &
done
wait
