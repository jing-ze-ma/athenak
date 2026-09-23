# sourced by ab*.sub / time*.sub (ck-fast GPU jobs)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
P=/viper/ptmp2/jinma/ckfast_0923
BIN=${BIN:-$P/athena.gpu.fast}
INP=$P/prod_fast.athinput     # = production input + default-valued ck_impl keys
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true \
problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true \
problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true"
TX="problem/ck_impl_tol=1e-10 problem/ck_impl_dtol=1e-10 problem/ck_impl_maxit=20"
CMB="problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
WS="problem/ck_impl_warm=true problem/ck_impl_warm_step=true"
args() {
  case $1 in
    h)    echo "problem/user_srcs=false" ;;
    s)    echo "" ;;
    t4)   echo "$T4" ;;
    t4x)  echo "$T4 $TX" ;;
    o)    echo "$T4 problem/ck_impl_once=true" ;;
    x2|x4|x8) echo "$T4 problem/ck_impl_xstep=${1#x}" ;;
    j)    echo "$T4 problem/ck_impl_jreuse=0.2" ;;
    j5)   echo "$T4 problem/ck_impl_jreuse=0.5" ;;
    j9)   echo "$T4 problem/ck_impl_jreuse=0.9" ;;
    ws)   echo "$T4 $WS" ;;
    x4t)  echo "$T4 problem/ck_impl_xstep=4 problem/ck_impl_xstep_thr=0.01" ;;
    x8t)  echo "$T4 problem/ck_impl_xstep=8 problem/ck_impl_xstep_thr=0.01" ;;
    w)    echo "$T4 problem/ck_impl_warm=true" ;;
    p)    echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.1" ;;
    p1)   echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=1.0" ;;
    p5)   echo "$T4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5" ;;
    dn)   echo "$T4 problem/ck_impl_once=true problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true" ;;
    d2|d4|d8) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=${1#d} problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true" ;;
    ox4)  echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=4" ;;
    ox4p) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=4 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5" ;;
    ox4j) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=4 problem/ck_impl_jreuse=0.2" ;;
    k)    echo "$T4 problem/ck_impl_cvkeep=true" ;;
    ck4|ck8) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=${1#ck} $CMB problem/ck_impl_cvkeep=true" ;;
    ck4t|ck8t) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=${1:2:1} problem/ck_impl_xstep_thr=0.01 $CMB problem/ck_impl_cvkeep=true" ;;
    c8t)  echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=8 problem/ck_impl_xstep_thr=0.01 $CMB" ;;
    n)    echo "$T4 problem/ck_impl_nosync=true" ;;
    c2|c4|c8) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=${1#c} $CMB" ;;
    cn4)  echo "$T4 problem/ck_impl_xstep=4 $CMB" ;;       # the combination without once
    cw4|cw8) echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=${1#cw} $CMB $WS" ;;
    c4t)  echo "$T4 problem/ck_impl_once=true problem/ck_impl_xstep=4 problem/ck_impl_xstep_thr=0.01 $CMB" ;;
    *)    echo "UNKNOWN_ARM_$1" ;;
  esac
}
