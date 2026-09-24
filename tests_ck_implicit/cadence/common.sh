# sourced by ab.sub / time.sub (ck-cadence GPU jobs)
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
export HSA_XNACK=1
export HSA_NO_SCRATCH_RECLAIM=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
P=/viper/ptmp2/jinma/ckcad_0924
BIN=${BIN:-$P/athena.gpu.cad}
INP=/viper/ptmp2/jinma/wt_ckcad/tests_ck_implicit/cadence/prod_cad.athinput  # + ck_impl_every keys
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true \
problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true \
problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true"
CMB="problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
C2="$T4 problem/ck_impl_once=true problem/ck_impl_xstep=2 $CMB"
# arms: h (hydro only), c2 (reference), eN (every N, no guard), eNgX (guard thr X/100)
args() {
  case $1 in
    h)  echo "problem/user_srcs=false" ;;
    c2) echo "$C2" ;;
    e[0-9]*g[0-9]*) n=${1#e}; n=${n%g*}; g=${1#*g}
        echo "$C2 problem/ck_impl_every=$n problem/ck_impl_every_thr=$(awk "BEGIN{print $g/100}")" ;;
    e[0-9]*) echo "$C2 problem/ck_impl_every=${1#e}" ;;
    *)  echo "UNKNOWN_ARM_$1" ;;
  esac
}
