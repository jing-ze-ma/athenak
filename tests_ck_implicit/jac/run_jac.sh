#!/bin/bash -l
# ck_impl_esc / ck_impl_aa gates on the well-posed suite (tests_ck_implicit/wellposed/),
# CPU, login node nice 10 (correctness only, no timing).  Binaries: athena.cpu.base (HEAD
# be02c647 + wellposed_hooks.patch) and athena.cpu.jac (ck-jac + the same patch).
# Arms: t4 (as wellposed/arms.sh, tol 1e-8, maxit 8), sub (ls_sub, sub_max 32 = HEAD),
# esc (sub + ck_impl_esc 1), aa (t4 + ck_impl_aa 4), aas (sub + aa 4), aae (esc + aa 4),
# jn (t4 + ck_impl_jneg).  usage: run_jac.sh G1|A1|A2|D|E [arms...]
set -u
W=/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed
source $W/arms.sh
P=/viper/ptmp2/jinma/ckjac_0923
WP=/viper/ptmp2/jinma/wellposed_0923
BJ=$P/athena.cpu.jac
BB=$P/athena.cpu.base
IN=$P/wpj.athinput          # wp.athinput with the glob / esc / aa keys declared
KEYS=$P/jac_keys.athinput   # copy of jac/jac_keys.athinput
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
T4E="$T4 problem/ck_impl_tol=1e-8 problem/ck_impl_dtol=1e-8 problem/ck_impl_maxit=8 \
problem/ck_impl_debug=-2"
SUB="problem/ck_impl_glob=ls_sub problem/ck_impl_sub_max=32"
jarm() {
  case $1 in
    t4)  echo "$T4E" ;;
    sub) echo "$T4E $SUB" ;;
    esc) echo "$T4E $SUB problem/ck_impl_esc=1" ;;
    aa)  echo "$T4E problem/ck_impl_aa=4" ;;
    aas) echo "$T4E $SUB problem/ck_impl_aa=4" ;;
    aae) echo "$T4E $SUB problem/ck_impl_esc=1 problem/ck_impl_aa=4" ;;
    jn)  echo "$T4E problem/ck_impl_jneg=true" ;;
    jna) echo "$T4E problem/ck_impl_jneg=true problem/ck_impl_aa=4" ;;
    semi) echo "problem/ck_impl_debug=-2" ;;
  esac
}
jrun() {
  local b=$1; shift; local d=$1; shift; local a=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && nice -n 10 $b -i $IN $GEO $BASE $(jarm $a) "$@" > run.log 2>&1 )
  echo "$? $d $(grep -c NOT-CONV $d/run.log) nonconv"
}
C=$1; shift
ARMS=${*:-"t4 sub esc aa aas aae jn"}
case $C in
G1)
  # the new switches off are bitwise: base vs jac binary, t4 / sub / semi, 20 and 2000 s,
  # plus D (H2) at 2000 s
  mkdir -p $P/G1; cd $P/G1
  for b in base jac; do
    B=$P/athena.cpu.$b
    for a in t4 sub semi; do
      jrun $B g_${a}_20_$b $a problem/rt_test_mu0=0.5 problem/rt_test_dt=20 time/nlim=20 &
      jrun $B g_${a}_2000_$b $a problem/rt_test_mu0=0.5 problem/rt_test_dt=2000 time/nlim=5 &
      jrun $B gd_${a}_2000_$b $a problem/rt_test_mu0=1.0 hydro/eos_h2=true \
        problem/rt_test_dt=2000 time/nlim=10 &
    done
    wait
  done
  for d in g_t4_20 g_t4_2000 gd_t4_2000 g_sub_20 g_sub_2000 gd_sub_2000 g_semi_20 \
           g_semi_2000 gd_semi_2000; do
    python3 $W/cmpdir.py ${d}_base ${d}_jac
    cmp ${d}_base/rec.txt ${d}_jac/rec.txt && echo "rec.txt identical $d"
  done
  ;;
A1)
  # from T** = the final restart of wellposed_0923/A/a1_t4x_20 (as README_glob.md)
  R=$WP/A/a1_t4x_20/rst/dhj.00003.rst
  mkdir -p $P/A; cd $P/A
  for a in $ARMS; do
    for dt in 20 200 2000; do
      d=a1_${a}_$dt; rm -rf $d; mkdir $d
      ( cd $d && nice -n 10 $BJ -r $R -i $KEYS $(jarm $a) problem/rt_test_dt=$dt \
          time/nlim=$((7200+200)) problem/rt_test_out=rec.txt \
          problem/rt_test_every=10 > run.log 2>&1 ) &
    done
    wait
  done
  for d in a1_*; do echo "$d $(grep -c NOT-CONV $d/run.log) nonconv"; done
  ;;
A2)
  mkdir -p $P/A; cd $P/A
  ln -sfn $WP/A/a2_ref a2_ref
  MU="problem/rt_test_mu0=0.5"
  for a in $ARMS; do
    jrun $BJ a2_${a}_20 $a $MU problem/rt_test_dt=20 time/nlim=100 problem/rt_test_every=100 &
    jrun $BJ a2_${a}_200 $a $MU problem/rt_test_dt=200 time/nlim=10 problem/rt_test_every=10 &
    jrun $BJ a2_${a}_2000 $a $MU problem/rt_test_dt=2000 time/nlim=1 problem/rt_test_every=1 &
    wait
  done
  ;;
D)
  mkdir -p $P/D; cd $P/D
  for a in $ARMS; do
    for dt in 20 200 2000; do
      jrun $BJ d_true_${dt}_$a $a problem/rt_test_mu0=1.0 hydro/eos_h2=true \
        problem/rt_test_dt=$dt time/nlim=100 problem/rt_test_every=10 &
    done
    wait
  done
  ;;
E)
  mkdir -p $P/E; cd $P/E
  for a in $ARMS; do
    for g in tm tmb; do
      case $g in tm) X="problem/ck_beam_sph=false" ;; tmb) X="" ;; esac
      for dt in 20 2000; do
        jrun $BJ e_${g}_${dt}_$a $a $X problem/rt_test_mu0=0.5 problem/rt_test_dt=$dt \
          time/nlim=20 problem/rt_test_every=1000 &
      done
    done
    wait
  done
  ;;
esac
