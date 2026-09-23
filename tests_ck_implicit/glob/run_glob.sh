#!/bin/bash -l
# ck_impl_glob gates on the well-posed suite (tests_ck_implicit/wellposed/), CPU, login
# node nice 10.  Binaries: athena.cpu.base (HEAD d80c84ca + wellposed_hooks.patch) and
# athena.cpu.glob (ck-glob + the same patch).  Arms: t4 (as wellposed/arms.sh), ls, sub
# (= t4 + ck_impl_glob = ls / ls_sub).  usage: run_glob.sh G1|A1|A2|D|Dref|E
set -u
W=/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed
source $W/arms.sh
P=/viper/ptmp2/jinma/ckglob_0924
WP=/viper/ptmp2/jinma/wellposed_0923
BG=$P/athena.cpu.glob
BB=$P/athena.cpu.base
# wp.athinput with the ck_impl_glob keys declared (command-line overrides need the key);
# built by: awk (see README_glob.md) from wellposed/wp.athinput + glob_keys.athinput
IN=$P/wpg.athinput
KEYS=$P/glob_keys.athinput   # copy of glob/glob_keys.athinput
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
# every T4 switch and tolerance set explicitly (restarts embed the spin-up arm's switches)
T4E="$T4 problem/ck_impl_tol=1e-8 problem/ck_impl_dtol=1e-8 problem/ck_impl_maxit=8"
garm() {
  case $1 in
    t4)  echo "$T4E" ;;
    ls)  echo "$T4E problem/ck_impl_glob=ls" ;;
    sub) echo "$T4E problem/ck_impl_glob=ls_sub" ;;
    s32) echo "$T4E problem/ck_impl_glob=ls_sub problem/ck_impl_sub_max=32" ;;
    semi) echo "" ;;
  esac
}
# grun <bin> <dir> <arm> <args...>: fresh run from the input file
grun() {
  local b=$1; shift; local d=$1; shift; local a=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && nice -n 10 $b -i $IN $GEO $BASE $(garm $a) "$@" > run.log 2>&1 )
  echo "$? $d $(grep -c NOT-CONV $d/run.log) nonconv"
}
case $1 in
G1)
  # switch none is bitwise: base vs glob binary, t4 and semi, three dt, plus D (H2)
  mkdir -p $P/G1; cd $P/G1
  for b in base glob; do
    B=$P/athena.cpu.$b
    for a in t4 semi; do
      grun $B g_${a}_20_$b $a problem/rt_test_mu0=0.5 problem/rt_test_dt=20 time/nlim=20 &
      grun $B g_${a}_2000_$b $a problem/rt_test_mu0=0.5 problem/rt_test_dt=2000 time/nlim=5 &
      grun $B gd_${a}_2000_$b $a problem/rt_test_mu0=1.0 hydro/eos_h2=true \
        problem/rt_test_dt=2000 time/nlim=10 &
    done
    wait
  done
  for d in g_t4_20 g_t4_2000 gd_t4_2000 g_semi_20 g_semi_2000 gd_semi_2000; do
    python3 $W/cmpdir.py ${d}_base ${d}_glob
    cmp ${d}_base/rec.txt ${d}_glob/rec.txt && echo "rec.txt identical $d"
  done
  ;;
A1)
  # from T** = the final restart of wellposed_0923/A/a1_t4x_20 (A0b + 200 calls t4x at
  # 20 s; the A0b restart itself is gone), 200 calls per arm and dt
  R=$WP/A/a1_t4x_20/rst/dhj.00003.rst
  mkdir -p $P/A; cd $P/A
  for a in t4 ls sub s32; do
    for dt in 20 200 2000; do
      d=a1_${a}_$dt; rm -rf $d; mkdir $d
      ( cd $d && nice -n 10 $BG -r $R -i $KEYS $(garm $a) problem/rt_test_dt=$dt \
          time/nlim=$((7200+200)) problem/rt_test_out=rec.txt \
          problem/rt_test_every=10 > run.log 2>&1 ) &
    done
  done
  wait
  for d in a1_*; do echo "$d $(grep -c NOT-CONV $d/run.log) nonconv"; done
  ;;
A2)
  mkdir -p $P/A; cd $P/A
  ln -sfn $WP/A/a2_ref a2_ref
  MU="problem/rt_test_mu0=0.5"
  for a in t4 ls sub s32; do
    grun $BG a2_${a}_20 $a $MU problem/rt_test_dt=20 time/nlim=100 problem/rt_test_every=100 &
    grun $BG a2_${a}_200 $a $MU problem/rt_test_dt=200 time/nlim=10 problem/rt_test_every=10 &
    grun $BG a2_${a}_2000 $a $MU problem/rt_test_dt=2000 time/nlim=1 problem/rt_test_every=1 &
  done
  wait
  ;;
D)
  mkdir -p $P/D; cd $P/D
  for dt in 20 200 2000; do
    for a in t4 ls sub s32; do
      grun $BG d_true_${dt}_$a $a problem/rt_test_mu0=1.0 hydro/eos_h2=true \
        problem/rt_test_dt=$dt time/nlim=100 problem/rt_test_every=10 &
    done
  done
  wait
  ;;
Dref)
  # small-dt reference for D (eos_h2 true): t4 at 20 s to t = 2e5 s (call 10000)
  mkdir -p $P/D; cd $P/D
  grun $BG d_true_ref20 t4 problem/rt_test_mu0=1.0 hydro/eos_h2=true \
    problem/rt_test_dt=20 time/nlim=10000 problem/rt_test_every=100
  ;;
E)
  # symmetry (identical columns): tm and tmb geometries, 20 calls
  mkdir -p $P/E; cd $P/E
  for g in tm tmb; do
    case $g in tm) X="problem/ck_beam_sph=false" ;; tmb) X="" ;; esac
    for dt in 20 2000; do
      for a in t4 ls sub s32; do
        grun $BG e_${g}_${dt}_$a $a $X problem/rt_test_mu0=0.5 problem/rt_test_dt=$dt \
          time/nlim=20 problem/rt_test_every=1000 &
      done
    done
  done
  wait
  ;;
esac
