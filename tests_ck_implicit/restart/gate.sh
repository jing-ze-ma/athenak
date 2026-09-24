#!/bin/bash
# ck-restart CPU gates (login node, serial, nice 10).  See README_restart.md.
#  R: N cycles straight vs N/2 + restart + N/2 (new binary), arms t4 c2 c2e4 (+ semi),
#     on the well-posed transient (wpc, N = 8, and 3 + 5) and the full-hydro 6-cycle
#     case (wpf, 3 + 3; the c2 restart there lands on an xstep re-apply -> rebuild).
#  D: base vs new, straight runs (hst + bin every cycle): the new code does not change
#     a run.  O: base and new restarted from the SAME old-format (base) restart file.
Q=/viper/ptmp2/jinma/ckrst_0924
NEW=$Q/athena.cpu.${NEWTAG:-r2}
BASE=$Q/athena.cpu.base
T=/viper/ptmp2/jinma/wt_ckrst/tests_ck_implicit
CMP=$T/restart/cmp.py
GEO="problem/rt_ck=true problem/rt_use_cons=true problem/ck_spherical=true problem/ck_beam_sph=true problem/ck_sweep_form=1"
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true problem/ck_impl_verbose=true"
C2="$T4 problem/ck_impl_once=true problem/ck_impl_xstep=2 problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
O="output1/dt=1e-6 output2/dt=1e30 output3/dt=1e-6 time/tlim=1e30"
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
args() { case $1 in semi) echo "" ;; t4) echo "$T4" ;; c2) echo "$C2" ;;
  c2e4) echo "$C2 problem/ck_impl_every=4" ;; c2e4g) echo "$C2 problem/ck_impl_every=4 problem/ck_impl_every_thr=0.02" ;; esac; }
inp() { echo $T/restart/wp_rst.athinput; }
G=$Q/gate
run() {  # run <dir> <bin> <inp> <arm> <nlim>
  rm -rf $1; mkdir -p $1
  (cd $1 && nice -n 10 $2 -i $(inp $3) $GEO $O $(args $4) time/nlim=$5 > run.log 2>&1; echo "rc=$? $1" >> $G/rc.txt)
}
rrun() { # rrun <dir> <bin> <from-dir> <nlim>
  rm -rf $1; mkdir -p $1
  f=$(ls $3/rst/*.rst | tail -1)
  (cd $1 && nice -n 10 $2 -r $f $O time/nlim=$4 > run.log 2>&1; echo "rc=$? $1" >> $G/rc.txt)
}
mkdir -p $G; : > $G/rc.txt
ARMS="semi t4 c2 c2e4 c2e4g"
# phase 1: straight and first halves
for a in $ARMS; do
  run $G/wpc/$a/s8 $NEW wpc $a 8 &
  run $G/wpc/$a/a4 $NEW wpc $a 4 &
  run $G/wpc/$a/a3 $NEW wpc $a 3 &
  run $G/wpf/$a/s6 $NEW wpf $a 6 &
  run $G/wpf/$a/a3 $NEW wpf $a 3 &
  run $G/wpf/$a/b_s6 $BASE wpf $a 6 &
  run $G/wpf/$a/b_a3 $BASE wpf $a 3 &
done
wait
# phase 2: restarts
for a in $ARMS; do
  rrun $G/wpc/$a/r4 $NEW $G/wpc/$a/a4 8 &
  rrun $G/wpc/$a/r3 $NEW $G/wpc/$a/a3 8 &
  rrun $G/wpf/$a/r3 $NEW $G/wpf/$a/a3 6 &
  rrun $G/wpf/$a/o_new $NEW $G/wpf/$a/b_a3 6 &
  rrun $G/wpf/$a/o_base $BASE $G/wpf/$a/b_a3 6 &
done
wait
{
echo "== R: restart continuation (new binary): straight vs half + restart + half"
for a in $ARMS; do
  python3 $CMP rst $G/wpc/$a/s8 $G/wpc/$a/r4
  python3 $CMP rst $G/wpc/$a/s8 $G/wpc/$a/r3
  python3 $CMP rst $G/wpf/$a/s6 $G/wpf/$a/r3
done
echo "== D: base vs new, straight 6 cycles (hst, bin every cycle; rst for semi)"
for a in $ARMS; do
  if [ $a = semi ]; then python3 $CMP run $G/wpf/$a/b_s6 $G/wpf/$a/s6 rst
  else python3 $CMP run $G/wpf/$a/b_s6 $G/wpf/$a/s6; fi
done
echo "== O: an OLD (base-written) restart read by base and by new"
for a in $ARMS; do python3 $CMP run $G/wpf/$a/o_base $G/wpf/$a/o_new; done
echo "== rebuild / warning lines"
grep -h "ck restart\|implicit-ck state\|xstep snapshot" $G/*/*/r*/run.log $G/wpf/*/o_new/run.log | sort | uniq -c
grep -v "rc=0" $G/rc.txt
} > $G/gate.out 2>&1
echo GATE_DONE
