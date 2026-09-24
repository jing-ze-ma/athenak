#!/bin/bash
# restart gate, history form: hst every cycle; 12 cycles straight vs 8 + restart + 4.
# c2 itself (e1) is not a bitwise restart (the Newton's cross-call state, e.g. the pass-0
# thin/thick classification ck_thk, is not in the restart); this measures whether the
# cadence adds to that difference.
Q=/viper/ptmp2/jinma/ckcad_0924
IN=/viper/ptmp2/jinma/wt_ckcad/tests_ck_implicit/cadence/wp_cad.athinput
GEO="problem/rt_ck=true problem/rt_use_cons=true problem/ck_spherical=true problem/ck_beam_sph=true problem/ck_sweep_form=1"
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true problem/ck_impl_verbose=true"
C2="$T4 problem/ck_impl_once=true problem/ck_impl_xstep=2 problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
O="output1/dt=1e-6 output2/dt=1e30 output3/dt=1e-6 time/tlim=1e30"
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
B=$Q/athena.cpu.cad1
R=$Q/rstgate3
for arm in e1 e4 e4g2; do
  case $arm in e1) X="$C2" ;; e4) X="$C2 problem/ck_impl_every=4" ;;
    e4g2) X="$C2 problem/ck_impl_every=4 problem/ck_impl_every_thr=0.02" ;; esac
  for run in s8 a4; do d=$R/$arm/$run; rm -rf $d; mkdir -p $d
    (cd $d && nice -n 10 $B -i $IN $GEO $O $X time/nlim=${run:1} > run.log 2>&1) &
  done
done
wait
for arm in e1 e4 e4g2; do d=$R/$arm/r4; rm -rf $d; mkdir -p $d
  f=$(ls $R/$arm/a4/rst/*.rst | tail -1)
  (cd $d && nice -n 10 $B -r $f time/nlim=8 > run.log 2>&1) &
done
wait
python3 - <<PY
import glob, os, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert
R="$R"
for arm in ['e1','e4','e4g2']:
    out=[]
    for fb in sorted(glob.glob(R+'/'+arm+'/r4/bin/*.bin')):
        fa=R+'/'+arm+'/s8/bin/'+os.path.basename(fb)
        A=bin_convert.read_binary(fa); B=bin_convert.read_binary(fb)
        d=[]
        for v in ('dens','eint'):
            x=np.asarray(A['mb_data'][v]); y=np.asarray(B['mb_data'][v])
            d.append(np.max(np.abs(y/x-1.0)))
        out.append('%s d %.2e e %.2e' % (os.path.basename(fb).split('.')[-2], d[0], d[1]))
    print('RSTGATE3', arm, '; '.join(out))
PY
echo RSTGATE3_DONE
