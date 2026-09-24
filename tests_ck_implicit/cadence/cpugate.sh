#!/bin/bash
# ck-cadence CPU gates (login node, serial, nice 10).
#  A-off: base (rt-integration 90b01f2f) vs cad with every = 1: well-posed A2 (t4, semi)
#         rec.txt + rst, and the full hydro 6 cycles (T4, and c2) rst + hst.
#  A-rst: cad, full hydro with c2 + every = 4 (no guard, and thr 0.02): 12 cycles straight
#         against 8 + restart + 4 (restart on a multiple of N), and 6 + 6 (not a multiple);
#         the same for every = 1 as the baseline restart behaviour of c2.
Q=/viper/ptmp2/jinma/ckcad_0924
F=/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/fast
C=/viper/u2/jinma/ATHENAK/athenak/tests_ck_implicit/wellposed/cmpdir.py
IN=$F/wp_fast.athinput
GEO="problem/rt_ck=true problem/rt_use_cons=true problem/ck_spherical=true problem/ck_beam_sph=true problem/ck_sweep_form=1"
T4="problem/ck_implicit=true problem/ck_impl_arat=1e30 problem/ck_impl_frozen_op=true problem/ck_impl_lin=true problem/ck_impl_lin_thr=1 problem/ck_impl_fuse=true problem/ck_impl_jac_lin=true problem/ck_impl_cvsec=true problem/ck_impl_verbose=true"
C2="$T4 problem/ck_impl_once=true problem/ck_impl_xstep=2 problem/ck_impl_jreuse=0.2 problem/ck_impl_pred=true problem/ck_impl_pred_fac=0.5 problem/ck_impl_nosync=true"
O="output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 time/tlim=1e30"
module purge; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
export OMP_NUM_THREADS=1
if [ "${1:-all}" != rst ]; then
for b in base cad1; do bash $F/wp_run.sh $Q/athena.cpu.$b $Q/wp_$b A2 t4 semi & done
IN=$F/wp_fast.athinput
for b in base cad1; do for a in t4 c2; do d=$Q/hyd/${a}_$b; rm -rf $d; mkdir -p $d
  if [ $a = t4 ]; then X="$T4"; else X="$C2"; fi
  (cd $d && nice -n 10 $Q/athena.cpu.$b -i $IN $GEO $O $X time/nlim=6 > run.log 2>&1; echo "rc=$? hyd $a $b") &
done; done
wait
for a in t4 semi; do cmp $Q/wp_base/A2/$a/rec.txt $Q/wp_cad1/A2/$a/rec.txt && echo "rec $a BITWISE"
  python3 $C $Q/wp_base/A2/$a $Q/wp_cad1/A2/$a; done
for a in t4 c2; do python3 $C $Q/hyd/${a}_base $Q/hyd/${a}_cad1; done
fi
# ---- restart gate (cad binary; the input copy carries the new keys)
IN=/viper/ptmp2/jinma/wt_ckcad/tests_ck_implicit/cadence/wp_cad.athinput
B=$Q/athena.cpu.${RB:-cad1}
R=$Q/rstgate
for arm in e1 e4 e4g2; do
  case $arm in e1) X="$C2" ;; e4) X="$C2 problem/ck_impl_every=4" ;;
    e4g2) X="$C2 problem/ck_impl_every=4 problem/ck_impl_every_thr=0.02" ;; esac
  for run in s12 a8 a6; do d=$R/$arm/$run; rm -rf $d; mkdir -p $d
    n=${run:1}
    (cd $d && nice -n 10 $B -i $IN $GEO $O $X time/nlim=$n > run.log 2>&1; echo "rc=$? $arm $run") &
  done
done
wait
for arm in e1 e4 e4g2; do for run in a8 a6; do
  d=$R/$arm/r${run:1}; rm -rf $d; mkdir -p $d
  f=$(ls $R/$arm/$run/rst/*.rst | tail -1)
  (cd $d && nice -n 10 $B -r $f time/nlim=12 > run.log 2>&1; echo "rc=$? $arm restart ${run}") &
done; done
wait
python3 - <<PY
import glob
import numpy as np
R="$R"
def data(f):
    d=open(f,'rb').read(); i=d.find(b'<par_end>'); return d[i:] if i>=0 else d
def rel(a,b):
    n=min(len(a),len(b))//8*8
    x=np.frombuffer(a[len(a)-n:],dtype='<f8'); y=np.frombuffer(b[len(b)-n:],dtype='<f8')
    ok=np.isfinite(x)&np.isfinite(y)&(np.abs(x)>1e-300)
    return np.max(np.abs(y[ok]/x[ok]-1.0))
for arm in ['e1','e4','e4g2']:
    s=sorted(glob.glob(R+'/'+arm+'/s12/rst/*.rst'))[-1]
    for r in ['r8','r6']:
        f=sorted(glob.glob(R+'/'+arm+'/'+r+'/rst/*.rst'))[-1]
        a,b=data(s),data(f)
        print('RSTGATE', arm, r, 'BITWISE' if a==b else 'DIFFER max rel %.3e' % rel(a,b))
PY
echo CPUGATE_DONE
