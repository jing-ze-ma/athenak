#!/bin/bash
. ./common.sh
N=100
tr_() { local d=$1; shift; local b=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && /usr/bin/time -f "%e" -o t.txt $b -i $IN $BASE time/nlim=$N \
      problem/ck_spherical=true problem/ck_beam_sph=true "$@" > run.log 2>&1 )
  local sw=$(grep -o "nsweep=[0-9]*" $d/run.log | tail -1 | cut -d= -f2)
  local gp=$(grep -o "ckdesum=[-0-9.eE+]*" $d/run.log | tail -60 | cut -d= -f2 | \
             python3 -c "import sys;v=[abs(float(x)) for x in sys.stdin];print('%.2e'%(max(v) if v else 0))")
  echo "COST $d wall=$(cat $d/t.txt) sweeps=${sw:-200} maxgap=$gp"
}
tr_ c_off $REF
tr_ c_p1     $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_arat=1e30 problem/ck_impl_colskip=false
tr_ c_split  $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_colskip=false
tr_ c_skip   $NEW problem/ck_implicit=true problem/ck_impl_verbose=true
tr_ c_once   $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true
tr_ c_once_m6 $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true problem/ck_impl_maxit=6
tr_ c_once_m5 $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true problem/ck_impl_maxit=5
tr_ c_once_m4 $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true problem/ck_impl_maxit=4
tr_ c_once_m6_nosk $NEW problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true problem/ck_impl_maxit=6 problem/ck_impl_colskip=false
echo "COSTDONE"
# ---- once-per-step relaxed arm, for gate (f)
run d_relax_once $NEW time/nlim=200 problem/ck_spherical=true problem/ck_beam_sph=true \
  problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_once=true \
  problem/ck_impl_maxit=6 output3/dt=1e30 > /dev/null 2>&1
echo "RELAXONCE done"
./bigdt.sh > bigdt.out 2>&1
echo "BIGDT done"
for arm in off on once; do
  RST=$(ls -1 d_relax_$arm/rst/*.rst | tail -1)
  for col in "night 0 2" "day92 1 2" "day38 0 4"; do
    set -- $col
    D=col_${1}_$arm; rm -rf $D; mkdir -p $D
    ( cd $D && $NEW -r $PWD/../$RST -t 00:00:30 time/nlim=999999 \
        problem/ck_dump_file=col.txt problem/ck_dump_m=$2 problem/ck_dump_k=$3 \
        output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 \
        > run.log 2>&1 )
    rm -rf $D/bin $D/rst
    echo "COL $D $(wc -l < $D/col.txt 2>/dev/null) lines"
  done
done
echo ALLGATES_DONE
