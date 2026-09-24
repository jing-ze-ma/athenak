#!/bin/bash -l
# runs_5c_thinstab gate: implicit_closure_thin_relax = 1.5 ON vs 0 (named, = off) on the
# thick / wave problems, serial CPU.  usage: gate_thick.sh <athena> <run dir>
# Every output file of the pair is compared byte for byte; t10_radwave.py / t6 numbers are
# printed for whatever differs.
EXE=$1; R=$2; H=$(dirname $(readlink -f $0))/../../..
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
mkdir -p $R/inp
for f in tests_m1/runs_3b8/rad_m1_radwave.athinput \
         tests_m1/runs_5a_sp_s1/inp/marshak_cart.athinput \
         tests_m1/runs_3b5/pulse_md.athinput; do
  b=$(basename $f)
  # name the keys the arms override (at their defaults, so the input is unchanged)
  python3 - $H/$f $R/inp/$b <<'PY'
import re, sys
s = open(sys.argv[1]).read()
add = ''
for k, v in (('implicit_closure_thin_relax', '0.0'), ('implicit_closure_lag', 'pass'),
             ('implicit_solver', 'line_jacobi'), ('implicit_offdiag', 'auto')):
    if not re.search(r'^\s*' + k + r'\s*=', s, re.M):
        add += '%s = %s\n' % (k, v)
open(sys.argv[2], 'w').write(s.replace('<rad_m1>\n', '<rad_m1>\n' + add, 1))
PY
done
IMPL="rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"
declare -A C
C[rw_x2]="rad_m1_radwave mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 meshblock/nx2=64 problem/radwave_dir=x2 $IMPL time/tlim=1.184313050927584 output1/dt=0.0370097828"
C[rw_xy]="rad_m1_radwave mesh/nx1=64 mesh/nx2=64 meshblock/nx1=64 meshblock/nx2=64 problem/radwave_dir=xy $IMPL time/tlim=0.8374357893586236 output1/dt=0.0261698684"
C[mk_m1]="marshak_cart rad_m1/closure=m1 rad_m1/implicit_closure_lag=step"
C[pulse2d]="pulse_md mesh/nx2=64 meshblock/nx2=64 $IMPL rad_m1/implicit_cfl=20.0"
for c in "${!C[@]}"; do
  set -- ${C[$c]}; inp=$1; shift
  for w in off on; do
    v=0.0; [ $w = on ] && v=1.5
    d=$R/${c}_$w; rm -rf $d; mkdir -p $d
    (cd $d && OMP_NUM_THREADS=1 nice $EXE -i $R/inp/$inp.athinput -d $d "$@" \
       rad_m1/implicit_closure_thin_relax=$v > log.txt 2>&1; echo "rc=$?" >> log.txt) &
  done
done
wait
for c in $(echo "${!C[@]}" | tr ' ' '\n' | sort); do
  echo "$c: $(python3 $H/tests_m1/runs_5c_thinstab/cmpdata.py $R/${c}_off $R/${c}_on); $(grep -h 'rc=' $R/${c}_off/log.txt $R/${c}_on/log.txt | tr '\n' ' ')"
done
