#!/bin/bash
# bitwise A/B of problem/ck_sweep_cache at both settings of problem/ck_spherical.
# 20 cycles of inputs/tests/dhj_ck_spherical.athinput, hst + a bin dump every cycle.
set -e
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cksw/src/athena
cd $R/tests_ck_sweep_cost/cpu_ab
for sph in false true; do
  for cch in 0 1 2; do
    d=s${sph}_c${cch}
    rm -rf $d && mkdir -p $d
    (cd $d && $X -i ../base.athinput -d . \
        time/nlim=20 \
        problem/ck_spherical=$sph problem/ck_sweep_cache=$cch \
        output1/dt=1.0 output3/dt=1.0 > run.log 2>&1) || { echo "FAIL $d"; tail -20 $d/run.log; exit 1; }
    echo "done $d"
  done
done
for sph in false true; do
  a=s${sph}_c0
  for b in s${sph}_c1 s${sph}_c2; do
    echo "== ck_spherical=$sph : $a vs $b =="
    diff -q <(grep -v '^#' $a/*.hst) <(grep -v '^#' $b/*.hst) \
      && echo "  hst BITWISE identical"
    nd=0; nb=0
    for f in $a/bin/*.bin; do
      g=$b/bin/$(basename $f)
      n=$(cmp -l $f $g | wc -l)
      nb=$((nb+1)); [ "$n" -gt 5 ] && nd=$((nd+1))
      o=$(cmp -l $f $g | awk '{print $1}' | tr '\n' ' ')
      [ "$n" -gt 5 ] && echo "  $(basename $f): $n bytes differ at [$o]"
    done
    echo "  $nb dumps, $nd with a difference outside the 5-byte parameter echo"
  done
done
