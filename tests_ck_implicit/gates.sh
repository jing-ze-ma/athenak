#!/bin/bash
# Gates (b) .. (f) for problem/ck_implicit.  See README.md.  Run from this directory.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
REF=$R/build_ckimp_ref/src/athena
NEW=$R/build_ckimp_new/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
BASE="output3/dt=1e30 output4/dt=1e30 problem/rt_use_cons=true"

run() {  # run <dir> <binary> <extra args...>
  local d=$1; shift
  local b=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && $b -i $IN $BASE "$@" > run.log 2>&1 )
  echo "$? $d"
}

echo "=== (b) the sweep-to-gas gap, four combinations, 20 cycles"
for sp in false true; do
  for bs in false true; do
    # today: the rt_desum line, which needs the per-cell report
    run b_off_sp${sp}_bs${bs} $REF time/nlim=20 problem/ck_spherical=$sp \
        problem/ck_beam_sph=$bs problem/rt_cell_report=true problem/rt_report_every=5
    # ck_implicit: the ckdesum line, the same ratio formed from the converged state
    run b_on_sp${sp}_bs${bs} $NEW time/nlim=20 problem/ck_spherical=$sp \
        problem/ck_beam_sph=$bs problem/ck_implicit=true problem/ck_impl_verbose=true
  done
done

echo "=== (e) robustness at large dt: cfl 0.3 / 3 / 30, ck_spherical on and off"
for sp in false true; do
  for cfl in 0.3 3.0 30.0; do
    run e_sp${sp}_cfl${cfl} $NEW time/nlim=4 time/cfl_number=$cfl \
        problem/ck_spherical=$sp problem/ck_implicit=true problem/ck_impl_verbose=true \
        problem/ck_impl_maxit=20
  done
done

echo "=== (f) cost, 100 cycles"
for sp in false true; do
  for bs in false true; do
    d=f_off_sp${sp}_bs${bs}; rm -rf $d; mkdir -p $d
    ( cd $d && /usr/bin/time -f "%e" -o t.txt $REF -i $IN $BASE time/nlim=100 \
        problem/ck_spherical=$sp problem/ck_beam_sph=$bs > run.log 2>&1 )
    d=f_on_sp${sp}_bs${bs}; rm -rf $d; mkdir -p $d
    ( cd $d && /usr/bin/time -f "%e" -o t.txt $NEW -i $IN $BASE time/nlim=100 \
        problem/ck_spherical=$sp problem/ck_beam_sph=$bs problem/ck_implicit=true \
        problem/ck_impl_verbose=true > run.log 2>&1 )
  done
done
echo "=== done"
