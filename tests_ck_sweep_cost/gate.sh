#!/bin/bash
# The tests_ck_sph thermal gate (README.md section 0.2: per-cell flux identity and the
# column budget) re-run with problem/ck_sweep_cache off and on, on the spherical form.
# One-shot column dump at t = 0 from inputs/tests/dhj_ck_spherical.athinput.
# The two dumps must be byte-identical and the residuals must stay at 1e-10 / 1e-9.
set -e
R=/viper/u2/jinma/ATHENAK/athenak
BIN=$R/build_cksw/src/athena
cd $R/tests_ck_sweep_cost
for cch in 0 1 2; do
  D=gate_c$cch
  rm -rf $D; mkdir -p $D
  ( cd $D && $BIN -i ../cpu_ab/base.athinput -d . time/nlim=1 \
      problem/ck_spherical=true problem/ck_sweep_cache=$cch \
      problem/ck_dump_file=col.txt problem/ck_dump_m=0 problem/ck_dump_k=5 \
      output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 \
      > run.log 2>&1 ) \
    || { echo "FAIL $D"; tail -15 $D/run.log; exit 1; }
  rm -rf $D/bin $D/rst
done
echo "== column dumps byte-identical? =="
cmp gate_c0/col.txt gate_c1/col.txt && cmp gate_c0/col.txt gate_c2/col.txt && echo "  IDENTICAL"
echo "== budget2.py, ck_sweep_cache = 2 =="
python3 $R/tests_ck_sph/budget2.py gate_c2/col.txt
