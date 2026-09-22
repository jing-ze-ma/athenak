#!/bin/bash
. ./common.sh
RST=$PWD/d_relax_off/rst/dhj.00001.rst
for cfl in 0.3 3.0 30.0; do for ar in 2.0 1e30; do
  d=k_cfl${cfl}_ar${ar}; rm -rf $d; mkdir -p $d
  ( cd $d && $NEW -r $RST time/nlim=202 time/cfl_number=$cfl \
      problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_maxit=20 \
      problem/ck_impl_arat=$ar problem/ck_impl_demax=0.0 \
      output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 > run.log 2>&1 )
  echo "$d exit=$?"
done; done
for d in k_*; do
  echo "== $d"
  grep "### ck_implicit" $d/run.log | head -4 | sed 's/### ck_implicit //'
done
