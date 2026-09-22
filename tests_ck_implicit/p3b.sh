#!/bin/bash
# phase 3, the WARM START measured where it can be measured: maxit raised so that the
# calls actually reach ck_impl_tol, otherwise the pass count is just the cap.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
NEW=$R/build_ckfast/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
BASE="output1/dt=1e30 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 \
problem/rt_use_cons=true problem/ck_spherical=true problem/ck_beam_sph=true \
problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true \
problem/ck_impl_maxit=20"
run() { local d=$1; shift; rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && /usr/bin/time -f "%e" -o t.txt $NEW -i $IN $BASE "$@" > run.log 2>&1 )
  echo "$? $d $(cat $d/t.txt)"; }
run p3w_cold time/nlim=20
run p3w_warm time/nlim=20 problem/ck_impl_warm=true
# and the same pair with a dump every cycle, to compare the converged states
run p3w_cold_d time/nlim=20 output3/dt=1e-6
run p3w_warm_d time/nlim=20 output3/dt=1e-6 problem/ck_impl_warm=true
echo done
