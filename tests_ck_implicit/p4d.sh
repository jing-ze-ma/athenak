#!/bin/bash
# phase 4 gate (b'): the pass count as a CONVERGENCE measure -- ck_impl_maxit = 20, so
# that the cap is not what ends the iteration, and the converged states can be compared.
set -u
R=/viper/u2/jinma/ATHENAK/athenak
NEW=$R/build_ckq/src/athena
IN=$R/inputs/tests/dhj_ck_implicit.athinput
BASE="output1/dt=1e30 output2/dt=1e30 output4/dt=1e30 problem/rt_use_cons=true"
SPH="problem/ck_spherical=true problem/ck_beam_sph=true"
IMP="problem/ck_implicit=true problem/ck_impl_verbose=true problem/ck_impl_frozen_op=true
     problem/ck_impl_maxit=20"
run() { local d=$1; shift; rm -rf "$d"; mkdir -p "$d"
  ( cd "$d" && /usr/bin/time -f "%e" -o t.txt $NEW -i $IN "$@" > run.log 2>&1 )
  echo "$? $d $(cat $d/t.txt 2>/dev/null)"; }
run p4d_base time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP
run p4d_r1   time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP problem/ck_impl_reuse_jac=1
run p4d_s2   time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP problem/ck_impl_seed=2
run p4d_b1   time/nlim=20 output3/dt=1e-6 $BASE $SPH $IMP problem/ck_impl_reuse_jac=1 \
    problem/ck_impl_seed=2
echo "=== done"
