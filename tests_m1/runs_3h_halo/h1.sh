#!/bin/bash
# H1: bitwise gates, reference binary (HEAD) vs the optimised one.
set -u
H=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
A=$H/../..
REFB=$H/athena_ref_box;  NEWB=$A/build_cpu_box/src/athena
REFM=$H/athena_ref_m1;   NEWM=$A/build_cpu_m1/src/athena
run () {  # run <dir> <exe> <inp> <overrides...>
  d=$H/$1; shift; exe=$1; shift; inp=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  (cd "$d" && timeout 1500 "$exe" -i "$inp" -d . "$@" > log.txt 2>&1)
  echo "rc=$?"
}
EDD=$H/../runs_3g_newton_T/he_slab_m1_2d_V3edd.athinput
COL=$H/../runs_3g_newton_T/he_box_m1_1d_impl.athinput
LEV=$H/../runs_3b8/he_slab_m1_2d.athinput
PUL=$H/../runs_3b5/pulse_md.athinput

O_STEP="problem/vpert=1.0e-3 time/tlim=10 rad_m1/implicit_gas_newton=true"
run a_step_ref $REFB $EDD $O_STEP rad_m1/implicit_closure_lag=step
run a_step_new $NEWB $EDD $O_STEP rad_m1/implicit_closure_lag=step
run a_pass_ref $REFB $EDD $O_STEP rad_m1/implicit_closure_lag=pass
run a_pass_new $NEWB $EDD $O_STEP rad_m1/implicit_closure_lag=pass

O_LEV="problem/vpert=1.0e-3 rad_m1/implicit_offdiag=operator \
rad_m1/implicit_closure_lag=step rad_m1/implicit_trans_limit=lp time/tlim=3"
run b_lev_ref $REFB $LEV $O_LEV
run b_lev_new $NEWB $LEV $O_LEV

run c_col_ref $REFB $COL time/tlim=10
run c_col_new $NEWB $COL time/tlim=10

# multi-D pulse, runs_3b5 G4 style: 1 vs 2x2 MeshBlocks, both binaries
P="mesh/nx1=64 mesh/nx2=64 mesh/x2min=0.0 mesh/x2max=1.0 problem/pulse_y0=0.5 \
time/cfl_number=1e6 rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab \
rad_m1/kappa_s=64000 rad_m1/implicit_cfl=1e4 time/tlim=1875.0 output1/dt=1e30 \
output3/dt=1875.0 rad_m1/implicit_offdiag=operator"
B1="meshblock/nx1=64 meshblock/nx2=64"
B2="meshblock/nx1=64 meshblock/nx2=32"
B4="meshblock/nx1=32 meshblock/nx2=32 mesh/ix1_bc=reflect mesh/ox1_bc=reflect rad_m1/implicit_partition=gather"
run d_p1_ref $REFM $PUL $P $B1
run d_p1_new $NEWM $PUL $P $B1
run d_p2_ref $REFM $PUL $P $B2
run d_p2_new $NEWM $PUL $P $B2
run e_g1_ref $REFM $PUL $P $B1 mesh/ix1_bc=reflect mesh/ox1_bc=reflect rad_m1/implicit_partition=gather
run e_g1_new $NEWM $PUL $P $B1 mesh/ix1_bc=reflect mesh/ox1_bc=reflect rad_m1/implicit_partition=gather
run e_g4_ref $REFM $PUL $P $B4
run e_g4_new $NEWM $PUL $P $B4
