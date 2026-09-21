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
PUL=$H/../runs_3b5/pulse_md.athinput
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
