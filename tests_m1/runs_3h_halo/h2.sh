#!/bin/bash
set -u
H=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
A=$H/../..
EDD=$H/../runs_3g_newton_T/he_slab_m1_2d_V3edd.athinput
run () { d=$H/$1; shift; exe=$1; shift
  rm -rf "$d"; mkdir -p "$d"
  (cd "$d" && timeout 1500 "$exe" -i "$EDD" -d . problem/vpert=1.0e-3 time/tlim=40 \
     output5/dt=1000 output2/dt=1000 "$@" > log.txt 2>&1)
  printf "%-14s rc=%d  " "$(basename $d)" $?
  grep -E "zone-cycles/cpu_second|cpu time used" "$d/log.txt" | tr '\n' ' '; echo
}
run t_newt_ref $H/athena_ref_box               rad_m1/implicit_gas_newton=true
run t_newt_new $A/build_cpu_box/src/athena     rad_m1/implicit_gas_newton=true
run t_base_ref $H/athena_ref_box               rad_m1/implicit_gas_newton=false
run t_base_new $A/build_cpu_box/src/athena     rad_m1/implicit_gas_newton=false
