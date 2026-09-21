#!/bin/bash
# tests_m1/runs_3g_newton_T -- G2 unit gates of milestone 3g.  The four arms (reference,
# A = implicit_gas_newton, B = implicit_eos_cache, A+B) of every implicit gate of the
# module that exercises the GAS coupling, plus two that do not (T3, T3b: pure scattering,
# where the options are unreachable and the answer must be bitwise).
# Serial CPU, build_cpu_m1.  Usage: bash run_gates.sh
set -u
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R=$(cd "$HERE/../.." && pwd)
X=$R/build_cpu_m1/src/athena
I=$HERE/inp

run() { local name=$1; shift; local dir=$HERE/gates/$name
  rm -rf "$dir"; mkdir -p "$dir"
  (cd "$dir" && timeout 1500 "$X" -i "$1" -d . "${@:2}" > log.txt 2>&1)
  echo "=== $name rc=$?"; }

IMP="rad_m1/transport=implicit_x1 time/cfl_number=1e6"
IMPJ="$IMP rad_m1/implicit_bc_x1min=efix rad_m1/implicit_bc_x1max=efix"
IMPM="$IMP rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 \
rad_m1/implicit_bc_x1max=marshak"
IMPT="rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab \
rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step \
rad_m1/implicit_allow_multid=true"

for arm in "ref:" "A:rad_m1/implicit_gas_newton=true" \
           "B:rad_m1/implicit_eos_cache=true" \
           "AB:rad_m1/implicit_gas_newton=true rad_m1/implicit_eos_cache=true"; do
  a=${arm%%:*}; o=${arm#*:}
  run t5i_$a  $I/rad_m1_equil.athinput       $IMP rad_m1/implicit_cfl=4.4e-7 $o
  run t5t_$a  $I/rad_m1_equil_table.athinput $IMP rad_m1/implicit_cfl=1.0e-4 $o
  run t6_$a   $I/rad_m1_marshak.athinput     $IMPM rad_m1/implicit_cfl=10 \
              mesh/nx1=128 meshblock/nx1=128 $o
  run t3_$a   $I/rad_m1_thick_pulse.athinput $IMP rad_m1/implicit_cfl=10 $o
  run t3b_$a  $I/rad_m1_jump.athinput        $IMPJ rad_m1/implicit_cfl=10 $o
  run t10a_$a $I/rad_m1_radwave.athinput     mesh/nx2=1 meshblock/nx2=1 \
              problem/radwave_dir=x1 rad_m1/transport=implicit_x1 \
              time/tlim=1.184313050927584 output1/dt=0.0370097828 $o
  run t10c_$a $I/rad_m1_radwave.athinput     mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 \
              meshblock/nx2=64 problem/radwave_dir=x2 $IMPT \
              time/tlim=1.184313050927584 output1/dt=0.0370097828 $o
done
