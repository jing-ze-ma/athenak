#!/bin/bash
# MILESTONE 3b PHASE E -- gate T10, the radiation-modified acoustic wave.
#
# The same linear wave (design sect. 14 / inputs/tests/rad_m1_radwave.athinput) is
# propagated along x1 and along x2 under the validated column path (implicit_x1), under
# the full 7-point implicit transport, and under the explicit scheme, so that the
# TRANSVERSE gas-radiation coupling is measured against an x1 path that milestones 3a/3b
# already gate.  Serial CPU, build_cpu_m1 (no -D PROBLEM: radwave is a BUILT-IN test).
#
#   usage:  bash tests_m1/runs_3b6/run_gates_3b6.sh [arm ...]      (no arm = all)
set -u
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
EXE=$ROOT/build_cpu_m1/src/athena
INP=$ROOT/inputs/tests/rad_m1_radwave.athinput
OUT=$ROOT/tests_m1/runs_3b6
TLIM2=1.184313050927584        # 2 periods, lambda = 1 (x1, x2, x3 arms)
TLIMD=0.8374357893586236       # 2 periods, lambda = 1/sqrt(2) (xy arm)
ODT2=0.0370097828             # TLIM2/32
ODTD=0.0261698684             # TLIMD/32

# the recommended multi-D production settings (runs_3b5)
IMPL="rad_m1/transport=implicit rad_m1/implicit_solver=bicgstab \
rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"

run() {          # run <name> <extra athinput overrides...>
  local name=$1; shift
  local dir=$OUT/$name
  rm -rf "$dir"; mkdir -p "$dir"
  echo "=== $name"
  (cd "$dir" && "$EXE" -i "$INP" -d . "$@" > log.txt 2>&1)
  tail -1 "$dir/log.txt"
}


ARMS="${*:-a b c d1 d2 e f1 f2}"
has() { case " $ARMS " in *" $1 "*) return 0;; esac; return 1; }

# (a) x1, 1-D mesh, implicit_x1 -- the validated path
has a && run a_x1_implx1 \
  mesh/nx2=1 meshblock/nx2=1 problem/radwave_dir=x1 \
  rad_m1/transport=implicit_x1 time/tlim=$TLIM2 output1/dt=$ODT2

# (b) x1 on a 2-D mesh, full implicit transport
has b && run b_x1_impl2d \
  mesh/nx2=4 meshblock/nx1=64 meshblock/nx2=4 problem/radwave_dir=x1 \
  $IMPL time/tlim=$TLIM2 output1/dt=$ODT2

# (c) x2 on a 2-D mesh, full implicit transport -- the wave rotated by 90 degrees
has c && run c_x2_impl2d \
  mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 meshblock/nx2=64 problem/radwave_dir=x2 \
  $IMPL time/tlim=$TLIM2 output1/dt=$ODT2

# (d) the EXPLICIT scheme along x1 and along x2 (itself only ever gated in 1-D)
has d1 && run d1_x1_expl2d \
  mesh/nx2=4 meshblock/nx1=64 meshblock/nx2=4 problem/radwave_dir=x1 \
  rad_m1/transport=explicit time/tlim=$TLIM2 output1/dt=$ODT2

has d2 && run d2_x2_expl2d \
  mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 meshblock/nx2=64 problem/radwave_dir=x2 \
  rad_m1/transport=explicit time/tlim=$TLIM2 output1/dt=$ODT2

# (f) the SAME wave in a MARGINALLY THICK medium (tau per wavelength 1e2, i.e. 1.6 per
# cell), where the radiation is no longer in equilibrium diffusion and the transverse
# face flux is no longer a pure diffusion: there is no closed-form phase speed, so the
# gate is x2 against x1 alone.
has f1 && run f1_x1_thin \
  mesh/nx2=4 meshblock/nx1=64 meshblock/nx2=4 problem/radwave_dir=x1 \
  rad_m1/kappa_p=1.0e2 rad_m1/kappa_e=1.0e2 rad_m1/kappa_f=1.0e2 \
  $IMPL time/tlim=$TLIM2 output1/dt=$ODT2

has f2 && run f2_x2_thin \
  mesh/nx1=4 mesh/nx2=64 meshblock/nx1=4 meshblock/nx2=64 problem/radwave_dir=x2 \
  rad_m1/kappa_p=1.0e2 rad_m1/kappa_e=1.0e2 rad_m1/kappa_f=1.0e2 \
  $IMPL time/tlim=$TLIM2 output1/dt=$ODT2

# (e) the 45-degree wave on a square mesh, full implicit transport
has e && run e_xy_impl2d \
  mesh/nx1=64 mesh/nx2=64 meshblock/nx1=64 meshblock/nx2=64 problem/radwave_dir=xy \
  $IMPL time/tlim=$TLIMD output1/dt=$ODTD

if has a || has b || has c || has d1 || has d2 || has e || has f1 || has f2; then
echo "=== analysis"
# the fourth field is what the arm is GATED against: empty = the analytic c_s,
# a label = that arm, `none` = reported only (the explicit scheme is an accuracy
# cross-check, not a gate).
ARGS=""
for a in "a_x1_implx1:x1:" "b_x1_impl2d:x1:a_x1_implx1" \
         "c_x2_impl2d:x2:a_x1_implx1" "d1_x1_expl2d:x1:none" \
         "d2_x2_expl2d:x2:d1_x1_expl2d" "e_xy_impl2d:xy:" \
         "f1_x1_thin:x1:none" "f2_x2_thin:x2:f1_x1_thin"; do
  n=${a%%:*}; r=${a#*:}; d=${r%%:*}; pr=${r#*:}
  [ -d "$OUT/$n/bin" ] && ARGS="$ARGS --arm $OUT/$n,$d,$n,$pr"
done
python3 "$ROOT/tests_m1/t10_radwave.py" $ARGS \
  --json "$ROOT/tests_m1/plots/impl3b6_radwave.json"
fi

# ======================================================================================
# G-SEED (phase E): the SEEDED 2-D He slab, 200 s, with the gas coupling split into its
# pieces.  build_cpu_box (-D PROBLEM=box_convection).  Run as
#   bash tests_m1/runs_3b6/run_gates_3b6.sh heslab
# ======================================================================================
if has heslab; then
  BOX=$ROOT/build_cpu_box/src/athena
  S=$OUT/he_slab_m1_2d.athinput
  SEED="time/tlim=200.0 problem/vpert=1.0e-3 time/ndiag=200 \
rad_m1/implicit_offdiag=operator rad_m1/implicit_closure_lag=step"
  for arm in "s_base:" "s_notrans:rad_m1/dbg_gas_force_trans=false" \
             "s_noforce:rad_m1/dbg_gas_force=false" \
             "s_noheat:rad_m1/dbg_gas_heat=false"; do
    n=${arm%%:*}; o=${arm#*:}
    d=$OUT/$n; rm -rf "$d"; mkdir -p "$d"
    ( cd "$d" && "$BOX" -i "$S" -d . $SEED $o > log.txt 2>&1 ) &
  done
  wait
  for n in s_base s_notrans s_noforce s_noheat; do
    echo "=== $n"; tail -3 "$OUT/$n/log.txt"
  done
fi
