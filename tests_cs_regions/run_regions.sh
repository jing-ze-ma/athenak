#!/bin/bash
# =====================================================================================
# PER-REGION ERROR GATE for the cubed sphere (cs).
#
#   ./run_regions.sh gate     baseline (plm, wb off) x {strat, rot, loop} x n=16,32,64
#   ./run_regions.sh matrix   reconstruction x well-balanced-source matrix, n=32
#   ./run_regions.sh wide     the improved configs at n=16 (and n=64 for strat)
#   ./run_regions.sh km       the K&M dynamic radial well-balancing + stale-cache arm
#   ./run_regions.sh long     n=16, one radial sound-crossing time (spurious v at rest)
#   ./run_regions.sh all      all of the above, then the tables
#   ./run_regions.sh table    (re)build the tables from the logs already present
#
# Everything is serial CPU, Release, no MPI.  Logs land in logs/, tables in tables/,
# scratch output (history files, *-cs-errs.dat) in out/ and is disposable.
#
# THE THREE TESTS (all existing cs_test problems; no new pgen code):
#   strat = inputs/tests/cubed_sphere_mhd_strat.athinput   iprob 13
#           isothermal hydrostatic atmosphere AT REST (constant radial g as a user source
#           term) carrying a uniform, force-free Cartesian field.  Exact solution = the
#           initial condition forever and v == 0, so |v| IS the error.
#   rot   = inputs/tests/cubed_sphere_mhd_conv.athinput    iprob 9
#           rigid rotation v = Omega zhat x r (pressure supplies the centripetal force)
#           carrying a uniform field that precesses with it.  Exact at every time.
#           This is the advection-across-seams-and-a-vertex test.
#   loop  = inputs/tests/cubed_sphere_toroidal.athinput    iprob 11
#           the purely azimuthal field B = b0c(-y,x,0) seeded from its VECTOR POTENTIAL
#           (faces_from_potential), balanced by p0 - b0c^2 R^2.  Static; carries a real
#           current, so it is the field-loop analogue.  Error = evolved face field vs the
#           exact static one.
#
# REGIONS.  Classified by ANGLE, not by block index, so the answer does not depend on how
# many MeshBlocks span a panel.  A cell/face is "near an edge" when it lies within NBAND
# cells of |xi| = pi/4 (or |eta| = pi/4):
#     CUBE VERTEX    near in BOTH xi and eta   (codim 2)
#     panel SEAM     near in exactly one       (codim 1)
#     panel INTERIOR neither
# NBAND is <problem>/conv_nband and is scaled with the resolution (NBAND32 at n=32, so
# 2/4/8 at n=16/32/64) so that the band has a FIXED PHYSICAL WIDTH.  A fixed CELL count
# is a region that shrinks as the grid refines, and an L1 over a shrinking region is not
# a convergence rate -- that trap once read the seam as first order.
# =====================================================================================
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HERE="$ROOT/tests_cs_regions"
ATH="$ROOT/build_cs_reg/src/athena"
LOGS="$HERE/logs"; OUT="$HERE/out"; TAB="$HERE/tables"
mkdir -p "$LOGS" "$OUT" "$TAB"

NBAND32=${NBAND32:-4}      # band half-width in cells AT n=32; scaled below
NGHOST=${NGHOST:-3}        # 3 for every run: ppm4/ppmx/wenoz need >=3, and holding it
                           # fixed keeps the halo width out of the comparison
NPROC=${NPROC:-2}          # at most 2 concurrent runs (login-node memory cap)

# THE ORDER TABLES COMPARE A FIXED PHYSICAL TIME, not a fixed cycle count, so every
# resolution has integrated the same problem.  The time is PER TEST, chosen so that the
# COARSEST grid still takes O(20) steps -- at 3 or 4 steps the answer is the start-up
# transient and the order comes out ~1.5 whatever the scheme does.  Measured dt at n=16:
#   strat 2.4e-3 (Alfven-limited at the top of the atmosphere, rho = e^-5 d0) -> 17 steps
#   rot   1.1e-2 -> 45 steps ;  loop 1.5e-2 -> 17 steps
tlim_for() {
  case "$1" in
    strat) echo "${TLIM_STRAT:-0.04}" ;;
    rot)   echo "${TLIM_ROT:-0.25}" ;;
    loop)  echo "${TLIM_LOOP:-0.25}" ;;
  esac
}
# ~1 radial sound crossing of the strat atmosphere (dr = 0.5, c_s = 0.408)
TLIM_LONG=${TLIM_LONG:-1.20}

# NOTE AthenaK REFUSES a command-line override of a parameter the input file never
# mentions, so the rot and loop inputs are copies of the stock ones with the gate's three
# knobs declared.  strat uses the stock file, which already declares all three.
inp_for() {
  case "$1" in
    strat) echo "$ROOT/inputs/tests/cs_regions_strat.athinput" ;;
    rot)   echo "$ROOT/inputs/tests/cs_regions_rot.athinput" ;;
    loop)  echo "$ROOT/inputs/tests/cs_regions_loop.athinput" ;;
  esac
}

# iprob 11 has its OWN finaliser (the per-region FACE-field table); conv_errors would
# replace it with the rigid-rotation comparison, which does not apply to it.
cerr_for() { [ "$1" = loop ] && echo 0 || echo 1; }

# one run: test recon wb n tlim [extra...]   (tlim = "-" takes the per-test default)
run_one() {
  local t=$1 recon=$2 wb=$3 n=$4 tlim=$5; shift 5
  local tag="${t}_${recon}_wb${wb}_n${n}"
  if [ "$tlim" = "-" ]; then tlim=$(tlim_for "$t"); else tag="${tag}_t${tlim}"; fi
  for x in "$@"; do tag="${tag}_$(echo "$x" | tr '/=' '__')"; done
  local log="$LOGS/$tag.log"
  if [ -s "$log" ] && grep -q "REGION" "$log"; then echo "  skip $tag"; return 0; fi
  local nb=$(( NBAND32 * n / 32 )); [ $nb -lt 1 ] && nb=1
  echo "  run  $tag  (nband=$nb)"
  ( cd "$OUT" && /usr/bin/time -f "WALLCLOCK %e s" "$ATH" -i "$(inp_for "$t")" -d . \
      job/basename="$tag" \
      mesh/nghost=$NGHOST \
      mesh/nx1=$n mesh/nx2=$n mesh/nx3=$n \
      meshblock/nx1=$n meshblock/nx2=$n meshblock/nx3=$n \
      time/tlim=$tlim time/nlim=-1 \
      mhd/reconstruct=$recon \
      mhd/cs_wellbalanced_src=$( [ "$wb" = 1 ] && echo true || echo false ) \
      problem/conv_errors=$(cerr_for "$t") problem/conv_nband=$nb \
      "$@" ) > "$log" 2>&1
  local rc=$?
  [ $rc -ne 0 ] && echo "  *** FAILED rc=$rc : $tag"
  return 0
}

# run a list of "test recon wb n tlim [extra]" lines, NPROC at a time.  The slot count is
# enforced by counting LIVE children, not by `wait -n`: a bare `wait -n` on this bash
# returns as soon as ANY child is reaped, which let the batch drift up to five concurrent
# solvers on the login node the first time this ran.
run_batch() {
  while read -r line; do
    [ -z "$line" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do sleep 2; done
    run_one $line &
  done
  wait
}

do_gate() {
  echo "== GATE: baseline plm, wb off, three tests, n = 16 32 64"
  run_batch <<EOF
strat plm 0 16 -
rot plm 0 16 -
loop plm 0 16 -
strat plm 0 32 -
rot plm 0 32 -
loop plm 0 32 -
strat plm 0 64 -
rot plm 0 64 -
loop plm 0 64 -
EOF
}

do_matrix() {
  echo "== MATRIX: reconstruct x cs_wellbalanced_src at n = 32"
  local lines=""
  for t in strat rot loop; do
    for r in plm ppm4 ppmx wenoz; do
      for w in 0 1; do
        lines="$lines$t $r $w 32 -
"
      done
    done
  done
  echo "$lines" | run_batch
}

do_wide() {
  echo "== WIDE: the improved configs at n = 16 and 64 (n=32 comes from the matrix)"
  run_batch <<EOF
strat plm 1 16 -
strat wenoz 0 16 -
strat wenoz 1 16 -
rot plm 1 16 -
rot wenoz 0 16 -
rot wenoz 1 16 -
loop wenoz 1 16 -
loop plm 1 16 -
strat plm 1 64 -
strat wenoz 1 64 -
EOF
}

# The K&M (dynamic) well-balanced RADIAL reconstruction, and what a STALE background cache
# does to it.  Only the strat test has a potential (cs_test fills it for iprob 13).
do_km() {
  echo "== K&M: <mhd>/wellbalance_dynamic + wb_x1, and wb_cache_every = 0 / 1 / 10"
  run_batch <<EOF
strat plm 0 16 - mhd/wellbalance_dynamic=true mhd/wb_x1=true
strat plm 0 32 - mhd/wellbalance_dynamic=true mhd/wb_x1=true
strat plm 0 64 - mhd/wellbalance_dynamic=true mhd/wb_x1=true
strat plm 0 16 $TLIM_LONG mhd/wellbalance_dynamic=true mhd/wb_x1=true mhd/wb_cache_every=1
strat plm 0 16 $TLIM_LONG mhd/wellbalance_dynamic=true mhd/wb_x1=true mhd/wb_cache_every=10
strat plm 0 16 $TLIM_LONG mhd/wellbalance_dynamic=true mhd/wb_x1=true
EOF
}

do_long() {
  echo "== LONG: strat, n = 16, tlim = $TLIM_LONG (~1 radial sound crossing)"
  run_batch <<EOF
strat plm 0 16 $TLIM_LONG
strat plm 1 16 $TLIM_LONG
strat wenoz 0 16 $TLIM_LONG
strat wenoz 1 16 $TLIM_LONG
EOF
}

case "${1:-all}" in
  gate)   do_gate ;;
  matrix) do_matrix ;;
  wide)   do_wide ;;
  long)   do_long ;;
  table)  ;;
  km)     do_km ;;
  all)    do_gate; do_matrix; do_wide; do_km; do_long ;;
  *) echo "usage: $0 [gate|matrix|wide|km|long|table|all]"; exit 1 ;;
esac

python3 "$HERE/regions.py" "$LOGS" > "$TAB/regions.txt"
echo "== tables written to $TAB/regions.txt"
cat "$TAB/regions.txt"
