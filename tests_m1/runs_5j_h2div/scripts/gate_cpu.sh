#!/bin/bash -l
# m1-h2div CPU gates, base (7c16ff86) vs fix.  usage: gate_cpu.sh [A|B|F|D ...]
# inputs: the runs_5g_defaults2 copies, read in place (/viper/ptmp2/jinma/m1def2_0924/inp)
W=/viper/ptmp2/jinma/h2div_0924; I=/viper/ptmp2/jinma/m1def2_0924/inp; K=$W/inp
S=/viper/ptmp2/jinma/wt_h2div/tests_m1/runs_5j_h2div/scripts
CMP=/viper/ptmp2/jinma/wt_h2div/tests_m1/gates/cmp.py
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
exe() { local b=${1%T}; [ $2 = box ] && echo $W/bin/athena_${b}_box_cpu || echo $W/bin/athena_${b}_cpu; }
run() {  # run <name> <base|fix> <box|none> <np> <input> [args]
  local n=$1 t=$2 p=$3 np=$4 inp=$5; shift 5
  local d=$W/cpu/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $(exe $t $p) -i $inp -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
rrun() {  # rrun <name> <tag> <box|none> <np> <rstfile> [args]
  local n=$1 t=$2 p=$3 np=$4 r=$5; shift 5
  local d=$W/cpu/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $(exe ${t%%_*} $p) -r $r -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
info() { grep -h "stage fallbacks\|NON-CONVERGED=" $1/log.txt | sed 's/.*\(NON-CONVERGED=[^ ]*\).*/\1/; s/.*\(stage fallbacks=[^ ]*\).*/\1/' | tr '\n' ' '; tail -n1 $1/log.txt; }
cmpd() { echo "== $1: $(python3 $S/gcmp.py $W/cpu/$1/base $W/cpu/$1/fix | tail -n1) | base: $(info $W/cpu/$1/base) | fix: $(info $W/cpu/$1/fix)"; }
BOX="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16"
LEV="rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2"
VP="problem/vpert=1.0e-2"
TIGHT="rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-11 rad_m1/time2_lin_tol_fac=1"
for g in "$@"; do case $g in
A)  # be named: bitwise
  for t in base fix; do
    run A_box1 $t box 1 $I/box3d_be_x.athinput $BOX time/nlim=30 $VP $LEV
    run A_slabv2 $t box 2 $I/slab2d_plm_vimp_be_x.athinput meshblock/nx2=16 time/nlim=100 $VP $LEV rad_m1/closure=vet_sc rad_m1/vet_tensor=full
    run A_nd2 $t box 2 $I/slab2d_nd_be.athinput meshblock/nx2=16 time/tlim=200 $VP
    for f in cart_sym marshak_cart rw_cart; do run A_$f $t none 1 $I/${f}_be.athinput; done
  done; wait
  for c in A_box1 A_slabv2 A_nd2 A_cart_sym A_marshak_cart A_rw_cart; do cmpd $c; done;;
B)  # hesdirk2 by default (key absent): bitwise where no stage fails
  for t in base fix; do
    run B_nd2 $t box 2 $I/slab2d_nd.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP
    run B_boxnd2 $t box 2 $I/box3d_nd.athinput $BOX time/nlim=30 $VP
    run B_boxv2 $t box 2 $I/box3d_nd.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc
    for f in cart_sym marshak_cart rw_cart; do run B_$f $t none 1 $I/$f.athinput; done
  done; wait
  for c in B_nd2 B_boxnd2 B_boxv2 B_cart_sym B_marshak_cart B_rw_cart; do cmpd $c; done;;
F)  # a FORCED stage-1 failure (time2_dbg_fail = cycle): base, fix and a tight fix reference
  for t in base fix; do
    run F_nd2 $t box 2 $K/slab2d_nd_k.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP rad_m1/time2_dbg_fail=40
    run F_boxv2 $t box 2 $K/box3d_nd_k.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc rad_m1/time2_dbg_fail=8
    run F_box2 $t box 2 $K/box3d_nd_k.athinput $BOX time/nlim=30 $VP rad_m1/time2_dbg_fail=8
  done
  run F_nd2 fixT box 2 $K/slab2d_nd_k.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP rad_m1/time2_dbg_fail=40 $TIGHT
  run F_boxv2 fixT box 2 $K/box3d_nd_k.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc rad_m1/time2_dbg_fail=8 $TIGHT
  run F_box2 fixT box 2 $K/box3d_nd_k.athinput $BOX time/nlim=30 $VP rad_m1/time2_dbg_fail=8 $TIGHT
  run F_box2 baseT box 2 $K/box3d_nd_k.athinput $BOX time/nlim=30 $VP rad_m1/time2_dbg_fail=8 $TIGHT
  wait
  for c in F_nd2 F_boxv2 F_box2; do
    echo "== $c | base: $(info $W/cpu/$c/base) | fix: $(info $W/cpu/$c/fix) | fixT: $(info $W/cpu/$c/fixT)"
    echo "   fix vs base: $(python3 $CMP $W/cpu/$c/base $W/cpu/$c/fix)"
    echo "   base vs fixT: $(python3 $CMP $W/cpu/$c/fixT $W/cpu/$c/base)"
    echo "   fix  vs fixT: $(python3 $CMP $W/cpu/$c/fixT $W/cpu/$c/fix)"
  done
  echo "   F_box2 baseT vs fixT: $(python3 $CMP $W/cpu/F_box2/fixT $W/cpu/F_box2/baseT)";;
D)  # restart bitwise across a forced fallback (fix): straight F_nd2/fix vs restart at t=100,
    # with the forced failure after the restart point
  run D_src fix box 2 $K/slab2d_nd_k.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP rad_m1/time2_dbg_fail=900
  wait
  rrun D_h2 fix box 2 $W/cpu/D_src/fix/rst/m1slab.00001.rst time/tlim=200
  wait
  echo "== D_h2: fallbacks straight $(info $W/cpu/D_src/fix) cont $(info $W/cpu/D_h2/fix)"
  python3 $S/rstcmp.py $W/cpu/D_src/fix $W/cpu/D_h2/fix;;
esac; done
echo "GATES $* DONE"
