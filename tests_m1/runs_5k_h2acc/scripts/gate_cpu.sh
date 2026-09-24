#!/bin/bash -l
# m1-h2acc CPU gates, base (d2572b4f) vs new.  usage: gate_cpu.sh [A|B|D|S ...]
W=/viper/ptmp2/jinma/h2acc_0924; I=/viper/ptmp2/jinma/m1def2_0924/inp; K=/viper/ptmp2/jinma/h2div_0924/inp
Q=/viper/ptmp2/jinma/sph2_0924/inp; S=$W/scripts
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
exe() { echo $W/bin/athena_${1%%_*}_$2_cpu; }
run() {  # run <name> <tag> <box|none> <np> <input> [args]
  local n=$1 t=$2 p=$3 np=$4 inp=$5; shift 5
  local d=$W/cpu/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $(exe $t $p) -i $inp -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
rrun() { local n=$1 t=$2 p=$3 np=$4 r=$5; shift 5
  local d=$W/cpu/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $(exe $t $p) -r $r -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
info() { grep -h "stage fallbacks\|NON-CONVERGED=" $1/log.txt | sed 's/.*\(NON-CONVERGED=[^ ]*\).*/\1/; s/.*\(stage fallbacks=[^ ]*\).*/\1/' | tr '\n' ' '; tail -n1 $1/log.txt; }
cmpd() { echo "== $1: $(python3 $S/gcmp.py $W/cpu/$1/${2:-base} $W/cpu/$1/${3:-new} | tail -n1) | ${2:-base}: $(info $W/cpu/$1/${2:-base}) | ${3:-new}: $(info $W/cpu/$1/${3:-new})"; }
BOX="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16"
LEV="rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2"
VP="problem/vpert=1.0e-2"
for g in "$@"; do case $g in
A)  # be named: bitwise
  for t in base new; do
    run A_box1 $t box 1 $I/box3d_be_x.athinput $BOX time/nlim=30 $VP $LEV
    run A_slabv2 $t box 2 $I/slab2d_plm_vimp_be_x.athinput meshblock/nx2=16 time/nlim=100 $VP $LEV rad_m1/closure=vet_sc rad_m1/vet_tensor=full
    run A_nd2 $t box 2 $I/slab2d_nd_be.athinput meshblock/nx2=16 time/tlim=200 $VP
    for f in cart_sym marshak_cart rw_cart; do run A_$f $t none 1 $I/${f}_be.athinput; done
    run A_vatm $t none 1 $Q/sp_sph_atm_vc_be.athinput time/nlim=300
  done; wait
  for c in A_box1 A_slabv2 A_nd2 A_cart_sym A_marshak_cart A_rw_cart A_vatm; do cmpd $c; done;;
B)  # hesdirk2 by default: bitwise where no VET closure reads a gas state; VET arms with old
  for t in base new; do
    run B_nd2 $t box 2 $I/slab2d_nd.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP
    run B_boxnd2 $t box 2 $I/box3d_nd.athinput $BOX time/nlim=30 $VP
    for f in cart_sym marshak_cart rw_cart; do run B_$f $t none 1 $I/$f.athinput; done
  done
  run B_boxv2 base box 2 $I/box3d_nd.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc
  run B_boxv2 new box 2 $I/box3d_nd.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc rad_m1/time2_vet_gas=old
  run B_boxv2 newS box 2 $I/box3d_nd.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc
  run B_vatm base none 1 $Q/sp_sph_atm_vc.athinput time/nlim=300
  run B_vatm new none 1 $Q/sp_sph_atm_vc.athinput time/nlim=300 rad_m1/time2_vet_gas=old
  run B_vatm newS none 1 $Q/sp_sph_atm_vc.athinput time/nlim=300
  wait
  for c in B_nd2 B_boxnd2 B_cart_sym B_marshak_cart B_rw_cart B_boxv2 B_vatm; do cmpd $c; done
  cmpd B_boxv2 base newS; cmpd B_vatm base newS;;
D)  # restart bitwise (new): slab2d Eddington and slab2d vet_sc, rst at t = 100
  run D_src new box 2 $K/slab2d_nd_k.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP
  run D_srcv new box 2 $K/slab2d_nd_k.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP rad_m1/closure=vet_sc
  wait
  rrun D_h2 new box 2 $W/cpu/D_src/new/rst/m1slab.00001.rst time/tlim=200
  rrun D_h2v new box 2 $W/cpu/D_srcv/new/rst/m1slab.00001.rst time/tlim=200
  wait
  echo "== D_h2 Eddington: straight $(info $W/cpu/D_src/new) cont $(info $W/cpu/D_h2/new)"
  python3 $S/rstcmp.py $W/cpu/D_src/new $W/cpu/D_h2/new
  echo "== D_h2v vet_sc: straight $(info $W/cpu/D_srcv/new) cont $(info $W/cpu/D_h2v/new)"
  grep -h "time2_vet_gas" $W/cpu/D_h2v/new/log.txt | head -2
  python3 $S/rstcmp.py $W/cpu/D_srcv/new $W/cpu/D_h2v/new;;
esac; done
echo "GATES $* DONE"
