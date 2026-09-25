#!/bin/bash -l
# m1-sp-order2b CPU gates: ref (fe12a518) vs new, the inputs of runs_5h_sph2
# (/viper/ptmp2/jinma/sph2_0924/inp, read only) plus copies with the new keys.
# usage: gate_cpu.sh [A|B|C|CO|F|FO|O|R ...]
# A, B: Cartesian (bitwise).  C, F: sp with the keys absent (the new sp defaults: a
# DIFF is legitimate where a Marshak face or vet_col is involved).  CO, FO: the same
# inputs with the OLD values named (bitwise).  O: an OLD restart (ref-written, keys
# absent) continued by ref and new (bitwise).  R: restart with the new defaults.
W=/viper/ptmp2/jinma/sporder2b_0925; I=/viper/ptmp2/jinma/sph2_0924/inp; J=$W/inp
S=$(dirname $(readlink -f $0))
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() {  # run <name> <ref|new> <box|none> <np> <input> [args]
  local n=$1 t=$2 p=$3 np=$4 inp=$5; shift 5
  local d=$W/gate/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $W/bin/athena_${t}_${p}_cpu -i $inp -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
rrun() {  # rrun <name> <ref|new> <box|none> <np> <rstfile> [args]
  local n=$1 t=$2 p=$3 np=$4 r=$5; shift 5
  local d=$W/gate/$n/$t; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $W/bin/athena_${t}_${p}_cpu -r $r -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
cmpd() { echo "== $1: $(python3 $S/gcmp.py $W/gate/$1/ref $W/gate/$1/new | tail -n1) $(tail -qn1 $W/gate/$1/ref/log.txt $W/gate/$1/new/log.txt | tr '\n' ' ')"; }
BOX="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16"
LEV="rad_m1/implicit_vimp_fold=true rad_m1/implicit_fast_kernels=true rad_m1/implicit_one_pass=4 rad_m1/implicit_predictor_order=2"
OVL="rad_m1/implicit_halo_mpi=true rad_m1/implicit_halo_overlap=true rad_m1/implicit_halo_ovl_faces=true"
VP="problem/vpert=1.0e-2"
mkdir -p $J
for g in "$@"; do case $g in
A)  # Cartesian, be named: bitwise
  for t in ref new; do
    run A_box1 $t box 1 $I/box3d_be_x.athinput $BOX time/nlim=30 $VP $LEV
    run A_box2o $t box 2 $I/box3d_be_x.athinput $BOX time/nlim=30 $VP $LEV $OVL
    run A_slab1 $t box 1 $I/slab2d_plm_vimp_be_x.athinput time/nlim=100 $VP $LEV
    run A_slabv2 $t box 2 $I/slab2d_plm_vimp_be_x.athinput meshblock/nx2=16 time/nlim=100 $VP $LEV rad_m1/closure=vet_sc rad_m1/vet_tensor=full
    run A_nd2 $t box 2 $I/slab2d_nd_be.athinput meshblock/nx2=16 time/tlim=200 $VP
    run A_boxnd2 $t box 2 $I/box3d_nd_be.athinput $BOX time/nlim=30 $VP
    for f in cart_sym marshak_cart rw_cart pp_np; do run A_$f $t none 1 $I/${f}_be.athinput; done
  done; wait
  for c in A_box1 A_box2o A_slab1 A_slabv2 A_nd2 A_boxnd2 A_cart_sym A_marshak_cart A_rw_cart A_pp_np; do cmpd $c; done;;
B)  # Cartesian, key absent (hesdirk2 + vimp defaults): bitwise
  for t in ref new; do
    run B_nd2 $t box 2 $I/slab2d_nd.athinput meshblock/nx2=16 time/tlim=200 output3/dt=100 $VP
    run B_boxnd2 $t box 2 $I/box3d_nd.athinput $BOX time/nlim=30 $VP
    run B_boxv2 $t box 2 $I/box3d_nd.athinput $BOX time/nlim=20 $VP rad_m1/closure=vet_sc
    for f in cart_sym marshak_cart rw_cart milne_vc; do run B_$f $t none 1 $I/$f.athinput; done
  done; wait
  for c in B_nd2 B_boxnd2 B_boxv2 B_cart_sym B_marshak_cart B_rw_cart B_milne_vc; do cmpd $c; done;;
C)  # spherical-polar S1/S2/S5 + Cartesian vet_col, time_scheme = be named: bitwise
  for t in ref new; do
    run C_sd $t none 1 $I/sp_sph_diff_be.athinput
    run C_sds $t none 1 $I/sp_sph_diff_str_be.athinput
    run C_sym $t none 4 $I/sp_sph_sym_be.athinput time/nlim=200
    run C_symg $t none 4 $I/sp_sph_sym_gas_be.athinput time/nlim=200
    run C_mk $t none 1 $I/sp_marshak_sph_be.athinput
    run C_fs $t none 1 $I/sp_sph_fs_be.athinput
    run C_atm $t none 1 $I/sp_sph_atm_be.athinput mesh/nx1=32 meshblock/nx1=32
    run C_rw $t none 1 $I/sp_rw_sph_R100_be.athinput
    run C_msym $t none 4 $I/sp_sph_sym_be.athinput time/nlim=200 rad_m1/closure=m1 rad_m1/implicit_precond=line
    run C_vatm $t none 1 $I/sp_sph_atm_vc_be.athinput mesh/nx1=64 meshblock/nx1=64
    run C_vsym $t none 4 $I/sp_sph_sym_be.athinput time/nlim=200 rad_m1/closure=vet_col
    run C_vsymg $t none 4 $I/sp_sph_sym_gas_be.athinput time/nlim=200 rad_m1/closure=vet_col rad_m1/c_light=100.0
    run C_vstr $t none 1 $I/sp_sph_atm_str_vc_be.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=200
    run C_vpp $t none 1 $I/milne_vc_be.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=200
    run C_vbin $t none 1 $I/sp_sph_atm_bin_vc_be.athinput time/nlim=100
  done; wait
  for c in C_sd C_sds C_sym C_symg C_mk C_fs C_atm C_rw C_msym C_vatm C_vsym C_vsymg C_vstr C_vpp C_vbin; do cmpd $c; done;;
F)  # spherical-polar, key absent (the hesdirk2 default on the wedge): bitwise
  for t in ref new; do
    run F_sym $t none 4 $I/sp_sph_sym.athinput time/nlim=200
    run F_symg $t none 4 $I/sp_sph_sym_gas.athinput time/nlim=200
    run F_vsym $t none 4 $I/sp_sph_sym.athinput time/nlim=200 rad_m1/closure=vet_col
    run F_vatm $t none 1 $I/sp_sph_atm_vc.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=300
    run F_fs $t none 1 $I/sp_sph_fs.athinput
  done; wait
  for c in F_sym F_symg F_vsym F_vatm F_fs; do cmpd $c; done;;
CO|FO)  # the old sp behaviour named: bitwise
  OLDK="implicit_marshak_face=cell vet_col_order2=false vet_col_fk_min=0.3333333333333333 vet_col_reflect_top=false vet_col_surface_face=false time2_vet_col=lag time2_vstage=false"
  if [ $g = CO ]; then L="sp_sph_diff_be sp_sph_diff_str_be sp_sph_sym_be sp_sph_sym_gas_be sp_marshak_sph_be sp_sph_fs_be sp_sph_atm_be sp_rw_sph_R100_be sp_sph_atm_vc_be sp_sph_atm_str_vc_be sp_sph_atm_bin_vc_be"
  else L="sp_sph_sym sp_sph_sym_gas sp_sph_atm_vc sp_sph_fs"; fi
  for f in $L; do python3 $S/addkeys.py $I/$f.athinput $J/o_$f.athinput rad_m1 $OLDK; done
  for t in ref new; do
    for f in $L; do
      a=""; np=1
      case $f in *sym*) a="time/nlim=200"; np=4;; *atm_vc*) a="mesh/nx1=64 meshblock/nx1=64 time/nlim=300";; *atm_be*|*str_vc*) a="mesh/nx1=32 meshblock/nx1=32 time/nlim=200";; *bin_vc*) a="time/nlim=100";; esac
      run ${g}_$f $t none $np $J/o_$f.athinput $a
      case $f in *sym*) run ${g}_vc_$f $t none $np $J/o_$f.athinput $a rad_m1/closure=vet_col rad_m1/c_light=100.0;;
                 *sph_fs*) run ${g}_vc_$f $t none $np $J/o_$f.athinput rad_m1/closure=vet_col;; esac
    done
  done; wait
  for d in $W/gate/${g}_*; do cmpd $(basename $d); done;;
O)  # an OLD restart (written by ref, keys absent) continued by ref and by new: bitwise
  for c in vc vcr ed; do
    a="rad_m1/c_light=100.0"
    [ $c = vc ] && a="$a rad_m1/closure=vet_col"
    [ $c = vcr ] && a="$a rad_m1/closure=vet_col rad_m1/implicit_bc_x1max=reflect"
    run O_src_$c ref none 4 $I/sph_sym_gas_rst.athinput $a
  done; wait
  for c in vc vcr ed; do
    for t in ref new; do rrun O_$c $t none 4 $W/gate/O_src_$c/ref/rst/sd.00001.rst; done
  done; wait
  for c in vc vcr ed; do cmpd O_$c; grep -h "NON-CONV" $W/gate/O_$c/new/log.txt | cut -c1-140; done;;
R)  # restart with the NEW defaults (keys absent): straight vs restart at cycle 250
  for c in vc vcr ed; do
    a="rad_m1/c_light=100.0"
    [ $c = vc ] && a="$a rad_m1/closure=vet_col"
    [ $c = vcr ] && a="$a rad_m1/closure=vet_col rad_m1/implicit_bc_x1max=reflect"
    run R_src_$c new none 4 $I/sph_sym_gas_rst.athinput $a
  done; wait
  for c in vc vcr ed; do rrun R_rst_$c new none 4 $W/gate/R_src_$c/new/rst/sd.00001.rst; done; wait
  for c in vc vcr ed; do
    python3 $S/rstcmp.py $W/gate/R_src_$c/new $W/gate/R_rst_$c/new
    grep -h "NON-CONV\|time2_vet_col=\|reflecting top" $W/gate/R_src_$c/new/log.txt | cut -c1-160 | head -3
    grep -h "time2_vet_col\|vet_col_order2\|marshak_face\|fk_min\|reflect_top" $W/gate/R_rst_$c/new/log.txt | head -3
  done;;
esac; done
