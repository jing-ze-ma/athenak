#!/bin/bash -l
# m1-sph2 validation runs (new binary only).  usage: val_cpu.sh <group ...>
W=/viper/ptmp2/jinma/sph2_0924; I=$W/inp; X=$W/bin/athena_new_none_cpu
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
run() {  # run <name> <np> <input> [args]
  local n=$1 np=$2 inp=$3; shift 3
  local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $X -i $inp -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
rrun() {  # rrun <name> <np> <rst> [args]
  local n=$1 np=$2 r=$3; shift 3
  local d=$W/cpu/$n; rm -rf $d; mkdir -p $d
  (cd $d && OMP_NUM_THREADS=1 nice mpirun -np $np --oversubscribe --bind-to none \
     $X -r $r -d $d "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) &
}
A=$W/scripts/ana.py
for g in "$@"; do case $g in
T5)  # radwave along theta in a shell at R = 100, error vs dt at fixed resolution
  O="time/nlim=-1 output1/dt=0.592156525463792"
  for c in 0.4 0.2 0.1 0.05; do
    run t5_be_$c 1 $I/rw_be.athinput $O time/cfl_number=$c
    run t5_h2_$c 1 $I/rw_h2.athinput $O time/cfl_number=$c
    run t5_h2nv_$c 1 $I/rw_h2nv.athinput $O time/cfl_number=$c
    run t5m_be_$c 1 $I/rwm_be.athinput $O time/cfl_number=$c
    run t5m_h2_$c 1 $I/rwm_h2.athinput $O time/cfl_number=$c
  done
  run t5_ref 1 $I/rw_h2.athinput $O time/cfl_number=0.0125
  run t5m_ref 1 $I/rwm_h2.athinput $O time/cfl_number=0.0125
  wait
  for a in be h2 h2nv; do echo "T5 eddington $a: $(python3 $A order 'tab/*.tab' dens $W/cpu/t5_ref $W/cpu/t5_${a}_{0.4,0.2,0.1,0.05} --pert 1.0)"; done
  for a in be h2; do echo "T5 m1 $a: $(python3 $A order 'tab/*.tab' dens $W/cpu/t5m_ref $W/cpu/t5m_${a}_{0.4,0.2,0.1,0.05} --pert 1.0)"; done;;
T1)  # transient spherical diffusion (T-S1 shell, stretched r, tau 2000), t = 100
  O="time/nlim=-1 time/tlim=100.0"
  for cl in eddington m1 vet_col; do
    for c in 100.0 50.0 25.0 12.5; do
      run t1_${cl}_be_$c 4 $I/sym_be.athinput $O rad_m1/implicit_cfl=$c rad_m1/closure=$cl
      run t1_${cl}_h2_$c 4 $I/sym_h2.athinput $O rad_m1/implicit_cfl=$c rad_m1/closure=$cl
    done
    run t1_${cl}_ref 4 $I/sym_h2.athinput $O rad_m1/implicit_cfl=3.125 rad_m1/closure=$cl
    wait
  done
  for cl in eddington m1 vet_col; do for a in be h2; do
    echo "T1 $cl $a: $(python3 $A order 'tab/sd.j0k0.*.tab' m1_e $W/cpu/t1_${cl}_ref $W/cpu/t1_${cl}_${a}_{100.0,50.0,25.0,12.5} --pert 0.0)"; done; done;;
TV)  # the T-S4 atmosphere relaxing from the Eddington start (tensor changes), t = 6.4
  O="time/nlim=-1 time/tlim=6.4"
  for c in 0.3 0.15 0.075 0.0375; do
    run tv_be_$c 1 $I/atmvc_be.athinput $O time/cfl_number=$c
    run tv_h2_$c 1 $I/atmvc_h2.athinput $O time/cfl_number=$c
    run tvm_be_$c 1 $I/atm_be.athinput $O time/cfl_number=$c mesh/nx1=32 meshblock/nx1=32
    run tvm_h2_$c 1 $I/atm_h2.athinput $O time/cfl_number=$c mesh/nx1=32 meshblock/nx1=32
  done
  run tv_ref 1 $I/atmvc_h2.athinput $O time/cfl_number=0.009375
  run tvm_ref 1 $I/atm_h2.athinput $O time/cfl_number=0.009375 mesh/nx1=32 meshblock/nx1=32
  wait
  for a in be h2; do echo "TV vet_col $a: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/tv_ref $W/cpu/tv_${a}_{0.3,0.15,0.075,0.0375})"; done
  for a in be h2; do echo "TV m1 $a: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/tvm_ref $W/cpu/tvm_${a}_{0.3,0.15,0.075,0.0375})"; done;;
TV2)  # the same atmosphere pure scattering (no gas coupling): E and the vet_col tensor
  # relax from the Eddington start, t = 6.4
  O="time/nlim=-1 time/tlim=6.4 rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=100.0 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_lin_tol=1.0e-13"
  for c in 0.3 0.15 0.075 0.0375; do
    run tw_be_$c 1 $I/atmvc_be.athinput $O time/cfl_number=$c
    run tw_h2_$c 1 $I/atmvc_h2.athinput $O time/cfl_number=$c
  done
  run tw_ref 1 $I/atmvc_h2.athinput $O time/cfl_number=0.009375
  wait
  for a in be h2; do echo "TV2 vet_col $a: $(python3 $A order 'tab/sa.m1.*.tab' m1_e $W/cpu/tw_ref $W/cpu/tw_${a}_{0.3,0.15,0.075,0.0375})"; done
  grep -h "NON-CONVERGED\|time_scheme=hesdirk2:" $W/cpu/tw_h2_0.3/log.txt;;
SS)  # steady states under the new default (key absent: hesdirk2) vs the be runs of gate C
  run ss_sd 1 $I/sp_sph_diff.athinput
  run ss_sds 1 $I/sp_sph_diff_str.athinput
  run ss_fs 1 $I/sp_sph_fs.athinput
  run ss_atm 1 $I/sp_sph_atm.athinput mesh/nx1=32 meshblock/nx1=32
  run ss_vatm 1 $I/sp_sph_atm_vc.athinput mesh/nx1=64 meshblock/nx1=64
  run ss_vstr 1 $I/sp_sph_atm_str_vc.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=200
  run ss_vpp 1 $I/milne_vc.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=200
  wait
  for p in "sd C_sd" "sds C_sds" "fs C_fs" "atm C_atm" "vatm C_vatm" "vstr C_vstr" "vpp C_vpp"; do set -- $p
    echo "SS $1: $(python3 $A steady 'tab/*m1*.tab' $W/cpu/ss_$1 $W/cpu/$2/new) | $(grep -h 'time_scheme=hesdirk2:' $W/cpu/ss_$1/log.txt) $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/ss_$1/log.txt)"; done;;
SYM)  # T-sym under hesdirk2 (the default), 4 ranks, 500 steps
  run sym_edd 4 $I/sp_sph_sym.athinput rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_precond=line
  run sym_vc 4 $I/sp_sph_sym.athinput rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_precond=line rad_m1/closure=vet_col
  run sym_m1 4 $I/sp_sph_sym.athinput rad_m1/implicit_lin_tol=1.0e-14 rad_m1/implicit_tol=1.0e-12 rad_m1/implicit_precond=line rad_m1/closure=m1
  run symg_edd 4 $I/sp_sph_sym_gas.athinput rad_m1/c_light=100.0
  run symg_vc 4 $I/sp_sph_sym_gas.athinput rad_m1/c_light=100.0 rad_m1/closure=vet_col
  wait
  for d in sym_edd sym_vc sym_m1 symg_edd symg_vc; do echo "SYM $d: $(python3 $W/scripts/sym.py $W/cpu/$d) $(grep -h 'time_scheme=hesdirk2:\|implicit_vimp=true' $W/cpu/$d/log.txt | tr '\n' ' ') $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/$d/log.txt)"; done;;
THIN)  # the runs_5b Finding-1 reproducer with vet_col under hesdirk2 (seed 1e-6)
  N2="mesh/nx1=32 meshblock/nx1=32 time/nlim=400 rad_m1/implicit_lin_tol=1.0e-10 rad_m1/implicit_tol=1.0e-10 rad_m1/implicit_maxit=100 output3/dcycle=50 problem/atm_seed=1.0e-6"
  S="rad_m1/kappa_p=0.0 rad_m1/kappa_e=0.0 rad_m1/kappa_f=0.0 rad_m1/kappa_s=25.0"
  C="mesh/use_spherical_polar=false mesh/x2min=0 mesh/x2max=0.2 output1/slice_x2=0.1 output2/slice_x2=0.1"
  run zs_vc 1 $I/sp_sph_atm_bin_vc.athinput $N2 $S rad_m1/closure=vet_col
  run zs_vcF 1 $I/sp_sph_atm_bin_vc.athinput $N2 $S rad_m1/closure=vet_col rad_m1/vet_col_axis=flux
  run zc_vc 1 $I/sp_sph_atm_bin_vc.athinput $N2 $S $C rad_m1/closure=vet_col
  wait
  for d in zs_vc zs_vcF zc_vc; do echo "THIN $d: $(python3 $W/scripts/spread.py $W/cpu/$d) $(grep -h 'time_scheme=hesdirk2:' $W/cpu/$d/log.txt) $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/$d/log.txt)"; done;;
R6)  # risk 6: static pure-scattering atmosphere, gas momentum on; be vs hesdirk2 (+vimp)
  run r6A 1 $I/r6_be.athinput problem/atm_temp=1.0 time/cfl_number=100.0 output4/dcycle=1000
  wait
  rrun r6A3 1 $W/cpu/r6A/rst/r6.00001.rst time/nlim=1010 time/cfl_number=0.3 output3/dcycle=10
  wait
  R=$(ls $W/cpu/r6A3/rst/r6.*.rst | tail -1)
  B="time/nlim=1110 time/cfl_number=0.3 problem/atm_hold=false output1/dcycle=50 output2/dcycle=50 output4/dcycle=10 rad_m1/dbg_gas_force=true rad_m1/force_reference=wb_arad"
  rrun r6B_be 1 $R $B
  rrun r6B_h2 1 $R $B rad_m1/time_scheme=hesdirk2
  rrun r6B_h2v 1 $R $B rad_m1/time_scheme=hesdirk2 rad_m1/implicit_vimp=true
  rrun r6B_bev 1 $R $B rad_m1/implicit_vimp=true
  wait
  for d in r6B_be r6B_h2 r6B_h2v r6B_bev; do echo "R6 $d: $(python3 $W/scripts/vrms.py $W/cpu/$d) $(grep -h 'time_scheme=hesdirk2:\|implicit_vimp=true\|vimp positivity' $W/cpu/$d/log.txt | tr '\n' ' ') $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/$d/log.txt) $(tail -n1 $W/cpu/$d/log.txt)"; done;;
HE)  # the He wedge grid (hewedge, stretched r), reduced to 96 x 32 x 32, gas free, 60 steps
  H="mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60"
  run he_be 4 $I/hw_be.athinput $H
  run he_h2 4 $I/hw_h2.athinput $H
  run he_h2e 4 $I/hw_h2.athinput $H rad_m1/closure=eddington
  wait
  for d in he_be he_h2 he_h2e; do echo "HE $d: $(grep -ho 'NON-CONVERGED=[^ ]*' $W/cpu/$d/log.txt) $(grep -h 'time_scheme=hesdirk2:\|implicit_vimp=true\|vimp positivity\|NOT ADMISS' $W/cpu/$d/log.txt | sort | uniq -c | tr '\n' ' ') $(tail -n1 $W/cpu/$d/log.txt)"; done;;
esac; done
echo "VAL $* DONE"
