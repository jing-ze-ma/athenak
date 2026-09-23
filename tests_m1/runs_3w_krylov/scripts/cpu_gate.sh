#!/bin/bash
# CPU gates of m1-krylov (switches off bitwise vs 89dc28d7; new paths; restart)
W=/viper/ptmp2/jinma/krylov_0923; R=$W/scripts/run.sh
md5sum $W/bin/athena_ref_cpu $W/bin/athena_new_cpu
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=200"
C="hydro/rad_flux_inner=2475202000000000.5"
P="rad_m1/implicit_krylov_pipe=true"
H="rad_m1/implicit_halo_mpi=true"
X2="meshblock/nx2=16"; X4="meshblock/nx2=8"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60"
# --- 2-D slab, Eddington
$R e_ref   ref 1 slab2d_def $S &
$R e_off   new 1 slab2d_def $S &
$R e_ctl   new 1 slab2d_def $S $C &
$R e_pipe  new 1 slab2d_def $S $P &
$R e_ref2  ref 2 slab2d_def $S $X2 &
$R e_off2  new 2 slab2d_def $S $X2 &
$R e_hm2   new 2 slab2d_def $S $X2 $H &
$R e_ph2   new 2 slab2d_def $S $X2 $H $P &
$R e_ref4  ref 4 slab2d_def $S $X4 &
$R e_off4  new 4 slab2d_def $S $X4 &
$R e_hm4   new 4 slab2d_def $S $X4 $H &
$R e_ph4   new 4 slab2d_def $S $X4 $H $P &
# --- 2-D slab, vet_sc full
$R v_ref   ref 1 slab2d_def $S $V &
$R v_off   new 1 slab2d_def $S $V &
$R v_pipe  new 1 slab2d_def $S $V $P &
$R v_ref2  ref 2 slab2d_def $S $V $X2 &
$R v_off2  new 2 slab2d_def $S $V $X2 &
$R v_hm2   new 2 slab2d_def $S $V $X2 $H &
$R v_ph2   new 2 slab2d_def $S $V $X2 $H $P &
# --- 3-D box 84x32x32, 4 blocks, 60 cycles
$R e3_ref  ref 1 box3d_def $B3 &
$R e3_off  new 1 box3d_def $B3 &
$R e3_ctl  new 1 box3d_def $B3 $C &
$R e3_pipe new 1 box3d_def $B3 $P &
$R e3_ref2 ref 2 box3d_def $B3 &
$R e3_off2 new 2 box3d_def $B3 &
$R e3_hm2  new 2 box3d_def $B3 $H &
$R e3_ref4 ref 4 box3d_def $B3 &
$R e3_off4 new 4 box3d_def $B3 &
$R e3_hm4  new 4 box3d_def $B3 $H &
$R e3_ph4  new 4 box3d_def $B3 $H $P &
$R v3_ref4 ref 4 box3d_def $B3 $V &
$R v3_off4 new 4 box3d_def $B3 $V &
$R v3_hm4  new 4 box3d_def $B3 $V $H &
$R v3_ph4  new 4 box3d_def $B3 $V $H $P &
# --- restart (new paths on, 2 ranks): full run to t=40 with rst at t=20
$R rs_full new 2 slab2d_def time/tlim=40 output3/dt=20 output2/dt=40 output5/dt=40 $X2 $H $P &
wait
# restart from t=20
d=$W/cpu/rs_rst; rm -rf $d; mkdir -p $d; cd $d
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
rf=$(ls $W/cpu/rs_full/rst/*.00001.rst)
OMP_NUM_THREADS=1 mpirun -np 2 --oversubscribe --bind-to none $W/bin/athena_new_cpu -r $rf -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
echo CPU GATE DONE
