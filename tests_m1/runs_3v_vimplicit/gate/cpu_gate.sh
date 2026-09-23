#!/bin/bash
# CPU gates of implicit_vimp (runs_3v_vimplicit): switch-off bitwise vs 89dc28d7, the 2-D
# He slab (200 s) and 3-D box (60 cycles) with vimp, 1 block / 2 blocks / 2 ranks, and
# the restart (slab to 100 s, restart to 200 s, against the continuous run).
W=/viper/ptmp2/jinma/vimp_3v/gate; R=$W/cpu_run.sh
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=200 output3/dt=100"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60"
$R base e_def_base 1 slab2d_def $S &
$R new  e_def_new  1 slab2d_def $S &
$R base e_plm_base 1 slab2d_plm $S &
$R new  e_plm_new  1 slab2d_plm $S &
$R new  e_vim  1 slab2d_plm_vimp $S &
$R new  e_vim_mb 1 slab2d_plm_vimp $S meshblock/nx2=16 &
$R new  e_vim_m2 2 slab2d_plm_vimp $S meshblock/nx2=16 &
$R new  e_vim_nost 1 slab2d_plm_vimp_nost $S &
$R new  v_plm  1 slab2d_plm $S $V &
$R new  v_vim  1 slab2d_plm_vimp $S $V &
wait
$R new  e3_plm 1 box3d_plm $B3 &
$R new  e3_vim 1 box3d_plm_vimp $B3 &
$R new  e3_vim_m2 2 box3d_plm_vimp $B3 &
wait
echo CPU GATE DONE
