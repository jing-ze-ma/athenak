#!/bin/bash -l
# T-bit CPU gates of m1-vetcol (ref = m1-sp2 b5f3afea binaries, new = m1-vetcol)
W=/viper/ptmp2/jinma/vetcol_0924; R=$W/run.sh; D=/viper/ptmp2/jinma/defaults_0923/cpu
BX=box_convection
S="time/tlim=200 problem/vpert=1.0e-2"
X2="meshblock/nx2=16"
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
IMPM="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak time/cfl_number=1e6 rad_m1/implicit_cfl=10"
for c in ref new; do
  $R t_slab1_$c $c $BX 1 $D/slab2d_old.athinput $S &
  $R t_slab2_$c $c $BX 2 $D/slab2d_old.athinput $S $X2 &
  $R t_slabnd2_$c $c $BX 2 $D/slab2d_nd.athinput $S $X2 &
  $R t_vet2_$c $c $BX 2 $D/slab2d_old.athinput $S $V $X2 &
  $R t_vimp2_$c $c $BX 2 /viper/ptmp2/jinma/s0_0923/inp/slab2d_vimp.athinput time/tlim=50 problem/vpert=1.0e-2 $X2 &
  $R t_box4_$c $c $BX 4 $D/box3d_old.athinput $B3 &
  $R t_boxd4_$c $c $BX 4 $D/box3d_nd.athinput $B3 &
  $R t_expl2_$c $c $BX 2 $D/slab2d_old.athinput time/nlim=40 problem/vpert=1.0e-2 rad_m1/transport=explicit output1/dt=1.0e-6 $X2 &
  $R t_mkx1_$c $c none 1 $W/inp/marshak_x1.athinput $IMPM &
  $R t_mkcart_$c $c none 1 $W/inp/marshak_cart.athinput &
  $R t_spblast_$c $c sp_test 2 $W/inp/sp_blast.athinput time/nlim=100 &
done
wait
for t in slab1 slab2 slabnd2 vet2 vimp2 box4 boxd4 expl2 mkx1 mkcart spblast; do
  echo "== $t: $(python3 $W/cmp.py $W/cpu/t_${t}_ref $W/cpu/t_${t}_new | tail -n1) rc $(tail -n1 $W/cpu/t_${t}_ref/log.txt) $(tail -n1 $W/cpu/t_${t}_new/log.txt)"
done
echo TBIT DONE
