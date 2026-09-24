#!/bin/bash
# Cartesian T-bit CPU set (runs_5a_sp_s1 table minus the sp hydro case), ref vs new
W=/viper/ptmp2/jinma/sphhalo_0924; R=$W/scripts/run.sh; D=/viper/ptmp2/jinma/defaults_0923/cpu
S1=/viper/ptmp2/jinma/s1_0924/inp
S="time/tlim=200 problem/vpert=1.0e-2"
X2="meshblock/nx2=16"
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
B3="mesh/nx2=32 mesh/nx3=32 mesh/x2max=3.39480832e8 mesh/x3max=3.39480832e8 meshblock/nx2=16 meshblock/nx3=16 time/nlim=60 problem/vpert=1.0e-2"
IMPM="rad_m1/transport=implicit_x1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak time/cfl_number=1e6 rad_m1/implicit_cfl=10"
for c in ref new; do
  BX=$W/bin/athena_${c}_box_convection_cpu; BN=$W/bin/athena_${c}_none_cpu
  $R c_slab1_$c $BX 1 $D/slab2d_old.athinput $S &
  $R c_slab2_$c $BX 2 $D/slab2d_old.athinput $S $X2 &
  $R c_slabnd2_$c $BX 2 $D/slab2d_nd.athinput $S $X2 &
  $R c_vet2_$c $BX 2 $D/slab2d_old.athinput $S $V $X2 &
  $R c_vimp2_$c $BX 2 /viper/ptmp2/jinma/s0_0923/inp/slab2d_vimp.athinput time/tlim=50 problem/vpert=1.0e-2 $X2 &
  $R c_box4_$c $BX 4 $D/box3d_old.athinput $B3 &
  $R c_boxd4_$c $BX 4 $D/box3d_nd.athinput $B3 &
  $R c_expl2_$c $BX 2 $D/slab2d_old.athinput time/nlim=40 problem/vpert=1.0e-2 rad_m1/transport=explicit output1/dt=1.0e-6 $X2 &
  $R c_mkx1_$c $BN 1 $S1/marshak_x1.athinput $IMPM &
  $R c_mkcart_$c $BN 1 $S1/marshak_cart.athinput &
done
wait
for t in slab1 slab2 slabnd2 vet2 vimp2 box4 boxd4 expl2 mkx1 mkcart; do
  echo "== $t: $(python3 $W/scripts/cmpdir.py $W/runs/c_${t}_ref $W/runs/c_${t}_new) rc $(tail -n1 $W/runs/c_${t}_ref/log.txt) $(tail -n1 $W/runs/c_${t}_new/log.txt)"
done
echo CART_DONE
