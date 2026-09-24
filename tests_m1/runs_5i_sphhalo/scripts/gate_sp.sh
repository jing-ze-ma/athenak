#!/bin/bash
# sp gates (CPU): halo_mpi on/off, restart, T-sym / T-S1 / T-S4 ref vs new
W=/viper/ptmp2/jinma/sphhalo_0924; R=$W/scripts/run.sh; I=$W/inp
N=$W/bin/athena_new_none_cpu; F=$W/bin/athena_ref_none_cpu
# 1. on/off pairs, 4 ranks, 10 cycles
for c in vc edd m1; do for tl in L T; do
  $R w_${c}_on_$tl $N 4 $I/w_${c}_on_$tl.athinput &
  $R w_${c}_off_$tl $N 4 $I/w_${c}_off_$tl.athinput &
done; wait; done
# 2. restart
$R rsA $N 4 $I/rs_vc.athinput; wait
d=$W/rst0; rm -rf $d; cp -r $W/runs/rsA/rst $d
$R rsB $N 4 $d/hw.00001.rst
# 3. T-sym (4 ranks), T-S1 (1 rank), T-S4 1 rank and 4 ranks (8x8 angular)
$R tsym_ref $F 4 $I/sph_sym.athinput &
$R tsym_new $N 4 $I/sph_sym.athinput &
$R tsym_off $N 4 $I/sph_sym_off.athinput &
$R ts1_ref $F 1 $I/sph_diff.athinput mesh/nx1=64 meshblock/nx1=64 &
$R ts1_new $N 1 $I/sph_diff.athinput mesh/nx1=64 meshblock/nx1=64 &
$R ts1s_ref $F 1 $I/sph_diff_str.athinput mesh/nx1=64 meshblock/nx1=64 &
$R ts1s_new $N 1 $I/sph_diff_str.athinput mesh/nx1=64 meshblock/nx1=64 &
wait
Q="mesh/nx2=8 mesh/nx3=8"
for f in sph_atm_vc sph_atm; do
  $R ${f}_ref1 $F 1 $I/$f.athinput &
  $R ${f}_new1 $N 1 $I/$f.athinput &
  $R ${f}_ref4 $F 4 $I/$f.athinput $Q &
  $R ${f}_new4 $N 4 $I/$f.athinput $Q &
  $R ${f}_off4 $N 4 $I/${f}_off.athinput $Q &
  wait
done
echo GATESP_DONE
