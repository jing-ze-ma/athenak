#!/bin/bash -l
# the S1/S2/S5 spherical-polar gates (+ vet_col Cartesian), ref (9c1a12d7) vs new: bitwise
W=/viper/ptmp2/jinma/vetcol2_0924; R=$W/run.sh
for c in ref new; do
  $R e_sd_$c $c none 1 $W/inp/sph_diff.athinput &
  $R e_sds_$c $c none 1 $W/inp/sph_diff_str.athinput &
  $R e_sym_$c $c none 4 $W/inp/sph_sym.athinput time/nlim=200 &
  $R e_symg_$c $c none 4 $W/inp/sph_sym_gas.athinput time/nlim=200 &
  $R e_mk_$c $c none 1 $W/inp/marshak_sph.athinput &
  $R m_fs_$c $c none 1 $W/inp/sph_fs.athinput &
  $R m_atm_$c $c none 1 $W/inp/sph_atm.athinput mesh/nx1=32 meshblock/nx1=32 &
  $R m_rw_$c $c none 1 $W/inp/rw_sph_R100.athinput &
  $R m_sym_$c $c none 4 $W/inp/sph_sym.athinput time/nlim=200 rad_m1/closure=m1 rad_m1/implicit_precond=line &
  $R v_atm_$c $c none 1 $W/inp/sph_atm_vc.athinput mesh/nx1=64 meshblock/nx1=64 &
  $R v_sym_$c $c none 4 $W/inp/sph_sym.athinput time/nlim=200 rad_m1/closure=vet_col &
  $R v_symg_$c $c none 4 $W/inp/sph_sym_gas.athinput time/nlim=200 rad_m1/closure=vet_col rad_m1/c_light=100.0 &
  $R v_str_$c $c none 1 $W/inp/sph_atm_str_vc.athinput mesh/nx1=32 meshblock/nx1=32 time/nlim=200 &
  $R v_pp_$c $c none 1 $W/inp2/milne_vc.athinput mesh/nx1=64 meshblock/nx1=64 time/nlim=200 &
  $R v_bin_$c $c none 1 $W/inp/sph_atm_bin_vc.athinput time/nlim=100 &
done
wait
for t in e_sd e_sds e_sym e_symg e_mk m_fs m_atm m_rw m_sym v_atm v_sym v_symg v_str v_pp v_bin; do
  echo "== $t: $(python3 $W/cmp.py $W/cpu/${t}_ref $W/cpu/${t}_new | tail -n1) rc $(tail -n1 $W/cpu/${t}_ref/log.txt) $(tail -n1 $W/cpu/${t}_new/log.txt)"
done
echo SP GATES DONE
