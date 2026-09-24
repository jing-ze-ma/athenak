#!/bin/bash -l
# the S1 (Eddington) and S2 (M1 closure) spherical-polar gates, ref (m1-sp2) vs new: bitwise
W=/viper/ptmp2/jinma/vetcol_0924; R=$W/run.sh
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
done
wait
for t in e_sd e_sds e_sym e_symg e_mk m_fs m_atm m_rw m_sym; do
  echo "== $t: $(python3 $W/cmp.py $W/cpu/${t}_ref $W/cpu/${t}_new | tail -n1) rc $(tail -n1 $W/cpu/${t}_ref/log.txt) $(tail -n1 $W/cpu/${t}_new/log.txt)"
done
echo SP GATES DONE
