#!/bin/bash -l
# the S1 spherical-polar Eddington gates, ref (4e428a30) vs new (m1-sp2): bitwise
W=/viper/ptmp2/jinma/s2_0924; R=$W/run.sh
for c in ref new; do
  $R e_sd_$c $c none 1 $W/inp/sph_diff.athinput &
  $R e_sds_$c $c none 1 $W/inp/sph_diff_str.athinput &
  $R e_sym_$c $c none 4 $W/inp/sph_sym.athinput time/nlim=200 &
  $R e_symg_$c $c none 4 $W/inp/sph_sym_gas.athinput time/nlim=200 &
  $R e_mk_$c $c none 1 $W/inp/marshak_sph.athinput &
done
wait
for t in sd sds sym symg mk; do
  echo "== $t: $(python3 $W/cmp.py $W/cpu/e_${t}_ref $W/cpu/e_${t}_new | tail -n1) rc $(tail -n1 $W/cpu/e_${t}_ref/log.txt) $(tail -n1 $W/cpu/e_${t}_new/log.txt)"
done
echo S1SP DONE
