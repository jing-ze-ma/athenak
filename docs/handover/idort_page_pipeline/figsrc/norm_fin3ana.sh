#!/bin/bash
# FIN3 page: normalised copies of IN's evaluation files (ana/out_281dd686) with the old arm names the page scripts parse:
#   sp arm HFfZ+<FIN3 keys> -> HFfZ (FIN3), F2O's HFfZ+<keys> -> HFfW (FIN2 control, open markers); Cartesian FIN3_x -> FIN_x
O=/raven/ptmp/jinma/fin3_1012/ana/out_281dd686; D=/raven/ptmp/jinma/idort_page/cache/fin3ana
for f in xord tord exact rswtq lr cons_lz picard; do F2=0; case $f in xord|tord) F2=1;; esac
  { sed -E 's/(HF[A-Za-z0-9]*|HCF[A-Za-z0-9]*)\+[^ ]*/\1/g' $O/FIN3/sp/$f.txt
    [ $F2 = 1 ] && [ -f $O/F2O/sp/$f.txt ] && sed -E 's/HFfZ\+[^ ]*/HFfW/g; s/(HF[A-Za-z0-9]*|HCF[A-Za-z0-9]*)\+[^ ]*/\1_F2O/g' $O/F2O/sp/$f.txt; } > $D/$f.txt
done
sed -E 's#FIN3_#FIN_#g' $O/FIN3/cart.txt > $D/cart_fin.txt
# Cartesian ladders: FIN3 (arm FIN / FINs) + FIN2 control F2O renamed to arm PRODd (open symbols; label rewritten to FIN2)
{ cat $O/FIN3/lad.txt; sed -E 's/^(\S+\s+)FIN(s?)(\s)/\1PRODd\3/' $O/F2O/lad.txt; } > $D/lad_FIN.txt
