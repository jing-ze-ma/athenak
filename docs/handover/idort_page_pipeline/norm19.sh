#!/bin/bash
# normalised eval tables for the v10 plot scripts: PROD (eval/out_p) arm names as FIN3 page names; v10 FIN3 (idort_page cache/fin3ana) as the
# open-symbol control (HFfW / PRODd)
O=/raven/ptmp/jinma/idort1020/eval/out_p; D=/raven/ptmp/jinma/idort1020/page/cache/p19ana; V=/raven/ptmp/jinma/idort_page/cache/fin3ana; mkdir -p $D
for f in xord tord exact rswtq lr cons_lz picard; do
  { sed -E 's/(HF[A-Za-z0-9]*|HCF[A-Za-z0-9]*)\+[^ ]*/\1/g' $O/sp/$f.txt
    case $f in xord|tord) grep -E '^\S+\s+HFfZ\s' $V/$f.txt | sed -E 's/^(\S+\s+)HFfZ(\s)/\1HFfW\2/';; esac; } > $D/$f.txt
done
sed -E 's#FIN3_#FIN_#g; s#NP_#FIN_#g' $O/cart.txt > $D/cart_fin.txt
{ cat $O/lad.txt; sed -E 's/^(\S+\s+)FIN(s?)(\s)/\1PRODd\3/' <(grep -v '^==' $V/lad_FIN.txt | grep -E '^\S+\s+FINs?\s'); } > $D/lad_FIN.txt
ls -la $D
