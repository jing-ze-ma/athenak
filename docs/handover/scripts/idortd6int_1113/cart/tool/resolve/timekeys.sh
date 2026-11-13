#!/bin/bash
# idort-x1-1112: p arm = PROD deck + TIME keys (as idort_time_1111: blend_aniso 1.0 + KT on Cartesian, KTS on sp / cs)
# + X1 keys (keys_time.sh X1) + FC keys (11-12: FCA all decks, FCC sp/cs) + vet_gd_rotate_xfade=true where tool/xfade.py says the rotation is on the physical clock.
# Only decks with a <rad_m1> block. Log: timekeys.txt in the arm dir (geometry, xfade decision, every key applied).
f=$1; arm=${2:-p}; [ "$arm" = q ] && export CAND=fmm; T=/viper/ptmp2/jinma/d6int_1113/orion_task/tool; source $T/keys_time.sh
grep -q "^<rad_m1>" $f || { echo "timekeys: no rad_m1 block" >> timekeys.txt; exit 0; }
if grep -qiE "^ *(use_spherical_polar|use_cubed_sphere|cubed_sphere) *= *true" $f; then K="$KTS"; FK="$FCA $FCC"; g=curv; else K="$KT"; FK="$FCA"; g=cart; fi
python3 $T/setkeys_inplace.py $f rad_m1/implicit_blend_aniso=1.0 $K $X1 $FK
xf=$(python3 $T/xfade.py $f); A="rad_m1/implicit_blend_aniso=1.0 $K $X1 $FK"
case "$xf" in yes*) python3 $T/setkeys_inplace.py $f $XF; A="$A $XF";; esac
{ echo "timekeys: $g arm=$arm cand=$CAND"; echo "xfade: $xf"; echo "keys: $A"; } > timekeys.txt
