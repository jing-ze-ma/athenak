#!/bin/bash
# runs_3z_validate: all tables.  bash analyse_all.sh > RESULTS_all.txt
V=/viper/ptmp2/jinma/validate_0923; S=$V/scripts; T=/viper/u2/jinma/ATHENAK/athenak/tests_m1
echo "### T7 radiative shocks (idx 2 = t 1e-10 for M2/M5, 1e-9 for M7; last = 2e-10 / 4e-9)"
python3 $S/t7_table.py $V/t7 2 last
echo "### packet (beam width)"
python3 $S/packet_anal.py $V/gpu2d/pk_*
echo "### shadow"
for d in $V/gpu2d/sh_*; do
  n=$(basename $d); a=${n/sh_/arr_}; arr=""
  [ -d $V/gpu2d/$a ] && arr="--arrival-dumps $V/gpu2d/$a/bin/*.bin"
  echo "== $n $(grep -o 'NON-CONVERGED=[0-9.e+-]*' $d/run.log | tail -1) $(grep -o 'stage fallbacks=[0-9.e+-]*' $d/run.log | tail -1)"
  python3 $T/t2_shadow.py $d/bin/*.bin $arr 2>&1 | grep -i "arrival\|depth\|speed\|PASS\|FAIL" | head -6
done
echo "### He 3-D vet_sc (gate3d.py; first dir = reference)"
python3 $T/runs_3o_vet3d_cfl/gate3d.py $V/he3d/h2_c0.15 $V/he3d/h2_c0.3 $V/he3d/be_c0.15 $V/he3d/be_c0.3
python3 $S/he_pairs.py
