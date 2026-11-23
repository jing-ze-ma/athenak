#!/bin/bash
# gate_ia11.sh NODE: 1-GPU tiny BSG smoke gate for IA11 (orion 23i reference: pass 1 3.686e-11, converged after 3 passes, F/c 2.9006e+42); IA9 gate (gate2.sh, 4.020e-06 ref) unchanged
J=/resnick/groups/carnegie_poc/jingze; S=$J/orion_runs_1121; C=$S/ctl; G=$J/wt_bsgnote1101; FR=$J/fresh_prod_1120; F=$FR/fg_v2; OUT=$S/h200chk; R=$G/docs/handover/caltech_results/H200CHK
B=$J/wt_orion1121/bin/athena_ia11_bd9553537; ND=$1; HJ=7333072; TG=GATE_ia11; X="mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32"
env RT=0 NR=1 KEEPBIN=0 TL=15 HOLDER=$HJ NODE=$ND OUT=$OUT FRESH=$FR MEMWRAP=$C/memwrap.sh $F/armg_cal.sh $TG $B $F/decks/bsg_wngl_smoke.athinput 1 $X > $OUT/armg_$TG.out 2>&1
Lg=$OUT/$TG/a/run.log; p1=$(grep -m1 'clean_init steady pass 1:' $Lg | grep -oE '[0-9.]+e[-+][0-9]+$'); v=FAIL
[ "$p1" = 3.686e-11 ] && grep -q 'converged after 3 passes' $Lg && grep 'clean_init steady: converged' $Lg | grep -q 'F/c [^>]*-> 2.9006e+42' && v=PASS
{ echo "GATE $v: $(basename $B) md5 $(cut -c1-8 $B.md5) node $ND pass1=$p1 (ref 3.686e-11, 3 passes, F/c 2.9006e+42)"; grep -E 'clean_init steady (pass [0-3]:|:)' $Lg; } > $R/gate_ia11.txt; cat $R/gate_ia11.txt; grep -c FATAL $Lg
