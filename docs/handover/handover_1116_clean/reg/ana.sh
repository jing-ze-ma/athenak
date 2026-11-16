#!/bin/bash
# ana.sh ARM [ARM ...]: scores, iterations, CLEAN WARNING counts and conservation of each arm
: "${REG:?set REG}"; export ATHENA_VIS=${ATHENA_VIS:-$REG/cc/ana}; cd $REG
its() { grep -h "^CLEAN cycle" "$1" | grep -o "(sum [0-9]*)" | tr -d '()' | awk '{s+=$2}END{printf "%d", s}'; }
row() { printf "%-10s it %7s warn %2s cpu %8.2f | %s\n" "$2" "$(its $1)" "$(grep -c 'CLEAN WARNING' $1)" "$(grep -h 'cpu time used' $1 | awk '{print $NF}')" "$3"; }
for n in "$@"; do echo "== $n"
 row runs/thin/${n}_p_32/run.log thin_p "$(THPER=1 python3 gtg.py score thin ${n}_p 32 | tail -1 | grep -o 'L1(E) [^ ]* L1(F) [^ ]*')"
 row runs/thin/${n}_r_32/run.log thin_r "$(python3 gtg.py score thin ${n}_r 32 | tail -1 | grep -o 'L1(E) [^ ]* L1(F) [^ ]*')"
 row runs/thin/${n}_f_32/run.log thin_f "$(FULL=1 python3 gtg.py score thin ${n}_f 32 | tail -1 | grep -o 'L1(E) [^ ]* L1(F) [^ ]*')"
 row runs/dif/${n}_32/run.log dif "$(python3 gtg.py score dif $n 32 | tail -1 | grep -o 'L1(E) [^ ]*')"
 row runs/lat/${n}_32/run.log lat "$(python3 gtg.py score lat $n 32 | tail -1 | grep -o 'L1(E) [^ ]*')"
 row runs/pulse0/${n}_64/run.log pulse0 "energy $(python3 energy.py runs/pulse0/${n}_64 | awk '{print $5}' | sort -g | sed -n '$p')"
 row runs/thin/${n}_s88_32/run.log s88 "$(grep -h 'steady pass' runs/thin/${n}_s88_32/run.log | tail -1 | cut -c20-110)"
 row runs/thin/${n}_s81_32/run.log s81 "$(grep -h 'steady pass' runs/thin/${n}_s81_32/run.log | tail -1 | cut -c20-110)"
 row rhse/$n/log.txt rhse "$(grep -h 'clean_init steady:' rhse/$n/log.txt | cut -c1-70) vcs_max $(awk 'NR>2{print $5}' rhse/$n/rh.user.hst | sort -g | tail -1)"
 row cs/uni_$n/run.log cs_uni "$(python3 $REG/cs_tests/ct.py uni cs/uni_$n | tail -1)"
 row cs/closed_$n/run.log cs_closed "$(python3 $REG/cs_tests/ct.py cons cs/closed_$n | tail -1 | cut -c1-110)"
 row cs/thin1_$n/run.log cs_thin1 "$(python3 $REG/cs_tests/ct.py score thin cs/thin1_$n | tail -1 | cut -c1-80)"
 row cs/feed1_$n/run.log cs_feed1 "$(python3 -c "
import sys, os; sys.path.insert(0,'$REG/cs_tests'); sys.path.insert(0,os.environ['ATHENA_VIS']); import ct; ct.budget('cs/feed1_$n', 10.0)" 2>&1 | head -2 | tr '\n' ' ' | cut -c1-330)"
 for q in q Q; do row fb/runsf/flz/${q}${n}_32/run.log flz_$q "$(cd fb && python3 fb.py hst flz ${q}${n} 32 2>&1 | tail -1 | cut -c1-120)"; done
  row cs/thin1mb_$n/run.log cs_thin1mb "$(python3 $REG/cs_tests/ct.py score thin cs/thin1mb_$n | tail -1 | cut -c1-80)"
  (cd cc && python3 score7.py $n | sed "s/^/  cc /")
done
