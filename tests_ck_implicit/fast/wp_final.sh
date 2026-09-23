#!/bin/bash
# final binary (fast6 = 2b3e6917): off gate, and the final combination + jreuse alone
F=/viper/ptmp2/jinma/wt_ckfast/tests_ck_implicit/fast
P=/viper/ptmp2/jinma/ckfast_0923
B=$P/athena.cpu.fast6
bash $F/wp_run.sh $B $P/wp9 A2 t4 semi n
for t in A2 D A1; do
  bash $F/wp_run.sh $B $P/wp9 $t c2 c4 c8 j p5
done
echo CHAIN7_DONE
