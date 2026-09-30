#!/bin/bash
# He box accuracy gate (m1-positivity): C (off), P (off, 1-ulp cfl kick = noise), N (on),
# window T1..T2 (default 60000..74000, equal physical time).  Same tools as hebox_cfl2_0927.
W=/viper/ptmp2/jinma/m1pos_0930/hebox; cd $W; T1=${1:-60000}; T2=${2:-74000}
export PYTHONPATH=/viper/ptmp2/jinma/hebox_cfl2_0927/res:$PYTHONPATH
ANA_REF=C nice python3 /viper/ptmp2/jinma/hebox_cfl2_0927/scripts/ana.py runs $T1 $T2 C,P,N > ana_refC_${T1}_${T2}.txt 2>&1
nice python3 conv.py $T1 $T2 C P N > conv_${T1}_${T2}.txt 2>&1
python3 - <<PY > kehst_${T1}_${T2}.txt
import numpy as np
for a in ['C','P','N']:
    d=np.loadtxt('runs/%s/m1slab.hydro.hst'%a); t=d[:,0]; w=(t>=$T1)&(t<=$T2)
    k1=d[w,7]; kh=d[w,8]+d[w,9]
    print(a, 'KE1 mean %.4e std %.2e  KEh mean %.4e std %.2e  n %d' % (k1.mean(), k1.std(), kh.mean(), kh.std(), w.sum()))
PY
grep -h "m1-positivity\|Picard iter\|floor clips" runs/*/run.log > counters.txt
