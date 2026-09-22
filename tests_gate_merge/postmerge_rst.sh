#!/bin/bash -l
# restart gate of the merged tree: box mode 3 with the production warm start, 50 cycles
# straight vs 25 + restart + 25 (history every cycle); uses build_cpu_box as built
module purge; module load gcc/14 cmake/4.0
A=/viper/u2/jinma/ATHENAK/athenak; G=$A/tests_gate_merge; X=$A/build_cpu_box/src/athena
IN=/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/he_box_w8.athinput
O="mesh/nx2=16 mesh/nx3=16 meshblock/nx1=134 meshblock/nx2=8 meshblock/nx3=8"
H="output1/dt=1e-9 output2/dt=1e-9"
D=$G/rst_gate; rm -rf $D; mkdir -p $D/s $D/c; 
cd $D/s; $X -i $IN $O time/nlim=50 > run.log 2>&1
cd $D/c; $X -i $IN $O time/nlim=25 output5/dt=1e9 > run1.log 2>&1
R=$(ls rst/*.rst | tail -1); echo "restart from $R"
$X -r $R time/nlim=50 > run2.log 2>&1; echo "restart exit=$?"
grep -i -E 'warn|RTWARM|warm' run2.log | head -5
python3 - <<'PY'
import numpy as np
for f in ['feczrt.hydro.hst','feczrt.user.hst']:
    a=[l for l in open('../s/'+f) if not l.startswith('#')]
    b=[l for l in open(f) if not l.startswith('#')]
    sb=set(b); miss=[i for i,l in enumerate(a) if l not in sb]
    print(f,'straight rows',len(a),'chain rows',len(b),'straight rows missing bit-for-bit in chain:',len(miss), miss[:5])
PY
for f in rt_surface.bin rt_profile.bin; do cmp -s $f ../s/$f && echo "$f IDENTICAL" || echo "$f DIFFER"; done
rm -rf $D/s/bin $D/s/cbin* $D/s/rst $D/c/bin $D/c/cbin* $D/c/rst
echo RST_GATE_DONE
