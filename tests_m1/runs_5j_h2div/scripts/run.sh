#!/bin/bash -l
# usage: run.sh <bin-tag> <name> [block/key=val ...]   CPU, 1 rank, correctness only; n=64 shadow to t=20
# overrides are written into a copy of the input (so keys absent from the input work too)
b=$1; name=$2; shift 2
module purge >/dev/null 2>&1; module load gcc/14 openmpi/5.0 >/dev/null 2>&1
W=/viper/ptmp2/jinma/h2div_0924
INP=${INP:-$W/inp/m1_shadow.athinput}
d=$W/runs/$name; rm -rf $d; mkdir -p $d; cd $d
cp $INP in.athinput
for o in mesh/nx1=64 meshblock/nx1=64 mesh/nx2=40 meshblock/nx2=40 time/tlim=20 output1/dt=0.5 "$@"; do
  blk=${o%%/*}; kv=${o#*/}; k=${kv%%=*}; v=${kv#*=}
  python3 - in.athinput "$blk" "$k" "$v" <<'PY'
import sys, re
f, blk, k, v = sys.argv[1:]
L = open(f).read().split('\n'); cur = None; done = False
for n, l in enumerate(L):
    m = re.match(r'<(\S+)>', l)
    if m: cur = m.group(1); continue
    if cur == blk and re.match(r'\s*' + re.escape(k) + r'\s*=', l):
        L[n] = '%s = %s' % (k, v); done = True
if not done:
    n = [i for i, l in enumerate(L) if l.strip() == '<%s>' % blk][0]
    L.insert(n + 1, '%s = %s' % (k, v))
open(f, 'w').write('\n'.join(L))
PY
done
OMP_NUM_THREADS=1 nice $W/bin/athena_${b}_cpu -i in.athinput -d $d > log.txt 2>&1
echo "rc=$?" >> log.txt
