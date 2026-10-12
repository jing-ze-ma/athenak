#!/bin/bash
# analyse.sh <runs/JOBID>: report numbers for TASK-2026-10-11-deltaai-bsgrec-smoke
S=$1; G=/work/nvme/bivj/jma20/bsgrec_1011; PY=/u/jma20/ATHENAK/venv/bin/python
grep -E '^(RUN|END)' $G/logs/*.$(basename $S).out
for n in R4 RC4 R2; do
  D=$S/$n; [ -f $D/out.txt ] || continue
  echo; echo "=== $n"
  echo "FATAL $(grep -c FATAL $D/out.txt)  NaN-lines $(grep -c -i -w nan $D/out.txt)  NON-CONVERGED: $(grep -o 'NON-CONVERGED=[^ ]*' $D/out.txt | tail -1)"
  grep -h 'implicit transport: solves' $D/out.txt | tail -1
  grep -h -i -E 'NOT ADMISSIBLE|fallback|time2_pred_cap|limited|implicit_gas_src|gas_src' $D/out.txt | sort | uniq -c | head -12
  $PY -I - $D <<'PY'
import sys, re
o = open(sys.argv[1] + '/out.txt', errors='replace').read()
el = [(int(m.group(2)), float(m.group(1)), float(m.group(3))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+) time=\S+ dt=(\S+)', o)]
if not el: print('no cycles'); sys.exit()
d = {c1: t1 - t0 for (c0, t0, _), (c1, t1, _) in zip(el[:-1], el[1:]) if c1 == c0 + 1}
last = el[-1][0]; rot = 10
st = [d[c] for c in d if 3 <= c and c % rot and c != last]; ro = [(c, round(d[c], 3)) for c in d if c % rot == 0]
print(f'setup (elapsed at cycle 0) {el[0][1]:.1f} s; cycles {last}; dt {el[1][2] if len(el)>1 else el[0][2]:.6e} s (cycle 0 dt {el[0][2]:.6e})')
print('per-cycle wall s: ' + ' '.join(f'{c}:{d[c]:.3f}' for c in sorted(d)))
if st: print(f'steady s/cycle (cycles 3..{last} excl. rotation and last): mean {sum(st)/len(st):.3f} median {sorted(st)[len(st)//2]:.3f} (n {len(st)}); rotation cycles {ro}')
PY
done
for n in R4 RC4; do
  D=$S/$n; [ -f $D/agcar3d.user.hst ] || [ -f $D/bsg3d.user.hst ] || continue
  echo; echo "=== $n user.hst: t, L_top, L_in, Etot, Picard"
  $PY -I - $D/bsg3d.user.hst <<'PY'
import sys, re, numpy as np
f = sys.argv[1]; hdr = [l for l in open(f) if l.startswith('#')][-1]
c = {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', hdr)}
h = np.loadtxt(f, ndmin=2); ks = ['time', 'L_top', 'L_in', 'Etot', 'Picard']
print('  ' + '  '.join(f'{k:>15s}' for k in ks)); [print('  ' + '  '.join(f'{r[c[k]]:15.8e}' for k in ks)) for r in h]
PY
  echo "--- ic0check $n"
  PYTHONPATH=$G/bundle/docs/handover/bsg-hrdet-1009/page $PY $G/bundle/docs/handover/bsgrec-1011/ic0check.py check $D $G/bundle/docs/handover/bsgrec-1011/ic0_ref_1009.json 2>&1 | tail -12
done
