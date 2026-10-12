#!/bin/bash
# analyse.sh <runs/JOBID>: TASK-2026-10-12-deltaai-bsgaccel-fullbsg report numbers
S=$1; G=/work/nvme/bivj/jma20/bsgaccel_1012; PY=/u/jma20/ATHENAK/venv/bin/python
grep -E '^(RUN|END|SKIP|C4|N|B)' $G/logs/*.$(basename $S).out | grep -v '^Note'
for n in C4 C4g N4 N4n N2 B4 B2; do
  D=$S/$n; [ -f $D/out.txt ] || continue
  echo; echo "=== $n"
  echo "FATAL $(grep -c FATAL $D/out.txt)  NaN-lines $(grep -c -i -w nan $D/out.txt)  $(grep -o 'NON-CONVERGED=[^ ]*' $D/out.txt | tail -1)"
  grep -h -E 'implicit transport: solves|inner iterations mean|vet_gd_halo_ipc: ON|time_scheme=hesdirk2: stage steps|positivity fallbacks|gas coupling: newton|out of memory|cudaErrorMemory|bad_alloc' $D/out.txt | sed 's/^/  /' | cut -c1-330
  echo "  NOT ADMISSIBLE lines: $(grep -c 'NOT ADMISSIBLE' $D/out.txt) $(grep -o 'NOT ADMISSIBLE at cycle [0-9]*' $D/out.txt | grep -o '[0-9]*$' | tr '\n' ',')"
  $PY -I - $D <<'PY'
import sys, re
o = open(sys.argv[1] + '/out.txt', errors='replace').read()
el = [(int(m.group(2)), float(m.group(1)), float(m.group(3))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+) time=\S+ dt=(\S+)', o)]
if len(el) < 2: print('  cycles:', len(el)); sys.exit()
d = {c1: t1 - t0 for (c0, t0, _), (c1, t1, _) in zip(el[:-1], el[1:]) if c1 == c0 + 1}
st = [d[c] for c in d if 3 <= c <= 22 and c % 10]; al = list(d.values())
print(f'  setup {el[0][1]:.1f} s; last cycle {el[-1][0]}; steady (3..22 not x10) mean {sum(st)/max(len(st),1):.3f} median {sorted(st)[len(st)//2] if st else float("nan"):.3f} (n {len(st)}); '
      f'rotation c10 {d.get(10, float("nan")):.3f} c20 {d.get(20, float("nan")):.3f}; average {sum(al)/len(al):.3f} s/cycle')
print('  dt per cycle: ' + ' '.join(f'{c}:{dt:.4g}' for c, _, dt in el))
PY
  pk=$(awk '{for(i=2;i<=NF;i+=2){g=$i; sub(",","",g); v=$(i+1); if(v+0>m[g]) m[g]=v+0}} END{for(g in m) printf "GPU%s %d  ", g, m[g]}' $D/gpumem.log 2>/dev/null); echo "  peak MiB: $pk"
done
$PY -I - $S <<'PY'
import sys, re, os, numpy as np
S = sys.argv[1]
def hst(n):
    f = f'{S}/{n}/bsg3d.user.hst'
    if not os.path.exists(f): return None
    hdr = [l for l in open(f) if l.startswith('#')][-1]
    return np.loadtxt(f, ndmin=2), {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', hdr)}
b, n = hst('B4'), hst('N4')
if b and n:
    m = min(len(b[0]), len(n[0])); print(f'\n=== B4 vs N4 user.hst ({m} rows): max relative difference over rows')
    for k in ('L_int', 'L_top', 'E_rad', 'Etot'):
        x, y = b[0][:m, b[1][k]], n[0][:m, n[1][k]]; den = np.maximum(np.abs(x), 1e-300)
        print(f'  {k:6s} {np.max(np.abs(x - y) / den):.3e}  (last row: B4 {x[-1]:.10e} N4 {y[-1]:.10e})')
for nm in ('C4', 'C4g'):
    h = hst(nm)
    if h is None: continue
    print(f'\n=== {nm} user.hst: t, dt, L_top, L_in, Etot, Picard')
    for r in h[0]: print('  ' + '  '.join(f'{r[h[1][k]]:.8e}' for k in ('time', 'dt', 'L_top', 'L_in', 'Etot', 'Picard')))
PY
