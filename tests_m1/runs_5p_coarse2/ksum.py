"""ksum.py <prof dir> [c0=10] [c1=40]: kernel time per cycle (ms) and launches per cycle
of the mg / preconditioner / Krylov kernels between the NewTimeStep kernels of cycles
c0 and c1, plus the GPU idle time before each group (gap to the previous kernel end)."""
import collections
import csv
import glob
import re
import sys

d = sys.argv[1]
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 10
c1 = int(sys.argv[3]) if len(sys.argv) > 3 else 40
kt = glob.glob(d + '/**/*kernel_trace.csv', recursive=True)[0]
rows = []
for r in csv.DictReader(open(kt)):
    rows.append((int(r['Start_Timestamp']), int(r['End_Timestamp']), r['Kernel_Name']))
rows.sort()


def name(n):
    for key in ('M1PCRMG', 'M1PCRX'):
        if key in n:
            return key + ('<float>' if '<float>' in n else '')
    m = re.search(r'RadiationM1::(\w+)\(', n)
    s = m.group(1) if m else n[:40]
    lam = re.search(r'\{lambda\([^)]*\)#(\d+)\}', n)
    tag = ('#' + lam.group(1)) if lam else ''
    return s + tag + ('R' if 'ParallelReduce' in n else '')


nts = [i for i, r in enumerate(rows) if 'NewTimeStep' in r[2]]
per = max(1, round(len(nts)/41))
ends = nts[per - 1::per]
i0, i1 = ends[c0], ends[c1]
tk = collections.defaultdict(float)
ti = collections.defaultdict(float)
nk = collections.Counter()
prev = rows[i0][1]
for s, e, n in rows[i0 + 1:i1 + 1]:
    k = name(n)
    tk[k] += (e - s)*1e-6
    ti[k] += max(0, s - prev)*1e-6
    nk[k] += 1
    prev = max(prev, e)
nc = c1 - c0
tot = (rows[i1][1] - rows[i0][1])*1e-6/nc
print(f'window cycles {c0}-{c1}: {tot:.2f} ms/cycle (GPU timeline)')
print(f'{"kernel":40s} {"ms/cyc":>7s} {"idle":>6s} {"n/cyc":>6s} {"us/launch":>9s}')
for k in sorted(tk, key=lambda x: -(tk[x] + ti[x]))[:30]:
    print(f'{k:40s} {tk[k]/nc:7.3f} {ti[k]/nc:6.3f} {nk[k]/nc:6.1f} '
          f'{tk[k]/nk[k]*1e3:9.1f}')
g = collections.defaultdict(lambda: [0.0, 0.0, 0])
for k in tk:
    grp = ('mg' if ('MG' in k or 'M1PCRMG' in k) else
           'fine PCR' if 'M1PCRX' in k else 'other')
    g[grp][0] += tk[k]/nc
    g[grp][1] += ti[k]/nc
    g[grp][2] += nk[k]/nc
for k, v in g.items():
    print(f'GROUP {k:10s} kernel {v[0]:7.3f} idle {v[1]:6.3f} launches {v[2]:6.1f}')
