#!/usr/bin/env python3
"""tests_r19 part 3: measured F_2s/F_req (run's own face table) vs the porosity
prediction <F_col>/F_req and the shell-mean-state flux F_mean/F_req."""
import numpy as np
import poros as P

TS = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge9d/mltfaces_wedge9d.txt.ts'
TURN = 4705.0


def read_blocks(path):
    blocks, t, rows = [], None, []
    for line in open(path):
        if line.startswith('# t '):
            if t is not None:
                blocks.append((t, np.array(rows)))
            t = float(line.split()[2])
            rows = []
        elif line.strip() and not line.startswith('#'):
            rows.append([float(x) for x in line.split()])
    blocks.append((t, np.array(rows)))
    return blocks


B = read_blocks(TS)
ids = [5, 14, 16, 20]
res = []
r0 = None
for i in ids:
    A = P.analyse(i, None, [])
    res.append(A)
    if r0 is None:
        r0 = A['rc']/P.RS
rows = [k for k in range(0, len(r0), 4) if 0.55 <= r0[k] <= 1.001]

print('%4s %7s | ' % ('i', 'r/R') + ' | '.join(
    '%-26s' % ('t=%.2f  2s / <Fcol> / Fmean' % (A['t']/TURN)) for A in res))
for k in rows:
    ln = '%4d %7.4f | ' % (k, r0[k])
    for A in res:
        j = min(range(len(B)), key=lambda q: abs(B[q][0]-A['t']))
        a = B[j][1]
        kf = int(np.argmin(np.abs(a[:, 1] - A['rc'][k])))
        f2s = a[kf, 17]/a[kf, 10]
        ln += '%8.3f %8.3f %8.3f | ' % (f2s, A['Fcm'][k]/A['Freq'][k],
                                        A['Fmean'][k]/A['Freq'][k])
    print(ln)

print('\nresidual  (F_2s/Freq) - (<Fcol>/Freq)   and   (F_2s/Freq) - (Fmean/Freq)')
print('%4s %7s | ' % ('i', 'r/R') + ' | '.join('%-17s' % ('t=%.2f' % (A['t']/TURN))
                                               for A in res))
for k in rows:
    ln = '%4d %7.4f | ' % (k, r0[k])
    for A in res:
        j = min(range(len(B)), key=lambda q: abs(B[q][0]-A['t']))
        a = B[j][1]
        kf = int(np.argmin(np.abs(a[:, 1] - A['rc'][k])))
        f2s = a[kf, 17]/a[kf, 10]
        ln += '%8.3f %8.3f | ' % (f2s - A['Fcm'][k]/A['Freq'][k],
                                  f2s - A['Fmean'][k]/A['Freq'][k])
    print(ln)
