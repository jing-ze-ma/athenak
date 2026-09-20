#!/usr/bin/env python3
"""Part 1: MLT shell-mean face table -> F_2s/F_req, F_dep/F_req, F_res/F_req, SUM.

Columns of mltfaces_*.txt.ts (1-based):
 1 i  2 r  3 T_f  4 p_f  5 grad  6 grad_ad  7 x  8 F_mlt  9 F_rad 10 F_used
11 F_req 12 F_raddiff 13 w_blend 14 F_rad_col0 15 F_conv_res 16 grad_rad
17 D 18 F_2s 19 F_deposited 20 F_cap
"""
import numpy as np, sys

F = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge9d/mltfaces_wedge9d.txt.ts'
TURN = 4705.0
R = 2.3717e11
TARG = [0.60, 0.65, 0.72, 0.80, 0.88, 0.93, 0.96]


def read_blocks(path):
    blocks = []
    t = None
    rows = []
    for line in open(path):
        if line.startswith('# t '):
            if t is not None:
                blocks.append((t, np.array(rows)))
            t = float(line.split()[2])
            rows = []
        elif line.strip() and not line.startswith('#'):
            rows.append([float(x) for x in line.split()])
    if t is not None:
        blocks.append((t, np.array(rows)))
    return blocks


blocks = read_blocks(F)
print('# blocks %d  t %.1f .. %.1f' % (len(blocks), blocks[0][0], blocks[-1][0]))

# pick ~every 0.1 turnover
want = [2.40 + 0.1 * k for k in range(12)]
sel = []
for w in want:
    tt = w * TURN
    j = min(range(len(blocks)), key=lambda k: abs(blocks[k][0] - tt))
    if j not in [s for s, _ in sel]:
        sel.append((j, w))

r0 = blocks[0][1][:, 1]
idx = [int(np.argmin(np.abs(r0 / R - x))) for x in TARG]
print('# faces: ' + '  '.join('r/R=%.3f(i=%d)' % (r0[k] / R, blocks[0][1][k, 0])
                              for k in idx))

hdr = '%6s' % 'turn'
for x in TARG:
    hdr += ' | %-27s' % ('r/R=%.2f  2s/dep/res/SUM' % x)
print(hdr)
for j, w in sel:
    t, a = blocks[j]
    line = '%6.3f' % (t / TURN)
    for k in idx:
        freq = a[k, 10]
        f2s = a[k, 17] / freq
        fdep = a[k, 18] / freq
        fres = a[k, 14] / freq
        line += ' | %6.3f %6.3f %6.3f %6.3f' % (f2s, fdep, fres, f2s + fdep + fres)
    print(line)

# also: shell-integrated total carried flux vs F_req at each face, full radial profile
print('\n# radial profile of SUM/F_req at selected times')
prof_t = [2.40, 2.80, 3.00, 3.20, 3.50]
cols = []
for w in prof_t:
    tt = w * TURN
    j = min(range(len(blocks)), key=lambda k: abs(blocks[k][0] - tt))
    cols.append((blocks[j][0] / TURN, blocks[j][1]))
print('%8s ' % 'r/R' + ' '.join('%8.3f' % c[0] for c in cols))
for k in range(0, len(r0), 6):
    print('%8.4f ' % (r0[k] / R) + ' '.join(
        '%8.3f' % ((c[1][k, 17] + c[1][k, 18] + c[1][k, 14]) / c[1][k, 10])
        for c in cols))
