"""raycheap15: WHERE the self-convergence difference of m1_e sits on the cubed sphere.
python3 conv_loc.py [--time T] RUN_c0.15 RUN_c0.075 ... (coarse first, finest last)
Prints, for every successive pair, d = |E_k - E_k+1|_1/|E_finest|_1 split by region
(seam band, cube-vertex patch, inner / outer radial shells, the panel-centre patch of the
clump shell r in [2.6, 3.4], the rest), the region p, the 5 MeshBlocks with the largest
share, and (3+ runs) a Richardson estimate: p_R from the 3 finest runs, the extrapolated
E_R = E_f + (E_f - E_f-1)/(2^p_R - 1) and the error |E_k - E_R|_1/|E_R|_1 of every run.
Positions come from the LOGICAL mb_geometry (x2, x3 in [-1, 1] per panel; x1 uniform),
which is what a v1.2 cs file keeps there.  Needs vis/python/bin_convert.py."""
import sys, glob, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'vis', 'python'))
import bin_convert as bc

TSEL = None
if len(sys.argv) > 2 and sys.argv[1] == '--time':
    TSEL = float(sys.argv[2]); sys.argv = [sys.argv[0]] + sys.argv[3:]
NB = 2   # band width in cells


def load(run):
    fs = sorted(glob.glob(f'{run}/bin/*.m1.*.bin'))
    f = fs[-1] if TSEL is None else min(fs, key=lambda x: abs(bc.read_binary(x)['time'] - TSEL))
    d = bc.read_binary(f)
    g = np.asarray(d['mb_geometry'])
    o = sorted(range(len(g)), key=lambda m: tuple(np.round(d['mb_logical'][m], 0)) + tuple(np.round(g[m], 12)))
    E, reg, mbi = [], [], []
    x1min, x1max = d['x1min'], d['x1max']
    for m in o:
        e = np.asarray(d['mb_data']['m1_e'][m])          # (nk, nj, ni)
        nk, nj, ni = e.shape
        x1 = g[m][0] + g[m][3]*np.arange(ni)
        x2 = g[m][1] + g[m][4]*np.arange(nj)
        x3 = g[m][2] + g[m][5]*np.arange(nk)
        X3, X2, X1 = np.meshgrid(x3, x2, x1, indexing='ij')
        b2 = (1.0 - np.abs(X2)) < NB*g[m][4]
        b3 = (1.0 - np.abs(X3)) < NB*g[m][5]
        r = np.zeros(e.shape, dtype=int)                  # 0 rest
        r[(X1 - x1min) < NB*g[m][3]] = 3                  # inner shell
        r[(x1max - X1) < NB*g[m][3]] = 4                  # outer shell
        r[(X1 > 2.6) & (X1 < 3.4) & (np.abs(X2) < 0.35) & (np.abs(X3) < 0.35)] = 5
        r[b2 ^ b3] = 1                                    # seam band
        r[b2 & b3] = 2                                    # vertex patch
        E.append(e.ravel()); reg.append(r.ravel()); mbi.append(np.full(e.size, m))
    g2 = [g[m][1] for m in o] + [g[m][1] + g[m][4]*(d['mb_data']['m1_e'][m].shape[1] - 1) for m in o]
    g3 = [g[m][2] for m in o] + [g[m][2] + g[m][5]*(d['mb_data']['m1_e'][m].shape[0] - 1) for m in o]
    print(run, 'n_mb', len(o), 'logical x2 centres [%.4f, %.4f] x3 [%.4f, %.4f]'
          % (min(g2), max(g2), min(g3), max(g3)), '(the region split assumes a panel spans [-1, 1])')
    return d['time'], np.concatenate(E), np.concatenate(reg), np.concatenate(mbi)


L = [load(r) for r in sys.argv[1:]]
print('t_end', [round(u[0], 4) for u in L])
reg, mbi = L[-1][2], L[-1][3]
nrm = np.abs(L[-1][1]).sum()
names = ['rest', 'seam', 'vertex', 'inner', 'outer', 'clumppatch']
print('cells per region', {names[q]: int((reg == q).sum()) for q in range(6)})
D = []
for k in range(len(L) - 1):
    a = np.abs(L[k][1] - L[k+1][1])
    D.append([a[reg == q].sum()/nrm for q in range(6)] + [a.sum()/nrm])
    top = sorted(set(mbi), key=lambda m: -a[mbi == m].sum())[:5]
    print(f'pair {k}: d', '%.3e' % D[-1][-1], ' '.join(f'{names[q]}=%.2e' % D[-1][q] for q in range(6)),
          '| top MBs', [(int(m), '%.2f' % (a[mbi == m].sum()/a.sum())) for m in top])
for k in range(len(D) - 1):
    print(f'p pair {k}/{k+1}:', ' '.join(
        f'{(names + ["all"])[q]}=%.2f' % (math.log2(D[k][q]/D[k+1][q]) if D[k+1][q] > 0 and D[k][q] > 0
                                           else float('nan')) for q in range(7)))
if len(L) >= 3:
    d1 = np.abs(L[-3][1] - L[-2][1]).sum(); d2 = np.abs(L[-2][1] - L[-1][1]).sum()
    pr = math.log2(d1/d2)
    pe = max(pr, 0.5)
    ER = L[-1][1] + (L[-1][1] - L[-2][1])/(2.0**pe - 1.0)
    nR = np.abs(ER).sum()
    print('Richardson p_R (3 finest) %.2f (used %.2f); err vs E_R:' % (pr, pe),
          ['%.3e' % (np.abs(u[1] - ER).sum()/nR) for u in L])
