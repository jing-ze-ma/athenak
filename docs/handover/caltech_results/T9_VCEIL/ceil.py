#!/usr/bin/env python3
"""ceil.py OUT.txt DUMP   (orion zzzzc4 items 1-3, READ ONLY).  One t9l62 hydro_w dump (float32 primitives after c2p, real cells only),
block by block.  Wedge values (theta pi/4-3pi/4, phi 0-pi), NOT scaled.  cgs.  T = (2/3) eint/rho * 8.1063567271e-9 K; c_g = sqrt(5/3 P/rho);
tau_cell = kappa_R(T, rho) rho dr (run's Rosseland table)."""
import os, sys, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rsgE_page/vispy_4de8168d')
import bin_convert as bc
VC = 3.0e8; TU = 8.1063567271e-9
TAB = '/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/tables/rosseland_mesakap_gs98_x0.70_z0.02.txt'
L_ = open(TAB).readlines()
g = [x.split() for x in L_ if x.startswith('# ') and len(x.split()) == 7 and x.split()[1].isdigit()][-1]
nT, nD, lT0, dlT, lD0, dlD = int(g[1]), int(g[2]), *map(float, g[3:7])
K = np.array([float(x) for x in L_ if not x.startswith('#')]).reshape(nT, nD)
def kap(T, r):
    x = np.clip((np.log10(np.maximum(T, 1)) - lT0) / dlT, 0, nT - 1.001); y = np.clip((np.log10(np.maximum(r, 1e-30)) - lD0) / dlD, 0, nD - 1.001)
    i, j = x.astype(int), y.astype(int); fx, fy = x - i, y - j
    return 10 ** (K[i, j] * (1 - fx) * (1 - fy) + K[i + 1, j] * fx * (1 - fy) + K[i, j + 1] * (1 - fx) * fy + K[i + 1, j + 1] * fx * fy)
out, f = sys.argv[1], sys.argv[2]
h = bc.read_binary(f); nb = h['n_mbs']
TH = (0.5, 0.8, 0.9, 0.99, 0.999)
S = {t: np.zeros(4) for t in TH}                       # n, volume, mass, KE of |v| >= t vceil
ref = dict(all=np.zeros(4), r2e12=np.zeros(4), rho10=np.zeros(4))
clip = np.zeros(3, int)                                # within 1e-6 of vceil; above (1 + 1e-6); all >= (1 - 1e-6)
tail = np.zeros(10, int)                               # |v|/vceil bins 0.99 .. 1.00 (width 1e-3)
edges = np.logspace(7, np.log10(3e8), 16); H = np.zeros((15, 3)); below = np.zeros(3); above = np.zeros(3)
best = None; keep = []
for b in range(nb):
    rf = np.asarray(h['mb_x1f'][b], float); tf = np.asarray(h['mb_x2f'][b], float); pf = np.asarray(h['mb_x3f'][b], float)
    dV = np.diff(pf)[:, None, None] * (np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :, None] * ((rf[1:] ** 3 - rf[:-1] ** 3) / 3.)[None, None, :]
    rc = np.broadcast_to(np.asarray(h['mb_x1v'][b], float)[None, None, :], dV.shape)
    D = {k: np.asarray(h['mb_data'][k][b], float) for k in ('dens', 'velx', 'vely', 'velz', 'eint')}
    rho = D['dens']; v = np.sqrt(D['velx'] ** 2 + D['vely'] ** 2 + D['velz'] ** 2); m = rho * dV; ke = 0.5 * m * v * v
    acc = lambda q: np.array([q.sum(), dV[q].sum(), m[q].sum(), ke[q].sum()])
    ref['all'] += acc(np.ones_like(v, bool)); ref['r2e12'] += acc(rc > 2e12); ref['rho10'] += acc(rho > 1e-10)
    for t in TH:
        S[t] += acc(v >= t * VC)
    u = v / VC
    clip += (int((abs(u - 1) <= 1e-6).sum()), int((u > 1 + 1e-6).sum()), int((u >= 1 - 1e-6).sum()))
    q = (u >= 0.99) & (u < 1 - 1e-6); tail += np.histogram(u[q], bins=np.linspace(0.99, 1.0, 11))[0]
    lo = rho < 1e-11
    for c, w in enumerate((np.ones_like(m), dV, m)):
        H[:, c] += np.histogram(v[lo], bins=edges, weights=w[lo])[0]
        below[c] += w[lo & (v < edges[0])].sum(); above[c] += w[lo & (v > edges[-1])].sum()
    k = np.unravel_index(np.argmax(v), v.shape)
    if best is None or v[k] > best[0]:
        best = (v[k], rc[k], float(h['mb_x2v'][b][k[1]]), float(h['mb_x3v'][b][k[0]]), rho[k], 2 / 3. * D['eint'][k] / rho[k] * TU, D['velx'][k], D['vely'][k], D['velz'][k])
    q = u >= 0.9
    if q.any():
        kk = np.argwhere(q); T = 2 / 3. * D['eint'][q] / rho[q] * TU; cg = np.sqrt(5 / 3. * 2 / 3. * D['eint'][q] / rho[q])
        dr = np.broadcast_to(np.diff(rf)[None, None, :], dV.shape)[q]
        keep.append(np.column_stack([rc[q], rho[q], T, kap(T, rho[q]) * rho[q] * dr, cg / v[q], v[q], dV[q], m[q], D['velx'][q], np.asarray(h['mb_x2v'][b])[kk[:, 1]], np.asarray(h['mb_x3v'][b])[kk[:, 0]]]))
P = []
p = lambda *a: P.append(' '.join(str(x) for x in a))
p('FILE %s' % f); p('t = %.1f s, cycle %d, mtime %s, %d blocks, wedge cells %d' % (h['time'], h['cycle'], __import__('time').strftime('%m-%d %H:%M', __import__('time').localtime(os.path.getmtime(f))), nb, int(ref['all'][0])))
p('wedge totals: volume %.4e mass %.4e g KE %.4e erg | r > 2e12: mass %.4e KE %.4e | rho > 1e-10: mass %.4e KE %.4e' % (ref['all'][1], ref['all'][2], ref['all'][3], ref['r2e12'][2], ref['r2e12'][3], ref['rho10'][2], ref['rho10'][3]))
p('|v| >=      N        volume[cm3]  mass[g]     KE[erg]     M/M(r>2e12) KE/KE(r>2e12) M/M(rho>1e-10) KE/KE(rho>1e-10)')
for t in TH:
    s = S[t]; p('%.3f vceil %8d  %.3e  %.3e  %.3e   %.2e    %.2e      %.2e       %.2e' % (t, s[0], s[1], s[2], s[3], s[2] / ref['r2e12'][2], s[3] / ref['r2e12'][3], s[2] / ref['rho10'][2], s[3] / ref['rho10'][3]))
p('|v|max = %.6e (%.6f vceil) at r %.4e cm, theta %.4f, phi %.4f rad: rho %.3e, T %.3e K, v_r %.4e v_theta %.4e v_phi %.4e' % (best[0], best[0] / VC, *best[1:]))
p('clips: |v| within 1e-6 of vceil: %d; above vceil (1 + 1e-6): %d; smooth tail 0.99 <= |v|/vceil < 1 - 1e-6, bins of 1e-3: %s (sum %d)' % (clip[0], clip[1], tail.tolist(), tail.sum()))
if keep:
    q = np.vstack(keep); pc = lambda x: ' / '.join('%.3e' % y for y in np.percentile(x, [0, 1, 10, 50, 90, 99, 100]))
    p('>= 0.9 vceil set (%d cells), min / p1 / p10 / p50 / p90 / p99 / max:' % len(q))
    for nm, c in (('r [cm]', 0), ('rho', 1), ('T [K]', 2), ('tau_cell', 3), ('c_g/|v|', 4)):
        p('  %-9s %s' % (nm, pc(q[:, c])))
    p('  cells with v_r < 0: %d; theta %.4f-%.4f, phi %.4f-%.4f; mass-weighted <r> %.4e' % ((q[:, 8] < 0).sum(), q[:, 9].min(), q[:, 9].max(), q[:, 10].min(), q[:, 10].max(), np.average(q[:, 0], weights=q[:, 7])))
else:
    p('>= 0.9 vceil set: empty')
p('histogram of |v| in rho < 1e-11 (15 log bins 1e7-3e8 cm/s): lo hi N volume[cm3] mass[g]   (below 1e7: N %d vol %.3e mass %.3e; above 3e8: N %d vol %.3e mass %.3e)' % (below[0], below[1], below[2], above[0], above[1], above[2]))
for i in range(15):
    p('  %.3e %.3e %9d %.3e %.3e' % (edges[i], edges[i + 1], H[i, 0], H[i, 1], H[i, 2]))
open(out, 'w').write('\n'.join(P) + '\n'); print('done', f, flush=True)
