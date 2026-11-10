"""collect.py: Fe I Kp-marg blueshift (= -(Kp-marg (a) median), Seidel+26 observable as on dhj_r1) for the MHD runs
and their hydro twins; points = mean of 2 consecutive dumps (rot k, k+0.5), spread = half the pair difference.
W189: ESPRESSO/MAROON-X mixed by wE = 0.59 as fe_series.py.  Writes fe_table.txt, fe_table.json, fe_vs_rot.png."""
import glob, json, os, re
import numpy as np
O = os.path.dirname(os.path.abspath(__file__))
RUNS = [('r1m7', 'W121 MHD r1m7 (v7, 16x16)', 'w121', None, '#c0392b', 's'),
        ('r1m8', 'W121 MHD r1m8 (v8, 32x32)', 'w121', None, '#e67e22', 'D'),
        ('w121p1', 'W121 hydro p1', 'w121', None, '#1f4e79', 'o'),
        ('w189m1', 'W189 MHD m1 (v7)', 'w189', 0.59, '#c0392b', 's'),
        ('w189m2', 'W189 MHD m2 (v8)', 'w189', 0.59, '#e67e22', 'D'),
        ('w189p1', 'W189 hydro p1', 'w189', 0.59, '#1f4e79', 'o')]
SEI = {'w121': (4.75, 0.30), 'w189': (1.05, 1.59)}

def blue(run, n, c, wE):
    f = '%s/out/ana_%s_%05d_%s.npz' % (O, run, n, c)
    if not os.path.exists(f):
        return None
    z = np.load(f, allow_pickle=True)
    q = -np.array(z['kpm_full'].item()['q'])[::-1]
    if wE is not None:
        m = '%s/out_mx/ana_%s_%05d_%s.npz' % (O, run, n, c)
        if not os.path.exists(m):
            return None
        q = wE*q + (1 - wE)*(-np.array(np.load(m, allow_pickle=True)['kpm_full'].item()['q'])[::-1])
    return float(z['rot']), q

T = {}
for run, lab, pl, wE, col, mk in RUNS:
    ns = sorted(int(re.search(r'_(\d{5})_pc6', f).group(1)) for f in glob.glob('%s/out/ana_%s_*_pc6.npz' % (O, run)))
    pts = {}
    ok = [n for n in ns if blue(run, n, 'pc6', wE) is not None]
    used = set()
    for n in sorted(ok, reverse=True):          # pair each dump with the previous one (newest first)
        if n in used:
            continue
        grp = [n] + ([n - 1] if (n - 1) in ok and (n - 1) not in used else [])
        used.update(grp)
        for m in grp:
            a, b = blue(run, m, 'pc6', wE), blue(run, m, 'pc0', wE)
            pts.setdefault(min(grp), []).append((m, a[0], a[1], b[1] if b is not None else np.full(3, np.nan)))
    rows = []
    for k in sorted(pts):
        P = pts[k]
        v6 = np.array([p[2][1] for p in P]); v0 = np.array([p[3][1] for p in P])
        w = np.mean([p[2][2] - p[2][0] for p in P])/2
        rows.append(dict(rot=float(np.mean([p[1] for p in P])), dumps=[p[0] for p in P], fe=float(v6.mean()),
                         spread=float(np.ptp(v6)/2) if len(P) > 1 else float('nan'), kpm_halfwidth=float(w),
                         fe_nocut=float(v0.mean()), per_dump=[[p[1], float(p[2][1])] for p in P]))
    T[run] = dict(label=lab, planet=pl, rows=rows)
json.dump(T, open(O + '/fe_table.json', 'w'), indent=1)
L = ['run  rot  dumps  Fe blueshift pc6 [km/s]  +-pair spread  (Kp-marg 16-84 half-width)  no-cut']
for run in T:
    for r in T[run]['rows']:
        L.append('%-7s %5.2f %-10s %6.2f  +-%4.2f  (%4.2f)  %6.2f' % (run, r['rot'], ','.join(map(str, r['dumps'])),
                 r['fe'], r['spread'], r['kpm_halfwidth'], r['fe_nocut']))
open(O + '/fe_table.txt', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for ax, pl, ttl in ((axs[0], 'w121', 'WASP-121b'), (axs[1], 'w189', 'WASP-189b')):
    s, e = SEI[pl]
    ax.axhspan(s - e, s + e, color='#999999', alpha=0.25, lw=0, label='Seidel+2026 1 sigma')
    ax.axhline(s, color='#555555', lw=1)
    for run, lab, p, wE, col, mk in RUNS:
        if p != pl or not T[run]['rows']:
            continue
        R = T[run]['rows']
        x = [r['rot'] for r in R]; y = [r['fe'] for r in R]; ye = [0 if np.isnan(r['spread']) else r['spread'] for r in R]
        ax.errorbar(x, y, yerr=ye, fmt=mk + '-', color=col, ms=5, lw=1.2, capsize=2, label=lab + ' (pair mean)')
        pd = np.array([d for r in R for d in r['per_dump']])
        ax.plot(pd[:, 0], pd[:, 1], '.', color=col, ms=4, alpha=0.6)
    ax.set_title(ttl + ': Fe I Kp-marg blueshift (Fe cut p<1e-6 bar)', fontsize=9)
    ax.set_xlabel('rotations (dots: single dumps)'); ax.grid(alpha=0.25, lw=0.5); ax.legend(fontsize=7)
axs[0].set_ylabel('Fe I net blueshift (km/s)')
fig.tight_layout(); fig.savefig(O + '/fe_vs_rot.png', dpi=130)
print(O + '/fe_vs_rot.png')
