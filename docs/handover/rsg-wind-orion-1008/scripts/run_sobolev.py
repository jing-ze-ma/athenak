"""Driver: sanity checks, beta-law grid (Table A), repro winds + ft chromospheres (Table B), plots.
usage: python run_sobolev.py   (CK_DATA set, python-waterboa)"""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import rsglib as L  # noqa: E402
import sobolev as S  # noqa: E402

KMS = 1e5
STARS = ('golden16', 'betelgeuse')
MODES = ('clamp', 'extrap')
BETAS = (0.5, 1., 2.)
VINF = (10., 20., 40.)
MDOT = (1e-7, 1e-6, 1e-5, 1e-4)
XIS = (0., 2., 5.)
MUS = (20., 56.)
IR = slice(0, 7)       # bands longward of 1.32 um (molecular IR: H2O/CO/...)
REPRO, FT = '../repro/', '../ft/'
lines = []


def P(s=''):
    print(s)
    lines.append(s)


def evaluate(st, r, rho, T, dvdr, K, tauR, mode, xi, mu, flux, sel):
    rph = S.r_photo(r, tauR)
    G, c = S.gamma_sob(st, r, rho, T, dvdr, mode, xi*KMS, mu, flux, K=K, r_ph=rph)
    Gs = np.where(sel, G, -1.)
    i = int(np.argmax(Gs))
    fr = c[i]/max(G[i], 1e-300)
    return G, dict(G=G[i], i=i, x=r[i]/st['R'], T=T[i], rho=rho[i],
                   b10=fr[10], b8=fr[8], ir=fr[IR].sum())


# ---------------------------------------------------------------- sanity
P('## Sanity')
for s in STARS:
    st = L.star(s)
    for mode in MODES:
        T = np.array([1500., 2000., 2500., 3000., 3500.])
        rho = np.full(5, 1e-14)
        r = st['R']*np.linspace(1.5, 2., 5)
        K = S.opacities(T, rho, mode)[0]
        th = S.gamma_thin(st, K)
        hi = S.gamma_sob(st, r, rho, T, 1e30, mode, 2e5, 20., 'A', K=K)[0]
        lo = S.gamma_sob(st, r, rho, T, 1e-30, mode, 2e5, 20., 'A', K=K)[0]
        z = S.gamma_sob(st, r, rho, T, 0., mode, 2e5, 20., 'A', K=K)[0]
        P(f'{s} {mode} rho=1e-14 T={T.astype(int).tolist()}: thin {np.round(th, 4).tolist()} '
          f'| Sob(dv/dr=1e30) {np.round(hi, 4).tolist()} | Sob(1e-30) max {lo.max():.2e} '
          f'| Sob(0) max {z.max():.1e}')
P('gmap_table.md golden16 clamp T=2500 rho=1e-14: 0.38 (compare thin entry 3 above)')

# ---------------------------------------------------------------- Table A
rowsA = []
store = {}
hdr = ('| star | mode | beta | v_inf | Mdot | mu_abs | xi | flux | Gmax | r/R | T | rho | '
       'thin there | Gmax r>=1.05R | f(0.26-0.42um) | f(0.61-0.85um) | f(IR>1.3um) | '
       'tau_R(1.001R) |')
rowsA.append(hdr)
rowsA.append('|---'*18 + '|')
viol = 0.
for s, mode in itertools.product(STARS, MODES):
    st = L.star(s)
    for beta, vinf, md in itertools.product(BETAS, VINF, MDOT):
        w = S.beta_wind(st, mode, beta, vinf*KMS, md*L.MSUN/L.YR)
        r, rho, T, K, tauR = w['r'], w['rho'], w['T'], w['K'], w['tauR']
        th = S.gamma_thin(st, K)
        sel = tauR < 1
        sel2 = sel & (r >= 1.05*st['R'])
        for mu, xi, fl in itertools.product(MUS, XIS, 'AB'):
            G, e = evaluate(st, r, rho, T, w['dvdr'], K, tauR, mode, xi, mu, fl, sel)
            if fl == 'A':
                viol = max(viol, np.max(G - th))
            g2 = np.max(np.where(sel2, G, 0.))
            store[(s, mode, beta, vinf, md, mu, xi, fl)] = (e, G, th, w)
            rowsA.append(f'| {s} | {mode} | {beta} | {vinf:.0f} | {md:.0e} | {mu:.0f} | '
                         f'{xi:.0f} | {fl} | {e["G"]:.3f} | {e["x"]:.3f} | {e["T"]:.0f} | '
                         f'{e["rho"]:.1e} | {th[e["i"]]:.3f} | {g2:.3f} | {e["b10"]:.2f} | '
                         f'{e["b8"]:.2f} | {e["ir"]:.2f} | {tauR[0]:.2e} |')
P(f'max(Gamma_Sob(A) - Gamma_thin) over the whole beta grid = {viol:.2e} (must be <= 0)')
open('tableA_full.md', 'w').write('\n'.join(rowsA) + '\n')

P('\n## Table A summary: max over beta of Gamma_Sob max (tau_R<1), mu_abs=20; '
  'cells "A / B (beta at max of A)"')
for s, mode in itertools.product(STARS, MODES):
    P(f'\n### {s} {mode}')
    P('| Mdot \\ v_inf | xi | ' + ' | '.join(f'{v:.0f} km/s' for v in VINF) + ' |')
    P('|---'*(2 + len(VINF)) + '|')
    for md, xi in itertools.product(MDOT, XIS):
        cells = []
        for v in VINF:
            a = [store[(s, mode, b, v, md, 20., xi, 'A')][0]['G'] for b in BETAS]
            bb = [store[(s, mode, b, v, md, 20., xi, 'B')][0]['G'] for b in BETAS]
            cells.append(f'{max(a):.3f} / {max(bb):.3f} ({BETAS[int(np.argmax(a))]})')
        P(f'| {md:.0e} | {xi:.0f} | ' + ' | '.join(cells) + ' |')
for fl in 'AB':
    for mu in MUS:
        k = max((kk for kk in store if kk[7] == fl and kk[5] == mu),
                key=lambda kk: store[kk][0]['G'])
        e = store[k][0]
        P(f'overall max flux {fl} mu {mu:.0f}: {e["G"]:.3f} at {k[:5]} xi={k[6]:.0f} '
          f'r/R={e["x"]:.3f} T={e["T"]:.0f} rho={e["rho"]:.1e} '
          f'band fractions 0.26-0.42um {e["b10"]:.2f} 0.61-0.85um {e["b8"]:.2f} IR {e["ir"]:.2f}')
over = [(k, store[k][0]) for k in store if store[k][0]['G'] > 0.5]
P(f'cases with Gamma_Sob > 0.5: {len(over)} of {len(store)}')
for fl in 'AB':
    sub = [(k, e) for k, e in over if k[7] == fl]
    if sub:
        xs = np.array([e['x'] for k, e in sub])
        Ts = np.array([e['T'] for k, e in sub])
        mds = sorted(set(k[4] for k, e in sub))
        vs = sorted(set(k[3] for k, e in sub))
        P(f'  flux {fl}: {len(sub)} cases, Mdot {mds}, v_inf {vs}, r/R {xs.min():.3f}-'
          f'{xs.max():.3f}, T {Ts.min():.0f}-{Ts.max():.0f} K, '
          f'modes {sorted(set(k[1] for k, e in sub))}, stars {sorted(set(k[0] for k, e in sub))}, '
          f'mu {sorted(set(k[5] for k, e in sub))}, xi {sorted(set(k[6] for k, e in sub))}')
sub = [(k, e) for k, e in over if store[k][0]['x'] >= 1.05]
P(f'  of these with the max at r >= 1.05 R: {len(sub)}')

# ---------------------------------------------------------------- Table B
P('\n## Table B: static Gamma_F max vs Gamma_Sob max (tau_R<1, r>=R), xi = 2 km/s')
P('| star | case | mode | dv/dr | static GF max (r/R, T) | thin max | '
  'Sob A mu20 (r/R, T) | Sob B mu20 (r/R, T) | Sob A mu56 | Sob B mu56 | '
  'f(0.26-0.42) / f(0.61-0.85) / f(IR) at A mu20 max |')
P('|---'*11 + '|')
colB = {}
for s in STARS:
    st = L.star(s)
    R = st['R']
    cases = [(c, None, 'beta1') for c in ('w5v10', 'w5v30', 'w6v10', 'w6v30')] + \
        [(c, None, g) for c in ('ft8', 'ft10', 'ft12', 'ft15') for g in ('v/r', 'v/H')]
    for c, _, how in cases:
        for mode in MODES:
            d = np.load(f'{REPRO if c[0] == "w" else FT}col_{s}_{c}_{mode}.npz')
            r, rho, T, tauR = d['r']*R, d['rho'], d['T'], d['tauR']
            if how == 'beta1':
                md = {'5': 1e-5, '6': 1e-6}[c[1]]
                vinf = float(c.split('v')[1])*KMS
                x = np.maximum(1 - R/r, 0.)
                dvdr = (vinf - 3e5)*R/r**2
            else:
                vc = float(d['v'])
                if how == 'v/r':
                    dvdr = vc/r
                else:
                    H = np.abs(1./np.gradient(np.log(rho), r))
                    dvdr = vc/H
            key = (s, c, mode)
            if key not in colB:
                colB[key] = S.opacities(T, rho, mode)[0]
            K = colB[key]
            sel = (tauR < 1) & (r >= R)
            gs = np.where(sel, d['GF'], -1)
            js = int(np.argmax(gs))
            th = S.gamma_thin(st, K)
            out = {}
            for mu, fl in itertools.product(MUS, 'AB'):
                G, e = evaluate(st, r, rho, T, dvdr, K, tauR, mode, 2., mu, fl, sel)
                out[(mu, fl)] = (e, G)
            store[('B', s, c, mode, how)] = (out, d, th, dvdr)
            a, b = out[(20., 'A')][0], out[(20., 'B')][0]
            P(f'| {s} | {c} | {mode} | {how} | {d["GF"][js]:.3f} ({d["r"][js]:.2f}, '
              f'{T[js]:.0f}) | {np.max(np.where(sel, th, 0)):.3f} | {a["G"]:.3f} '
              f'({a["x"]:.2f}, {a["T"]:.0f}) | {b["G"]:.3f} ({b["x"]:.2f}, {b["T"]:.0f}) | '
              f'{out[(56., "A")][0]["G"]:.3f} | {out[(56., "B")][0]["G"]:.3f} | '
              f'{a["b10"]:.2f} / {a["b8"]:.2f} / {a["ir"]:.2f} |')

open('sobolev_tables.md', 'w').write('\n'.join(lines) + '\n')

# ---------------------------------------------------------------- plots
fig, ax = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
reps = [('golden16', 'clamp', 1., 20., 1e-6), ('golden16', 'clamp', 1., 20., 1e-5),
        ('betelgeuse', 'clamp', 1., 20., 1e-5)]
for a, (s, mode, b, v, md) in zip(ax[0], reps):
    for xi, ls in ((0., ':'), (2., '-'), (5., '--')):
        for fl, col in (('A', 'C0'), ('B', 'C1')):
            e, G, th, w = store[(s, mode, b, v, md, 20., xi, fl)]
            a.plot(w['r']/L.star(s)['R'], G, ls, color=col, label=f'Sob {fl} xi={xi:.0f}')
    a.plot(w['r']/L.star(s)['R'], th, 'k-', lw=0.8, label='thin')
    a.set(xscale='log', yscale='log', ylim=(1e-4, 2), xlabel='r/R', ylabel='Gamma',
          title=f'{s} {mode} beta={b} v_inf={v:.0f} Mdot={md:.0e} (mu 20)')
    a.legend(fontsize=7)
repsB = [('golden16', 'w5v10', 'clamp', 'beta1'), ('golden16', 'ft12', 'clamp', 'v/r'),
         ('betelgeuse', 'ft15', 'clamp', 'v/r')]
for a, (s, c, mode, how) in zip(ax[1], repsB):
    out, d, th, dvdr = store[('B', s, c, mode, how)]
    for fl, col in (('A', 'C0'), ('B', 'C1')):
        a.plot(d['r'], out[(20., fl)][1], '-', color=col, label=f'Sob {fl} ({how})')
        a.plot(d['r'], out[(56., fl)][1], ':', color=col, label=f'Sob {fl} mu56')
    if how == 'v/r':
        o2 = store[('B', s, c, mode, 'v/H')][0]
        a.plot(d['r'], o2[(20., 'A')][1], '-', color='C2', label='Sob A (v/H)')
    a.plot(d['r'], d['GF'], 'r-', label='static Gamma_F')
    a.plot(d['r'], th, 'k-', lw=0.8, label='thin')
    a.set(xscale='log', yscale='log', ylim=(1e-4, 2), xlim=(1, 10), xlabel='r/R',
          title=f'{s} {c} {mode} xi=2')
    a.legend(fontsize=7)
fig.savefig('gamma_sob_r.png', dpi=110)

fig, ax = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True)
for a, (s, mode) in zip(ax, itertools.product(STARS, MODES)):
    for v, col in zip(VINF, ('C0', 'C1', 'C2')):
        for fl, ls in (('A', '-'), ('B', '--')):
            y = [max(store[(s, mode, b, v, md, 20., 2., fl)][0]['G'] for b in BETAS)
                 for md in MDOT]
            a.plot(MDOT, y, ls, marker='o', color=col, label=f'v_inf {v:.0f} {fl}')
    a.axhline(0.5, color='k', lw=0.5)
    a.set(xscale='log', yscale='log', ylim=(1e-3, 1.5), xlabel='Mdot [Msun/yr]',
          ylabel='max_beta max_r Gamma_Sob', title=f'{s} {mode}, xi=2, mu=20')
    a.legend(fontsize=7)
fig.savefig('gamma_sob_mdot.png', dpi=110)
