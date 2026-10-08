"""Gamma_Sob on the MESA FT columns (../mesa/col_<model>_<mft|mft12>_<mode>.npz, read-only).
usage: python mesa_sob.py  -> mesa_sob.md, mesa_sob.png"""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import rsglib as L  # noqa: E402
import sobolev as S  # noqa: E402

D = '../mesa/'
MODELS = ['m10lgl4.4a2', 'm10lgl4.6', 'm10lgl4.8', 'm10lgl5.0', 'm15lgl5.1', 'm20lgl5.3',
          'm20lgl5.4', 'm20lgl5.5', 'm26_p455']
WL = L.WL
BNAME = [f'{WL[b+1]:.2f}-{WL[b]:.2f}um' for b in range(L.NB)]
static_tab = {}
for ln in open(D + 'mesa_table.md'):
    p = [x.strip() for x in ln.split('|')]
    if len(p) > 5 and p[1] in MODELS and p[2] in ('mft', 'mft12'):
        static_tab[(p[1], p[2], p[3])] = float(p[4])

rows = ['| model | log L | Teff | case | mode | static GF max (table; npz r>=R) | '
        'Sob A v/r ; v/H ; beta1 | Sob B v/r ; v/H ; beta1 | largest (r/R, T, dominant band, '
        'its fraction) | Sob A v/H mu56 | Sob B v/H mu56 |', '|---'*11 + '|']
res = {}
for m, case, mode in itertools.product(MODELS, ('mft', 'mft12'), ('clamp', 'extrap')):
    d = np.load(f'{D}col_{m}_{case}_{mode}.npz')
    st = L.star(m)
    R = st['R']
    st = dict(st, Teff=float(d['Teff']), kE=float(d['kE']))   # same Teff, kE as the column
    r, rho, T, tauR = d['r']*R, d['rho'], d['T'], d['tauR']
    vc = float(d['v'])
    K = S.opacities(T, rho, mode)[0]
    rph = S.r_photo(r, tauR)
    sel = (tauR < 1) & (r >= R)
    H = np.abs(1./np.gradient(np.log(rho), r))
    grads = {'v/r': vc/r, 'v/H': vc/H}
    if case == 'mft':
        grads['beta1'] = np.where(r >= R, (30e5 - vc)*R/r**2, vc/r)
    out = {}
    for (g, dv), mu, fl in itertools.product(grads.items(), (20., 56.), 'AB'):
        G, c = S.gamma_sob(st, r, rho, T, dv, mode, 2e5, mu, fl, K=K, r_ph=rph)
        i = int(np.argmax(np.where(sel, G, -1.)))
        out[(g, mu, fl)] = (G[i], i, c[i]/G[i])
    gs = np.max(np.where(sel, d['GF'], 0.))
    res[(m, case, mode)] = (out, st, gs, rph/R)
    best = max(((k, v) for k, v in out.items() if k[1] == 20.), key=lambda kv: kv[1][0])
    (g, mu, fl), (Gb, ib, fr) = best
    bb = int(np.argmax(fr))

    def f3(fl):
        return ' ; '.join(f'{out[(g, 20., fl)][0]:.3f}' if (g, 20., fl) in out else '-'
                          for g in ('v/r', 'v/H', 'beta1'))
    rows.append(f'| {m} | {np.log10(st["L"]/L.LSUN):.2f} | {st["Teff"]:.0f} | {case} | {mode} | '
                f'{static_tab.get((m, case, mode), np.nan):.3f}; {gs:.3f} | {f3("A")} | '
                f'{f3("B")} | {fl} {g}: {r[ib]/R:.2f}, {T[ib]:.0f} K, {BNAME[bb]}, '
                f'{fr[bb]:.2f} | {out[("v/H", 56., "A")][0]:.3f} | '
                f'{out[("v/H", 56., "B")][0]:.3f} |')
allb = [(k, kk, v) for k, (o, *_) in res.items() for kk, v in o.items()]
for mu in (20., 56.):
    k, kk, v = max((a for a in allb if a[1][1] == mu), key=lambda a: a[2][0])
    d = np.load(f'{D}col_{k[0]}_{k[1]}_{k[2]}.npz')
    fr = v[2]
    rows.append(f'\noverall max mu {mu:.0f}: {v[0]:.3f} at {k} grad {kk[0]} flux {kk[2]} '
                f'r/R {d["r"][v[1]]:.2f} T {d["T"][v[1]]:.0f} rho {d["rho"][v[1]]:.1e}; '
                f'bands: ' + ', '.join(f'{BNAME[b]} {fr[b]:.2f}' for b in np.argsort(fr)[::-1][:3]))
for g in ('v/r', 'v/H', 'beta1'):
    k, kk, v = max((a for a in allb if a[1][0] == g and a[1][1] == 20.),
                   key=lambda a: a[2][0])
    rows.append(f'max over models, grad {g}, mu 20: {v[0]:.3f} at {k} flux {kk[2]}')
rows.append('r(tau_R=2/3)/R used for flux B: ' + ', '.join(
    f'{k[0]}/{k[1]}/{k[2]} {v[3]:.4f}' for k, v in res.items() if k[2] == 'clamp'))
s = '\n'.join(rows) + '\n'
open('mesa_sob.md', 'w').write(s)
print(s)

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
for a, case in zip(ax, ('mft', 'mft12')):
    for mode, mk in (('clamp', 'o'), ('extrap', 's')):
        ks = [k for k in res if k[1] == case and k[2] == mode]
        lL = [np.log10(res[k][1]['L']/L.LSUN) for k in ks]
        o = np.argsort(lL)
        lL = np.array(lL)[o]
        ks = [ks[i] for i in o]
        lab = mode
        a.plot(lL, [res[k][0][('v/H', 20., 'A')][0] for k in ks], '-', marker=mk,
               color='C0', label=f'Sob A v/H {lab}')
        a.plot(lL, [res[k][0][('v/H', 20., 'B')][0] for k in ks], '-', marker=mk,
               color='C1', label=f'Sob B v/H {lab}')
        a.plot(lL, [res[k][2] for k in ks], '-', marker=mk, color='r',
               label=f'static GF {lab}')
    a.set(yscale='log', ylim=(5e-3, 1), xlabel='log L/Lsun', ylabel='max Gamma (tau_R<1, r>=R)',
          title=f'MESA {case}, xi=2 km/s, mu_abs=20')
    a.axhline(0.5, color='k', lw=0.5)
    a.legend(fontsize=7)
fig.savefig('mesa_sob.png', dpi=110)
