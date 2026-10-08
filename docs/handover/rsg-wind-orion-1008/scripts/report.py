"""Tables (tables.md) and plots for the two-phase test; reads part1_K_*.npz, rows2/3.pkl."""
import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import twophase as TP             # noqa: E402
import rsglib as L                # noqa: E402

C3 = ['#2a78d6', '#eb6834', '#555550']
out = []
P = out.append

# ---------- Part 1 ----------
P('## Part 1: optically thin RE roots (s = stable, u = unstable), T grid 300-4000 K\n')
P('Fixed P (isobaric, Field). Columns: P = 1e-6, 1e-4, 1e-2 dyn/cm^2 '
  '(rho ~ 1e-17..1e-13 at these T); fixed-rho roots are listed in part1_roots_rho.txt.\n')
P('| star | mode | r/R | W | P=1e-6 | P=1e-4 | P=1e-2 | two stable (Tw/Tc) at 1e-4 | rho_c/rho_w |')
P('|---|---|---|---|---|---|---|---|---|')
txt = []
for star in ('golden16', 'm20lgl5.5'):
    for mode in ('clamp', 'extrap'):
        st, tab, res = TP.part1(star, mode)
        for x in TP.XR:
            cells = []
            for Pv in (1e-6, 1e-4, 1e-2):
                i = int(np.argmin(abs(np.log(TP.PG/Pv))))
                cells.append(','.join('%d%s' % (t, 's' if q else 'u')
                                      for t, q in res['P'][x][i]) or 'none')
            i = int(np.argmin(abs(np.log(TP.PG/1e-4))))
            pr = TP.two_stable(res['P'][x][i])
            if pr:
                mw = np.interp(pr[0], TP.TG, tab['mup'][i])
                mc = np.interp(pr[1], TP.TG, tab['mup'][i])
                ratio = '%.1f' % (pr[0]*mc/(pr[1]*mw))
                prs = '%d / %d' % pr
            else:
                ratio, prs = '-', 'no'
            P(f'| {star} | {mode} | {x} | {TP.dilution(x):.4f} | ' + ' | '.join(cells)
              + f' | {prs} | {ratio} |')
            for j, rho in enumerate(TP.RHOG):
                txt.append(f'{star} {mode} x={x} rho={rho:.2e} ' + ','.join(
                    '%d%s' % (t, 's' if q else 'u') for t, q in res['rho'][x][j]))
open(TP.OUT + 'part1_roots_rho.txt', 'w').write('\n'.join(txt) + '\n')

# plot: T_RE vs rho per r
fig, ax = plt.subplots(1, 5, figsize=(17, 4), sharey=True)
for k, x in enumerate(TP.XR):
    a = ax[k]
    for c, star in zip(C3, ('golden16', 'm20lgl5.5')):
        for mode, mk in (('clamp', 'o'), ('extrap', 'D')):
            st, tab, res = TP.part1(star, mode)
            for j, rho in enumerate(TP.RHOG):
                for t, q in res['rho'][x][j]:
                    a.plot(rho*(1.06 if mode == 'extrap' else 1), t, mk, ms=5 if q else 4,
                           mfc=c if q else 'none', mec=c, alpha=0.9 if mode == 'clamp'
                           else 0.5)
    a.axhspan(300, 800, color='0.9', zorder=0)
    a.set_xscale('log'); a.set_yscale('log')
    a.set_title(f'r = {x} R  (W = {TP.dilution(x):.3f})')
    a.set_xlabel(r'$\rho$ [g cm$^{-3}$]')
    a.grid(alpha=0.25, lw=0.5)
ax[0].set_ylabel(r'$T_{\rm RE}$ [K]')
for c, star in zip(C3, ('golden16 (Teff 4106)', 'm20lgl5.5 (Teff 3776)')):
    ax[0].plot([], [], 'o', color=c, label=star)
ax[0].plot([], [], 'o', color='k', label='stable (filled)')
ax[0].plot([], [], 'o', mfc='none', mec='k', label='unstable (open)')
ax[0].plot([], [], 'D', color='k', alpha=0.5, label='extrap mode (diamond)')
ax[0].legend(fontsize=7, loc='upper right')
fig.suptitle('Optically thin radiative-equilibrium gas temperatures, ck opacities '
             '(grey band: below the requested 800 K range)', fontsize=10)
fig.tight_layout(); fig.savefig(TP.OUT + 'part1_TRE_rho.png', dpi=130)

# plot: heating curves golden16 clamp at P=1e-4
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for a, star in zip(ax, ('golden16', 'm20lgl5.5')):
    st, tab, res = TP.part1(star, 'clamp')
    i = int(np.argmin(abs(np.log(TP.PG/1e-4))))
    for k, x in enumerate(TP.XR):
        W = TP.dilution(x)
        kb = (L.GW*tab['Kp'][i]).sum(-1)
        ab = (kb*W*TP.Bb(st['Teff'])).sum(-1)
        em = (kb*TP.Bb(TP.TG)).sum(-1)
        a.plot(TP.TG, ab/em, lw=1.8, label=f'{x} R', color=plt.cm.viridis(k/4.5))
    a.axhline(1, color='k', lw=0.8)
    a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(0.05, 50)
    a.set_title(f'{star}, clamp, P = 1e-4 dyn/cm$^2$')
    a.set_xlabel('gas T [K]'); a.grid(alpha=0.25, lw=0.5)
ax[0].set_ylabel('absorbed / emitted (heating > 1)')
ax[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(TP.OUT + 'part1_heating.png', dpi=130)

# ---------- Part 2/3 ----------
R2 = pickle.load(open(TP.OUT + 'rows2.pkl', 'rb'))
R3 = pickle.load(open(TP.OUT + 'rows3.pkl', 'rb'))

P('\n## Part 2: gas only, max over (case, mode, grad, f, T_c choice)\n')
P('| model | r/R | Tbar | rho_bar | Gamma_w(Sob) max | Gamma_c gas(Sob) max | Gamma_bar gas max '
  '| Gamma_c thin max | RE pair (mft clamp, f=0.1) |')
P('|---|---|---|---|---|---|---|---|---|')
for m in TP.MODELS:
    for x in TP.XR:
        rr = [r for r in R2 if r['model'] == m and r['x'] == x]
        if not rr:
            continue
        gb = max(r['fm']*r['Gcgas'] + (1-r['fm'])*r['Gw'] for r in rr)
        ref = [r for r in rr if r['case'] == 'mft' and r['mode'] == 'clamp'
               and r['src'] == 'RE' and r['f'] == 0.1]
        pair = '%d/%d K, f_m %.2f' % (ref[0]['Tw'], ref[0]['Tc'], ref[0]['fm']) if ref \
            else 'none'
        r0 = [r for r in rr if r['case'] == 'mft' and r['mode'] == 'clamp'][0]
        P(f'| {m} | {x} | {r0["Tbar"]:.0f} | {r0["rhobar"]:.1e} | '
          f'{max(r["Gw"] for r in rr):.3f} | {max(r["Gcgas"] for r in rr):.3f} | {gb:.3f} | '
          f'{max(r["Gcthin"] for r in rr):.3f} | {pair} |')

# condensation map
P('\n## Part 3a: where the condensation criteria are met (T_cond optimistic / nominal)\n')
P('| model | r/R | T_d p=1 | T_d p=0 | T_d p=-1 | Al2O3 (1600/1400) grain p=1, 0, -1 '
  '| silicate (1200/1000) grain p=1, 0, -1 |')
P('|---|---|---|---|---|---|---|')
for m in TP.MODELS:
    st = L.star(m)
    for x in TP.XR:
        W = TP.dilution(x)
        td = [TP.tgrain(st['Teff'], W, p) for p in (1, 0, -1)]
        cells = []
        for sp in ('Al2O3', 'silicate'):
            lo, hi = TP.SPECIES[sp]['Tc']
            cells.append(', '.join(('yes' if t < lo else 'opt' if t < hi else 'no')
                                   for t in td))
        P(f'| {m} | {x} | {td[0]:.0f} | {td[1]:.0f} | {td[2]:.0f} | ' + ' | '.join(cells)
          + ' |')
P('\n(yes = below the nominal (lower) T_cond; opt = only below the optimistic (upper) '
  'T_cond.)  Gas criterion: T_c (cool-phase gas T) < T_cond; with the RE pair T_c = '
  '700-940 K at 3-5 R, so both species pass the gas criterion wherever a two-phase RE '
  'solution exists; in the T_c scan, 1000/1200/1400 K.\n')


def where(r):
    return (f"r={r['x']} {r['case']}/{r['mode']} {r['grad']} f={r['f']} {r['src']} "
            f"Tc={r['Tc']:.0f} fc={r['fcond']} kd={r['kd']:.0f} l={r['l']}")


P('\n## Part 3b: max Gamma_c and Gamma_bar (T_cond optimistic; porosity on, l = 0.01 r '
  'or 0.1 r; l=0 column = no porosity)\n')
P('| model | species | criterion | r range | max Gamma_c | where | max Gamma_bar | where '
  '| max Gamma_c, l=0 | max Gamma_bar, l=0 |')
P('|---|---|---|---|---|---|---|---|---|---|')
for m in TP.MODELS:
    for sp in TP.SPECIES:
        for crit in ('gas', 'grain_p1', 'grain_p0', 'grain_pm1'):
            for lab, xr in (('2-4', (2., 4.)), ('1.5-5', (1.5, 5.))):
                rr = [r for r in R3 if r['model'] == m and r['sp'] == sp and
                      r['crit'] == crit and r['Tcond'] == 'opt' and
                      xr[0] <= r['x'] <= xr[1]]
                rp = [r for r in rr if r['l'] > 0]
                r0 = [r for r in rr if r['l'] == 0]
                a = max(rp, key=lambda r: r['Gc'])
                b = max(rp, key=lambda r: r['Gbar'])
                P(f'| {m} | {sp} | {crit} | {lab} | {a["Gc"]:.3f} | {where(a)} | '
                  f'{b["Gbar"]:.3f} | {where(b)} | {max(r["Gc"] for r in r0):.3f} | '
                  f'{max(r["Gbar"] for r in r0):.3f} |')

P('\n## Part 3c: minimum delta*f_cond*kappa_d [cm^2/g gas] for Gamma_c > 1 '
  '(case mft, clamp, grad v/H, f = 0.1; cool phase = RE pair if it exists, else '
  'T_c = 1200 K scan); inf = impossible: max effective opacity 1/(rho_c l) < kappa_Edd\n')
P('| model | r/R | kappa_Edd | T_c | rho_c | kappa_gas,c(Sob) | X_req l=0 | X_req l=0.01r '
  '| 1/(rho_c l) l=0.01r | X_req l=0.1r | 1/(rho_c l) l=0.1r |')
P('|---|---|---|---|---|---|---|---|---|---|---|')
for m in TP.MODELS:
    for x in TP.XR:
        rr = [r for r in R2 if r['model'] == m and r['x'] == x and r['case'] == 'mft'
              and r['mode'] == 'clamp' and r['grad'] == 'v/H' and r['f'] == 0.1]
        sel = [r for r in rr if r['src'] == 'RE'] or [r for r in rr if r['Tc'] == 1200.]
        if not sel:
            continue
        r = sel[0]
        xs = [TP.x_required(r, lc) for lc in (0., 0.01, 0.1)]
        mx = [1/(r['rc']*lc*r['x']*r['R']) for lc in (0.01, 0.1)]
        P(f"| {m} | {x} | {r['kE']:.2f} | {r['Tc']:.0f} | {r['rc']:.1e} | "
          f"{r['Gcgas']*r['kE']:.3f} | {xs[0]:.2f} | {xs[1]:.2f} | {mx[0]:.0f} | "
          f"{xs[2]:.2f} | {mx[1]:.1f} |")
P('\nSupply (delta_max * f_cond=1 * kappa_d=3000): Al2O3 1e-4*3000 = 0.30; silicate '
  '4e-3*3000 = 12 cm^2/g gas (kappa_d = 300: 0.03 / 1.2).\n')

open(TP.OUT + 'tables.md', 'w').write('\n'.join(out) + '\n')

# plot: Gamma_c, Gamma_bar vs r
fig, ax = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
for j, m in enumerate(TP.MODELS):
    for i, sp in enumerate(('silicate', 'Al2O3')):
        a = ax[i, j]
        for c, crit in zip(C3, ('gas', 'grain_pm1', 'grain_p1')):
            for q, ls in (('Gc', '-'), ('Gbar', '--')):
                ys = []
                for x in TP.XR:
                    rr = [r for r in R3 if r['model'] == m and r['x'] == x and
                          r['case'] == 'mft' and r['mode'] == 'clamp' and
                          r['grad'] == 'v/H' and r['f'] == 0.1 and r['sp'] == sp and
                          r['crit'] == crit and r['Tcond'] == 'opt' and
                          r['fcond'] == 1. and r['kd'] == 3000. and r['l'] == 0.01]
                    sel = [r for r in rr if r['src'] == 'RE'] or \
                        [r for r in rr if r['Tc'] == 1200.]
                    ys.append(sel[0][q] if sel else np.nan)
                a.plot(TP.XR, ys, ls, color=c, lw=2, marker='o', ms=5,
                       label=f'{crit} {"Gamma_c" if q == "Gc" else "Gamma_bar"}')
        a.axhline(1, color='k', lw=0.8); a.axhline(0.5, color='k', lw=0.5, ls=':')
        a.set_yscale('log'); a.set_ylim(1e-3, 30); a.grid(alpha=0.25, lw=0.5)
        a.set_title(f'{m}: {sp}, f_cond=1, kappa_d=3000, f=0.1, l=0.01 r', fontsize=9)
        if i == 1:
            a.set_xlabel('r / R')
ax[0, 0].set_ylabel(r'$\Gamma$'); ax[1, 0].set_ylabel(r'$\Gamma$')
ax[0, 0].legend(fontsize=7)
fig.suptitle('Cool-phase (solid) and mass-weighted (dashed) Eddington factor; mft clamp, '
             'grad v/H; cool phase = RE pair (3-5 R) or T_c = 1200 K scan', fontsize=10)
fig.tight_layout(); fig.savefig(TP.OUT + 'gamma_r.png', dpi=130)
