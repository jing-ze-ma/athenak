#!/usr/bin/env python3
"""lowbeta-tests-1102: collect metrics -> results.json, draw the page figures."""
import glob, json, math, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

B = '/resnick/groups/carnegie_poc/jingze/lowbeta_1102'
IMG = '/resnick/home/jingze/ATHENAK/pages/lowbeta_cart/img'
os.makedirs(IMG, exist_ok=True)

R = {}
for f in sorted(glob.glob(f'{B}/runs/*/*/metrics.json')):
    m = json.load(open(f))
    R.setdefault(m['test'], {})[m['case']] = m

# colours by Riemann solver, dash by corner EMF
COL = {'hlle': '#2a78d6', 'hlld': '#eb6834', 'lhlld': '#1baf7a'}
DASH = {'gs': '-', 'sm': '--', 'uct': ':'}
SCH = ['hlle_gs', 'hlle_sm', 'hlld_gs', 'hlld_sm', 'hlld_uct', 'lhlld_gs', 'lhlld_sm', 'lhlld_uct']
NAME = {'hlle': 'HLLE', 'hlld': 'HLLD', 'lhlld': 'LHLLD', 'gs': 'GS07', 'sm': 'emf_smooth',
        'uct': 'UCT-HLLD'}
def lab(s):
    a = s.split('_')
    t = NAME[a[0]] + ' + ' + NAME[a[1]]
    if len(a) > 2: t += ' + FOFC'
    return t
def sty(s):
    a = s.split('_')
    if len(a) > 2:        # FOFC arm
        return dict(color=COL[a[0]], ls='-.', lw=1.4, marker='s', ms=5)
    if a[0] == 'lhlld':   # hollow + thinner so an identical HLLD curve shows beneath
        return dict(color=COL[a[0]], ls=DASH[a[1]], lw=1.3, marker='o', ms=7, mfc='none')
    return dict(color=COL[a[0]], ls=DASH[a[1]], lw=2.2, marker='o', ms=5)

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e6e6e6', 'grid.linewidth': 0.6,
                     'axes.edgecolor': '#888', 'savefig.dpi': 130, 'figure.facecolor': 'white'})
FLOOR = 1e-30


def save(fig, name):
    for ax in fig.axes:
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.savefig(f'{IMG}/{name}.png', bbox_inches='tight')
    plt.close(fig)


def legend_all(fig, ax, **kw):
    seen = {}
    for a in fig.axes:
        for h, l in zip(*a.get_legend_handles_labels()):
            seen.setdefault(l, h)
    order_ = [lab(x) for x in SCH + ['hlld_gs_fofc', 'hlle_gs_fofc']]
    ls_ = sorted(seen, key=lambda l: order_.index(l) if l in order_ else 99)
    ax.legend([seen[l] for l in ls_], ls_, **kw)


def slope_guides(ax, Ns, ref):
    Ns = np.array(Ns, float)
    for o, c in ((1, '#aaa'), (2, '#666')):
        y = ref*(Ns/Ns[0])**(-o)
        ax.plot(Ns, y, color=c, lw=1, ls='-.', zorder=0)
        ax.text(Ns[-1]*1.05, y[-1], f'{o}st' if o == 1 else '2nd', fontsize=8, color=c,
                va='center')


def order(Ns, E):
    Ns = np.array(Ns, float); E = np.array(E, float)
    return np.log(E[:-1]/E[1:])/np.log(Ns[1:]/Ns[:-1])


summary = {}

# ------------------------------------------------------------------ CPAW
BETAS = [0.2, 1e-2, 1e-3, 1e-4, 1e-5]
def cpaw_series(test, rec, sc, beta, key):
    Ns, E = [], []
    for N in (32, 64, 128, 256):
        m = R.get(test, {}).get(f'b{beta:.0e}_{rec}_{sc}_N{N}')
        if m and key in m and m.get('completed') and np.isfinite(m[key]):
            Ns.append(N); E.append(m[key])
    return Ns, E

for test, schemes in (('cpaw2d', SCH), ('cpaw1d', ['hlle_gs', 'hlld_gs', 'lhlld_gs'])):
    if test not in R: continue
    summary[test] = {}
    for rec in ('plm', 'ppm4'):
        for key, kl in (('L1_Bperp', 'L1 error of B_perp'), ('L1_p', 'L1 of p - p0'),
                        ('L1_vperp', 'L1 error of v_perp')):
            fig, axs = plt.subplots(1, 5, figsize=(17, 3.6), sharey=False)
            for ax, beta in zip(axs, BETAS):
                for sc in schemes:
                    Ns, E = cpaw_series(test, rec, sc, beta, key)
                    if Ns:
                        ax.loglog(Ns, E, label=lab(sc), **sty(sc))
                        summary[test].setdefault(rec, {}).setdefault(key, {}).setdefault(sc, {})[f'{beta:g}'] = dict(N=Ns, E=E, order=order(Ns, E).tolist() if len(Ns) > 1 else [])
                Ns0, E0 = cpaw_series(test, rec, 'hlld_gs', beta, key)
                if Ns0: slope_guides(ax, [32, 64, 128, 256], E0[0]*1.5)
                ax.set_title(f'beta = {beta:g}')
                ax.set_xlabel('N')
                ax.set_xticks([32, 64, 128, 256]); ax.set_xticklabels(['32', '64', '128', '256'])
            axs[0].set_ylabel(kl)
            legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.12, 0.5),
                           frameon=False)
            save(fig, f'{test}_{key}_{rec}')
    # order vs beta (slope between N=128 and 256)
    for key in ('L1_Bperp', 'L1_p'):
        fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
        for ax, rec in zip(axs, ('plm', 'ppm4')):
            for sc in schemes:
                bs, os_ = [], []
                for beta in BETAS:
                    Ns, E = cpaw_series(test, rec, sc, beta, key)
                    if len(Ns) >= 2 and Ns[-1] == 256 and Ns[-2] == 128:
                        bs.append(beta); os_.append(order(Ns, E)[-1])
                if bs: ax.semilogx(bs, os_, label=lab(sc), **sty(sc))
            ax.axhline(2, color='#666', lw=1, ls='-.'); ax.axhline(0, color='#999', lw=0.8)
            ax.set_title(f'{test}, {rec.upper()}: order of {key[3:]} (N = 128 -> 256)')
            ax.set_xlabel('plasma beta'); ax.invert_xaxis()
        axs[0].set_ylabel('measured order')
        axs[1].legend(fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
        save(fig, f'{test}_order_{key}')
    # min p vs time, N=128 PLM, beta 1e-3 and 1e-4 and 1e-5
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.6))
    for ax, beta in zip(axs, (1e-3, 1e-4, 1e-5)):
        p0 = beta*0.505
        for sc in schemes:
            m = R[test].get(f'b{beta:.0e}_plm_{sc}_N128')
            if not m or 'pmin_t' not in m: continue
            t = np.array(m['p_t'])/m['T']; pm = np.array(m['pmin_t'])/p0
            ax.plot(t, np.maximum(pm, 1e-3), label=lab(sc), **dict(sty(sc), marker=None))
        ax.set_yscale('log'); ax.set_ylim(1e-3*0.8, 2)
        ax.axhline(1, color='#666', lw=1, ls='-.')
        ax.set_title(f'beta = {beta:g}, N = 128 (PLM)'); ax.set_xlabel('t / wave period')
    axs[0].set_ylabel('min p / exact p  (floored at 1e-3)')
    legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, f'{test}_pmin_t')
    # floor hits map: total efloor counts at each (beta, N), PLM
    tab = {}
    for sc in schemes:
        for beta in BETAS:
            for N in (32, 64, 128, 256):
                for rec in ('plm', 'ppm4'):
                    m = R[test].get(f'b{beta:.0e}_{rec}_{sc}_N{N}')
                    if m:
                        tab[f'{rec}|{sc}|{beta:g}|{N}'] = dict(efloor=m['counts']['efloor'],
                            tfloor=m['counts']['tfloor'], dfloor=m['counts']['dfloor'],
                            pmin_over_p0=min(m.get('pmin_t', [np.nan]))/(beta*0.505),
                            completed=m.get('completed'), wall=m.get('wall_s'))
    summary[test]['floors'] = tab

# ------------------------------------------------------------------ linear waves
if 'lw' in R:
    summary['lw'] = {}
    LB = [0.369, 1e-2, 1e-3, 1e-4]
    for mode, mt in (('fast', 'fast'), ('alf', 'Alfven'), ('slow', 'slow')):
        fig, axs = plt.subplots(1, 4, figsize=(15, 3.5))
        for ax, beta in zip(axs, LB):
            for sc in SCH:
                Ns, E = [], []
                for N in (32, 64, 128, 256):
                    m = R['lw'].get(f'b{beta:.0e}_{mode}_{sc}_N{N}')
                    if m and 'err_rms' in m:
                        Ns.append(N); E.append(m['err_rms'])
                if Ns:
                    ax.loglog(Ns, E, label=lab(sc), **sty(sc))
                    summary['lw'].setdefault(mode, {}).setdefault(sc, {})[f'{beta:g}'] = dict(N=Ns, E=E, order=order(Ns, E).tolist() if len(Ns) > 1 else [])
            ax.set_title(f'{mt}, beta = {beta:g}'); ax.set_xlabel('N')
            ax.set_xticks([32, 64, 128, 256]); ax.set_xticklabels(['32', '64', '128', '256'])
            Ns0 = [32, 64, 128, 256]
            m0 = R['lw'].get(f'b{beta:.0e}_{mode}_hlld_gs_N32')
            if m0 and 'err_rms' in m0: slope_guides(ax, Ns0, m0['err_rms']*1.5)
        if mode == 'slow':
            for ax, beta in zip(axs, LB):
                if beta <= 1e-3:
                    ax.text(0.97, 0.97, 'HLLD / LHLLD + GS07:\nblew up at N >= %d\n(p floor every step, dt collapse)' % (64 if beta > 5e-4 else 32), transform=ax.transAxes,
                            ha='right', va='top', fontsize=8, color='#c62828')
        axs[0].set_ylabel('RMS L1 error (amp 1e-6)')
        legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.12, 0.5), frameon=False)
        save(fig, f'lw_{mode}')
    summary['lw']['status'] = {c: dict(completed=bool(m.get('completed')), rc=m.get('rc'),
        floor_hits=m['counts']['efloor'] + m['counts']['tfloor'],
        first_floor_cycle=m.get('first_floor_cycle'), end_time=m.get('end_time'))
        for c, m in R['lw'].items()}
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.6), sharey=True)
    for ax, mode in zip(axs, ('fast', 'alf', 'slow')):
        for sc in SCH:
            d = summary['lw'].get(mode, {}).get(sc, {})
            bs = [b for b in LB if f'{b:g}' in d and d[f'{b:g}']['order']]
            if bs: ax.semilogx(bs, [d[f'{b:g}']['order'][-1] for b in bs], label=lab(sc), **sty(sc))
        ax.axhline(2, color='#666', lw=1, ls='-.'); ax.set_title(f'{mode}: order at the finest pair')
        ax.set_xlabel('plasma beta'); ax.invert_xaxis()
    axs[0].set_ylabel('measured order')
    legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, 'lw_order')

# ------------------------------------------------------------------ vortex
if 'vortex' in R:
    summary['vortex'] = {}
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.6))
    NV = [32, 64, 128, 256, 512, 1024]
    for sc in SCH:
        Ns, Ep, Ev, Pm, Fl = [], [], [], [], []
        for N in NV:
            m = R['vortex'].get(f'{sc}_N{N}')
            if m and 'L1_p' in m:
                Ns.append(N); Ep.append(m['L1_p']); Ev.append(m['L1_v'])
                Pm.append(min(m.get('pmin_t', [m.get('pmin_final', np.nan)])))
                Fl.append(m['counts']['efloor'] + m['counts']['tfloor'])
        if not Ns: continue
        axs[0].loglog(Ns, Ev, label=lab(sc), **sty(sc))
        axs[1].loglog(Ns, Ep, label=lab(sc), **sty(sc))
        axs[2].loglog(Ns, np.maximum(Pm, 1e-16), label=lab(sc), **sty(sc))
        summary['vortex'][sc] = dict(N=Ns, L1_p=Ep, L1_v=Ev, pmin=Pm, floor_hits=Fl,
                                     order_p=order(Ns, Ep).tolist(), order_v=order(Ns, Ev).tolist())
    m0 = R['vortex'].get('hlld_gs_N32')
    if m0:
        slope_guides(axs[0], NV, m0['L1_v']*1.5); slope_guides(axs[1], NV, m0['L1_p']*1.5)
    axs[0].set_ylabel('L1 error of v'); axs[1].set_ylabel('L1 error of p')
    axs[2].set_ylabel('min p over the run (floored 1e-16)')
    axs[2].axhline(5.3e-12, color='#666', lw=1, ls='-.'); axs[2].text(40, 8e-12, 'exact centre p = 5.3e-12', fontsize=8)
    for ax in axs: ax.set_xlabel('N'); ax.set_xticks(NV); ax.set_xticklabels([str(n) for n in NV])
    legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, 'vortex')

# ------------------------------------------------------------------ RD
if 'rd' in R:
    summary['rd'] = {}
    fig, axs = plt.subplots(2, 3, figsize=(15, 6.4), sharex=True)
    shades = {200: '#9ec5f4', 800: '#3987e5', 3200: '#184f95', 25600: '#111111'}
    for j, sc in enumerate(['hlle_gs', 'hlld_gs', 'lhlld_gs']):
        for N in (200, 800, 3200, 25600):
            f = f'{B}/runs/rd/{sc}_N{N}/profile.npz'
            if not os.path.exists(f): continue
            z = np.load(f)
            axs[0, j].plot(z['x'], z['rho'], color=shades[N], lw=1.5, label=f'N = {N}')
            axs[1, j].plot(z['x'], z['eint'], color=shades[N], lw=1.5, label=f'N = {N}')
            m = R['rd'][f'{sc}_N{N}']
            summary['rd'].setdefault(sc, {})[N] = dict(rho_min=m.get('rho_min'), eint_max=m.get('eint_max'),
                                                     eint_min=m.get('eint_min'))
        axs[0, j].set_title(NAME[sc.split('_')[0]]); axs[1, j].set_xlabel('x')
        axs[0, j].set_xlim(0.3, 0.7)
    axs[0, 0].set_ylabel('density'); axs[1, 0].set_ylabel('internal energy e = p/(gamma-1)')
    axs[0, 2].legend(fontsize=8, frameon=False)
    save(fig, 'rd_profiles')

# ------------------------------------------------------------------ rotor / blast
ARMS = SCH + ['hlld_gs_fofc', 'hlle_gs_fofc']
if 'rotor' in R:
    summary['rotor'] = {}
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.6), sharey=True)
    for ax, N in zip(axs, (200, 400, 800)):
        for sc in ARMS:
            m = R['rotor'].get(f'{sc}_N{N}')
            if not m or 'hist' not in m: continue
            t = np.array(m['hist']['time']); e = np.array(m['hist']['emin'])*(2/3)
            ax.plot(t, np.where(e > 0, e, np.nan), label=lab(sc), **dict(sty(sc), marker=None, lw=1.6))
            summary['rotor'].setdefault(sc, {})[N] = dict(pmin=m.get('pmin_all'), emin=m.get('emin_all'),
                counts=m['counts'], completed=m.get('completed'), t_end=m.get('hist_tmax'),
                first_floor_cycle=m.get('first_floor_cycle'), wall=m.get('wall_s'))
        ax.set_yscale('log'); ax.set_title(f'low-beta rotor, N = {N}'); ax.set_xlabel('t')
    axs[0].set_ylabel('min (E - KE - ME)(gamma-1) over mesh\n(gaps = negative)')
    legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, 'rotor_pmin_t')
    # pressure maps N=400
    sel = ['hlld_gs', 'hlld_sm', 'hlld_uct', 'lhlld_gs']
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.9))
    for ax, sc in zip(axs, sel):
        f = f'{B}/runs/rotor/{sc}_N400/final.npz'
        if os.path.exists(f):
            p = np.load(f)['p']
            im = ax.imshow(np.log10(np.maximum(p, 1e-6)), origin='lower', extent=(0, 1, 0, 1),
                           cmap='Blues', vmin=-3, vmax=-0.5)
        ax.set_title(lab(sc), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axs, shrink=0.8, label='log10 p (t = 0.12)')
    save(fig, 'rotor_maps')

if 'blast' in R:
    summary['blast'] = {}
    groups = [('b2e-2', 2, 256, 'beta 0.016 (2-D, N=256)'), ('b2.5e-4', 2, 256, 'beta 2.5e-4 (2-D, N=256)'),
              ('b2.5e-4', 3, 64, 'beta 2.5e-4 (3-D, N=64)'), ('b2.5e-4', 3, 128, 'beta 2.5e-4 (3-D, N=128)'),
              ('b3e-6', 2, 256, 'beta 3e-6 (2-D, N=256)')]
    fig, axs = plt.subplots(1, len(groups), figsize=(18, 4.2), sharey=True)
    for ax, (lb, dim, N, title) in zip(axs, groups):
        ys, vals, fl, cols, labels = [], [], [], [], []
        for i, sc in enumerate(ARMS):
            m = R['blast'].get(f'{lb}_{dim}d_N{N}_{sc}')
            if not m: continue
            hits = m['counts']['efloor'] + m['counts']['tfloor']
            ok = bool(m.get('completed'))
            summary['blast'][f'{lb}_{dim}d_N{N}_{sc}'] = dict(completed=ok, floor_hits=hits,
                dfloor=m['counts']['dfloor'], fofc=m['counts']['fofc'], pmin=m.get('pmin_all'),
                emin=m.get('emin_all'), t_end=m.get('hist_tmax'), tend=m['tend'], wall=m.get('wall_s'),
                rc=m.get('rc'))
            ncyc = m.get('end_cycle')
            if not ncyc:
                ls_ = [l for l in open(f'{B}/runs/blast/{lb}_{dim}d_N{N}_{sc}/run.log') if not l.startswith('#')]
                ncyc = int(ls_[-1].split()[0]) if ls_ else 1
            frac = m['counts']['efloor']/(2.0*ncyc*N**dim)
            summary['blast'][f'{lb}_{dim}d_N{N}_{sc}']['floor_fraction'] = frac
            summary['blast'][f'{lb}_{dim}d_N{N}_{sc}']['cycles'] = ncyc
            ax.barh(i, max(frac, 1e-7), color=COL[sc.split('_')[0]], height=0.6,
                    hatch=None if ok else '///', edgecolor='white', lw=0)
        ax.set_xscale('log'); ax.set_title(title, fontsize=9.5)
        ax.set_xlabel('fraction of cell updates at the p floor'); ax.set_xlim(1e-7, 1)
    axs[0].set_yticks(range(len(ARMS))); axs[0].set_yticklabels([lab(s) for s in ARMS], fontsize=8)
    axs[0].invert_yaxis()
    fig.suptitle('hatched bar = run did not reach t_end (time step collapsed, stopped at the 90-min cap)', fontsize=9, y=1.02)
    save(fig, 'blast_floors')
    # min p vs time for the beta 3e-6 and 2.5e-4 2-D
    fig, ax = plt.subplots(1, 1, figsize=(8, 3.6))
    for sc in ARMS:
        f = f'{B}/runs/blast/b3e-6_2d_N256_{sc}/run.mhd.hst'
        if not os.path.exists(f): continue
        a = np.loadtxt(f, ndmin=2)
        ax.semilogy(a[:, 0]/0.02, a[:, 1], label=lab(sc), **dict(sty(sc), marker=None, lw=1.6))
    ax.set_xlabel('t / t_end'); ax.set_ylabel('time step dt')
    ax.set_title('beta 3e-6 blast (2-D, N = 256): time step', fontsize=9.5)
    ax.legend(fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, 'blast_dt')
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.6))
    for ax, (lb, dim, N, title) in zip(axs, [groups[0], groups[1], groups[4]]):
        for sc in ARMS:
            m = R['blast'].get(f'{lb}_{dim}d_N{N}_{sc}')
            if not m or 'hist' not in m: continue
            t = np.array(m['hist']['time'])/m['tend']; e = np.array(m['hist']['pmin'])/0.1
            ax.plot(t, np.maximum(e, 1e-10), label=lab(sc), **dict(sty(sc), marker=None, lw=1.6))
        ax.set_yscale('log'); ax.set_title(title, fontsize=9.5); ax.set_xlabel('t / t_end')
        ax.set_ylim(5e-11, 2)
    axs[0].set_ylabel('min p / ambient p\n(1e-10 = at the floor)')
    legend_all(fig, axs[-1], fontsize=7.5, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    save(fig, 'blast_emin_t')

json.dump(dict(summary=summary, runs=R), open(f'{B}/results.json', 'w'), default=float)
print('ok', {k: len(v) for k, v in R.items()})
