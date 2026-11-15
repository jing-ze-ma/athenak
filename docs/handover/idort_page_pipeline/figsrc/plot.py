"""IDORT results page figures from the cached FIN3 final data (fin3-1012 281dd686; extract.py) and the evaluation tables.
FIN3 page (worker PG, 2026-10-13): every label FIN -> FIN3 and PRODd -> FIN2 (the FIN2 control F2O) is rewritten in save();
old version in figsrc.bak_fin3.
usage: plot.py [name ...]   (default: all).  Writes fig/<file>.png and merges fig/figures.json.
No simulation is run; every number in a caption is computed here from the cached run data or read from
the suite's own evaluation output files (named in the caption source)."""
import json
import math
import os
import pickle
import re
import sys
import traceback

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

P = '/raven/ptmp/jinma/idort1020/page'
FIG = P + '/fig'
CACHE = P + '/cache'
SF = '/raven/ptmp/jinma/suiteF_1012'
BIN_D = os.environ.get('BIN_D', '563ab474')   # idort-1020: rerun-row figures get the rerun bin; fin3-1012 FIN3 final (= 1fdbe56e/914e5433 binaries + implicit_predictor_frozen false)
BIN_C = '281dd686 v10'   # comparison arm FIN2 control (F2O: suiteF runs_js / same binary, fin3 features at old values)
F3 = '/raven/ptmp/jinma/fin3_1012'
F3SP = P + '/cache/p19ana'   # IN's evaluation of the final scheme (fin3_1012/ana/out_281dd686, eval_arm.sh), arm names normalised by norm_fin3ana.sh
BIN_O = '0197f020'

# ---- one style for every figure: white surface, dark ink, recessive grid
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#6250d6'   # FIN, comparison arm, third, fourth
INK, INK2, GRID = '#1a1a1a', '#52514e', '#e4e4e2'
REF = dict(color=INK, ls='--', lw=1.3)
plt.rcParams.update({
    'font.size': 14, 'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 12,
    'xtick.labelsize': 13, 'ytick.labelsize': 13, 'figure.titlesize': 15,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.edgecolor': '#9a9a96', 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7,
    'lines.linewidth': 1.8, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.facecolor': 'white', 'legend.frameon': True, 'legend.framealpha': 0.92,
    'legend.edgecolor': '#d0d0cc', 'image.cmap': 'magma', 'axes.titleweight': 'normal'})
LEVC = ['#9cc3ef', '#4f93e0', '#1f5fb0', '#0d3570']   # sequential blues: coarse -> fine
DPI = 150        # figsize width 12 in -> 1200 px

FIGS = {}


def load(name):
    with open(os.path.join(CACHE, name + '.pkl'), 'rb') as fh:
        return pickle.load(fh)


def _levelticks(fig):
    """log x axes of ladders: major ticks at the levels (plain numbers), no minor ticks"""
    from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
    for ax in fig.axes:
        if ax.get_xscale() != 'log':
            continue
        xs = set()
        for ln in ax.get_lines():
            if ln.get_marker() not in (None, 'None', '', ' '):
                xs.update(float(v) for v in np.atleast_1d(ln.get_xdata()))
        if 1 < len(xs) <= 10:
            ax.xaxis.set_major_locator(FixedLocator(sorted(xs)))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: '%g' % v))
            ax.xaxis.set_minor_locator(NullLocator())


def relabel(s):
    s = re.sub(r'\bFIN\b', 'PROD', s)
    s = re.sub(r'\bPRODd\b', 'v10 FIN3', s)
    return s.replace('final_1012', 'fin3-1012').replace('method-1012', 'fin3-1012')


def save(fig, fname, key, caption, status='ok', extra=None):
    import matplotlib.text as mtext
    for t in fig.findobj(mtext.Text):
        s0 = t.get_text()
        if s0:
            s1 = relabel(s0)
            if s1 != s0:
                t.set_text(s1)
    caption = relabel(caption)
    _levelticks(fig)
    fig.savefig(os.path.join(FIG, fname), dpi=DPI)
    plt.close(fig)
    d = {'file': fname, 'caption': caption, 'status': status}
    if extra:
        d.update(extra)
    FIGS[key] = d
    print('wrote', fname)


def nodata(key, reason, extra=None):
    d = {'file': None, 'caption': reason, 'status': 'no data'}
    if extra:
        d.update(extra)
    FIGS[key] = d
    print('no data', key)


def suptitle(fig, s, sub=None):
    fig.suptitle(s, x=0.01, ha='left', fontsize=13, color=INK)
    if sub:
        fig.text(0.01, 0.935, sub, ha='left', va='top', fontsize=10, color=INK2)


def slopes(ax, x0, y0, xs, orders=(1, 2), lab=True):
    """reference slopes through (x0, y0) for error vs N (decreasing)"""
    xs = np.array(xs, float)
    for o, ls in zip(orders, (':', '--')):
        ax.plot(xs, y0*(xs/x0)**(-o), color='#8a8a86', ls=ls, lw=1.1, label=('slope -%d' % o) if lab else None)


def noimgrid(*axs):
    for a in axs:
        a.grid(False)


# ============================================================ 01-04 crossing beams (IDORT Fig. 2)
def pool(x1, x2, z, f):
    """f x f max-pool of z (nx2, nx1) for display of 1-cell-wide beams on large grids"""
    if f == 1:
        return x1, x2, z
    n2, n1 = z.shape
    z = z[:n2//f*f, :n1//f*f].reshape(n2//f, f, n1//f, f).max(axis=(1, 3))
    return x1[:n1//f*f].reshape(-1, f).mean(1), x2[:n2//f*f].reshape(-1, f).mean(1), z


def cb_map(ax, r, z, norm, outline=True):
    f = 4 if len(r['x1']) > 512 else 1
    x1, x2, zz = pool(r['x1'], r['x2'], z, f)
    im = ax.pcolormesh(x2, x1, np.maximum(zz, 1e-8).T, norm=norm, cmap='magma', shading='nearest',
                       rasterized=True)
    if f > 1:
        ax.text(0.03, 0.01, '4x4 max-pooled', transform=ax.transAxes, color='white', fontsize=8)
    if outline:
        ax.contour(r['x2'], r['x1'], r['ex'].T, levels=[0.5*max(b['e'] for b in r['beams'])],
                   colors=['#7fe3ff'], linewidths=0.8)
    ax.set_aspect('equal')
    noimgrid(ax)
    return im


def cb_rowratio(r):
    rat = r['rowE']/np.maximum(r['rowex'], 1e-300)
    x1 = r['x1']
    return x1, rat


def cb_numbers(r):
    x1, rat = cb_rowratio(r)
    mid = (x1 > x1.min() + 0.1*(x1.max() - x1.min())) & (x1 < x1.max() - 0.1*(x1.max() - x1.min()))
    eb = max(b['e'] for b in r['beams'])
    i = int(np.argmin(np.abs(x1 - 0.65)))
    return float(np.median(rat[mid])), float(r['E'][:, i].max()/eb), float(np.min(r['E']))


def f01():
    d = load('cb')
    r = d['cb/FIN_r4']
    eb = max(b['e'] for b in r['beams'])
    norm = LogNorm(0.03*eb, eb)
    fig = plt.figure(figsize=(12, 7.2))
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 0.08, 0.32, 2.7], wspace=0.18, left=0.06, right=0.98, top=0.86, bottom=0.09)
    for k, (z, t) in enumerate(((r['ex'], 'exact (cell avg)'), (r['E'], 'FIN E'), (r['J'], 'J_sc (formal sol.)'))):
        a = fig.add_subplot(gs[0, k])
        im = cb_map(a, r, z, norm, outline=(k == 1))
        a.set_title(t, fontsize=10.5)
        a.set_xlabel('x')
        if k == 0:
            a.set_ylabel('y')
        else:
            a.set_yticklabels([])
    cax = fig.add_subplot(gs[0, 3])
    fig.colorbar(im, cax=cax).ax.set_title('E/E_b', fontsize=9.5)
    a = fig.add_subplot(gs[0, 5])
    lab = {'cb/FIN_r4': 'ring 4, c dt/dx 256'}
    x1, rat = cb_rowratio(r)
    a.plot(x1, rat, color=C1, label='FIN ' + lab['cb/FIN_r4'])
    if r['J'] is not None:
        a.plot(x1, r['J'].sum(axis=0)/np.maximum(r['rowex'], 1e-300), color=C3, ls=':', lw=2, label='J_sc')
    a.axhline(1.0, **REF, label='exact')
    a.axhline(0.9, color=C2, lw=1, ls='-.', label='criterion 0.9')
    a.set_ylim(0.6, 1.25)
    a.set_xlabel('y (height)')
    a.set_ylabel('beam energy per row / exact')
    a.legend(loc='lower left')
    row, pk, fl = cb_numbers(r)
    st = r['st']
    suptitle(fig, 'Crossing beams (IDORT 4.1.1, Fig. 2; J21 4.2): FIN, ring 4, c dt/dx = 256',
             'cyan line = analytic beam outline (half of E_beam); binary %s; Picard %.1f, NON-CONVERGED %g' % (
                 BIN_D, st.get('picard', float('nan')), st.get('nonconv', float('nan'))))
    save(fig, '01_cb.png', '01_cb',
         'E maps (exact / FIN / J_sc) and per-row beam energy of the two crossing beams, FIN on %s (ring 4, '
         '64x256, c dt/dx 256): median row energy %.3f of exact, peak %.3f E_b at y = 0.65, Picard %.1f with %g '
         'non-converged steps.%s (Raven: Picard 33.3, 6 NC, floor 4e-19).' % (BIN_D, row, pk, st.get('picard', float('nan')), st.get('nonconv', float('nan')), ' Plotted runs: Viper 1fdbe56e CPU build (bitwise to 281dd686 there); the Raven 281dd686 rerun reproduces the metrics to print precision, only Picard / NC counts and round-off floors differ'))
    nodata('01_cb_sp', 'Spherical counterpart T8 (vet_gd_beam injected beams) is parked and was not run in suiteF '
           '(see 27_pinc_sp for the production thermal channel beam).')


def f02():
    d = load('cb')
    runs = [('cb/FIN_r12', 'ring 12, c dt/dx 256'), ('cb/FIN_r4c0.5', 'ring 4, c dt/dx 0.5'),
            ('cb/FIN_r4c32', 'ring 4, c dt/dx 32'), ('cbh/FIN_r4', 'ring 4, cbh')]
    eb = max(b['e'] for b in d[runs[0][0]]['beams'])
    norm = LogNorm(0.03*eb, eb)
    fig = plt.figure(figsize=(12, 7.2))
    gs = fig.add_gridspec(1, 7, width_ratios=[1, 1, 1, 1, 0.08, 0.32, 2.3], wspace=0.18, left=0.06, right=0.98, top=0.86, bottom=0.09)
    cols = [C1, C2, C3, C4]
    a2 = fig.add_subplot(gs[0, 6])
    nums = []
    for k, (run, lab) in enumerate(runs):
        r = d[run]
        a = fig.add_subplot(gs[0, k])
        im = cb_map(a, r, r['E'], norm)
        a.set_title(lab.replace(', ', '\n'), fontsize=10)
        a.set_xlabel('x')
        if k:
            a.set_yticklabels([])
        else:
            a.set_ylabel('y')
        x1, rat = cb_rowratio(r)
        a2.plot(x1, rat, color=cols[k], lw=1.6, ls=['-', '--', '-.', ':'][k], label=lab)
        row, pk, _ = cb_numbers(r)
        nums.append('%s: row %.3f, peak %.3f, Picard %.1f' % (lab, row, pk, r['st'].get('picard', float('nan'))))
    fig.colorbar(im, cax=fig.add_subplot(gs[0, 4])).ax.set_title('E/E_b', fontsize=9.5)
    a2.axhline(1.0, **REF, label='exact')
    a2.set_ylim(0.6, 1.25)
    a2.set_xlabel('y (height)')
    a2.set_ylabel('beam energy per row / exact')
    a2.legend(loc='lower left', fontsize=8.5)
    suptitle(fig, 'Crossing beams, variants (IDORT 4.1.1): FIN', 'binary %s; cyan = analytic beam outline' % BIN_D)
    save(fig, '02_cbvar.png', '02_cbvar', 'Crossing-beam E maps and row energies for the ring-12 / c dt/dx 0.5 / 32 / '
         'cbh variants, FIN on %s: %s.%s (cbh was run on Raven with the 1fdbe56e build at 10 threads; Raven 281dd686 rerun: Picard 3.13 / 2.94 / 3.85, cbh 4.0).' % (BIN_D, '; '.join(nums), ' Plotted runs: Viper 1fdbe56e CPU build (bitwise to 281dd686 there); the Raven 281dd686 rerun reproduces the metrics to print precision, only Picard / NC counts and round-off floors differ'))


def f03():
    d = load('cb')
    r = d['cb/FIN_r4n4x']
    eb = max(b['e'] for b in r['beams'])
    norm = LogNorm(0.03*eb, eb)
    fig = plt.figure(figsize=(12, 7.2))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.08, 0.32, 2.8], wspace=0.2, left=0.06, right=0.98, top=0.86, bottom=0.09)
    for k, (z, t) in enumerate(((r['ex'], 'exact'), (r['E'], 'FIN E'))):
        a = fig.add_subplot(gs[0, k])
        im = cb_map(a, r, z, norm, outline=(k == 1))
        a.set_title(t)
        a.set_xlabel('x')
        if k:
            a.set_yticklabels([])
        else:
            a.set_ylabel('y')
    fig.colorbar(im, cax=fig.add_subplot(gs[0, 2])).ax.set_title('E/E_b', fontsize=9.5)
    a = fig.add_subplot(gs[0, 4])
    x1, rat = cb_rowratio(r)
    a.plot(x1, rat, color=C1, label='FIN')
    a.axhline(1.0, **REF, label='exact')
    a.set_ylim(0.6, 1.25)
    a.set_xlabel('y (height)')
    a.set_ylabel('beam energy per row / exact')
    a.legend(loc='lower left')
    row, pk, _ = cb_numbers(r)
    nx = '%dx%d' % (len(r['x2']), len(r['x1']))
    suptitle(fig, 'Crossing beams at 4x resolution (%s cells), FIN' % nx, 'binary %s; ring 4' % BIN_D)
    save(fig, '03_cbn4x.png', '03_cbn4x', 'Crossing beams at %s cells (r4n4x), FIN on %s: median row energy %.3f '
         'of exact, peak %.3f E_b at y = 0.65.' % (nx, BIN_D, row, pk))


def f04():
    d = load('cbg')
    r = d['FIN']
    st = r['st']
    ex = r['ex']
    eb = ex.max()
    norm = LogNorm(0.03*eb, eb)
    fig = plt.figure(figsize=(12, 7.2))
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 0.08, 0.32, 2.7], wspace=0.18, left=0.06, right=0.98, top=0.84, bottom=0.09)
    for k, (z, t) in enumerate(((ex, 'exact'), (r['E'], 'FIN E'), (r['J'], 'J_sc (FIN)'))):
        a = fig.add_subplot(gs[0, k])
        im = a.pcolormesh(r['x2'], r['x1'], np.maximum(z, 1e-8).T, norm=norm, cmap='magma', shading='nearest', rasterized=True)
        a.set_aspect('equal')
        noimgrid(a)
        a.set_title(t, fontsize=10)
        a.set_xlabel('x')
        if k:
            a.set_yticklabels([])
        else:
            a.set_ylabel('y')
    fig.colorbar(im, cax=fig.add_subplot(gs[0, 3])).ax.set_title('E/E_b', fontsize=9.5)
    a = fig.add_subplot(gs[0, 5])
    a.plot(r['x1'], r['rowE']/np.maximum(r['rowex'], 1e-300), color=C1, label='FIN (%s)' % BIN_D)
    a.plot(r['x1'], r['J'].sum(axis=0)/np.maximum(r['rowex'], 1e-300), color=C3, ls=':', lw=2, label='J_sc')
    a.axhline(1.0, **REF, label='exact')
    a.set_ylim(0.5, 1.3)
    a.set_xlabel('y (height)')
    a.set_ylabel('beam energy per row / exact')
    a.legend(loc='lower left', fontsize=9)
    fat = d['PRODd'].get('fatal', '')
    suptitle(fig, 'Crossing beams, Gauss 8x32 angles (cbg), c dt/dx 256: FIN completes',
             'binary %s; Picard %.1f, NON-CONVERGED %g; FIN2 control: Picard %.1f, NON-CONVERGED %g'
             % (BIN_D, st.get('picard', float('nan')), st.get('nonconv', float('nan')), d['PRODd'].get('picard', float('nan')),
                d['PRODd'].get('nonconv', float('nan'))))
    rat = r['rowE']/np.maximum(r['rowex'], 1e-300)
    i = int(np.argmin(np.abs(r['x1'] + 1.5)))
    j = int(np.argmin(np.abs(r['x1'] - 1.9)))
    save(fig, '04_cbg.png', '04_cbg', 'cbg (Gauss 8x32, c dt/dx 256), FIN on %s: completes, row energy %.2f / %.2f of exact at '
         'y = -1.5 / 1.9, Picard %.1f with %g non-converged steps (FIN2 control: Picard %.1f, %g non-converged; the earlier '
         'candidate 4253af30 was FATAL).' % (BIN_D, rat[i], rat[j], st.get('picard', float('nan')), st.get('nonconv', float('nan')),
                                     d['PRODd'].get('picard', float('nan')), d['PRODd'].get('nonconv', float('nan'))))


# ============================================================ 05 homogeneous sphere (IDORT Fig. 3)
def f05():
    d = load('hcyl')
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    a = ax[0]
    RUNS5 = [x for x in (('FIN_g8x32', C1, '-', '256^2'), ('FIN_n512', C3, '--', '512^2')) if x[0] in d]
    big = 'FIN_n512' if 'FIN_n512' in d else 'FIN_g8x32'
    for run, col, ls, lab in RUNS5:
        f = d[run]['prof']
        a.loglog(f['x'], f['E'], color=col, ls=ls, label='E  FIN ' + lab)
        a.loglog(f['x'], np.abs(f['F']), color=col, ls=ls, lw=1.2, alpha=0.7)
        if 'fxx' in d[run] and run == big:
            a.loglog(f['x'], d[run]['fxx'], color=C2, ls=ls, label='P_xx/E  FIN ' + lab)
    ex = d['exact512']
    a.loglog(ex['x'], ex['E'], ':', color=INK, lw=2, label='exact')
    a.loglog(ex['x'], ex['F'], ':', color=INK, lw=2)
    a.loglog(ex['x'], ex['Pxx']/ex['E'], ':', color=INK, lw=2)
    a.set_ylim(5e-3, 1.5)
    a.set_xlim(0.03, 4)
    a.set_ylabel('E/(aT^4), F_x/(c aT^4), P_xx/E')
    a.axvline(1.0, color='#aaaaaa', lw=0.8)
    a.text(0.05, 0.75, 'E/(aT^4)', color=C1, fontsize=10)
    a.text(0.05, 0.36, 'P_xx/E', color=C2, fontsize=10)
    a.text(0.4, 0.012, 'F_x/(c aT^4)', color=C1, fontsize=10)
    a.set_xlabel('r (y = 0 cut)')
    a.legend(loc='lower left', fontsize=8.5)
    a.set_title('profiles vs exact (dotted)')
    r = d[big]
    im = ax[1].pcolormesh(r['x1'], r['x2'], np.maximum(r['E2'], 1e-4), norm=LogNorm(1e-2, 1), shading='nearest', rasterized=True)
    ax[1].add_patch(plt.Circle((0, 0), 1.0, fill=False, ec='#7fe3ff', lw=1))
    ax[1].set_aspect('equal')
    noimgrid(ax[1])
    ax[1].set_title('FIN E, %s (cyan: R0)' % ('512^2' if big == 'FIN_n512' else '256^2'))
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('y')
    fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.02)
    for run, col, ls, lab in RUNS5:
        f = d[run]['prof']
        ax[2].plot(f['x'], f['E']/f['Eex'], color=col, ls=ls, label='E ' + lab)
        ax[2].plot(f['x'], f['F']/np.maximum(f['Fex'], 1e-30), color=col, ls=ls, lw=1.1, alpha=0.6, label='F ' + lab)
    ax[2].axhline(1, **REF)
    ax[2].axvline(1.0, color='#aaaaaa', lw=0.8)
    ax[2].set_ylim(0.8, 1.15)
    ax[2].set_xlabel('r')
    ax[2].set_ylabel('code / exact')
    ax[2].legend(fontsize=8.5)
    m1 = d['FIN_g8x32']['metrics']
    m2 = d[big]['metrics']
    suptitle(fig, 'Homogeneous sphere (IDORT 4.1.2, Fig. 3), Cartesian analogue: homogeneous cylinder, FIN',
             'rho kappa = 10, R0 = 1, held T = 1, vacuum outside; Gauss 8x32 angles; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, '05_hsphere.png', '05_hsphere', 'Homogeneous cylinder (Cartesian analogue of the IDORT sphere), FIN on %s: '
         'E, F and P_xx/E vs the exact ray solution; E L1 outside %.1f %%, F %.1f %% (256^2, Gauss 8x32), f_xx(0.5) %s; %s.' % (
             BIN_D, 100*m1['E_L1_out'], 100*m1['F_L1_out'], m2.get('fxx(0.5) code/ex', ''),
             ('512^2: E %.1f %%, F %.1f %%' % (100*d['FIN_n512']['metrics']['E_L1_out'], 100*d['FIN_n512']['metrics']['F_L1_out']))
             if 'FIN_n512' in d else '512^2 run: pending (final rerun running)'))


def f05sp():
    d = load('sp_core')
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    a = ax[0]
    o = d['core_L16_u_256']
    cf = d['core_fine']
    a.loglog(o['r'], o['E'], color=C1, label='E  (FIN, 16x16, N 256)')
    a.loglog(o['r'], np.abs(o['F']), color=C3, label='F_r')
    a.loglog(cf['r'], cf['J'], ':', color=INK, lw=2, label='exact')
    a.loglog(cf['r'], np.abs(cf['H']), ':', color=INK, lw=2)
    a.set_ylim(2e-3, 2)
    a.axvline(1, color='#aaaaaa', lw=0.8)
    a.set_xlabel('r')
    a.set_ylabel('E/(aT0^4), F_r/(c aT0^4)')
    a.legend(loc='lower left', fontsize=9)
    a.set_title('IDORT sphere (T6b) on the sp wedge')
    for k, n in enumerate((64, 128, 256)):
        o = d['core_L16_u_%d' % n]
        ax[1].plot(o['r'], 100*(o['E']/o['J'] - 1), color=LEVC[k + 1], label='N_r %d' % n)
    ax[1].axhline(0, **REF)
    ax[1].axvline(1, color='#aaaaaa', lw=0.8)
    ax[1].set_ylim(-15, 15)
    ax[1].set_xscale('log')
    ax[1].set_xlabel('r')
    ax[1].set_ylabel('E / J_exact - 1  [%]')
    ax[1].legend(fontsize=9)
    ax[1].set_title('T6b, 16x16 lateral cells')
    for k, n in enumerate((64, 128, 256)):
        o = d['thin_L16_u_%d' % n]
        r = o['r']
        exE = 2.0*(1.0 - np.sqrt(1.0 - 1.0/r**2))
        ax[2].plot(r, 100*(o['E']/exE - 1), color=LEVC[k + 1], label='E, N_r %d' % n)
        ax[2].plot(r, 100*(o['F']*r*r - 1), color=LEVC[k + 1], ls='--', lw=1.2, label='F r^2, N_r %d' % n)
    ax[2].axhline(0, **REF)
    ax[2].set_xscale('log')
    ax[2].set_ylim(-3, 3)
    ax[2].set_xlabel('r')
    ax[2].set_ylabel('code / exact - 1  [%]')
    ax[2].legend(fontsize=8, ncol=2)
    ax[2].set_title('T6a thin black-body sphere, 16x16')
    o = d['core_L16_u_256']
    out = (o['r'] > 3)
    ext = 100*np.mean(o['E'][out]/o['J'][out] - 1)
    inn = 100*np.mean((o['E']/o['J'] - 1)[o['r'] < 0.9])
    t6a = d['thin_L16_u_256']
    rr = t6a['r']
    t6aE = 100*np.mean(np.abs(t6a['E']/(2.0*(1.0 - np.sqrt(1.0 - 1.0/rr**2))) - 1))
    suptitle(fig, 'Homogeneous sphere, spherical counterpart (vet_col + vet_gd wedge), FIN (HFfZ)',
             'binary %s; uniform r; left panel N_r 256' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, '05_hsphere_sp.png', '05_hsphere_sp', 'IDORT sphere (T6b) and thin black-body sphere (T6a) on the sp wedge, '
         'FIN (HFfZ) on %s: T6b E/J-1 %.2f %% inside r < 0.9 and %+.1f %% at r > 3 (16x16, N_r 256; vet_seg jump + hcface tau: converging, about first order); '
         'T6a mean |E/E_ex-1| %.2f %% at N_r 256.' % (BIN_D, inn, ext, t6aE))


# ============================================================ 06 non-LTE (IDORT Fig. 4)
def f06():
    d = load('nlte')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    cols = plt.cm.viridis(np.linspace(0, 0.9, 8))
    for l in range(1, 9):
        m = d[l]
        ax[0].loglog(m['tau'], m['ex'], color=cols[l - 1], lw=3, alpha=0.3)
        ax[0].loglog(m['tau'], m['E'], color=cols[l - 1], lw=1.4, label='eps = 1e-%d' % l)
        ax[1].semilogx(m['tau'], 100*(m['E']/m['ex'] - 1), color=cols[l - 1], lw=1.3, label='1e-%d: L1 %.2f %%' % (l, 100*m['L1']))
    ax[0].set_xlim(1e-5, 5e5)
    ax[0].set_ylim(1e-5, 2)
    ax[0].set_xlabel('tau')
    ax[0].set_ylabel('E / (a T^4)')
    ax[0].legend(fontsize=8.5, ncol=2, loc='lower right')
    ax[0].set_title('FIN (thin) vs Eddington analytic (thick, faint)')
    ax[1].set_xlim(1e-5, 5e5)
    ax[1].set_ylim(-15, 15)
    ax[1].axhline(0, **REF)
    ax[1].set_xlabel('tau')
    ax[1].set_ylabel('E / E_analytic - 1  [%]')
    ax[1].legend(fontsize=8, ncol=2)
    suptitle(fig, 'Non-LTE scattering atmosphere (IDORT 4.1.3, Fig. 4; J21 4.5), FIN',
             '1280 cells, Gauss 8x32 rays, steady; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, '06_nlte.png', '06_nlte', 'Non-LTE atmosphere E(tau) for eps = 1e-1 ... 1e-8, FIN on %s, against the '
         'Eddington analytic solution: L1 %.2f %% (eps 1e-1) ... %.2f %% (1e-8), E_top/ex %.4f-%.4f.' % (
             BIN_D, 100*d[1]['L1'], 100*d[8]['L1'], min(d[l]['top'] for l in d), max(d[l]['top'] for l in d)))


def f06sp():
    d = load('sp_nlte')
    n = 320
    x = (np.arange(n) + 0.5)/n*20.0
    tau = 1e-3*(np.exp(20.0 - x) - 1.0)
    eps = 1e-2
    edd = 1.0 - np.exp(-np.sqrt(3*eps)*tau)/(1.0 + np.sqrt(eps))
    R = d['ref']
    jr = np.interp(x, R['bz'], R['bJ'])
    es, ec = d['nlte2_r1000_tb']['E'], d['nlte2_tb_cart']['E']
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    ax[0].loglog(tau, es, color=C1, label='sp wedge r0/w = 1000 (FIN, vet_gd)')
    ax[0].loglog(tau, ec, color=C2, ls='--', label='Cartesian vet_sc (same keys)')
    ax[0].loglog(tau, jr, ':', color=INK, lw=2, label='exact (1-D transport)')
    ax[0].loglog(tau, edd, color='#8a8a86', lw=1, label='Eddington analytic')
    ax[0].set_xlim(1e-4, 1e4)
    ax[0].set_ylim(1e-2, 2)
    ax[0].set_xlabel('tau')
    ax[0].set_ylabel('E / (a T^4)')
    ax[0].legend(fontsize=9)
    ax[1].semilogx(tau, 100*(es/jr - 1), color=C1, label='sp / exact - 1')
    ax[1].semilogx(tau, 100*(ec/jr - 1), color=C2, ls='--', label='Cartesian / exact - 1')
    ax[1].axhline(0, **REF)
    ax[1].set_xlim(1e-4, 1e4)
    ax[1].set_ylim(-3, 3)
    ax[1].set_xlabel('tau')
    ax[1].set_ylabel('[%]')
    ax[1].legend(fontsize=9)
    bands = [(0, 0.1), (0.1, 1), (1, 10)]
    bb = [100*np.mean((es/jr - 1)[(tau >= a) & (tau < b)]) for a, b in bands]
    suptitle(fig, 'Non-LTE atmosphere, spherical counterpart (thin shell r0/width = 1000), eps = 1e-2',
             'FIN sp arm HFfSgZ (vet_col_source gas) vs Cartesian HCFfZ, N 320; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    save(fig, '06_nlte_sp.png', '06_nlte_sp', 'Non-LTE atmosphere (eps 1e-2) on a thin sp shell (r0/w 1000, vet_gd) and '
         'in Cartesian, FIN on %s, vs the exact 1-D transport solution: sp mean error %+.2f / %+.2f / %+.2f %% in tau '
         '< 0.1 / 0.1-1 / 1-10.' % (BIN_D, bb[0], bb[1], bb[2]))


# ============================================================ 07/08 shadow (IDORT Fig. 5)
def shd_fig(run, key, fname, title, extra_runs=()):
    d = load('shd')
    r = d[run]
    x1, x2 = r['x1'], r['x2']
    lit = np.percentile(r['ex'], 99)
    panels = [(r['ex'], 'exact'), (r['E'], 'FIN E'), (r['J'], 'J_sc')] + [(d[q]['E'], lab) for q, lab in extra_runs]
    fig = plt.figure(figsize=(12, 7.6))
    gs = fig.add_gridspec(2, len(panels), height_ratios=[1, 0.9], hspace=0.35, wspace=0.08, top=0.86, bottom=0.08,
                          left=0.06, right=0.9)
    cuts = []
    for k, (z, t) in enumerate(panels):
        a = fig.add_subplot(gs[0, k])
        im = a.pcolormesh(x1, x2, z, vmin=0, vmax=1.05*r['ex'].max(), cmap='viridis', shading='nearest', rasterized=True)
        a.contour(x1, x2, r['ex'], levels=[0.25*lit, 0.75*lit], colors=['white'], linewidths=0.7, linestyles=['--', '-'])
        a.axvline(0.4, color='#ffd166', lw=0.8, ls=':')
        a.set_aspect('equal')
        noimgrid(a)
        a.set_title(t, fontsize=10.5)
        a.set_xlabel('x')
        if k:
            a.set_yticklabels([])
        else:
            a.set_ylabel('y')
        cuts.append((t, z[:, np.argmin(np.abs(x1 - 0.4))]))
    cax = fig.add_axes([0.915, 0.5, 0.012, 0.36])
    fig.colorbar(im, cax=cax).set_label('E')
    a = fig.add_subplot(gs[1, :])
    sty = [REF, dict(color=C1), dict(color=C3, ls=':', lw=2), dict(color=C2, ls='-.'), dict(color=C4, ls='--')]
    for (t, c), s in zip(cuts, sty):
        a.plot(x2, c, label=t, **s)
    a.set_xlabel('y  (cut at x = 0.4, dotted line above)')
    a.set_ylabel('E')
    a.legend(ncol=5, fontsize=9, loc='upper center')
    i = int(np.argmin(np.abs(x1 - 0.4)))
    ex = r['ex'][:, i]
    l1 = np.sum(np.abs(r['E'][:, i] - ex))/np.sum(ex)
    st = r['st']
    suptitle(fig, title, 'white contours = exact shadow/beam edges (25 %% / 75 %% of lit); binary %s; Picard %.1f' % (
        BIN_D, st.get('picard', float('nan'))))
    return fig, l1


def f07():
    fig, l1 = shd_fig('shd/FIN_r12g2', '07_shd', '07_shd.png',
                      'Static shadow (IDORT 4.1.4, Fig. 5; J21 4.3), FIN, ring 12x2', extra_runs=(('shd/FIN_r36g2', 'FIN ring 36x2'),))
    um = _cart_metric('shd/FIN_r12g2', 'umbra.15')
    save(fig, '07_shd.png', '07_shd', 'Static shadow E maps (exact, FIN ring 12x2, J_sc, ring 36x2) and a cut at '
         'x = 0.4, FIN on %s: suite L1 %s and umbra %s (fin3_1012 ana/out_281dd686 cart.txt); L1 of the plotted x = 0.4 cut %.3f.' % (
             BIN_D, _cart_metric('shd/FIN_r12g2', 'L1'), um, l1))


def f08():
    fig, l1 = shd_fig('shdg/FIN_g8x32', '08_shdg', '08_shdg.png', 'Static shadow, Gauss 8x32 angles (IDORT 4.1.4), FIN')
    um = _cart_metric('shdg/FIN_g8x32', 'umbra.15')
    save(fig, '08_shdg.png', '08_shdg', 'Static shadow with the production Gauss 8x32 angle set, FIN on %s: suite L1 %s '
         'and umbra %s (fin3_1012 ana/out_281dd686 cart.txt); L1 of the plotted x = 0.4 cut %.3f.' % (
             BIN_D, _cart_metric('shdg/FIN_g8x32', 'L1'), um, l1))


def _cart_metric(run, col):
    """a column of the shadow table in suiteF ana/out/cart_fin.txt (row 'E')"""
    hdr = None
    for ln in open(F3SP + '/cart_fin.txt'):
        p = ln.split()
        if p and p[0] == 'run' and 'umbra.15' in p:
            hdr = p
        elif hdr and p and p[0].endswith(run) and len(p) > 1 and p[1] == 'E':
            return p[hdr.index(col)]
    return '?'


def f07sp():
    """FIN3 page: FIN3 (twin_state trans, vet_rotate_avg 4) with 120 directions and the same deck with 432
    (reruns/fin3/sp_shd_FIN3_432); clump interior masked; metrics on the whole wedge (t7eval.py)"""
    d = load('sp_shd')
    tv = json.load(open(os.path.join(CACHE, 't7eval.json')))
    rc, tc, rf = d['rc'], d['tc'], d['rf']
    lat = tc - math.pi/2
    th_f = np.linspace(d['th'][0], d['th'][1], len(tc) + 1) - math.pi/2
    R, TH = np.meshgrid(rf, th_f)
    X, Y = R*np.cos(TH), -R*np.sin(TH)
    fr = np.maximum(d['free'], 1e-30)
    ins = d['ins']
    fig = plt.figure(figsize=(12, 9.6))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.05], height_ratios=[1, 0.95], wspace=0.22, hspace=0.3,
                          left=0.07, right=0.94, top=0.88, bottom=0.07)
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#bdbdbd')
    has4 = d.get('E432') is not None
    for k, (z, t) in enumerate(((d['E']/fr, 'FIN, 120 directions'), ((d['E432']/fr) if has4 else None, 'FIN, 432 directions'),
                                (d['ex']/fr, 'exact'))):
        a = fig.add_subplot(gs[0, k])
        if z is None:
            a.text(0.5, 0.5, 'pending\n(432-dir FIN3 rerun running)', ha='center', va='center', transform=a.transAxes, color=INK2)
            a.set_title(t, fontsize=10.5)
            a.set_xticks([])
            a.set_yticks([])
            noimgrid(a)
            continue
        im = a.pcolormesh(X, Y, np.ma.masked_where(ins, z), vmin=0, vmax=1.5, cmap=cmap, shading='flat', rasterized=True)
        a.contour(rc[None, :]*np.cos(lat[:, None]), -rc[None, :]*np.sin(lat[:, None]), d['ex']/fr, levels=[0.02, 0.98],
                  colors=['white'], linewidths=0.8, linestyles=['-', '--'])
        cr, ct, cp, ca = d['cen']
        a.add_patch(plt.Circle((cr*math.cos(ct - math.pi/2), -cr*math.sin(ct - math.pi/2)), ca, fill=False, ec='#ffd166', lw=1))
        a.set_aspect('equal')
        noimgrid(a)
        a.set_title(t + ': E / E_free', fontsize=10.5)
        a.set_xlabel('r cos(lat)')
        if k == 0:
            a.set_ylabel('r sin(lat)')
    fig.colorbar(im, cax=fig.add_subplot(gs[0, 3]))
    for k, rr in enumerate((4.0, 5.0, 5.8)):
        a = fig.add_subplot(gs[1, k])
        i = int(np.argmin(np.abs(rc - rr)))
        a.plot(np.degrees(lat), d['E'][:, i]/fr[:, i], color=C1, label='FIN 120 dirs')
        if has4:
            a.plot(np.degrees(lat), d['E432'][:, i]/fr[:, i], color=C2, label='FIN 432 dirs')
        a.plot(np.degrees(lat), d['ex'][:, i]/fr[:, i], **REF, label='exact')
        a.set_title('r = %.2f' % rc[i], fontsize=10.5)
        a.set_xlabel('latitude [deg]')
        if k == 0:
            a.set_ylabel('E / E_free')
            a.legend(fontsize=8.5)
    f1, f4 = tv['FIN120'], tv.get('FIN432', {'lit_L1': {'-': float('nan')}, 'umbra': float('nan')})
    f2 = tv.get('F2O120', {})
    st = d['st']
    L = lambda o: ' / '.join('%.3f' % v for v in o['lit_L1'].values()) if not ('-' in o['lit_L1']) else '(432-dir rerun pending)'
    suptitle(fig, 'Static shadow, spherical counterpart (T7): black-body star r < 1 behind an opaque cold clump (r = 3, a = 0.4)',
             'slice through the clump (phi = %.2f), 48x48x32; grey = inside the clump; white: exact umbra (solid) / lit (dashed) '
             'edges; binary %s,\nFIN3 keys (vet_gd_twin_state trans, vet_rotate_avg 4); Picard %.2f (120 dirs)' % (
                 d['ph'], BIN_D, st.get('pic', float('nan'))))
    save(fig, '07_shd_sp.png', '07_shd_sp', 'Spherical shadow (T7), FIN3 final on %s: whole-wedge lit L1 at r = 4.05 / 4.99 / '
         '5.77 = %s with 120 directions and %s with 432; lit max %s (120); umbra <E>/E_free %.3f / %s (120 / 432); front cells '
         'E/E_free %.2f. Before (FIN2 control, twin_state mean, 120 dirs): lit L1 %s, lit max %s.' % (
             BIN_D, L(f1), L(f4), ' / '.join('%.2f' % v for v in f1['lit_max'].values()), f1['umbra'], ('%.3f' % f4['umbra']) if f4['umbra'] == f4['umbra'] else 'pending',
             f1['front_mean'], L(f2) if f2 else '?', ' / '.join('%.2f' % v for v in f2.get('lit_max', {}).values())))


# ============================================================ 09 static diffusion (J21 4.4), 10/11 dynamic (IDORT Fig. 6)
def f09():
    d = load('dif')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, m in enumerate(d['sdif']):
        ax[0].plot(m['x'], m['E'], color=LEVC[k + 1], label='FIN t = %.0f' % m['t'])
        ax[0].plot(m['x'], m['ex'], ':', color=INK, lw=2, label='analytic' if k == 0 else None)
        ax[1].plot(m['x'], 100*(m['E'] - m['ex'])/m['ex'].max(), color=LEVC[k + 1], label='t = %.0f: L1 %.2f %%' % (m['t'], 100*m['L1']))
    ax[0].set_xlim(-1.2, 1.2)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('E')
    ax[0].legend()
    ax[1].set_xlim(-1.2, 1.2)
    ax[1].axhline(0, **REF)
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('(E - E_analytic) / max E_analytic  [%]')
    ax[1].legend()
    suptitle(fig, 'Static diffusion (J21 4.4), FIN', 'rho kappa_s = 4e4, c = 10; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    s = d['sdif']
    save(fig, '09_sdif.png', '09_sdif', 'Static diffusion of a Gaussian pulse, FIN on %s, vs the analytic solution: L1 %s %% '
         'at the three output times, peak %.4f of analytic.' % (BIN_D, ' / '.join('%.2f' % (100*m['L1']) for m in s), s[-1]['peak']))


def f09sp():
    d = load('sp_dif')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    o = d[256]
    for k, p in enumerate(o['prof']):
        ax[0].plot(o['r'], p['E'], color=LEVC[k], label='t = %.0f' % p['t'] if p['t'] > 0 else 't = 0')
        if p['ex'] is not None:
            ax[0].plot(o['r'], p['ex'], ':', color=INK, lw=2, label='exact' if k == 1 else None)
    ax[0].set_xlim(2, 4)
    ax[0].set_xlabel('r')
    ax[0].set_ylabel('E')
    ax[0].legend()
    ax[0].set_title('FIN (HFfZ), N_r 256')
    ns = sorted(k for k in d if isinstance(k, int))
    for j in range(3):
        errs = []
        for n in ns:
            p = [q for q in d[n]['prof'] if q['ex'] is not None][j]
            eo = p['ex'].min()
            errs.append(np.sum(np.abs(p['E'] - p['ex']))/np.sum(p['ex'] - eo))
        ax[1].loglog(ns, errs, 'o-', color=LEVC[j + 1], ms=7, label='sqrt(Dt) = %.2f' % p['sq'])
    slopes(ax[1], ns[0], errs[0]*1.6, ns)
    ax[1].set_xlabel('N_r')
    ax[1].set_ylabel('L1(E - E_exact) / L1(E_exact)')
    ax[1].legend()
    o64 = [o2 for o2 in d[64]['prof'] if o2['ex'] is not None]
    suptitle(fig, 'Static diffusion, spherical counterpart (T4): Gaussian shell at r0 = 3, FIN (HFfZ)',
             'exact: u = rE obeys the 1-D heat equation; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    e = []
    for n in ns:
        p = [q for q in d[n]['prof'] if q['ex'] is not None][0]
        e.append(np.sum(np.abs(p['E'] - p['ex']))/np.sum(p['ex'] - p['ex'].min()))
    save(fig, '09_sdif_sp.png', '09_sdif_sp', 'Spherical diffusion (T4) on the sp wedge, FIN (HFfZ) on %s, vs the exact '
         'spherical solution: L1 %s at N_r %s (orders %s).' % (BIN_D, ' / '.join('%.1e' % x for x in e),
                                                             '/'.join(map(str, ns)),
                                                             ' '.join('%.2f' % math.log2(e[i]/e[i + 1]) for i in range(len(e) - 1))))


def dd_fig(k, key, fname, title):
    d = load('dif')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.6, 1]))
    for j, m in enumerate(d[k]):
        ax[0].plot(m['x'], m['E'], color=LEVC[min(j, 3)], label='FIN t = %.0f' % m['t'])
        ax[0].plot(m['x'], m['ex'], ':', color=INK, lw=2, label='analytic' if j == 0 else None)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('E')
    ax[0].legend(fontsize=9)
    ts = d[k + '_ts']
    rows = np.array(ts['rows'])
    ax[1].plot(rows[:, 0], rows[:, 1], 'o-', color=C1, ms=6, label='peak / analytic')
    ax[1].plot(rows[:, 0], 1 + rows[:, 2], 's--', color=C2, ms=6, label='1 + (x_c - v t)')
    ax[1].axhline(1, **REF)
    ax[1].set_ylim(0.85, 1.05)
    ax[1].set_xlabel('t')
    ax[1].legend(fontsize=9)
    var = rows[:, 3]
    deff = (var[-1] - var[1])/(rows[-1, 0] - rows[1, 0])/(2*ts['D'])
    suptitle(fig, title, 'v = 1, c = 1e3, rho kappa_s = 4e4, periodic L = 20; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return fig, rows, deff


def f10():
    fig, rows, deff = dd_fig('ddf', '10_ddf', '10_ddf.png', 'Dynamic diffusion (IDORT 4.2.1, Fig. 6; J21 4.4), gas_feedback on, FIN')
    save(fig, '10_ddf.png', '10_ddf', 'Advected diffusing pulse (gas_feedback on), FIN on %s: profiles vs the analytic '
         'moving Gaussian; at t = %.0f peak %.4f of analytic, centroid offset %+.4f, D_eff/D %.3f.' % (
             BIN_D, rows[-1, 0], rows[-1, 1], rows[-1, 2], deff))


def f11():
    fig, rows, deff = dd_fig('ddn', '11_ddn', '11_ddn.png', 'Dynamic diffusion, gas_feedback off (ddn), FIN')
    save(fig, '11_ddn.png', '11_ddn', 'Advected diffusing pulse with kinematic gas (ddn), FIN on %s: at t = %.0f peak %.4f '
         'of analytic, centroid offset %+.4f, D_eff/D %.3f.' % (BIN_D, rows[-1, 0], rows[-1, 1], rows[-1, 2], deff))


# ============================================================ 12 bubbles (IDORT Fig. 7)
def f12():
    """FIN3 page: BB's periodic-vs-outflow layout (bubbles_1012 ana/figbub.py) on the FIN3 final runs (vet_sc rotation and
    vet_rotate_avg 4 are part of FIN3): rows = periodic ring 100 (old suite deck), outflow ring 100, outflow ring 200, and the
    retarded (light-cone) exact solution of the OPEN box (bubref.py, IDORT BCs).  Metrics as figbub.py: E on the 256^2
    reference sample, L1 = sum|E - E_ref| / sum E_ref, E_tot = sum E / sum E_ref."""
    sys.path.insert(0, P + '/figsrc')
    import bubeval as be
    RD = P + '/fin3runs/cart/'
    R = np.load(CACHE + '/bubref_2.npz')
    rows = [('FIN periodic, ring 100 (old suite deck)', RD + 'bub/FIN_r100k5', 'per'),
            ('FIN outflow, ring 100', RD + 'bubu/FIN_r100', 'open'),
            ('FIN outflow, ring 200', RD + 'bubu/FIN_r200', 'open')]
    TS = be.TS
    fig = plt.figure(figsize=(11, 15))
    gs = fig.add_gridspec(len(rows) + 1, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.08, hspace=0.32, left=0.04, right=0.92,
                          top=0.88, bottom=0.02)
    nrm = LogNorm(1e-1, 1e4)
    res = {}
    for r, (lab, d, bcn) in enumerate(rows):
        ok = os.path.isdir(d)
        for c, t in enumerate(TS):
            a = fig.add_subplot(gs[r, c])
            a.set_xticks([])
            a.set_yticks([])
            noimgrid(a)
            if not ok:
                a.text(0.5, 0.5, 'not evaluated\n(test EXCLUDED, user decision)', ha='center', va='center', transform=a.transAxes, color=INK2)
                a.set_title('%s, t = %.3f' % (lab, t), fontsize=8.5)
                continue
            x = be.load_at(d, t)
            E, F1, F2 = be.fields(x, 1)
            Er = R['E2_%s' % bcn][c].reshape(256, 256)
            Es = E[1::2, 1::2]
            l1, et = np.abs(Es - Er).sum()/Er.sum(), Es.sum()/Er.sum()
            res.setdefault(lab, []).append((l1, et))
            im = a.imshow(np.maximum(E, 1e-1), origin='lower', extent=[-.5, .5, -.5, .5], norm=nrm, cmap='viridis')
            a.set_title('%s\nt = %.3f vs %s exact: L1 %.2f, E_tot %.2f' % (lab, t, 'periodic' if bcn == 'per' else 'open', l1, et),
                        fontsize=8.5)
    for c, t in enumerate(TS):
        a = fig.add_subplot(gs[-1, c])
        Er = R['E2_open'][c].reshape(256, 256)
        im = a.imshow(np.maximum(Er, 1e-1), origin='lower', extent=[-.5, .5, -.5, .5], norm=nrm, cmap='viridis')
        a.set_xticks([])
        a.set_yticks([])
        noimgrid(a)
        a.set_title('exact retarded, OPEN box, t = %.3f' % t, fontsize=8.5)
    fig.colorbar(im, cax=fig.add_subplot(gs[1:3, 3])).set_label('E')
    suptitle(fig, 'Multi-source bubbles (IDORT 4.2.2, Fig. 7), FIN: EXCLUDED from the suite (user decision; information only)',
             'binary %s; times 0.05 / 0.25 / 1 t_c\nRT (retard-1012): the time-dependent (retarded) formal-solution reference '
             'reaches L1 0.213 / 0.214 / 0.160' % BIN_D)
    txt = '; '.join('%s: L1 %s, E_tot %s' % (k.split(' (')[0], ' / '.join('%.3f' % v[0] for v in vals),
                                              ' / '.join('%.2f' % v[1] for v in vals)) for k, vals in res.items())
    pend = [lab.split(' (')[0] for lab, d, _ in rows if not os.path.isdir(d)]
    save(fig, '12_bub.png', '12_bub', 'EXCLUDED from the test suite (user decision; numbers for information only). Multi-source bubbles, FIN on %s, vs the retarded (light-cone) free-streaming reference '
         'at IDORT\'s times 0.05 / 0.25 / 1 t_c: %s.%s Retardation-limited: the steady VET tensor carries light that has not '
         'arrived yet; RT\'s time-dependent reference reaches 0.213 / 0.214 / 0.160 (IDORT\'s own L1 0.27).' % (
             BIN_D, txt, (' Not evaluated (test excluded): %s.' % ', '.join(pend)) if pend else ''))


def eq_fig(c, logy):
    d = load('eq')
    r = d[c]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].plot(r['t'], r['E'], 'o-', color=C2, ms=5, label='E_r  FIN')
    ax[0].plot(r['t'], r['T']**4, 's-', color=C1, ms=5, label='a T_gas^4  FIN')
    ax[0].plot(r['t'], r['Er'], ':', color=INK, lw=2, label='analytic (ODE)')
    ax[0].plot(r['t'], r['Tr']**4, ':', color=INK, lw=2)
    if logy:
        ax[0].set_yscale('log')
    ax[0].set_xlabel('t')
    ax[0].set_ylabel('E_r, a T^4')
    ax[0].legend()
    s = r['t'] > 0
    ax[1].semilogy(r['t'][s], np.abs(r['E'][s]/r['Er'][s] - 1) + 1e-16, 'o-', color=C2, ms=5, label='|E/E_ana - 1|')
    ax[1].semilogy(r['t'][s], np.abs(r['T'][s]/r['Tr'][s] - 1) + 1e-16, 's-', color=C1, ms=5, label='|T/T_ana - 1|')
    ax[1].set_xlabel('t')
    ax[1].set_ylabel('relative error')
    ax[1].legend()
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return fig, r


def f13():
    fig, r = eq_fig('eq1', False)
    suptitle(fig, 'Thermal equilibrium, case 1: (T, E) = (1, 100) (IDORT 4.3.1, Fig. 8; J21 4.1), FIN', 'binary %s' % BIN_D)
    save(fig, '13_eq1.png', '13_eq1', 'Thermal relaxation from (T, E) = (1, 100), FIN on %s: E and aT^4 vs the analytic '
         'ODE, final E %.5f, T %.5f (target 96.795, 3.13663); max rel. error E %.1e.' % (BIN_D, r['E'][-1], r['T'][-1], r['errE']))


def f14():
    fig, r = eq_fig('eq2', True)
    suptitle(fig, 'Thermal equilibrium, case 2: (T, E) = (100, 1) (IDORT 4.3.1, Fig. 8; J21 4.1), FIN', 'binary %s' % BIN_D)
    save(fig, '14_eq2.png', '14_eq2', 'Thermal relaxation from (T, E) = (100, 1), FIN on %s: final E %.5f, T %.5f (target '
         '145.788, 3.4747), final mean v %.1e.' % (BIN_D, r['E'][-1], r['T'][-1], r['v']))


def f15():
    d = load('eq')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    tgt = dict(E=1.1326, v=2.9554)
    for arm, col, ls in (('FIN', C1, '-'), ('PRODd', C2, '--')):
        s = np.asarray(d.get('eq3t_' + arm, []))
        if s.ndim != 2: continue   # idort-1020: single PROD arm (no v10 overlay)
        t = s[:, 0]
        ax[0].semilogy(t[1:], np.abs(s[1:, 1]/tgt['E'] - 1) + 1e-16, color=col, ls=ls, label='|E - E_J21|/E  ' + arm)
        ax[0].semilogy(t[1:], np.abs(s[1:, 3]/tgt['v'] - 1) + 1e-16, color=col, ls=ls, lw=1.1, alpha=0.6,
                       label='|v - v_J21|/v  ' + arm)
        ax[1].semilogy(t[1:], np.abs(s[1:, 5]/s[0, 5] - 1) + 1e-17, color=col, ls=ls, label='|P_tot - P_0|/P_0  ' + arm)
        ax[1].semilogy(t[1:], np.abs(s[1:, 6]/s[0, 6] - 1) + 1e-17, color=col, ls=ls, lw=1.1, alpha=0.6,
                       label='|E_tot - E_0|/E_0  ' + arm)
    for a in ax:
        a.set_xlabel('t')
        a.legend(fontsize=9)
    ax[0].set_ylabel('relative error vs code target')
    ax[1].set_ylabel('conservation error')
    sd = d['eq3t_FIN']
    suptitle(fig, 'Moving thermal equilibrium (IDORT 4.3.1 case 3; J21), v = 3, c = 10: PROD',
             'targets E 1.1326, v 2.9554 (NOTE-moveq); binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, '15_eq3t.png', '15_eq3t', 'Moving equilibrium (v = 3, c = 10), PROD on %s: final E %.5f, '
         'v %.5f, P_tot %.6f (exact 3).' % (BIN_D, sd[-1, 1], sd[-1, 3], sd[-1, 5]))


def f13sp():
    _relax_sp('relax1_early_T', '13_eq1_sp', '(T, E) = (1, 100)')


def f14sp():
    _relax_sp('relax2_early_T', '14_eq2_sp', '(T, E) = (100, 1) scaled')


def _relax_sp(case, key, lab):
    d = load('sp_relax')
    ts = d[case + '_ts']
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].plot(ts['t'], ts['E'], 'o', color=C2, ms=7, label='E_r  FIN (dt/1)')
    ax[0].plot(ts['t'], ts['T']**4, 's', color=C1, ms=7, label='a T^4  FIN (dt/1)')
    ax[0].plot(ts['tt'], ts['Eref'], ':', color=INK, lw=2, label='ODE')
    ax[0].plot(ts['tt'], ts['Tref']**4, ':', color=INK, lw=2)
    ax[0].set_yscale('log')
    ax[0].set_xlabel('t')
    ax[0].legend()
    ax[0].set_title('interior mean (r < 1.8)')
    ax[0].ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    lv = (1, 2, 4, 8, 16)
    errs = {}
    for n in lv:
        txt = d['%s_%d' % (case, n)]
        errs[n] = [(float(a), float(b)) for a, b in re.findall(r't ([0-9.e+-]+) E ([0-9.e+-]+) T', txt)]
    times = [e[0] for e in errs[1]]
    for j, t in enumerate(times):
        ax[1].loglog(lv, [max(errs[n][j][1], 1e-16) for n in lv], 'o-', color=LEVC[min(j, 3)], ms=7, label='t = %.1e' % t)
    slopes(ax[1], 1, errs[1][0][1], lv)
    ax[1].axhline(1e-12, color='#b0b0ac', lw=0.8)
    ax[1].set_xlabel('dt divisor')
    ax[1].set_ylabel('max |E/E_ODE - 1| (r < 1.8)')
    ax[1].legend(fontsize=8.5)
    suptitle(fig, 'Thermal relaxation, spherical counterpart (T1) %s, FIN (HFfZ, hesdirk2)' % lab,
             'uniform gas on the sp wedge r in [1, 2]; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    e1 = errs[1][-1][1]
    e8 = errs[8][-1][1]
    save(fig, key + '.png', key, 'Thermal relaxation %s on the sp wedge, FIN (HFfZ) on %s: interior E, aT^4 vs the ODE and '
         'the error vs dt divisor; at t = %.1e the error is %.1e at dt/1 and %.1e at dt/8 (round-off floor).' % (
             lab, BIN_D, times[-1], e1, e8))


# ============================================================ 16/17 linear waves (IDORT 4.3.2 conv figure)
def rw_fig(P, key, fname):
    d = load('rw')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    taus = (0.01, 1.0, 100.0)
    cols = (INK, '#e34948', C1)
    out = []
    for tau, col in zip(taus, cols):
        for arm, mk, fill in (('FIN', 'o', True), ('PRODd', 's', False)):
            k = (P, tau, arm)
            if k not in d:
                continue
            Ns = sorted(d[k])
            l1 = [d[k][n]['l1_rho'] for n in Ns]
            ax[0].loglog(Ns, l1, mk + ('-' if fill else '--'), color=col, ms=8 if fill else 6, mfc=col if fill else 'white',
                         lw=1.2, label='tau = %g  %s' % (tau, arm))
            w = [abs(complex(d[k][n]['w_re'], d[k][n]['w_im']) - complex(*d[k][n]['refs'][d[k][n]['ref_key']])) /
                 abs(complex(*d[k][n]['refs'][d[k][n]['ref_key']])) for n in Ns]
            ax[1].loglog(Ns, w, mk + ('-' if fill else '--'), color=col, ms=8 if fill else 6, mfc=col if fill else 'white', lw=1.2,
                         label='tau = %g  %s (ref %s)' % (tau, arm, d[k][Ns[-1]]['ref_key']))
            if arm == 'FIN':
                out.append('tau %g: L1(%d) %.1e, order %.2f%s' % (tau, Ns[-1], l1[-1], math.log2(l1[-2]/l1[-1]), '' if Ns[-1] >= 256 else ' (N %d-256: pending, final rerun running)' % (2*Ns[-1])))
    k0 = (P, 0.01, 'FIN')
    Ns = sorted(d[k0])
    slopes(ax[0], Ns[0], d[k0][Ns[0]]['l1_rho']*1.5, Ns)
    ax[0].set_xlabel('N')
    ax[0].set_ylabel('L1 norm of rho error / amplitude')
    ax[0].legend(fontsize=8, ncol=2)
    ax[0].set_title('a T0^4 / P0 = %g' % P)
    ax[1].set_xlabel('N')
    ax[1].set_ylabel('|omega - omega_ref| / |omega_ref|')
    ax[1].legend(fontsize=7.5)
    ax[1].set_title('complex frequency vs the model dispersion root')
    suptitle(fig, 'Linear radiation-acoustic waves (IDORT 4.3.2), P = %g: FIN' % P,
             'c = 10, one period, nt = 2N; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, fname, key, 'Linear-wave convergence at P = %g vs the exact eigenmode, FIN on %s: %s.' % (
        P, BIN_D, '; '.join(out)))


def f16():
    rw_fig(0.01, '16_rw_p001', '16_rw_p001.png')


def f17():
    rw_fig(10.0, '17_rw_p10', '17_rw_p10.png')


def parse_ord(fn):
    """sp ana/out_fin xord/tord: (case, arm, var) -> errors list"""
    out = {}
    for ln in open(fn):
        m = re.match(r'(\S+)\s+(\S+)\s+(\S+)\s+d?L1/mean\s+([0-9.e+\- ]+?)\s+order', ln)
        if m:
            out[(m.group(1), m.group(2), m.group(3))] = [float(x) for x in m.group(4).split()]
    return out


def wave_sp(case, key, title, status_note):
    d = load('sp_wave')
    xo = parse_ord(F3SP + '/xord.txt')
    to = parse_ord(F3SP + '/tord.txt')
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.4))
    for k, n in enumerate((64, 128, 256)):
        o = d['%s_%d' % (case, n)]
        ax[0].plot(o['r'], o['dens'] - 1, color=LEVC[k + 1], label='rho - 1, N_r %d' % n)
        ax[0].plot(o['r'], (o['E'] - 0.03)/0.03*0.1, color=LEVC[k + 1], ls='--', lw=1.1, label='0.1 dE/E0' if k == 2 else None)
    ax[0].set_xlabel('r')
    ax[0].set_title('final state, FIN (HFfZ)', fontsize=10.5)
    ax[0].legend(fontsize=8)
    vs = (('dens', 'rho', C1), ('velx', 'v', C2), ('m1_e', 'E', C3), ('m1_f1', 'F1', C4))
    txt = []
    for a, src, lev, lab in ((ax[1], xo, (128, 256), 'X: N_r (pairs 64-128, 128-256)'), (ax[2], to, None, 'T: dt divisor')):
        for v, nm, col in vs:
            e = src.get((case if src is xo else case + '_T', 'HFfZ', v))
            if not e:
                continue
            x = lev if lev else [2**(i + 1) for i in range(len(e))]
            a.loglog(x, e, 'o-', color=col, ms=7, label=nm)
            em = src.get((case if src is xo else case + '_T', 'HFfW', v))
            if em:
                a.loglog(x, em, 's--', color=col, ms=5, mfc='white', lw=1, alpha=0.8)
            txt.append('%s %s %s' % ('X' if src is xo else 'T', nm, ' '.join('%.2f' % math.log2(e[i]/e[i + 1]) for i in range(len(e) - 1))))
        a.set_xlabel(lab)
        a.set_ylabel('self-convergence L1 diff / mean')
        a.legend(fontsize=8.5)
        xs = lev if lev else [2, 4, 8]
        e0 = src.get((case if src is xo else case + '_T', 'HFfZ', 'm1_e'))
        if e0:
            slopes(a, xs[0], e0[0]*2, xs, lab=False)
    ax[2].plot([], [], 's--', color='#8a8a86', mfc='white', label='PRODd (open)')
    ax[2].legend(fontsize=8.5)
    suptitle(fig, title, 'dotted / dashed grey: slopes -1 / -2; binary %s; orders from fin3_1012 ana/out_281dd686 (sp) xord/tord' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, key + '.png', key, '%s, FIN (HFfZ) on %s: final profiles and X / T self-convergence; finest-pair orders %s%s.' % (
        title.split(':')[0], BIN_D, '; '.join(txt), status_note))


def f16sp():
    wave_sp('rsw', '16_rw_p001_sp', 'Radiation-dominated linear wave on the sp wedge (T9 rsw): spherical counterpart', '')


def f17sp():
    wave_sp('rswv', '17_rw_p10_sp', 'Moving-shell wave on the sp wedge (T9v rswv): spherical counterpart',
            ' (KNOWN-OPEN: X floor from the outer-wall boundary layer, T 1.2-1.8 with implicit_f0_init auto, 1.8-2.0 with zero)')


# ============================================================ 18 shocks (IDORT Figs. 11-13)
def f18():
    d = load('shock')
    fig, ax = plt.subplots(2, 3, figsize=(12, 7.2), sharex='col')
    txt = []
    for c, m in enumerate(('1.2', '2', '3')):
        R = d['ref' + m]
        sel = np.abs(R['x'] - R['xs']) < 0.012
        for n, col, ls in (('n512', C1, '-'), ('n1024', C3, '--')):
            k = '%s_%s' % (m, n)
            if k not in d:
                continue
            o = d[k]
            xx = (o['x'] - o['xs'])*1e3
            ax[0, c].plot(xx, o['tg'], color=col, ls=ls, label='T_gas FIN %s' % n[1:])
            ax[0, c].plot(xx, o['tr'], color=C2 if n == 'n512' else C4, ls=ls, label='T_rad FIN %s' % n[1:])
            ax[1, c].plot(xx, o['rho'], color=col, ls=ls, label='FIN %s' % n[1:])
            for ln in open(F3SP + '/cart_fin.txt'):
                if ln.startswith('rs_m%s/FIN_%s ' % (m, n)):
                    q = re.search(r'L1 rho ([0-9.]+) Tg ([0-9.]+) Tr ([0-9.]+)', ln)
                    txt.append('M %s (%s) L1 rho / T_gas / T_rad %s / %s / %s' % (m, n[1:], q.group(1), q.group(2), q.group(3)))
        ax[0, c].plot((R['x'][sel] - R['xs'])*1e3, R['tg'][sel], ':', color=INK, lw=2, label='IDORT ref')
        ax[0, c].plot((R['x'][sel] - R['xs'])*1e3, R['tr'][sel], ':', color=INK, lw=2)
        ax[1, c].plot((R['x'][sel] - R['xs'])*1e3, R['rho'][sel], ':', color=INK, lw=2, label='IDORT ref')
        ax[0, c].set_title('M = %s' % m)
        ax[1, c].set_xlabel('(x - x_shock) x 1e3')
        ax[1, c].set_xlim(-11, 6)
    ax[0, 0].set_ylabel('T / T0')
    ax[1, 0].set_ylabel('rho / rho0')
    ax[0, 0].legend(fontsize=8)
    ax[1, 0].legend(fontsize=8)
    ax[0, 2].legend(fontsize=8, loc='lower right')
    suptitle(fig, 'Radiation shocks (IDORT 4.3.3; J21 4.8), FIN vs the IDORT transport reference',
             'after 3 flow crossings; aligned on the density jump; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    pl = ', '.join('M %s %s' % (m, n[1:]) for m in ('1.2', '2', '3') for n in ('n512', 'n1024')
                   if '%s_%s' % (m, n) not in d and not (m == '1.2' and n == 'n1024'))
    save(fig, '18_shock.png', '18_shock', 'Radiating shocks M = 1.2 / 2 / 3, FIN on %s (512 and 1024 cells) vs the IDORT '
         'Zenodo transport reference at 512 cells: %s (fin3_1012 ana/out_281dd686 cart.txt, shockcmp.py).%s' % (
             BIN_D, '; '.join(txt), (' Pending (final rerun running): %s.' % pl) if pl else ''))



def f19():
    """Q2: IDORT/J21 conventions (rho, log 1-40, jet, t = 0.06 / 0.10 / 0.16 t0); FIN with the IDORT beam
    (540 a T0^4 per beam = 4 pi (1031.3 c a) w / c, w = 1/24) next to the old 200-per-beam run and IDORT's own snapshots"""
    d = load('cloud')
    fig = plt.figure(figsize=(10, 8.6))
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1.12, 0.035], wspace=0.08, hspace=0.22, left=0.08, right=0.94,
                          top=0.88, bottom=0.06)
    keys = [3, 5, 8]
    cols = [('FINb540', 'FIN, IDORT beam (540 per beam)')]
    for r, k in enumerate(keys):
        for c, (a_, lab) in enumerate(cols):
            a = fig.add_subplot(gs[r, c])
            if '%s_%d' % (a_, k) not in d:
                a.text(0.5, 0.5, 'pending\n(final rerun running)', ha='center', va='center', transform=a.transAxes, color=INK2)
                a.set_xticks([])
                a.set_yticks([])
                noimgrid(a)
                a.set_title(lab, fontsize=9)
                continue
            o = d['%s_%d' % (a_, k)]
            im = a.pcolormesh(o['x1'], o['x2'], o['dens'], norm=LogNorm(1, 40), cmap='jet', shading='nearest',
                              rasterized=True)
            a.set_aspect('equal')
            noimgrid(a)
            a.set_title('%s\nt = %.2f: halo %.3f, rho>3 %.3f..%.3f, max %.0f' % (lab, o['t'], o['halo'], o['r3'][0],
                                                                              o['r3'][1], o['rmax']), fontsize=9)
            if c:
                a.set_yticklabels([])
            if r < 2:
                a.set_xticklabels([])
        a = fig.add_subplot(gs[r, 1])
        a.imshow(d['idort_png'][k])
        a.set_axis_off()
        a.set_title('IDORT (Zenodo shadow_2d, snapshot %03d)' % k, fontsize=9)
    fig.colorbar(im, cax=fig.add_subplot(gs[:, 2])).ax.set_title('rho', fontsize=9.5)
    o = [d['FINb540_%d' % k] for k in keys]
    st = d.get('st_FINb540', {})
    suptitle(fig, 'Irradiated cloud (IDORT 4.3.4, Fig. 14; J21 Fig. 5): density, FIN with the IDORT beam',
             '+-15 deg beams from the left, 512x256, ring 12x2 (36 rays), dt 3e-4; halo = left edge of rho > 1.08 on y = 0; '
             'binary %s; Picard %.2f' % (BIN_D, st.get('picard', float('nan'))))
    save(fig, '19_cloud.png', '19_cloud', 'Irradiated cloud density at t = 0.06 / 0.10 / 0.16, FIN on %s with the IDORT beam '
         '(I = 1031.3 c a = 540 a T0^4 per beam): halo edge %s, dense (rho > 3) cloud %s, peak rho %s; IDORT snapshots (read off '
         'the PNGs on y = 0): halo %s, rho > 3 %s. Reproduces IDORT (qualitative test).' % (
             BIN_D, ' / '.join('%.3f' % x['halo'] for x in o), ' / '.join('%.3f..%.3f' % x['r3'] for x in o),
             ' / '.join('%.0f' % x['rmax'] for x in o), ' / '.join('%.3f' % d['idort_meas'][k]['halo'] for k in keys),
             ' / '.join('%.3f..%.3f' % d['idort_meas'][k]['r3'] for k in keys)))


# ============================================================ 20 Marshak (JSD12 5.1)
def f20():
    d = load('marsh')
    ref = d['ref']
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    txt = []
    for a, run, lab in ((ax[0], 'FIN_c1', 'c dt/dx = 1'), (ax[1], 'FIN_c10', 'c dt/dx = 10')):
        o = d[run]['tabs'][-1]
        a.plot(o['x'], o['E'], color=C1, label='FIN, t = %.3g' % o['t'])
        a.plot(ref[:, 0], ref[:, 1], ':', color=INK, lw=2, label='S_N reference')
        er = np.interp(o['x'], ref[:, 0], ref[:, 1])
        l1 = np.abs(o['E'] - er).sum()/np.abs(er).sum()
        a.set_title('%s: L1(E) = %.4f' % (lab, l1))
        a.set_xlabel('x')
        a.set_ylabel('E')
        a.legend()
        txt.append('%s L1(E) %.4f' % (lab, l1))
    suptitle(fig, 'Marshak wave (JSD12 5.1), FIN vs the S_N transport reference', 'binary %s; criterion L1(E) <= 0.02' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, '20_marsh.png', '20_marsh', 'Marshak wave E profile, FIN on %s, vs the S_N reference: %s.' % (BIN_D, '; '.join(txt)))


def f20sp():
    d = load('sp_marsh')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, n in enumerate((64, 128, 256)):
        o = d['marsh_r1_u_%d' % n]
        ax[0].plot(o['r'] - 1, o['E'], color=LEVC[k + 1], label='r0/w = 1, uniform, N_r %d' % n)
    o = d['marsh_r1_s_256']
    ax[0].plot(o['r'] - 1, o['E'], color=C2, ls='--', label='r0/w = 1, stretched, N_r 256')
    ax[0].set_xlabel('r - r_in')
    ax[0].set_ylabel('E')
    ax[0].legend(fontsize=8.5)
    ax[0].set_title('spherical Marshak wave, t = %.3g' % o['t'])
    a = d['marsh_r1000_u_128']
    c = d['marsh_cart_128']
    ax[1].plot(a['r'] - a['r'][0] + (c['r'][0] - 1), a['E'], color=C1, label='sp r0/w = 1000 (FIN)')
    ax[1].plot(c['r'] - 1, c['E'], ':', color=INK, lw=2, label='Cartesian vet_sc (same keys)')
    ax[1].set_xlabel('distance from the inner face')
    ax[1].set_ylabel('E')
    ax[1].legend()
    l1 = np.sum(np.abs(a['E'] - c['E']))/np.sum(np.abs(c['E']))
    ax[1].set_title('large-radius control: L1|sp - cart| / L1(cart) = %.1e' % l1)
    suptitle(fig, 'Marshak wave, spherical counterpart (T5), FIN (HFfZ)', 'binary %s; front converges at first order '
             '(sp X ladder: ladder_sp_marsh_r1.png)' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, '20_marsh_sp.png', '20_marsh_sp', 'Spherical Marshak wave (T5) on the sp wedge, FIN (HFfZ) on %s: E profiles at '
         'N_r 64-256 and the r0/w = 1000 shell vs the Cartesian twin (L1 difference %.1e).' % (BIN_D, l1))


# ============================================================ 21 crooked pipe (JSD12 5.2)
def f21():
    d = load('pipe')
    fig = plt.figure(figsize=(12, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.22, left=0.06, right=0.97, top=0.84, bottom=0.12)
    a = fig.add_subplot(gs[0, 0])
    names = ['(0.25,0)', '(2.75,0)', '(3.5,1.25)', '(4.25,0)', '(6.75,0)']
    cols = plt.cm.viridis(np.linspace(0, 0.85, 5))
    for run, ls in (('FIN_c220', '-'), ('FIN_c2200', '--')):
        pr = d[run]['probes']
        for k in range(5):
            a.plot(pr[:, 0], pr[:, k + 1]*0.5, color=cols[k], ls=ls, lw=1.5, label=names[k] if run == 'FIN_c220' else None)
    a.plot([], [], '-', color=INK2, label='c dt/dx 220')
    a.plot([], [], '--', color=INK2, label='c dt/dx 2200')
    a.set_yscale('log')
    a.set_xscale('log')
    a.set_xlabel('t [1e-8 s]')
    a.set_ylabel('T_gas [keV]')
    a.legend(fontsize=8, ncol=2)
    a.set_title('probe histories (FIN)')
    m = d['map']
    b = fig.add_subplot(gs[0, 1])
    im = b.pcolormesh(m['x1'], m['x2'], m['T']*0.5, cmap='inferno', vmin=0.05, vmax=0.5, shading='nearest', rasterized=True)
    b.set_aspect('equal')
    noimgrid(b)
    for (px, py), c in zip(((0.25, 0), (2.75, 0), (3.5, 1.25), (4.25, 0), (6.75, 0)), cols):
        b.plot(px, py, 'o', mfc=c, mec='white', ms=7)
    fig.colorbar(im, ax=b, fraction=0.03, pad=0.02).set_label('T_gas [keV]')
    b.set_title('T_gas, c dt/dx 220, t = %.1f' % m['t'])
    b.set_xlabel('x')
    b.set_ylabel('y')
    p1 = d['FIN_c220']['probes']
    p2 = d['FIN_c2200']['probes']
    i1 = int(np.argmin(np.abs(p1[:, 0] - 10)))
    i2 = int(np.argmin(np.abs(p2[:, 0] - 10)))
    dev = np.max(np.abs(p2[i2, 1:]/p1[i1, 1:] - 1))
    suptitle(fig, 'Crooked pipe (JSD12 5.2, Gentile 2001 geometry, Cartesian 512x256), FIN', 'static gas; binary %s; '
             'no quantitative Cartesian reference: dt self-consistency' % BIN_D)
    save(fig, '21_pipe.png', '21_pipe', 'Crooked pipe, FIN on %s: T_gas map and probe histories; probes at t = %.1f: %s keV '
         '(c dt/dx 220), c dt/dx 2200 within %.1f %%. Caveat (IN): the probe-5 arrival time t(T >= 0.5) is 4.20 / 3.31 (c dt/dx 220 / 2200) vs FIN2 3.81 / 3.71, so FIN3 is less dt-consistent than FIN2 there; PASS is on the end state with our own criterion.' % (BIN_D, p1[i1, 0], ' / '.join('%.3f' % (0.5*x) for x in p1[i1, 1:]),
                                                         100*dev))


def f22_24():
    for o, k, nm in ((22, 'photon_bubble', 'Photon bubble (JSD12 6.4)'), (23, 'dusty_wind', 'Dusty wind (J21 4.9)'),
                     (24, 'sincell', 'Sincell VET shock (JSD12 6.3)')):
        nodata('%02d_%s' % (o, k), '%s: not run (parked, not feasible here; NOTE-suite 3.6).' % nm)


# ============================================================ spherical extras
def f25():
    d = load('sp_patm')
    for case, key, ttl in (('patm_L16_s', '25_patm_sp', 'P1 grey atmosphere through the photosphere (patm)'),
                           ('patmT_L16_s', '26_patmT_sp', 'P1T thermal radiative-equilibrium atmosphere (patmT, coupled gas)')):
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
        eph = []
        for k, n in enumerate((64, 128, 256)):
            o = d['%s_%d' % (case, n)]
            m = (o['tau'] < 20) & (o['tau'] > 0)
            ax[0].semilogx(o['tau'][m], 100*o['rel'][m], color=LEVC[k + 1], label='N %d' % n)
            eph.append(o['e_ph'])
        ax[0].axvline(2/3, color='#aaaaaa', ls=':', lw=1)
        ax[0].axhline(0, **REF)
        ax[0].invert_xaxis()
        ax[0].set_xlabel('radial optical depth from the top')
        ax[0].set_ylabel('E / J_ref - 1  [%]')
        ax[0].legend()
        ax[0].set_title('vs 1-D spherical VEF transport reference')
        names = ('tau<0.1', '0.1-1', '1-10', '10-20')
        for j, nm in enumerate(names):
            e = [d['%s_%d' % (case, n)]['bands'][nm][0] for n in (64, 128, 256)]
            ax[1].loglog((64, 128, 256), e, 'o-', color=(C1, C2, C3, C4)[j], ms=7, label='band ' + nm)
        slopes(ax[1], 64, d['%s_64' % case]['bands']['10-20'][0]*1.5, (64, 128, 256))
        ax[1].set_xlabel('N_r (poly-stretched)')
        ax[1].set_ylabel('volume-weighted L1 of E/J_ref - 1')
        ax[1].legend(fontsize=8.5)
        suptitle(fig, 'Spherical extra: %s, FIN (HFfZ)' % ttl, '16x16 lateral, c dt/dx_min 1e3; binary %s' % BIN_D)
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        save(fig, key + '.png', key, '%s on the sp wedge, FIN (HFfZ) on %s: E/J_ref - 1 at R_ph %s %% at N 64 / 128 / 256 '
             '(thin band tau < 0.1 not converging at 16x16 lateral, KNOWN lateral-resolution plateau).' % (
                 ttl, BIN_D, ' / '.join('%+.2f' % (100*x) for x in eph)))


def f27():
    """FIN3 page: P2 in its NEW definition (user approved): the wide channel pinc_ha20_sr30 (hole 0.2 rad, screen r 2-3),
    AV's spt2.py / p2eval.py metric: M1 luminosity through the shells r > 4.1 over the staircase Monte Carlo L_ex of the
    run's own lateral grid.  PASS rule: G2 (GL 20x24 = 480 directions, rotation off) within 5 %% for both direction sets.
    G1 (120) with vet_rotate_avg 4 (the FIN3 default) and without rotation (R0) shown with both sets; the narrow pinc
    (hole 0.06) is the angular stress case (no pass/fail).  Runs: fin3_1012/p2_FIN3_281dd686."""
    sys.path.insert(0, SF + '/sp/ana')
    import glob as _g
    import p2eval
    import spl
    W = F3 + '/p2_FIN3_281dd686/'

    def run(case, arm, stag):
        for d in sorted(_g.glob(W + case + '/' + arm + '+*_64')):
            if ('gl_stagger~false' in d) == (stag == 'B'):
                return d
        return None

    def prof(d):
        fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if '.m1.' in f)
        b = spl.binall(os.path.join(d, 'bin', fs[-1]))
        F = b['m1_f1']
        n3, n2, n1 = F.shape
        rf = np.linspace(1.0, 5.0, n1 + 1)
        rc = 0.5*(rf[1:] + rf[:-1])
        tf = np.linspace(0.5*math.pi - 0.35, 0.5*math.pi + 0.35, n2 + 1)
        dom = (np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :]*(math.pi/3.0)/n3*np.ones((n3, 1))
        L = np.array([np.sum(F[:, :, i]*dom)*rc[i]**2 for i in range(n1)])
        g = p2eval.geom(d)
        lex = p2eval.LEX[(round(g[0], 3), round(g[1], 3))][n2]
        return rc, rf, tf, L/lex, b['m1_e']
    arms = [('pinc_ha20_sr30_n40', 'HFfZR0G2', 'A', C1, '-', 'G2 480 dirs, set A'),
            ('pinc_ha20_sr30_n40', 'HFfZR0G2', 'B', C1, '--', 'G2 480 dirs, set B'),
            ('pinc_ha20_sr30', 'HFfZ', 'A', C2, '-', 'G1 120 dirs, rotate_avg 4 (FIN3), set A'),
            ('pinc_ha20_sr30', 'HFfZ', 'B', C2, '--', 'G1 120 dirs, rotate_avg 4 (FIN3), set B'),
            ('pinc_ha20_sr30', 'HFfZR0', 'A', C3, '-', 'G1 120 dirs, no rotation, set A'),
            ('pinc_ha20_sr30', 'HFfZR0', 'B', C3, '--', 'G1 120 dirs, no rotation, set B')]
    ev = {}
    fig = plt.figure(figsize=(12, 5.4))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 0.05, 0.2, 1.25, 0.75], wspace=0.3, left=0.065, right=0.985, top=0.8, bottom=0.12)
    a = fig.add_subplot(gs[0, 3])
    for case, arm, stg, col, ls, lab in arms:
        d = run(case, arm, stg)
        if d is None:
            continue
        rc, rf, tf, LL, E = prof(d)
        e = p2eval.ev(d)
        ev[(case, arm, stg)] = e
        m = rc > 3.05
        a.plot(rc[m], LL[m], color=col, ls=ls, label='%s: %.3f' % (lab, e['m1']))
        if (case, arm, stg) == ('pinc_ha20_sr30', 'HFfZ', 'A'):
            Emap, rf0, tf0, st0 = E, rf, tf, spl.logstats(d)
    a.axhline(1, **REF, label='Monte Carlo L_ex (staircase)')
    a.axhspan(0.95, 1.05, color='#d9d9d5', alpha=0.5, lw=0, label='PASS band (G2, +-5 %)')
    a.axvline(p2eval.RMIN, color=INK2, lw=0.8, ls=':')
    a.set_xlabel('r')
    a.set_ylabel('channel luminosity L(r) / L_ex')
    a.set_ylim(0.6, 1.4)
    a.legend(fontsize=7.6, loc='upper left')
    a.set_title('wide channel: shell luminosity (mean over r > 4.1 in legend)', fontsize=10)
    # meridional slice through the hole (phi = 30 deg) of the FIN3 G1 run
    ax0 = fig.add_subplot(gs[0, 0])
    n3 = Emap.shape[0]
    k = int(n3//2)
    R, TH = np.meshgrid(rf0, tf0 - math.pi/2)
    X, Y = R*np.cos(TH), -R*np.sin(TH)
    z = np.maximum(Emap[k], 1e-12)
    vmax = np.nanmax(z[:, rf0[1:] > 4.2])
    im = ax0.pcolormesh(X, Y, z, norm=LogNorm(vmax*1e-3, vmax*30), cmap='magma', shading='flat', rasterized=True)
    ax0.set_aspect('equal')
    noimgrid(ax0)
    ax0.set_title('FIN E, G1 rotate_avg 4, set A', fontsize=10)
    ax0.set_xlabel('r cos(lat)')
    ax0.set_ylabel('r sin(lat)')
    fig.colorbar(im, cax=fig.add_subplot(gs[0, 1]))
    # narrow stress case (no pass / fail)
    a2 = fig.add_subplot(gs[0, 4])
    nar = []
    for j, stg in enumerate(('A', 'B')):
        d = run('pinc', 'HFfZ', stg)
        if d is None:
            continue
        e = p2eval.ev(d)
        nar.append(e['m1'])
        a2.bar(j, e['m1'], color=C2, alpha=1 if stg == 'A' else 0.55)
        a2.text(j, e['m1'] + 0.03, '%.2f' % e['m1'], ha='center', fontsize=9)
    a2.axhline(1, **REF)
    a2.set_xticks([0, 1])
    a2.set_xticklabels(['set A', 'set B'])
    a2.set_ylim(0, max(nar + [1.3]) + 0.2)
    a2.set_title('narrow channel (stress case)\nG1 FIN3, L / L_ex', fontsize=10)
    g2 = [ev.get(('pinc_ha20_sr30_n40', 'HFfZR0G2', s), {}).get('m1', float('nan')) for s in 'AB']
    g1 = [ev.get(('pinc_ha20_sr30', 'HFfZ', s), {}).get('m1', float('nan')) for s in 'AB']
    g0 = [ev.get(('pinc_ha20_sr30', 'HFfZR0', s), {}).get('m1', float('nan')) for s in 'AB']
    ver = 'PASS' if all(abs(x - 1) <= 0.05 for x in g2) else 'FAIL'
    suptitle(fig, 'Spherical extra P2 (new definition): thermal channel beam, wide channel pinc_ha20_sr30, FIN: %s' % ver,
             'hot core under an opaque shell (r 2-3) pierced by a 0.2-rad channel; c dt/dx 300; binary %s; Picard %.2f (G1 FIN3)\n'
             'PASS rule: G2 (480 dirs) within 5 %% of the Monte Carlo L_ex for both direction sets' % (BIN_D, st0.get('pic', float('nan'))))
    save(fig, '27_pinc_sp.png', '27_pinc_sp', 'P2 thermal channel beam, new definition (wide channel pinc_ha20_sr30, hole 0.2 rad, '
         'screen r 2-3), FIN on %s, M1 luminosity beyond the screen over the staircase Monte Carlo L_ex: G2 480 dirs %.3f / %.3f '
         '(sets A / B; %s, rule within 5 %%); G1 120 dirs with vet_rotate_avg 4 (FIN3 default) %.3f / %.3f; G1 without rotation '
         '%.3f / %.3f (ray-effect spread). Narrow channel (0.06 rad, stress case, no pass/fail): %s.' % (
             BIN_D, g2[0], g2[1], ver, g1[0], g1[1], g0[0], g0[1], ' / '.join('%.2f' % x for x in nar)))


def f28():
    d = load('sp_wave')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    c = d['pulse_cart_128']
    x = c['r'] - 1
    for k, (case, col, ls, lab) in enumerate((('pulse_r1_128', C2, '-', 'sp r0/w = 1'), ('pulse_r1000_128', C1, '-', 'sp r0/w = 1000'))):
        o = d[case]
        ax[0].plot(x, o['E'], color=col, ls=ls, label=lab)
        ax[1].plot(x, (o['E'] - c['E'])/c['E'].max(), color=col, ls=ls, label=lab)
    ax[0].plot(x, c['E'], ':', color=INK, lw=2, label='Cartesian vet_sc')
    ax[0].set_xlabel('distance from the inner face')
    ax[0].set_ylabel('E at t = %.3g' % c['t'])
    ax[0].legend()
    ax[1].axhline(0, **REF)
    ax[1].set_xlabel('distance from the inner face')
    ax[1].set_ylabel('(E_sp - E_cart) / max E_cart')
    ax[1].set_yscale('symlog', linthresh=1e-3)
    ax[1].legend()
    l1 = [np.mean(np.abs(d[k]['E'] - c['E']))/c['E'].max() for k in ('pulse_r1_128', 'pulse_r1000_128')]
    suptitle(fig, 'Spherical extra L: radiating pulse on a shell vs the Cartesian slab (large-radius control), FIN',
             'rho kappa_s = 1, N 128; r0/w = 1 differs by true spherical focusing; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, '28_pulse_lr_sp.png', '28_pulse_lr_sp', 'Pulse on a spherical shell vs the Cartesian slab, FIN on %s: '
         'L1 difference %.1e at r0/w = 1 (spherical focusing) and %.1e at r0/w = 1000.' % (BIN_D, l1[0], l1[1]))


def f29():
    d = load('sp_wave')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    txt = []
    for a, base in ((ax[0], 'rsw'), (ax[1], 'rswt')):
        s, c = d['%s_r1000_128' % base], d['%s_cart_128' % base]
        x = c['r'] - c['r'][0]
        a.plot(x, s['dens'] - 1, color=C1, label='rho - 1, sp r0/w 1000')
        a.plot(x, c['dens'] - 1, ':', color=INK, lw=2, label='rho - 1, Cartesian')
        a.plot(x, s['E']/np.mean(c['E']) - 1, color=C2, label='E/<E> - 1, sp')
        a.plot(x, c['E']/np.mean(c['E']) - 1, ':', color=INK2, lw=2, label='E/<E> - 1, Cartesian')
        a.set_xlabel('distance from the inner face')
        a.legend(fontsize=8.5)
        l1 = np.mean(np.abs(s['dens'] - c['dens']))/np.mean(np.abs(c['dens'] - np.mean(c['dens'])))
        a.set_title('%s: L1|sp - cart| / L1(rho pert.) = %.1e' % (base, l1))
        txt.append('%s %.1e' % (base, l1))
    suptitle(fig, 'Spherical extra L: radiation waves on a thin shell (r0/w = 1000) vs Cartesian, FIN (HFfZ / HCFfZ)',
             'N 128; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, '29_rsw_lr_sp.png', '29_rsw_lr_sp', 'Large-radius control for the radiation waves, FIN on %s: shell (r0/w 1000) vs '
         'Cartesian, relative L1 density difference %s.' % (BIN_D, ', '.join(txt)))


def f30():
    d = load('g5')
    fig, ax = plt.subplots(1, 4, figsize=(12, 3.9))
    txt = []
    for j, k in enumerate(('12.8', '128', '1280', '12800')):
        kap = float(k)
        D = 1.0/(3*kap)
        tau = 1.0/kap
        a = ax[j]
        for arm, col, mk in (('FIN', C1, 'o'), ('PROD', C2, 's')):
            if '%s_%s' % (arm, k) not in d:
                continue
            o = d['%s_%s' % (arm, k)]
            t, s2, use = o['s'][:, 0], o['s'][:, 1], o['s'][:, 2] > 0.5
            ex = s2[0] + 2*D*(t - tau*(1 - np.exp(-t/tau)))
            rat = o['res'].split('ratio=')[1].split()[0]
            a.plot(t, (s2 - ex)/ex[-1]*100, mk + '-', color=col, ms=5, lw=1, label='%s ratio %s' % (arm, rat))
            if arm == 'FIN':
                txt.append('kappa %s: %s' % (k, rat))
        a.axhline(0, **REF)
        a.set_title('kappa = %s' % k)
        a.set_xlabel('t')
        a.legend(fontsize=7.5)
    ax[0].set_ylabel('(sigma^2 - exact)/sigma^2_end [%]')
    suptitle(fig, 'G5 2-D diffusing pulse (g5fix gate): d sigma^2 / dt vs 2D, FIN on %s (PROD %s open)' % (BIN_D, BIN_O),
             'exact moment line sigma0^2 + 2D [t - tau (1 - exp(-t/tau))]')
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, '30_g5.png', '30_g5', 'G5 2-D pulse variance growth vs the exact moment solution, FIN on %s (PROD on %s for comparison): '
         'd sigma^2/dt / 2D = %s.' % (BIN_D, BIN_O, ', '.join(txt)))


# ============================================================ ladders
def parse_lad():
    out = {}
    for ln in open(F3SP + '/lad_FIN.txt'):
        m = re.match(r'(\S+)\s+(\S+)\s+\S+\s+([XTC])\s+(\S+)\s+lev\s+(\S+)\s+L1\s+([0-9.e+\- ]+)\|', ln)
        if m:
            levs = [int(x) for x in m.group(5).split('-')]
            out[(m.group(1), m.group(2), m.group(3), m.group(4))] = (levs, [float(x) for x in m.group(6).split()])
    return out


VARC = {'e': C3, 'm1_e': C3, 'f1': C4, 'm1_f1': C4, 'dens': C1, 'velx': C2, 'velz': C2, 'eint': '#8a8a86'}
VARN = {'e': 'E', 'm1_e': 'E', 'f1': 'F1', 'm1_f1': 'F1', 'dens': 'rho', 'velx': 'v', 'velz': 'v_phi', 'eint': 'e_int'}


def lad_cart():
    L = parse_lad()
    for case, ttl in (('rw_t10a', 'radiation wave, tau 10 (rw_t10a)'), ('rw_t1000a', 'radiation wave, tau 1000 (rw_t1000a)'),
                      ('rw_p10t100', 'radiation-dominated wave P 10, tau 100 (rw_p10t100)'),
                      ('rw_p10t1', 'radiation-dominated wave P 10, tau 1 (rw_p10t1)')):
        fig, ax = plt.subplots(1, 3, figsize=(12, 4.3))
        txt = []
        for a, kind, xl in ((ax[0], 'X', 'N (finer of the pair)'), (ax[1], 'T', 'dt divisor (finer of the pair)'),
                            (ax[2], 'C', 'N (dt ~ dx, finer of the pair)')):
            for v in ('dens', 'velx', 'e', 'f1'):
                for arm, ls, mk, mf in (('FIN', '-', 'o', None), ('PRODd', '--', 's', 'white')):
                    k = (case, arm, kind, v)
                    if k not in L:
                        continue
                    levs, e = L[k]
                    x = levs[1:len(e) + 1]
                    a.loglog(x, e, mk + ls, color=VARC[v], ms=6 if mf is None else 5, mfc=mf or VARC[v], lw=1.5 if mf is None else 1,
                             alpha=1 if mf is None else 0.75, label=VARN[v] if arm == 'FIN' else None)
                    if arm == 'FIN':
                        txt.append('%s %s %.2f' % (kind, VARN[v], math.log2(e[-2]/e[-1])))
            k0 = (case, 'FIN', kind, 'dens')
            if k0 in L:
                levs, e = L[k0]
                slopes(a, levs[1], e[0]*2.5, levs[1:len(e) + 1], lab=(kind == 'X'))
            a.set_xlabel(xl)
            a.set_title({'X': 'space (X)', 'T': 'time (T, hesdirk2)', 'C': 'combined (C)'}[kind])
        ax[0].set_ylabel('self-convergence L1 difference')
        ax[0].legend(fontsize=8.5)
        ax[2].plot([], [], 's--', color='#8a8a86', mfc='white', label='PRODd (open)')
        ax[2].legend(fontsize=8.5)
        suptitle(fig, 'Ladder: %s, FIN (filled) / PRODd (open)' % ttl, 'binary %s; data fin3_1012 ana/out_281dd686 lad.txt (FIN3) and F2O/lad.txt (FIN2 control)' % BIN_D)
        fig.tight_layout(rect=(0, 0, 1, 0.86))
        fin = '; '.join(t for t in txt if t.split()[1] in ('rho', 'v', 'E', 'F1'))
        save(fig, 'ladder_%s.png' % case, 'ladder_' + case, 'Self-convergence ladder %s (X, T, C), FIN on %s with PRODd open '
             'symbols; finest-pair orders %s.' % (case, BIN_D, fin))


def lad_sp():
    xo = parse_ord(F3SP + '/xord.txt')
    to = parse_ord(F3SP + '/tord.txt')
    cases = [('rsw', 'rsw_T', 'radiation-dominated wave (T9 rsw)', ('dens', 'velx', 'm1_e', 'm1_f1')),
             ('rswt', 'rswt_T', 'wave with thin top + bath (T9t rswt)', ('dens', 'velx', 'm1_e', 'm1_f1')),
             ('rswv', 'rswv_T', 'moving-shell wave (T9v rswv)', ('dens', 'velx', 'm1_e', 'm1_f1')),
             ('rswtv', 'rswtv_T', 'moving thin-top wave (T9tv rswtv)', ('dens', 'velx', 'm1_e', 'm1_f1')),
             ('rlat', 'rlat_T', 'lateral wave, thick (rlat)', ('dens', 'velz', 'm1_e', 'm1_f1')),
             ('rlat_om5', 'rlat_om5_T', 'lateral wave, thick, om5 (rlat_om5)', ('dens', 'velz', 'm1_e', 'm1_f1')),
             ('rlatt', None, 'lateral wave, thin (rlatt)', ('dens', 'velz', 'velx', 'm1_e', 'm1_f1')),
             ('rlatt_om5', 'rlatt_om5_T', 'lateral wave, thin, om5 (rlatt_om5)', ('dens', 'velz', 'm1_e', 'm1_f1')),
             ('marsh_r1_u', 'marsh_r1_T', 'spherical Marshak wave (T5 marsh_r1)', ('dens', 'velx', 'm1_e', 'm1_f1')),
             (None, 'relax1_T', 'thermal relaxation (T1 relax1_T)', ('dens', 'eint', 'm1_e')),
             ('patmT_L16_s', None, 'thermal photosphere (P1T patmT)', ('dens', 'm1_e', 'm1_f1'))]
    lat = ('rlat', 'rlatt', 'rlat_om5', 'rlatt_om5')
    for xc, tc, ttl, vs in cases:
        if xc in ('rswv', 'rlatt_om5'):
            continue   # replotted from OT's runs in lad_sp_ot
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
        txt = []
        for a, src, case, kind in ((ax[0], xo, xc, 'X'), (ax[1], to, tc, 'T')):
            if case is None:
                a.text(0.5, 0.5, 'no %s ladder for this case' % kind, ha='center', va='center', transform=a.transAxes, color=INK2)
                a.set_xticks([])
                a.set_yticks([])
                continue
            for v in vs:
                for arm, ls, mk, mf in (('HFfZ', '-', 'o', None), ('HFfW', '--', 's', 'white')):
                    e = src.get((case, arm, v))
                    if not e:
                        continue
                    if kind == 'X':
                        base = (32 if xc in lat else 64)
                        x = [base*2**(i + 1) for i in range(len(e))]
                    else:
                        x = [2**(i + 1) for i in range(len(e))]
                    a.loglog(x, e, mk + ls, color=VARC.get(v, INK), ms=6 if mf is None else 5, mfc=mf or VARC.get(v, INK),
                             lw=1.5 if mf is None else 1, alpha=1 if mf is None else 0.75, label=VARN[v] if arm == 'HFfZ' else None)
                    if arm == 'HFfZ':
                        txt.append('%s %s %s' % (kind, VARN[v], ' '.join('%.2f' % math.log2(e[i]/e[i + 1]) for i in range(len(e) - 1))))
            e0 = src.get((case, 'HFfZ', 'm1_e')) or src.get((case, 'HFfZ', 'dens'))
            if e0:
                x = [((32 if xc in lat else 64) if kind == 'X' else 1)*2**(i + 1) for i in range(len(e0))]
                slopes(a, x[0], e0[0]*2.5, x, lab=True)
            a.set_xlabel('N_r or N_phi (finer of the pair)' if kind == 'X' else 'dt divisor (finer of the pair)')
            a.set_title('space (X)' if kind == 'X' else 'time (T, hesdirk2)')
            a.set_ylabel('self-convergence L1 diff / mean')
            a.legend(fontsize=8.5)
        key = (xc or tc).replace('_u', '').replace('_L16_s', '')
        suptitle(fig, 'Ladder (sp wedge, vet_gd): %s, FIN HFfZ (filled) / FIN2 control (open)' % ttl,
                 'binary %s; data fin3_1012 ana/out_281dd686 (sp) xord.txt / tord.txt' % BIN_D)
        fig.tight_layout(rect=(0, 0, 1, 0.86))
        note = ' (dt error below the bath floor: errors flat at ~1e-7, NOTE-suiteF T1)' if tc == 'relax1_T' else ''
        save(fig, 'ladder_sp_%s.png' % key, 'ladder_sp_' + key, 'Self-convergence ladder for the %s on the sp wedge, FIN (HFfZ) '
             'on %s; orders %s%s.' % (ttl, BIN_D, '; '.join(txt), note))


def lad_sp_ot():
    """OT (open-1013): rswv with the user-approved deck (FX; FIN3 deck open), rlatt_om5_T incl. level 16 (criterion pairs)"""
    xo = parse_ord(F3SP + '/xord.txt')
    to = parse_ord(F3SP + '/tord.txt')
    ot = parse_ord(F3SP + '/ot_ord.txt')
    def ords(e):
        return ' '.join('%.2f' % math.log2(e[i]/e[i + 1]) for i in range(len(e) - 1))
    def panel(a, kind, series, base):
        for v, e, filled in series:
            if not e:
                continue
            x = [(base if kind == 'X' else 1)*2**(i + 1) for i in range(len(e))]
            col = VARC.get(v, INK)
            a.loglog(x, e, ('o-' if filled else 's--'), color=col, ms=6 if filled else 5, mfc=col if filled else 'white',
                     lw=1.5 if filled else 1, alpha=1 if filled else 0.75, label=VARN[v] if filled else None)
        a.set_xlabel('N_r or N_phi (finer of the pair)' if kind == 'X' else 'dt divisor (finer of the pair)')
        a.set_title('space (X)' if kind == 'X' else 'time (T, hesdirk2)')
        a.set_ylabel('self-convergence L1 diff / mean')
    # ---- rswv
    vs = ('dens', 'velx', 'm1_e', 'm1_f1')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    txt = []
    panel(ax[0], 'X', [(v, ot.get(('rswv', 'FX', v)), True) for v in vs] + [(v, xo.get(('rswv', 'HFfZ', v)), False) for v in vs], 64)
    panel(ax[1], 'T', [(v, ot.get(('rswv_T', 'FX', v)), True) for v in vs] + [(v, to.get(('rswv_T', 'HFfZ', v)), False) for v in vs], 64)
    for kind, case in (('X', 'rswv'), ('T', 'rswv_T')):
        txt.append(kind + ' ' + ', '.join('%s %s' % (VARN[v], ords(ot[(case, 'FX', v)])) for v in vs))
    e0 = ot[('rswv', 'FX', 'm1_e')]
    slopes(ax[0], 128, e0[0]*2.5, [128*2**i for i in range(len(e0))], lab=True)
    e0 = ot[('rswv_T', 'FX', 'm1_e')]
    slopes(ax[1], 2, e0[0]*2.5, [2**(i + 1) for i in range(len(e0))], lab=True)
    for a in ax:
        a.plot([], [], 's--', color=INK2, mfc='white', label='FIN3 deck (before)')
        lo, hi = a.get_ylim()
        a.set_ylim(lo, hi*40)
        a.legend(fontsize=8, ncol=3, loc='upper right')
    suptitle(fig, 'Ladder (sp wedge, vet_gd): moving-shell wave (T9v rswv), FIN3 with the user-approved rswv deck (filled) / FIN3 deck (open)',
             'deck: implicit_f0_init = zero + ox1_bc = outflow (rswv decks only, code defaults unchanged); binary %s; OT open-1013 runs' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, 'ladder_sp_rswv.png', 'ladder_sp_rswv', 'Self-convergence ladder for the moving-shell wave (T9v rswv) on the sp wedge, '
         'FIN3 on %s with the user-approved rswv deck change (implicit_f0_init = zero + outflow gas top; open symbols: the FIN3 deck, '
         'X 0.13-0.70 from the wall layer): X TEST-fixed, T OPEN. Orders %s (X N 64-1024; T levels 1-32, the late v / E / F1 drop is the '
         'per-step gas-radiation O(dt) term of rw_p10t1).' % (BIN_D, '; '.join(txt)))
    # ---- rlatt_om5
    vx = ('dens', 'velz', 'm1_e', 'm1_f1')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    panel(ax[0], 'X', [(v, xo.get(('rlatt_om5', 'HFfZ', v)), True) for v in vx] + [(v, xo.get(('rlatt_om5', 'HFfW', v)), False) for v in vx], 32)
    panel(ax[1], 'T', [(v, ot.get(('rlatt_om5_T', 'B', v)), True) for v in vx] + [(v, to.get(('rlatt_om5_T', 'HFfW', v)), False) for v in vx], 32)
    ax[1].axvspan(6, 20, color=INK2, alpha=0.10, lw=0)
    ax[1].text(10.5, 0.03, 'criterion pairs\n(levels 4-16)', transform=ax[1].get_xaxis_transform(), ha='center', va='bottom', fontsize=8.5, color=INK2)
    e0 = xo.get(('rlatt_om5', 'HFfZ', 'm1_e'))
    if e0:
        slopes(ax[0], 64, e0[0]*2.5, [64*2**i for i in range(len(e0))], lab=True)
    e0 = ot[('rlatt_om5_T', 'B', 'm1_e')]
    slopes(ax[1], 2, e0[0]*2.5, [2**(i + 1) for i in range(len(e0))], lab=True)
    for a in ax:
        a.plot([], [], 's--', color=INK2, mfc='white', label='FIN2 control')
        lo, hi = a.get_ylim()
        a.set_ylim(lo, hi*40)
        a.legend(fontsize=8, ncol=3, loc='upper right')
    tx = 'X ' + ', '.join('%s %s' % (VARN[v], ords(xo[('rlatt_om5', 'HFfZ', v)])) for v in vx if xo.get(('rlatt_om5', 'HFfZ', v)))
    tt = 'T ' + ', '.join('%s %s' % (VARN.get(v, 'F2'), ords(ot[('rlatt_om5_T', 'B', v)])) for v in vx + ('m1_f2',))
    suptitle(fig, 'Ladder (sp wedge, vet_gd): lateral wave, thin, om5 (rlatt_om5), FIN3 HFfZ (filled) / FIN2 control (open)',
             'T to level 16 (OT open-1013, FIN3 baseline); E / F1 judged on the finer pairs (levels 4-16); OPEN: v_phi 1.61 / 1.49, rho 1.73; binary %s' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    save(fig, 'ladder_sp_rlatt_om5.png', 'ladder_sp_rlatt_om5', 'Self-convergence ladder for the thin lateral wave (rlatt_om5) on the '
         'sp wedge, FIN3 on %s, T to level 16: rlatt_om5_T OPEN (user decision): E and F1 count as converged on the finer pairs '
         '(levels 4-16, user-approved criterion; the 1.3 was the coarsest pair), but v_phi is 1.61 / 1.49 on the finer pairs and rho 1.73 '
         'at the finest pair (OT is checking whether v_phi is the same per-step conservative2 velocity term); rlatt X F1 / v_r floor '
         'stays OPEN. Orders %s; %s.' % (BIN_D, tx, tt))


def lad_sp_exact():
    """ladders against exact solutions: dif_r1, thin / core (T6), patm bands, rswt_q"""
    ex = open(F3SP + '/exact.txt').read()
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.3))
    # dif_r1
    rows = re.findall(r'dif_r1 HFfZ\s+(\d+) Picard [0-9.]+ \| (.*)', ex)
    ns = [int(r[0]) for r in rows]
    vals = [[float(v) for v in re.findall(r'L1 ([0-9.e+-]+)', r[1])] for r in rows]
    for j in range(3):
        ax[0].loglog(ns, [v[j] for v in vals], 'o-', color=LEVC[j + 1], ms=7, label='output %d' % (j + 1))
    slopes(ax[0], ns[0], vals[0][0]*1.5, ns)
    ax[0].set_title('T4 spherical diffusion vs exact')
    ax[0].set_xlabel('N_r')
    ax[0].set_ylabel('L1 relative error')
    ax[0].legend(fontsize=8.5)
    # thin / core
    for case, col, nm in (('thin_L16_u', C1, 'T6a thin sphere E (16x16)'), ('thin_u', C2, 'T6a thin sphere E (4x4)')):
        r = re.findall(r'%s HFfZ\s+(\d+) .*?E/E_ex-1 L1 ([0-9.e+-]+) .*?max\|Fr2-1\| ([0-9.e+-]+)' % case, ex)
        ax[1].loglog([int(x[0]) for x in r], [float(x[1]) for x in r], 'o-', color=col, ms=7, label=nm)
        if case == 'thin_L16_u':
            ax[1].loglog([int(x[0]) for x in r], [float(x[2]) for x in r], 's--', color=C3, ms=6, label='T6a max|F r^2 - 1|')
    slopes(ax[1], 64, 3e-2, (64, 128, 256))
    ax[1].set_title('T6a vs exact (FIN: E at the 4e-6 floor)')
    ax[1].set_xlabel('N_r')
    ax[1].legend(fontsize=8, loc='center left')
    # rswt_q equilibrium
    rq = re.findall(r'HFfZ (\d+) max\|E/0.03-1\| ([0-9.e+-]+) max\|rho-1\| ([0-9.e+-]+)', open(F3SP + '/rswtq.txt').read())
    ax[2].loglog([int(x[0]) for x in rq], [float(x[1]) for x in rq], 'o-', color=C3, ms=7, label='max|E/E_b - 1|')
    ax[2].loglog([int(x[0]) for x in rq], [float(x[2]) for x in rq], 's-', color=C1, ms=7, label='max|rho - 1|')
    slopes(ax[2], 64, float(rq[0][1])*1.5, (64, 128, 256))
    ax[2].set_title('rswt_q: equilibrium with a top bath')
    ax[2].set_xlabel('N_r')
    ax[2].legend(fontsize=8.5)
    suptitle(fig, 'Ladders vs exact solutions on the sp wedge, FIN (HFfZ)', 'binary %s; data fin3_1012 ana/out_281dd686 (sp) exact.txt, rswtq.txt' % BIN_D)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    o = [math.log2(vals[i][0]/vals[i + 1][0]) for i in range(len(vals) - 1)]
    save(fig, 'ladder_sp_exact.png', 'ladder_sp_exact', 'Convergence vs exact solutions on the sp wedge, FIN (HFfZ) on %s: '
         'T4 diffusion orders %s; T6a thin-sphere E does not converge at fixed lateral resolution (F r^2 does); rswt_q '
         'equilibrium error %s.' % (BIN_D, ' '.join('%.2f' % x for x in o), ' -> '.join(x[1] for x in rq)))


def annotate(F):
    """order / names from the parallel worker's idort_results.json (same test order)"""
    try:
        R = json.load(open(P + '/idort_results.json'))
    except Exception:
        return
    byord = {t['order']: t for t in R['tests']}
    extra = {'25_patm_sp': 'P1 patm L16 poly', '26_patmT_sp': 'P1T patmT (coupled gas)', '27_pinc_sp': 'P2 pinc thermal channel beam, c dt/dx 300',
             '28_pulse_lr_sp': 'L pulse r0/w 1 / 1000 vs Cartesian', '29_rsw_lr_sp': 'L rsw / rswt r0/w 1000 vs Cartesian',
             '30_g5': 'G5 2-D pulse, kappa 12.8 ... 12800'}
    for k, v in F.items():
        v['spherical'] = k.endswith('_sp') or k.startswith('ladder_sp')
        if k[:2].isdigit() and int(k[:2]) in byord:
            t = byord[int(k[:2])]
            v['order'] = t['order']
            v['idort_ref'] = t['idort_ref']
            v['name'] = t['name']
        elif k in extra:
            v['order'] = int(k[:2])
            v['results_json_test_id'] = extra[k]
            v['results_json_section'] = 'extra_tests' if k == '30_g5' else 'spherical_extra'
        elif k.startswith('ladder'):
            v['results_json_section'] = 'extra_tests'
        if v.get('file'):
            v['arm'] = 'FIN3 final / sp HFfZ + FIN3 keys (fin3-1012 281dd686)'


ALL = ['f01', 'f02', 'f03', 'f04', 'f05', 'f05sp', 'f06', 'f06sp', 'f07', 'f07sp', 'f08', 'f09', 'f09sp', 'f10', 'f11', 'f12',
       'f13', 'f13sp', 'f14', 'f14sp', 'f15', 'f16', 'f16sp', 'f17', 'f17sp', 'f18', 'f19', 'f20', 'f20sp', 'f21', 'f22_24',
       'f25', 'f27', 'f28', 'f29', 'f30', 'lad_cart', 'lad_sp', 'lad_sp_ot', 'lad_sp_exact']

if __name__ == '__main__':
    os.makedirs(FIG, exist_ok=True)
    jf = os.path.join(FIG, 'figures.json')
    old = json.load(open(jf)) if os.path.exists(jf) else {}
    fails = []
    for n in (sys.argv[1:] or ALL):
        try:
            globals()[n]()
        except Exception:
            traceback.print_exc()
            fails.append(n)
    old.update(FIGS)
    annotate(old)
    json.dump(dict(sorted(old.items())), open(jf, 'w'), indent=1)
    print('FAILED:', fails)
