"""docfig.py: documentation figures and movie frames per star (NOTE 26o item 3, 26p). usage: docfig.py STAR fig | frames | limits | strip
maps: v_r cmasher fusion, T_eff cmocean ice, delta rho/<rho> cmocean balance. Reduced dumps only (pages_c1_1124 caches); read only."""
import sys, os, glob, json, re
import numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/wt_note_a0dc/docs/handover/plot_style')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, idort_style
import cmasher as cmr, cmocean
from matplotlib.ticker import MaxNLocator
from copy import copy
idort_style.use(font=12); plt.rcParams.update({'axes.linewidth': 1.0, 'figure.dpi': 100, 'savefig.dpi': 100, 'savefig.bbox': None, 'savefig.format': 'png', 'axes.labelsize': 12, 'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5, 'xtick.major.size': 4, 'ytick.major.size': 4, 'xtick.minor.size': 2, 'ytick.minor.size': 2, 'figure.constrained_layout.use': False})
P = '/resnick/groups/carnegie_poc/jingze/pages_c1_1124/'; O = '/resnick/groups/carnegie_poc/jingze/docfig_1125/'; RS = 6.957e10; SIG = 5.670374419e-5; GREY = '0.78'
# star: Rph_IC [cm], v_r shell tag, slice range [R_ph,IC], density floor
C = dict(m35mms=dict(rph=9.554e11, tag='1000', rr=(0.86, 1.06), fl=1e-17), m35zams=dict(rph=5.4972e11, tag='1000', rr=(0.95, 1.02), fl=1e-17, z=(0.972, 0.998)), t9l62=dict(rph=6.5326e12, tag='3800', rr=(0.42, 7.0), fl=3.6e-19), t9mlt=dict(rph=3.9282e12, tag='098', rr=(0.62, 1.12), fl=3.6e-19), bsg=dict(rph=3.4e12, tag=None, rr=(34., 54.), fl=6.15245e-14))
def cm_(c):
    c = copy(c); c.set_bad(GREY); return c
BAL, ICE, FUS = cmocean.cm.balance, cmocean.cm.ice, cmr.fusion
def files(s):
    if s == 'bsg': return sorted(glob.glob(P + 'bsg/pipe/cache/view/red_*.npz'))
    return sorted(glob.glob(P + s + '/out/data/red_*.npz'))
def load(s, f):
    z = np.load(f); d = {}
    if s == 'bsg':
        n = re.search(r'red_(\d+)', f).group(1); q = np.load(P + 'bsg/pipe/cache/view/int2_%s.npz' % n)
        d['t'] = float(z['time']) / 86400; d['vr'] = z['s36_vr'].T / 1e5; d['te'] = np.where(q['fr_own'].T > 0, (np.maximum(q['fr_own'].T, 1e-99) / SIG) ** 0.25, np.nan)     # [theta, phi]
        d['x1'] = np.degrees(z['x3f']); d['x2'] = 90 - np.degrees(z['x2f']); d['rf'] = z['x1f'] / RS; d['ang'] = z['x2f']; rho = z['p0_rho']; d['sl'] = np.ma.masked_where(rho <= 1.5 * C[s]['fl'], rho / z['rho_sh'] - 1); d['svr'] = np.ma.masked_array(z['p0_vr'] / 1e5, mask=np.ma.getmaskarray(d['sl'])); d['lab'] = 'r = 40.7 R_sun'; d['slk'] = 'theta'
    else:
        c = C[s]; d['t'] = float(z['t']) / 86400; d['vr'] = z['sh_vr_' + c['tag']][0].T; d['te'] = z['sh_Teff'][0].T; d['x1'] = np.degrees(z['ph_f']); d['x2'] = 90 - np.degrees(z['th_f']); rho = 10 ** z['eqa_lrho']; d['sl'] = np.ma.masked_where((rho <= 1.5 * c['fl']) | z['eqa_Tfl'].astype(bool) * 0, rho / z['rho'][:rho.shape[1]] - 1); d['rf'] = z['rf'][:rho.shape[1] + 1] / RS; d['ang'] = z['eqa_phi']; d['slk'] = 'phi'; d['rph_med'] = float(z['Rph_med'])
    if s != 'bsg': d['svr'] = np.ma.masked_array(z['eqa_vr'].astype(float), mask=np.ma.getmaskarray(d['sl']))   # 27zb: v_r [km/s] on the same slice, floor gas masked as the density slice
    return ph23(s, f, d)
def ph23(s, f, d):
    # 26z: v_r at tau_R = 2/3 of each column (p23.py), one definition for all stars; bsg also takes its T_eff from the same point
    n = re.search(r'red_(\d+)', f).group(1); q = O + 'p23/%s_%05d.npz' % (s, int(n)); z = np.load(q); assert abs(float(z['t']) - d['t']) < 1e-3, (s, n, float(z['t']), d['t'])
    d['vr'] = z['v'].T.astype(float); d['nnone'] = int(z['nnone']); d['rph23'] = float(np.nanmedian(z['r'])) / RS
    if s == 'bsg': d['te'] = z['te'].T.astype(float)
    else: d['te_d'] = float(np.nanmax(abs(z['te'].T - d['te']) / d['te']))
    return d
def srange(s, d):
    # the radial cells drawn in the slice panel (as draw())
    rr = C[s]['rr']; sc = 1 if s == 'bsg' else C[s]['rph'] / RS; rf = np.asarray(d['rf']); i0 = int(np.searchsorted(rf, rr[0] * sc)); i1 = min(len(rf) - 1, int(np.searchsorted(rf, rr[1] * sc))); return slice(i0, i1)
def zone(s, d):
    """slice cells below 0.97 R_ph (median, this dump) and inside the slice range"""
    sc = 1 if s == 'bsg' else C[s]['rph'] / RS; rm = 0.5 * (np.asarray(d['rf'][1:]) + np.asarray(d['rf'][:-1])); top = 51.0 if s == 'bsg' else 0.97 * d['rph_med'] / RS
    if 'z' in C[s]: return (rm >= C[s]['z'][0] * sc) & (rm <= C[s]['z'][1] * sc)
    return (rm >= C[s]['rr'][0] * sc) & (rm <= min(top, C[s]['rr'][1] * sc))
TEFIX = dict(t9l62=(4900., 10000.), t9mlt=(9000., 29500.))   # USER 11-25: brighter surfaces; same as the 3-D movies
def lims(s, pick=None):
    f = files(s); idx = pick or list(range(0, len(f), max(1, len(f) // 40))) + [len(f) - 1]
    V, T, S, SV = [], [], [], []
    for i in idx:
        d = load(s, f[i])
        if d['t'] < 0.05: continue
        V.append(np.nanpercentile(abs(d['vr']), 99.5)); T.append(np.nanpercentile(d['te'], [2, 98])); ok = (d['rf'][1:] / (C[s]['rph'] / RS) > 0) if s != 'bsg' else True
        rr = zone(s, d); q = np.ma.masked_invalid(d['sl'][:, rr]); S.append(np.nanpercentile(abs(q.compressed()), 99)); SV.append(np.abs(np.ma.masked_invalid(d['svr'][:, srange(s, d)]).compressed()))
    L = dict(v=float(np.percentile(V, 80)), te=[float(np.percentile([x[0] for x in T], 10)), float(np.percentile([x[1] for x in T], 90))], s=float(min(1.0, np.percentile(S, 60) * 0.8)))
    L['sv'] = float(np.percentile(np.concatenate(SV), 99))
    if s in TEFIX: L['te'] = list(TEFIX[s])
    return L
BARS = dict(bsg=10., m35mms=1., m35zams=0.1, t9l62=100., t9mlt=10.)                       # R_sun; lateral scale bar at the median photosphere radius of the dump
import matplotlib.patheffects as pe
def sbar(a, s, d):
    # horizontal bar of a round length [R_sun] at the photosphere: length = R_ph x angle (exact at the equator and for R_ph = the dump's median; the angle is L / R_ph)
    L = BARS[s]; deg = np.degrees(L / d['rph23']); W = d['x1'][-1] - d['x1'][0]; x0 = d['x1'][0] + max(0.05 * W, 0.12 * W - deg / 2); y0 = d['x2'][-1] + 0.07 * (d['x2'][0] - d['x2'][-1]); w = [pe.withStroke(linewidth=4.5, foreground='black')]
    if deg > 0.9 * (d['x1'][-1] - d['x1'][0]): return
    a.plot([x0, x0 + deg], [y0, y0], color='white', lw=2.2, solid_capstyle='butt', path_effects=w, zorder=9)
    a.text(x0 + deg / 2, y0 + 0.025 * (d['x2'][0] - d['x2'][-1]), r'$\rm %g\ R_\odot$' % L, color='white', fontsize=9, ha='center', va='bottom', path_effects=[pe.withStroke(linewidth=2.5, foreground='black')], zorder=9)
def draw(s, d, L, fn, stamp=True):
    c = C[s]; t9 = s.startswith('t9'); wide = t9
    fig = plt.figure(figsize=(12, 3.9 if wide else 4.5)); gs = fig.add_gridspec(1, 3, left=0.062 if wide else 0.05, right=0.985, top=0.90, bottom=0.10, wspace=0.28, width_ratios=[1, 1, 1.1 if not wide else 1.0])
    for j, (f, cm, lo, hi, lab, ext) in enumerate(((d['vr'], FUS, -L['v'], L['v'], r'$\rm v_r(\tau_R=2/3)\ [km\ s^{-1}]$', 'both'), (d['te'] / 1e3, ICE, L['te'][0] / 1e3, L['te'][1] / 1e3, r'$\rm T_{eff}\ [10^3\ K]$', 'both'))):
        a = fig.add_subplot(gs[0, j]); pc = a.pcolormesh(d['x1'], d['x2'], np.ma.masked_invalid(f), cmap=cm_(cm), vmin=lo, vmax=hi, shading='flat', rasterized=True); a.set_aspect('equal'); a.set_xlabel(r'$\rm \phi\ [deg]$'); a.set_ylabel(r'$\rm latitude\ [deg]$'); a.xaxis.set_major_locator(MaxNLocator(5)); a.yaxis.set_major_locator(MaxNLocator(5))
        cb = fig.colorbar(pc, ax=a, orientation='horizontal', shrink=0.85, pad=0.16, extend=ext); cb.set_label(lab); cb.locator = MaxNLocator(5); cb.update_ticks()
        sbar(a, s, d)
        if j == 0 and stamp: a.text(0.03, 1.02, 't = %.2f d' % d['t'], transform=a.transAxes, fontsize=13, color='black')
    a = fig.add_subplot(gs[0, 2]); rr = C[s]['rr']; sc = 1 if s == 'bsg' else c['rph'] / RS; rf = np.asarray(d['rf']); i0 = int(np.searchsorted(rf, rr[0] * sc)); i1 = min(len(rf) - 1, int(np.searchsorted(rf, rr[1] * sc)))
    if s == 'bsg': A = d['ang']; R, TH = np.meshgrid(rf[i0:i1 + 1], A); X, Y = R * np.sin(TH), R * np.cos(TH); xl, yl = r'$\rm x\ [R_\odot]$', r'$\rm z\ [R_\odot]$'
    else: A = d['ang']; R, TH = np.meshgrid(rf[i0:i1 + 1], A); X, Y = R * np.cos(TH), R * np.sin(TH); xl, yl = r'$\rm x\ [R_\odot]$', r'$\rm y\ [R_\odot]$'
    pc = a.pcolormesh(Y if s == 'bsg' else X, X if s == 'bsg' else Y, d['svr'][:, i0:i1], cmap=cm_(FUS), vmin=-L['sv'], vmax=L['sv'], shading='flat', rasterized=True); a.set_aspect('equal'); a.set_xlabel(xl if s != 'bsg' else r'$\rm z\ [R_\odot]$'); a.set_ylabel(yl if s != 'bsg' else r'$\rm x\ [R_\odot]$'); a.xaxis.set_major_locator(MaxNLocator(5)); a.yaxis.set_major_locator(MaxNLocator(5))
    cb = fig.colorbar(pc, ax=a, orientation='horizontal', shrink=0.85, pad=0.16, extend='both'); cb.set_label(r'$\rm v_r\ on\ the\ slice\ [km\ s^{-1}]$'); cb.locator = MaxNLocator(5); cb.update_ticks()
    fig.savefig(fn); plt.close(fig)
def frames(s):
    L = json.load(open(O + 'work/lim_%s.json' % s)); F = files(s); D = O + 'frames23/' + s; os.makedirs(D, exist_ok=True); n = len(F)
    from multiprocessing import Pool
    with Pool(int(os.environ.get('NPROC', '6'))) as p: p.map(fr1, [(s, F[i], i, L, D) for i in range(n)], chunksize=2)
    print(s, 'frames', len(glob.glob(D + '/*.png')))
def fr1(a):
    s, f, i, L, D = a; fn = D + '/f_%05d.png' % i
    if not os.path.exists(fn): draw(s, load(s, f), L, fn)
def strip():
    S = ['bsg', 'm35mms', 'm35zams', 't9l62', 't9mlt']; N = dict(bsg='M20 BSG', m35mms='M35 MMS', m35zams='M35 ZAMS', t9l62='M56 LBV (radiative IC)', t9mlt='M56 LBV')
    fig, ax = plt.subplots(2, 5, figsize=(15, 6.4)); fig.subplots_adjust(left=0.045, right=0.985, top=0.95, bottom=0.10, wspace=0.34, hspace=0.34)
    for j, s in enumerate(S):
        d = load(s, files(s)[-1]); L = lims(s, [len(files(s)) - 1])
        for i, (f, cm, lo, hi, lab) in enumerate(((d['vr'], FUS, -L['v'], L['v'], r'$\rm v_r(\tau_R=2/3)\ [km\ s^{-1}]$'), (d['te'] / 1e3, ICE, L['te'][0] / 1e3, L['te'][1] / 1e3, r'$\rm T_{eff}\ [10^3\ K]$'))):
            a = ax[i, j]; pc = a.pcolormesh(d['x1'], d['x2'], np.ma.masked_invalid(f), cmap=cm_(cm), vmin=lo, vmax=hi, shading='flat', rasterized=True); a.set_box_aspect(0.5 if s.startswith('t9') else 1.0); a.xaxis.set_major_locator(MaxNLocator(3)); a.yaxis.set_major_locator(MaxNLocator(3)); a.tick_params(labelsize=9); a.set_xlabel(r'$\rm \phi\ [deg]$', fontsize=10, labelpad=1); a.set_ylabel(r'$\rm latitude\ [deg]$', fontsize=10, labelpad=1)
            sbar(a, s, d); cb = fig.colorbar(pc, ax=a, orientation='horizontal', pad=0.24, shrink=0.9, extend='both'); cb.locator = MaxNLocator(3); cb.update_ticks(); cb.ax.tick_params(labelsize=9); cb.set_label(lab, fontsize=10)
            if i == 0: a.set_title('%s, t = %.1f d' % (N[s], d['t']), fontsize=11)
    fig.savefig(O + 'out/doc_strip.png'); plt.close(fig); print('strip', os.path.getsize(O + 'out/doc_strip.png'))
if __name__ == '__main__':
    s, m = sys.argv[1], sys.argv[2]
    Lf = O + 'work/lim_%s.json' % s
    if m == 'limits':
        _o = json.load(open(O + 'work_prev26z/lim_%s.json' % s)); _l = lims(s); _l['s'] = _o['s']
        if s != 'bsg' and s not in TEFIX: _l['te'] = _o['te']      # 26z: T_eff and density limits unchanged (bsg T_eff moves to tau 2/3)
        if False: _l['te'] = [6500., 12500.]   # the star cools from 20 kK to 7 kK: fixed limits fitted to the later part; early frames are saturated
        json.dump(_l, open(Lf, 'w'))
        print(s, json.load(open(Lf)), len(files(s)))
    if m == 'fig':
        d = load(s, files(s)[-1]); L = lims(s, [len(files(s)) - 1]); draw(s, d, L, O + 'out/doc_%s.png' % s); print(s, 'fig', 't = %.2f d' % d['t'], os.path.getsize(O + 'out/doc_%s.png' % s))
    if m == 'frames': frames(s)
    if m == 'strip': strip()
