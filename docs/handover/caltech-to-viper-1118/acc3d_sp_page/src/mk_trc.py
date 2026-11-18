#!/usr/bin/env python3
"""accsp_prod_page trc-1102 (USER): STREAM-TRACER section of the sp production page (READ ONLY on the run).
Passive scalars (<hydro>/nscalars = 2, pgen ry_per_accretor acctrc-1019): s0 = gas written by the nozzle (the stream), s1 = gas
that entered through the other user boundaries.  ryper.trc.hst (static run, then the stream run; cumulative registers continue
across the restart): Ms0, Ms1 = domain integrals of rho s; Nz0 = s0 mass written by the nozzle reset, Nz1 = s1 mass removed by
it; Sk0, Sk1 = removed by the sinks; Fi0, Fi1 = out through r_in; Fo0, Fo1 = out through r_out + theta walls (outward
positive); Bin = mass in through the non-nozzle boundaries (the s1 source).  Budget (pgen comment):
    Ms0 - Ms0(t0) = dNz0 - dSk0 - dFi0 - dFo0,      Ms1 - Ms1(t0) = - dNz1 - dSk1 - dFi1 - dFo1  (the s1 inflow is negative Fo1;
    Bin = the total mass that entered through those boundaries, shown for comparison).
Figures (img/): sp_trc_budget.png, sp_trc_maps.jpg (latest dump: orbital + x-z plane and the photosphere map of s0),
sp_trc_prof.png (s0 shell mass vs r for the stream dumps); img/trc.json = numbers + one-line takeaways for build_page.py.
Inputs: the red_*.npz of reduce_sp.py VER >= 2 (eq_s0, me_s0, map_s0ph, prof_*), so only new dumps cost anything."""
import os
import glob
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

E = os.environ.get
J = '/resnick/groups/carnegie_poc/jingze'
OUT = E('OUT', J + '/accsp_prod_page')
ST = E('STATIC', J + '/ringsmooth_1102/long/L_gz')
SR = E('STREAM', J + '/accsp_1102/prod/fresh/stream')
IMG, DATA = OUT + '/img', OUT + '/data'
TS, DAYS, RIN, DL1 = 0.01, 8.0521, 6.3, 16.2266
C1, C2, C3, C4, C5 = '#4cc3ff', '#ffb347', '#7dff8a', '#ff5c7a', '#c77dff'
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8', 'ytick.color': '#a3adb8',
                     'axes.edgecolor': '#3a424c', 'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
                     'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'})
res = dict(present=False)


def readhst(f):
    if not os.path.exists(f):
        return None
    lab = None
    for ln in open(f):
        if ln.startswith('#  ['):
            lab = [x.split('=')[1] for x in ln[1:].split()]
    a = np.atleast_2d(np.genfromtxt(f, comments='#', invalid_raise=False))
    if lab is None or a.size == 0:
        return None
    return {k: a[:, i] for i, k in enumerate(lab) if i < a.shape[1]}


def save():
    json.dump(res, open(IMG + '/trc.json', 'w'), indent=1)


# ------------------------------------------------------------------ budget from trc.hst (stream run; static run before TS)
Ta, Tb = readhst(ST + '/ryper.trc.hst'), readhst(SR + '/ryper.trc.hst')
if Tb is None:
    res['why'] = 'no ryper.trc.hst in the stream run yet'; save(); raise SystemExit(0)
if Ta is not None:
    keep = Ta['time'] < Tb['time'][0] - 1e-12
    T = {k: np.concatenate([Ta[k][keep], Tb[k]]) for k in Tb}
else:
    T = Tb
o = np.argsort(T['time'], kind='stable'); T = {k: v[o] for k, v in T.items()}
t = T['time']; i0 = int(np.argmin(abs(t - max(TS, t[0]))))
d = {k: T[k] - T[k][i0] for k in ('Ms0', 'Ms1', 'Nz0', 'Nz1', 'Sk0', 'Sk1', 'Fi0', 'Fi1', 'Fo0', 'Fo1', 'Bin')}
pred0 = d['Nz0'] - d['Sk0'] - d['Fi0'] - d['Fo0']
pred1 = - d['Nz1'] - d['Sk1'] - d['Fi1'] - d['Fo1']
inj = np.maximum(d['Nz0'], 1e-300)
s = t >= t[i0]
th = (t - TS) * DAYS * 24
fig, ax = plt.subplots(1, 3, figsize=(18, 4.8))
a = ax[0]
for k, c, lab in (('Nz0', C2, 'injected by the nozzle (Nz0)'), ('Ms0', C1, 'in the domain (Ms0)'), ('Sk0', C4, 'removed by the sinks (Sk0)'),
                  ('Fi0', C5, 'out through r_in (Fi0)'), ('Fo0', C3, 'out through r_out + walls (Fo0)')):
    a.plot(th[s], d[k][s], color=c, lw=1.2, label=lab)
a.set_xlabel('hours since stream on'); a.set_ylabel(r'$\Delta$ mass since stream on [code]'); a.set_title('stream tracer s0 budget'); a.legend(); a.grid(alpha=.3)
a = ax[1]
a.plot(th[s], (d['Ms0'][s] - pred0[s]) / inj[s], color=C1, lw=1, label=r's0: ($\Delta$Ms0 - booked) / injected')
a.plot(th[s], (d['Ms1'][s] - pred1[s]) / np.maximum(np.abs(d['Bin'][s]), 1e-300), color=C2, lw=1, label=r's1: ($\Delta$Ms1 - booked) / Bin')
a.axhline(0, color='#a3adb8', lw=.6, ls=':')
a.set_xlabel('hours since stream on'); a.set_title('budget closure (0 = every gram booked)'); a.legend(); a.grid(alpha=.3)
a = ax[2]
for k, c, lab in (('Bin', C2, 'entered via other boundaries (Bin)'), ('Ms1', C1, 'in the domain (Ms1)'), ('Sk1', C4, 'sinks (Sk1)'),
                  ('Fi1', C5, 'out through r_in (Fi1)'), ('Fo1', C3, 'out through r_out + walls (Fo1; < 0 = inflow)'), ('Nz1', '#9aa4ae', 'removed by the nozzle reset (Nz1)')):
    a.plot(th[s], d[k][s], color=c, lw=1.2, label=lab)
a.set_xlabel('hours since stream on'); a.set_title('boundary-inflow tracer s1'); a.legend(); a.grid(alpha=.3)
fig.tight_layout(); fig.savefig(IMG + '/sp_trc_budget.png', dpi=80); plt.close(fig)
j = -1
res.update(present=True, t=float(t[j]), inj=float(d['Nz0'][j]), dom=float(d['Ms0'][j]), sink=float(d['Sk0'][j]), rin=float(d['Fi0'][j]),
           rout=float(d['Fo0'][j]), res0=float(np.nanmax(np.abs((d['Ms0'][s] - pred0[s]) / inj[s])[th[s] > 0.5])) if (th[s] > 0.5).any() else None,
           bin=float(d['Bin'][j]), ms1=float(d['Ms1'][j]), s0max=float(T['s0max'][j]), s0min=float(T['s0min'][j]))
f_ = lambda x: 100 * x / res['inj'] if res['inj'] > 0 else float('nan')  # noqa: E731
res['take_budget'] = ('Of the stream gas injected so far (%.3g code), %.0f %% is still in the domain, %.0f %% went into the sinks, '
                      '%.0f %% through r_in and %.0f %% out through r_out/walls; the budget closes to %s of the injected mass.'
                      % (res['inj'], f_(res['dom']), f_(res['sink']), f_(res['rin']), f_(res['rout']),
                         ('%.1g' % res['res0']) if res['res0'] is not None else '–'))

# ------------------------------------------------------------------ dumps (red_*.npz VER >= 2 with s0)
R = []
for f in sorted(glob.glob(DATA + '/red_*.npz')):
    z = np.load(f)
    if 'prof_ms0' in z.files and int(z['ver']) >= 2:
        R.append(z)
R.sort(key=lambda z: float(z['t']))
Rs = [z for z in R if float(z['t']) > TS + 1e-9]
if Rs:
    L = Rs[-1]; tl = float(L['t'])
    # maps: orbital plane, x-z plane, photosphere
    fig = plt.figure(figsize=(21, 6.4)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.55])
    nm = LogNorm(1e-4, 1); cm = plt.get_cmap('viridis').copy(); cm.set_under('#12161b'); cm.set_bad('#12161b')
    lim = 16.5
    a = fig.add_subplot(gs[0])
    PH, RR = np.meshgrid(L['phf'], L['rf'], indexing='ij')
    pc = a.pcolormesh(RR * np.cos(PH), RR * np.sin(PH), np.maximum(L['eq_s0'], 1e-6), cmap=cm, norm=nm, shading='flat', rasterized=True)
    a.set_title('orbital plane: stream-gas fraction s0', fontsize=10)
    a2 = fig.add_subplot(gs[1])
    for h in (0, 1):
        TH, RR = np.meshgrid(L['thf'], L['rf'], indexing='ij')
        a2.pcolormesh((1 if h == 0 else -1) * RR * np.sin(TH), RR * np.cos(TH), np.maximum(L['me_s0'][h], 1e-6), cmap=cm, norm=nm,
                      shading='flat', rasterized=True)
    a2.set_title('meridional x-z plane (y = 0): s0', fontsize=10)
    for a_, yl in ((a, 'y'), (a2, 'z')):
        a_.add_patch(plt.Circle((0, 0), RIN, color='#3a424c')); a_.add_patch(plt.Circle((0, 0), DL1, fill=False, ec=C2, ls='--', lw=.7))
        a_.set_xlim(-lim, lim); a_.set_ylim(-lim, lim); a_.set_aspect('equal'); a_.set_xlabel(r'x [R$_\odot$]'); a_.set_ylabel(yl + r' [R$_\odot$]')
    a3 = fig.add_subplot(gs[2])
    q = L['map_s0ph']; NP = q.shape[0]; qq = np.roll(q, -NP // 2, axis=0)
    pc3 = a3.pcolormesh(np.degrees(L['phf']) - 180.0, 90.0 - np.degrees(L['thf']), np.maximum(qq.T, 1e-6), cmap=cm, norm=nm, shading='flat', rasterized=True)
    a3.set_xlabel(r'longitude $\phi$ [deg] (0 = towards the donor / L1)'); a3.set_ylabel('latitude [deg]')
    a3.set_title('s0 at the Rosseland photosphere cell (tau = 2/3)', fontsize=10); a3.set_xticks(range(-180, 181, 60)); a3.set_yticks(range(-90, 91, 30))
    fig.colorbar(pc, ax=[a, a2], shrink=.8, label='s0 (log)', extend='min'); fig.colorbar(pc3, ax=a3, shrink=.8, label='s0 (log)', extend='min')
    fig.suptitle('stream tracer s0 at t = %.4f code (stream on for %.1f h), dump %s' % (tl, (tl - TS) * DAYS * 24, str(L['idx'])), color='#e8ecf1')
    fig.savefig(IMG + '/sp_trc_maps.jpg', dpi=72, pil_kwargs={'quality': 88}); plt.close(fig)
    # radial profiles
    r = np.asarray(L['r']); rf = np.asarray(L['rf'])
    fig, ax = plt.subplots(1, 2, figsize=(15, 4.8))
    sel = Rs if len(Rs) <= 10 else [Rs[i] for i in np.unique(np.linspace(0, len(Rs) - 1, 10).astype(int))]
    cmap = plt.get_cmap('viridis')
    for n, z in enumerate(sel):
        c = cmap(n / max(len(sel) - 1, 1)); ms = np.asarray(z['prof_ms0']); dr = np.diff(np.asarray(z['rf']))
        lab = 't %.3f (+%.1f h)' % (float(z['t']), (float(z['t']) - TS) * DAYS * 24)
        ax[0].semilogy(np.asarray(z['r']), np.maximum(ms / dr, 1e-30), color=c, lw=1, label=lab)
        ax[1].semilogy(np.asarray(z['r']), np.maximum(ms / np.maximum(np.asarray(z['prof_m']), 1e-300), 1e-12), color=c, lw=1, label=lab)
    rph = float(L['Rph_med'])
    for a in ax:
        a.axvline(rph, color=C2, ls=':', lw=.8); a.axvline(RIN, color='#a3adb8', ls=':', lw=.6); a.set_xlabel(r'r [R$_\odot$]'); a.grid(alpha=.3)
    ax[0].set_ylabel(r'dM$_{s0}$/dr [code]'); ax[0].set_title('stream-gas mass per unit radius (dotted: median R$_{phot}$, r$_{in}$)')
    ax[1].set_ylabel('shell mass fraction of stream gas'); ax[1].set_title('stream-gas fraction of each shell'); ax[1].legend(fontsize=7, ncol=2)
    lo = max(1e-30, np.nanmax(np.asarray(L['prof_ms0']) / np.diff(rf)) * 1e-8); ax[0].set_ylim(bottom=lo)
    fig.tight_layout(); fig.savefig(IMG + '/sp_trc_prof.png', dpi=80); plt.close(fig)
    ms = np.asarray(L['prof_ms0']); cum = np.cumsum(ms) / max(ms.sum(), 1e-300)
    r10 = float(r[np.searchsorted(cum, 0.10)]); r50 = float(r[np.searchsorted(cum, 0.50)])
    below = float(ms[r < rph].sum() / max(ms.sum(), 1e-300))
    f01 = np.asarray(L['prof_ms0']) / np.maximum(np.asarray(L['prof_m']), 1e-300)
    deep = r[f01 > 1e-3]; rdeep = float(deep.min()) if deep.size else None
    w = np.outer(np.diff(L['phf']), -np.diff(np.cos(L['thf'])))
    cov = float((w * (np.nan_to_num(L['map_s0ph']) > 0.1)).sum() / w.sum())
    res.update(maps_t=tl, maps_idx=str(L['idx']), r10=r10, r50=r50, below_rph=below, rdeep=rdeep, cov_ph=cov, rph=rph, nprof=len(sel))
    res['take_maps'] = ('At t = %.4f the photosphere is stream gas (s0 > 0.1) over %.0f %% of its area; the maps show where the accreted gas '
                        'spreads in the orbital and x-z planes.' % (tl, 100 * cov))
    res['take_prof'] = ('Half of the stream gas in the domain is inside r = %.2f R_sun (10 %% inside %.2f; median R_phot %.2f), %.0f %% lies below '
                        'the photosphere, and the stream fraction exceeds 1e-3 of the shell mass down to r = %s R_sun.'
                        % (r50, r10, rph, 100 * below, ('%.2f' % rdeep) if rdeep is not None else '–'))
else:
    res['why_maps'] = 'no stream dump with the VER 2 tracer reduction yet'
save()
print('mk_trc: budget to t %.4f, %d stream dumps with tracer reductions' % (res['t'], len(Rs)))
