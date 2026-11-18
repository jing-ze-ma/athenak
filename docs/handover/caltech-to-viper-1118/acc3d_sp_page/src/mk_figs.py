#!/usr/bin/env python3
"""accsp_prod_page (11-02): figures, movie frames and img/status.json for the sp full-sphere accretor PRODUCTION
(accsp_1102/prod/fresh/stream, static phase = ringsmooth_1102/long/L_gz). Adapted from acc3d_cs128_page/src/mk_figs.py. READ ONLY.
env: OUT, STATIC (static-phase run dir), STREAM (production stream dir), CS128 (comparison run), FORCE=1 redraws frames"""
import os
import re
import glob
import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

E = os.environ.get
J = '/resnick/groups/carnegie_poc/jingze'
OUT = E('OUT', J + '/accsp_prod_page')
ST = E('STATIC', J + '/ringsmooth_1102/long/L_gz')
SR = E('STREAM', J + '/accsp_1102/prod/fresh/stream')
CS = E('CS128', J + '/csacc_1031/test128wb')
FORCE = E('FORCE', '') == '1'
IMG, DATA = OUT + '/img', OUT + '/data'
PORB, DAYS, TS = 0.473169, 8.0521, 0.01
DL1, RIN = 16.2266, 6.3
DARK = {'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
        'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
        'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'axes.titlesize': 11,
        'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
        'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'}
plt.rcParams.update(DARK)
C1, C2, C3, C4, C5 = '#ffb347', '#4cc3ff', '#7dff8a', '#c77dff', '#ff5c7a'
os.makedirs(IMG, exist_ok=True)
for s_ in ('fr_eq', 'fr_ph'):
    os.makedirs(f'{OUT}/movie/{s_}', exist_ok=True)


def readhst(f):
    if not os.path.exists(f):
        return None
    lab = None
    for ln in open(f):
        if ln.startswith('#  ['):
            lab = [x.split('=')[1] for x in ln[1:].split()]
    try:
        a = np.genfromtxt(f, comments='#', invalid_raise=False)
    except Exception:
        return None
    a = np.atleast_2d(a)
    if a.size == 0 or lab is None:
        return None
    # 11-03: after a restart from an earlier rst the hst has rewound rows; keep the newest pass only
    t = a[:, 0]; keep = np.ones(len(t), bool); tmin = np.inf
    for j in range(len(t) - 1, -1, -1):
        if t[j] >= tmin: keep[j] = False
        else: tmin = t[j]
    a = a[keep]
    return {k: a[:, i] for i, k in enumerate(lab) if i < a.shape[1]}


def cat(*ds):
    ds = [d for d in ds if d is not None]
    if not ds:
        return None
    keys = []
    for d in ds:
        keys += [k for k in d if k not in keys]
    out = {}
    for d in ds:
        n = len(d['time'])
        if out:
            keep = out['time'] < d['time'][0] - 1e-12
            out = {k: out[k][keep] for k in out}
        m = len(out['time']) if out else 0
        for k in keys:
            prev = out.get(k, np.full(m, np.nan))
            out[k] = np.concatenate([prev, d[k] if k in d else np.full(n, np.nan)])
    return out


def upto(d, tmax):
    if d is None:
        return None
    s = d['time'] <= tmax + 1e-9
    return {k: v[s] for k, v in d.items()}


def load(stem, static_cut=True):
    a = readhst(f'{ST}/ryper.{stem}.hst'); b = readhst(f'{SR}/ryper.{stem}.hst')
    return cat(upto(a, TS) if static_cut else a, b)


U, W, TR, HY = load('user'), load('wbt'), load('trc'), load('hydro')
CU = cat(readhst(CS + '/static/ryper.user.hst'), readhst(CS + '/stream/ryper.user.hst'))
CH = cat(readhst(CS + '/static/ryper.hydro.hst'), readhst(CS + '/stream/ryper.hydro.hst'))


def logscan(f, tmax=None):
    r = dict(el=[], BE=0, posfix=0, notpos=0, nan=0, fatal=0, ev=[], start='', warn=[])
    if not os.path.exists(f):
        return r
    for ln in open(f, errors='replace'):
        m = re.match(r'elapsed=([0-9.e+-]+) cycle=(\d+) time=([0-9.e+-]+) dt=([0-9.e+-]+)', ln)
        if m:
            if tmax is None or float(m.group(3)) <= tmax + 1e-9:
                r['el'].append([float(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))])
            continue
        m = re.search(r'at cycle (\d+)', ln)
        if 'redone with backward Euler' in ln:
            r['BE'] += 1; m and r['ev'].append(('BE', int(m.group(1))))
        if 'time2_posfix: stage' in ln:
            r['posfix'] += 1; m and r['ev'].append(('posfix', int(m.group(1))))
        if 'NOT POSITIVE' in ln:
            r['notpos'] += 1
        if re.search(r'\bnan\b', ln, re.I) and 'nanc' not in ln.lower():
            r['nan'] += 1
        if re.search(r'FATAL|DIVERGED|Segmentation|Aborted', ln):
            r['fatal'] += 1
        if ln.startswith('###   '):
            r['warn'].append(ln[6:].strip())
    return r


LS = logscan(ST + '/run.log', TS)
LR = logscan(SR + '/run.log')


def mspc(el, nlast=1):
    if len(el) < 2:
        return None, None
    a, b = el[-1 - nlast], el[-1]
    last = (b[0] - a[0]) / max(b[1] - a[1], 1) * 1e3
    allr = (el[-1][0] - el[0][0]) / max(el[-1][1] - el[0][1], 1) * 1e3
    return last, allr


class _LazyRed(dict):
    """pagefast-1102: a red_*.npz as a dict whose members are read on first access (was: every member of every dump read
    on every refresh; only the scalars of old dumps are needed, full arrays only for the newest dump and new frames)"""
    def __init__(self, f):
        super().__init__(); self._f = f; self._z = np.load(f); self._k = list(self._z.files)

    def __missing__(self, k):
        if k not in self._k:
            raise KeyError(k)
        v = self._z[k]; dict.__setitem__(self, k, v); return v

    def __contains__(self, k):
        return k in self._k

    def get(self, k, default=None):
        return self[k] if k in self._k else default

    def keys(self):
        return list(self._k)

    def __iter__(self):
        return iter(self._k)

    def __len__(self):
        return len(self._k)


reds = []
for f in sorted(glob.glob(DATA + '/red_*.npz')):
    reds.append(_LazyRed(f))
reds.sort(key=lambda d: float(d['t']))
red0 = next((d for d in reds if str(d['idx']) == 'a_00000'), None)
KC = json.load(open(DATA + '/cs128_keband.json')) if os.path.exists(DATA + '/cs128_keband.json') else {}
KC = sorted([v for v in KC.values() if v.get('t') is not None and v.get('tot')], key=lambda v: v['t'])

st = {}


def vline(a, xl='t [code]'):
    a.axvline(TS, color='#5b6470', ls=':', lw=1); a.grid(alpha=.3); a.set_xlabel(xl)


# ------------------------------------------------------------------ history
if U is not None:
    t = U['time']
    fig, ax = plt.subplots(3, 3, figsize=(16, 11.5)); ax = ax.ravel()
    a = ax[0]; a.semilogy(t, U['dt'], color=C1, lw=1, label='sp production')
    if CU is not None:
        a.semilogy(CU['time'], CU['dt'], color=C2, lw=1, ls='--', label='cs N=128 test')
    a.set_title('dt [code]'); a.legend(); vline(a)
    el = np.array(LS['el'] + LR['el']) if (LS['el'] or LR['el']) else np.zeros((0, 4))
    a = ax[1]
    if len(el) > 1:
        d0, d1 = np.diff(el[:, 0]), np.diff(el[:, 1]); ok = (d0 > 0) & (d1 > 0)
        a.plot(el[1:, 2][ok], (d0 / np.maximum(d1, 1) * 1e3)[ok], 'o-', color=C2, ms=3, lw=1)
        a.set_ylim(0, None)
    a.set_title('wall ms / cycle (run.log, every 50 cycles)'); vline(a)
    a = ax[2]; a.plot(t, U['Lmeas'], color=C1, lw=1, label=r'L$_{\rm meas}$'); a.plot(t, U['Ltop'], color=C2, lw=1, label=r'L$_{\rm top}$')
    q_ = np.r_[U['Lmeas'][t > 0.001], U['Ltop'][t > 0.001]]
    if q_.size:
        a.set_ylim(np.nanmin(q_) * 0.98, np.nanmax(q_) * 1.02)
    a.set_title(r'luminosity [L$_\odot$]'); a.legend(); vline(a)
    a = ax[3]; a.plot(t, U['Min'], color=C1, lw=1, label=r'M$_{\rm in}$ (nozzle, cumulative)')
    if 'Mnzf' in U:
        a.plot(t, U['Mnzf'], color=C3, lw=1, label=r'M$_{\rm nzf}$')
    if CU is not None:
        a.plot(CU['time'], CU['Min'], color=C2, lw=1, ls='--', label=r'cs128 M$_{\rm in}$')
    a.set_title('stream mass [code]'); a.legend(); vline(a)
    a = ax[4]; a.plot(t, U['KE'], color=C1, lw=1, label='sp production')
    if CU is not None:
        a.plot(CU['time'], CU['KE'], color=C2, lw=1, ls='--', label='cs N=128 test')
    a.set_yscale('log'); a.set_title('total kinetic energy KE [code] (user.hst)'); a.legend(); vline(a)
    a = ax[5]
    if HY is not None:
        a.plot(HY['time'], HY['2-KE'] + HY['3-KE'], color=C1, lw=1, label=r'sp: $\rho(v_\theta^2+v_\phi^2)/2$')
        a.plot(HY['time'], HY['1-KE'], color=C1, lw=.8, ls=':', label='sp radial')
    if CH is not None:
        a.plot(CH['time'], CH['2-KE'] + CH['3-KE'], color=C2, lw=1, ls='--', label='cs128 lateral (panel comps.)')
        a.plot(CH['time'], CH['1-KE'], color=C2, lw=.8, ls=':', label='cs128 radial')
    a.set_yscale('log'); a.set_title('lateral vs radial KE [code] (hydro.hst)'); a.legend(fontsize=7); vline(a)
    a = ax[6]; a.plot(t, U['Mfl'], color=C5, lw=1, label=r'M$_{\rm fl}$ floor mass booked')
    a.set_title('floors and M1 fallbacks'); vline(a); a.set_ylabel(r'M$_{\rm fl}$ [code]')
    if len(el) > 1:
        a2 = a.twinx(); a2.tick_params(colors='#a3adb8')
        for (kind, c, mk) in (('posfix', C3, 'o'), ('BE', C4, 's')):
            cyc = np.array([c_ for k_, c_ in LR['ev'] if k_ == kind], float)
            if cyc.size:
                tt = np.interp(cyc, el[:, 1][el[:, 2] >= TS - 1e-9], el[:, 2][el[:, 2] >= TS - 1e-9]) if (el[:, 2] >= TS - 1e-9).sum() > 1 else np.full(cyc.size, TS)
                a2.step(np.r_[TS, np.sort(tt)], np.r_[0, np.arange(1, cyc.size + 1)], where='post', color=c, lw=1,
                        label={'posfix': 'time2_posfix re-solves', 'BE': 'backward-Euler redos'}[kind])
        a2.set_ylabel('stream events (cumulative)')
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=7)
    a = ax[7]; a.plot(t, U['Pic'], color=C2, lw=1, label='Picard mean'); a.plot(t, U['Picmax'], color=C4, lw=1, label='Picard max')
    a.set_title('M1 Picard iterations'); a.legend(); vline(a)
    a = ax[8]
    if W is not None:
        a.plot(W['time'], W['Twb'], color=C3, lw=1, label=r'WB torque T$_{\rm wb}$')
    if TR is not None:
        a2 = a.twinx(); a2.plot(TR['time'], TR['Ms0'], color=C1, lw=1, label=r'stream tracer M$_{s0}$'); a2.tick_params(colors='#a3adb8')
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels(); a.legend(h1 + h2, l1 + l2, fontsize=7)
    a.set_title('WB tangential tide (wbt.hst) and stream tracer (trc.hst)'); vline(a)
    fig.suptitle('sp full-sphere production: static phase (t < 0.01, ringsmooth_1102/long/L_gz) + stream (dotted line = stream on); dashed = cs N=128 test',
                 color='#e8ecf1')
    fig.tight_layout(); fig.savefig(IMG + '/sp_hist.png', dpi=80); plt.close(fig)
    st['t'] = float(t[-1]); st['dt'] = float(U['dt'][-1])
    for k in ('Lmeas', 'Ltop', 'Min', 'Mnzf', 'Mfl', 'KE', 'Erad', 'Pic', 'Picmax', 'Mdom'):
        if k in U:
            st[k] = float(U[k][-1])
    r_ = t > TS - 1e-9
    st['dt_min_stream'] = float(U['dt'][r_].min()) if r_.any() else None
    st['hst_nan'] = int(sum((~np.isfinite(v)).sum() for f_ in (ST + '/ryper.user.hst', SR + '/ryper.user.hst') for v in (readhst(f_) or {}).values()))  # raw files (static/stream column sets differ)
    if HY is not None:
        st['KElat'] = float(HY['2-KE'][-1] + HY['3-KE'][-1]); st['KErad'] = float(HY['1-KE'][-1])
        if CH is not None and t[-1] <= CH['time'][-1] + 1e-9:
            j = np.argmin(abs(CH['time'] - t[-1])); st['cs_KElat'] = float(CH['2-KE'][j] + CH['3-KE'][j]); st['cs_t'] = float(CH['time'][j])
    if CU is not None:
        st['cs_dt_min_stream'] = float(CU['dt'][CU['time'] > TS].min())
    if W is not None:
        st['Twb'] = float(W['Twb'][-1])

st['static'] = dict(BE=LS['BE'], posfix=LS['posfix'], nan=LS['nan'], fatal=LS['fatal'], ms=mspc(LS['el']),
                    cycle=(LS['el'] or [[0, 0, 0, 0]])[-1][1])
st['stream'] = dict(BE=LR['BE'], posfix=LR['posfix'], notpos=LR['notpos'], nan=LR['nan'], fatal=LR['fatal'], ms=mspc(LR['el']),
                    cycle=LR['el'][-1][1] if LR['el'] else 0, el=LR['el'][-1][0] if LR['el'] else 0, warn=LR['warn'])
mt = os.path.getmtime(SR + '/run.log') if os.path.exists(SR + '/run.log') else 0
st['log_age_min'] = (time.time() - mt) / 60
pl = open(os.path.dirname(SR) + '/prod.log').read() if os.path.exists(os.path.dirname(SR) + '/prod.log') else ''
m = re.search(r'stream rc (\d+)', pl)
st['rc'] = int(m.group(1)) if m else None
st['state'] = (f"stopped (rc {st['rc']})" if m else ('running' if st['log_age_min'] < 30 else 'stalled? (run.log not written for %.0f min)' % st['log_age_min']))
st['ndump'] = len([x for x in glob.glob(SR + '/bin/ryper.hydro_w.*.bin') if not os.path.islink(x)])
st['nrst'] = len(glob.glob(SR + '/rst/ryper.*.rst'))
st['reds'] = [dict(idx=str(d['idx']), t=float(d['t']), m1=bool(d['has_m1']), Rph_med=float(d['Rph_med']), Rph_p2=float(d['Rph_p2']),
                   Rph_p98=float(d['Rph_p98']), Tph_med=float(d['Tph_med']), n_floor=int(d['n_floor']), n_cells=int(d['n_cells']),
                   n_floor_hot=int(d['n_floor_hot']), nonfinite=int(d['nonfinite']), vmax=float(d['vmax'])) for d in reds]

# ------------------------------------------------------------------ slices
SPEC = [('lrho', r'log$_{10}\rho$ [code]', 'magma', -8, 2.5),
        ('lT', r'log$_{10}T_{\rm gas}$ [K]', 'inferno', 3.6, 6.0),
        ('lbeta', r'log$_{10}P_{\rm rad}/P_{\rm gas}$', 'RdBu_r', -2, 5),
        ('vr', r'v$_r$ [km/s]', 'PuOr_r', -300, 300),
        ('vphi', r'v$_\phi$ [km/s] (rotating frame)', 'coolwarm', -300, 300)]


def nrm(key, lo, hi):
    if key == 'lbeta':
        from matplotlib.colors import TwoSlopeNorm
        return dict(norm=TwoSlopeNorm(0, lo, hi))
    return dict(vmin=lo, vmax=hi)


def geq(d):
    PH, RR = np.meshgrid(d['phf'], d['rf'], indexing='ij')
    return RR * np.cos(PH), RR * np.sin(PH)


def gme(d, h):
    TH, RR = np.meshgrid(d['thf'], d['rf'], indexing='ij')
    return (1 if h == 0 else -1) * RR * np.sin(TH), RR * np.cos(TH)


def title(d):
    t = float(d['t'])
    ph = 'static phase (stream off)' if t < TS + 1e-6 else f'stream on for {t - TS:.4f} code = {(t - TS) * DAYS * 24:.1f} h'
    return f't = {t:.4f} code = {t * DAYS * 24:.1f} h = {t / PORB:.4f} P_orb; {ph}'


def draw(a, d, pre, key, cmap, kw, mask=False):
    q = d.get(pre + '_' + key)
    if q is None:
        return None
    if mask and key in ('lT', 'lbeta'):
        q = np.where(d[pre + '_lrho'] < -7.0, np.nan, q)
    if pre == 'eq':
        X, Y = geq(d); return a.pcolormesh(X, Y, q, cmap=cmap, shading='flat', rasterized=True, **kw)
    pc = None
    for h in (0, 1):
        X, Y = gme(d, h); pc = a.pcolormesh(X, Y, q[h], cmap=cmap, shading='flat', rasterized=True, **kw)
    return pc


def slices(d, fn, lim=16.5):
    fig, ax = plt.subplots(2, 5, figsize=(24, 9.6))
    for i, (key, lab, cm, lo, hi) in enumerate(SPEC):
        cmap = plt.get_cmap(cm).copy(); cmap.set_bad('#5b6470')
        for row, (pre, xl, yl) in enumerate((('eq', 'x', 'y'), ('me', 'x', 'z'))):
            a = ax[row, i]
            a.add_patch(plt.Circle((0, 0), RIN, color='#3a424c'))
            pc = draw(a, d, pre, key, cmap, nrm(key, lo, hi))
            a.set_xlim(-lim, lim); a.set_ylim(-lim, lim); a.set_aspect('equal')
            if pc is None:
                a.text(0, 0, 'no m1 dump\nat this time', ha='center', va='center', color='#a3adb8'); a.set_title(lab, fontsize=9)
                continue
            a.add_patch(plt.Circle((0, 0), DL1, fill=False, ec=C2, ls='--', lw=.7))
            a.set_title(('orbital plane: ' if row == 0 else 'meridional x-z (y = 0): ') + lab, fontsize=9)
            a.set_xlabel(xl + r' [R$_\odot$]'); a.set_ylabel(yl + r' [R$_\odot$]')
            if row == 0:
                a.annotate(r'L1 radius (dashed); donor $\rightarrow$', xy=(lim * .97, -lim * .93), ha='right', fontsize=8, color=C2)
            fig.colorbar(pc, ax=a, shrink=.8)
    fig.suptitle(f"RY Per accretor, sp full sphere 500x116x512 (production): {title(d)}; dump {d['idx']}", color='#e8ecf1', fontsize=13)
    fig.tight_layout(); fig.savefig(fn, dpi=72); plt.close(fig)


def latlon(a, d, q, cm, lo, hi, lab):
    NP = q.shape[0]
    qq = np.roll(q, -NP // 2, axis=0)
    lonf = np.degrees(d['phf']) - 180.0; latf = 90.0 - np.degrees(d['thf'])
    pc = a.pcolormesh(lonf, latf, qq.T, cmap=cm, vmin=lo, vmax=hi, shading='flat', rasterized=True)
    a.set_xlabel(r'longitude $\phi$ [deg] (0 = towards the donor / L1)'); a.set_ylabel('latitude [deg]'); a.set_title(lab, fontsize=9)
    a.set_xticks(range(-180, 181, 60)); a.set_yticks(range(-90, 91, 30)); a.set_xlim(-180, 180); a.set_ylim(-90, 90)
    return pc


def photfig(d, fn):
    fig = plt.figure(figsize=(22, 10.5))
    gs = fig.add_gridspec(2, 4)
    a = fig.add_subplot(gs[0, 0:2])
    plt.colorbar(latlon(a, d, d['map_Rph'], 'viridis', 8.8, 9.6, r'Rosseland photosphere R($\tau_R$ = 2/3) [R$_\odot$] (radial rays from x1max; '
                        r'stream columns saturate)'), ax=a, shrink=.85)
    a = fig.add_subplot(gs[0, 2:4])
    plt.colorbar(latlon(a, d, np.log10(np.maximum(d['map_Tph'], 1)), 'inferno', 4.0, 4.9, r'log$_{10}T_{\rm gas}$ [K] at $\tau_R$ = 2/3'), ax=a, shrink=.85)
    zx, zy = (8.2, 11.0), (-1.4, 1.4)
    for i, (pre, key, lab, cm, lo, hi) in enumerate((('eq', 'lrho', r'log$_{10}\rho$', 'magma', -8, 2.5), ('eq', 'vr', r'v$_r$ [km/s]', 'PuOr_r', -20, 20),
                                                    ('me', 'lrho', r'log$_{10}\rho$', 'magma', -8, 2.5), ('me', 'vr', r'v$_r$ [km/s]', 'PuOr_r', -20, 20))):
        a = fig.add_subplot(gs[1, i])
        pc = draw(a, d, pre, key, cm, dict(vmin=lo, vmax=hi), mask=False)
        a.add_patch(plt.Circle((0, 0), float(d['Rph_med']), fill=False, ec=C3, ls='--', lw=.8))
        a.set_xlim(*zx); a.set_ylim(*zy); a.set_aspect('equal')
        a.set_title(('L1 tip, orbital plane: ' if pre == 'eq' else 'L1 tip, x-z plane: ') + lab, fontsize=9)
        a.set_xlabel(r'x [R$_\odot$]'); a.set_ylabel(('y' if pre == 'eq' else 'z') + r' [R$_\odot$]')
        fig.colorbar(pc, ax=a, shrink=.8)
    fig.suptitle(f"photosphere and L1 tip (dashed green = median R_phot {float(d['Rph_med']):.3f}): {title(d)}; dump {d['idx']}",
                 color='#e8ecf1', fontsize=13)
    fig.tight_layout(); fig.savefig(fn, dpi=72); plt.close(fig)


paired = [d for d in reds if bool(d['has_m1'])]
if reds:
    last = paired[-1] if paired else reds[-1]
    slices(last, IMG + '/sp_slices_latest.jpg'); photfig(reds[-1], IMG + '/sp_phot_latest.jpg')
    st['slice_dump'] = dict(idx=str(last['idx']), t=float(last['t'])); st['phot_dump'] = dict(idx=str(reds[-1]['idx']), t=float(reds[-1]['t']))
    _todo = []   # pagefast-1102: new movie frames rendered in parallel (FRAME_PROCS, fork; same functions, one file each)
    for d in reds:
        for sub, fn_ in (('fr_eq', slices), ('fr_ph', photfig)):
            f = f"{OUT}/movie/{sub}/f_{float(d['t']):09.5f}_{d['idx']}.png"
            if FORCE or not os.path.exists(f) or os.path.getmtime(f) < os.path.getmtime(f"{DATA}/red_{d['idx']}.npz"):
                _todo.append((fn_, d, f))

    def _rf(i):
        fn_, d, f = _todo[i]
        if _nfp > 1:
            d._z = np.load(d._f)   # own file handle in the forked worker (a shared zip offset corrupts reads)
        fn_(d, f); return i
    _nfp = min(int(E('FRAME_PROCS', '4')), len(_todo))
    if _nfp > 1:
        import multiprocessing as _mp
        with _mp.get_context('fork').Pool(_nfp) as _pl:
            list(_pl.imap_unordered(_rf, range(len(_todo))))
    else:
        for i in range(len(_todo)):
            _rf(i)
    # tip and photosphere vs time
    fig, ax = plt.subplots(1, 3, figsize=(18, 4.6))
    cmap = plt.get_cmap('viridis')
    sel = reds if len(reds) <= 12 else [reds[i] for i in np.unique(np.linspace(0, len(reds) - 1, 12).astype(int))]
    for n, d in enumerate(sel):
        c = cmap(n / max(len(sel) - 1, 1)); r = d['r']; s = (r > 8.6) & (r < 11)
        if red0 is not None:
            ax[0].plot(r[s], d['tip_rho'][s] / red0['tip_rho'][s], color=c, lw=1, label=f"t {float(d['t']):.4f}")
        ax[1].plot(r[s], d['tip_vr'][s], color=c, lw=1, label=f"t {float(d['t']):.4f}")
    ax[0].set_title(r'L1-tip column ($\phi$ = 0, equator): $\rho/\rho_{t=0}$'); ax[0].set_yscale('log')
    ax[1].set_title(r'L1-tip column: v$_r$ [km/s]')
    for a in ax[:2]:
        a.axvline(float(reds[-1]['Rph_med']), color=C3, ls=':', lw=.8); a.set_xlabel(r'r [R$_\odot$]'); a.grid(alpha=.3); a.legend(fontsize=6, ncol=2)
    tt = [float(d['t']) for d in reds]
    ax[2].plot(tt, [float(d['Rph_med']) for d in reds], 'o-', color=C1, ms=3, label='median (solid-angle weighted)')
    ax[2].plot(tt, [float(d['Rph_p2']) for d in reds], 's--', color=C2, ms=3, label='2nd percentile')
    ax[2].set_title(r'R$_{\rm phot}$ ($\tau_R$ = 2/3) over all columns [R$_\odot$]'); ax[2].legend(); vline(ax[2])
    fig.tight_layout(); fig.savefig(IMG + '/sp_tip.png', dpi=80); plt.close(fig)
    if red0 is not None:
        r = reds[-1]['r']; s2 = (r > 9.45) & (r < 9.7); s1 = (r > 9.25) & (r < 9.47)
        st['tip'] = dict(t=float(reds[-1]['t']), vr_min=float(reds[-1]['tip_vr'][s1].min()), vr_max=float(reds[-1]['tip_vr'][s1].max()),
                         rho_ratio_max=float((reds[-1]['tip_rho'][s2] / red0['tip_rho'][s2]).max()))

# ------------------------------------------------------------------ KE by latitude band vs cs128
BANDS = ['0-20', '20-40', '40-60', '60-80', '80-90']
RBIN = ['r<8', '8-9.3', '9.3-11', '>11']
kr = [d for d in reds if 'keb_lat' in d and np.size(d['keb_lat']) == 20]
if kr and KC:
    tcs_end = KC[-1]['t']
    tl = float(kr[-1]['t'])
    tm = min(tl, tcs_end)
    dsp = min(kr, key=lambda d: abs(float(d['t']) - tm)); dcs = min(KC, key=lambda v: abs(v['t'] - float(dsp['t'])))
    fig, ax = plt.subplots(1, 4, figsize=(21, 4.8))
    x = np.arange(5); wb = 0.27
    sp_lat = np.asarray(dsp['keb_lat']); cs_lat = np.asarray(dcs['lat'])
    a = ax[0]
    a.bar(x - wb / 2, sp_lat.sum(1), wb, color=C1, label=f"sp production t {float(dsp['t']):.4f}")
    a.bar(x + wb / 2, cs_lat.sum(1), wb, color=C2, label=f"cs N=128 test t {dcs['t']:.4f}")
    if tl > tm + 1e-6:
        a.plot(x, np.asarray(kr[-1]['keb_lat']).sum(1), 'D', color=C3, label=f"sp latest t {tl:.4f}")
    a.set_xticks(x); a.set_xticklabels([b + '°' for b in BANDS]); a.set_yscale('log'); a.set_xlabel('|latitude| band')
    a.set_title('lateral KE per |latitude| band, all r [code]'); a.legend(fontsize=7); a.grid(alpha=.3, axis='y')
    a = ax[1]
    a.bar(x - wb / 2, sp_lat[:, 0], wb, color=C1, label='sp, r < 8'); a.bar(x + wb / 2, cs_lat[:, 0], wb, color=C2, label='cs128, r < 8')
    a.plot(x - wb / 2, sp_lat[:, 1], '_', ms=14, mew=2, color=C5, label='sp, r 8-9.3'); a.plot(x + wb / 2, cs_lat[:, 1], '_', ms=14, mew=2, color=C4, label='cs128, r 8-9.3')
    a.set_xticks(x); a.set_xticklabels([b + '°' for b in BANDS]); a.set_yscale('log'); a.set_xlabel('|latitude| band')
    a.set_title(f"same, deep (r < 8) and envelope (r 8-9.3), t {float(dsp['t']):.4f}"); a.legend(fontsize=7); a.grid(alpha=.3, axis='y')
    a = ax[2]
    ts = [float(d['t']) for d in kr]
    a.plot(ts, [float(d['keb_tot'][2]) for d in kr], 'o-', color=C1, ms=3, label='sp lateral, all latitudes')
    a.plot([v['t'] for v in KC], [v['tot'][2] for v in KC], 's--', color=C2, ms=3, label='cs128 lateral, all latitudes')
    a.plot(ts, [float(np.asarray(d['keb_lat'])[3:].sum()) for d in kr], 'o-', color=C5, ms=3, label='sp lateral, |lat| > 60°')
    a.plot([v['t'] for v in KC], [float(np.asarray(v['lat'])[3:].sum()) for v in KC], 's--', color=C4, ms=3, label='cs128 lateral, |lat| > 60°')
    a.set_yscale('log'); a.set_title('lateral KE vs time [code] (keband.py)'); a.legend(fontsize=7); vline(a)
    a = ax[3]
    for d, c, lab in ((red0, '#5b6470', 't 0'), (next((d for d in reds if abs(float(d['t']) - TS) < 1e-6), None), C2, 't 0.01 (stream on)'), (reds[-1], C1, f"latest t {float(reds[-1]['t']):.4f}")):
        if d is not None and 'pole_vphi_rms' in d:
            a.plot(d['pole_th'], d['pole_vphi_rms'], 'o-', color=c, ms=3, label=lab)
    a.set_xlabel('colatitude of the row [deg] (both poles averaged)'); a.set_ylabel('km/s')
    a.set_title(r'$\phi$-rms of v$_\phi$ in the 10 polar rows, r 6.8'); a.legend(fontsize=7); a.grid(alpha=.3)
    a.axvspan(0, float(kr[-1]['pole_th'][6] + 2.3) if 'pole_th' in kr[-1] else 30, color='#262c33', zorder=0)
    fig.suptitle('lateral-KE check: KE = 0.5 ρ (v_θ² + v_φ²) dV by |latitude| band (accsp_1102/tools/keband.py, sp orthonormal / cs exact metric); '
                 'shaded = the 7 ring-averaged polar rows', color='#e8ecf1')
    fig.tight_layout(); fig.savefig(IMG + '/sp_keband.png', dpi=80); plt.close(fig)
    st['keband'] = dict(t_sp=float(dsp['t']), t_cs=dcs['t'], sp_bands=sp_lat.sum(1).tolist(), cs_bands=cs_lat.sum(1).tolist(),
                        sp_tot=float(dsp['keb_tot'][2]), cs_tot=dcs['tot'][2], t_latest=tl, latest_tot=float(kr[-1]['keb_tot'][2]),
                        sp_pole=float(sp_lat[3:].sum()), cs_pole=float(cs_lat[3:].sum()))
json.dump(st, open(IMG + '/status.json', 'w'), indent=1)
print('mk_figs: t %.5f, %d reductions, state %s' % (st.get('t', -1), len(reds), st['state']))
