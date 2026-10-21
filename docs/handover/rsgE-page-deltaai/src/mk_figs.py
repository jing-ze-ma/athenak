#!/usr/bin/env python3
"""RSGPAGE 10-19: figures + status table for the RSG g25 cs production page (READ ONLY on the run).
usage: mk_figs.py <red_dir> <img_dir>"""
import os
import re
import sys
import glob
import json
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import TwoSlopeNorm, LogNorm  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

RUN = os.environ.get('RUN', '/raven/ptmp/jinma/rsg_g25_1019/prod/R')
ICF = os.environ.get('ICF', '/raven/ptmp/jinma/rsg_g25_1019/ic/pass3b/ic_rsg_g25.txt')
OLD = os.environ.get('OLDRUN', '/raven/ptmp/jinma/rsg_g25_1019/prod/R_oldseed')
RED, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
RS, DAY, LSUN = 6.957e10, 86400., 3.828e33
RPH0 = float(os.environ.get('RPH0', '5.80019e13'))   # rsgE IC (ic_rsg_g25_deep.txt header)
RPHS = '%.0f R_sun' % (RPH0 / 6.957e10)
ZR = 0.2   # zoom radius / R_ph for the deep slice row
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
                     'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'axes.titlesize': 11,
                     'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
                     'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'})
TC = ['#9aa4ae', '#4cc3ff', '#7dff8a', '#ffb347', '#ff5c7a']


def readhst(f):
    lab = None
    for ln in open(f):
        if ln.startswith('#  ['):
            lab = [x.split('=')[1] for x in ln[1:].split()]
            break
    rows = [ln.split() for ln in open(f) if not ln.startswith('#')]
    nc = len(lab)
    rows = [x for x in rows if len(x) == nc]          # drop a partially written last line
    d = np.array(rows, float)
    t = d[:, 0]
    keep = np.zeros(len(t), bool)
    tmin = np.inf
    for i in range(len(t) - 1, -1, -1):
        if t[i] < tmin:
            keep[i] = True
            tmin = t[i]
    d = d[keep]
    return {n: d[:, i] for i, n in enumerate(lab)}


def smooth(y, n):
    if len(y) < 2 * n or n < 2:
        return y
    c = np.cumsum(np.insert(y, 0, 0.))
    s = (c[n:] - c[:-n]) / n
    return np.concatenate([np.full(n - 1, np.nan), s])


H = readhst(os.path.join(RUN, 'rsg.user.hst'))
HH = readhst(os.path.join(RUN, 'rsg.hydro.hst'))
td = H['time'] / DAY
HO = None
if OLD and os.path.exists(os.path.join(OLD, 'rsg.user.hst')):
    try:
        HO = readhst(os.path.join(OLD, 'rsg.user.hst'))
    except Exception:
        HO = None
reds = sorted(glob.glob(os.path.join(RED, 'red_*.npz')))
R = [np.load(f) for f in reds]
tr = np.array([float(d['t']) for d in R]) / DAY

# ---------------- history ----------------
fig, ax = plt.subplots(2, 3, figsize=(15, 7.8), constrained_layout=True)
n = max(1, min(100, len(td) // 60))
a = ax[0, 0]
Lin = H['L_in']
for k, c, lab in (('L_top', '#ffb347', 'L_top'), ('L_mid', '#c77dff', 'L_mid'), ('L_bot', '#4cc3ff', 'L_bot'),
                  ('L_MLT', '#7dff8a', 'L_MLT (sub-grid)')):
    a.plot(td, H[k] / Lin, color=c, lw=1.0, label=lab + ' / L_in')
a.axhline(1, color='#a3adb8', lw=0.6, ls=':')
a.set_ylabel('luminosity / L_in'); a.legend(loc='best'); a.set_title('luminosities (L_in = %.3g L_sun imposed at r_in)' % (Lin[-1] / LSUN))
a = ax[0, 1]
a.semilogy(td, np.maximum(H['KE_int'], 1e30), color='#ff5c7a', lw=1.2, label='KE_int')
a.semilogy(td, np.maximum(H['KEr_int'], 1e30), color='#c77dff', lw=1.0, label='KE_r,int (radial)')
if HO is not None:
    a.semilogy(HO['time'] / DAY, np.maximum(HO['KE_int'], 1e30), color='#9aa4ae', lw=0.9, ls='--', label='KE_int, old seed (archived)')
a.set_ylabel('erg'); a.legend(); a.set_title('kinetic energy inside the star (seed growth)')
a = ax[0, 2]
E0 = H['Etot'][0]
a.plot(td, (H['Etot'] - E0) / abs(E0), color='#e8ecf1', lw=1.2, label='(E_tot − E_tot,0)/|E_tot,0| (user hst)')
ct = np.concatenate([[0], np.cumsum(0.5 * (H['L_in'][1:] + H['L_in'][:-1] - H['L_top'][1:] - H['L_top'][:-1]) * np.diff(H['time']))])
a.plot(td, ct / abs(E0), color='#ffb347', lw=1.0, ls='--', label='∫(L_in − L_top) dt / |E_tot,0|')
M0 = H['M_tot'][0]
a.plot(td, (H['M_tot'] - M0) / M0, color='#4cc3ff', lw=1.0, label='(M_tot − M_0)/M_0')
a.set_title('energy and mass budget'); a.legend(loc='best'); a.set_ylabel('relative change')
a = ax[1, 0]
a.plot(td, H['dt'], color='#7dff8a', lw=1.0)
a.set_ylabel('dt [s]'); a.set_title('time step')
a2 = a.twinx()
a2.plot(td, smooth(H['Picard'], n), color='#e8ecf1', lw=0.8, alpha=0.8)
a2.set_ylabel('Picard per cycle (running mean, white)')
a = ax[1, 1]
if len(tr):
    a.fill_between(tr, [float(d['Rph_p2']) / RPH0 for d in R], [float(d['Rph_p98']) / RPH0 for d in R],
                   color='#4cc3ff', alpha=0.25, lw=0, label='2–98 % over the sphere')
    a.plot(tr, [float(d['Rph_med']) / RPH0 for d in R], 'o-', ms=2.5, color='#4cc3ff', lw=1.2, label='median R_ph')
a.axhline(1, color='#a3adb8', lw=0.6, ls=':')
a.set_ylabel('R_ph / R_ph,IC (%s)' % RPHS); a.legend(loc='lower left'); a.set_title('photosphere (from the 1 d dumps)')
a2 = a.twinx()
if len(tr):
    a2.plot(tr, [float(d['Teff_med']) for d in R], color='#ffb347', lw=0.9, ls='--')
a2.set_ylabel('median T_eff [K] (dashed)', color='#ffb347')
a = ax[1, 2]
a.semilogy(td, np.maximum(np.abs(H['Mdot_top']), 1e8), color='#66d9c6', lw=1.0, label='|Mdot_top| [g/s]')
a.semilogy(td, np.maximum(np.abs(H['P_sponge']), 1e20), color='#c77dff', lw=1.0, label='|P_sponge| [erg/s]')
a.semilogy(td, np.maximum(np.abs(H['P_esrc']), 1e20), color='#ffb347', lw=1.0, label='|P_esrc| [erg/s]')
a.set_title('outer boundary / sources'); a.legend()
for a in ax.flat:
    a.set_xlabel('t [d]'); a.grid(True, lw=0.5)
fig.savefig(os.path.join(OUT, 'fig_hist.png'), dpi=90)
plt.close(fig)

# ---------------- radial profiles ----------------
ic = np.loadtxt(ICF, comments='#')
icx, icrho, icT, icfm = ic[:, 0] / RPH0, ic[:, 1], ic[:, 5], ic[:, 6]
if len(R):
    nt = len(R)
    want = [0, 5, 10, tr[-1] / 2, tr[-1]]
    pick = sorted(set(int(np.argmin(abs(tr - w))) for w in want))
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for c, p in enumerate(pick):
        d = R[p]
        x = np.asarray(d['r']) / RPH0
        col = TC[c % len(TC)] if c < len(TC) else None
        lab = 't = %.1f d' % tr[p]
        ax[0, 0].loglog(x, d['rho'], color=col, lw=1, label=lab)
        if 'Tgas' in d.files:
            ax[0, 1].loglog(x, d['Tgas'], color=col, lw=1.1, label=lab)
            ax[0, 1].loglog(x, d['Trad'], color=col, lw=0.6, ls='--')
        else:
            ax[0, 1].loglog(x, d['Trad'], color=col, lw=1, label=lab + ' (T_rad)')
        ax[0, 2].semilogx(x, np.asarray(d['vrrms']) / 1e5, color=col, lw=1.1, label=lab + ' v_r')
        ax[0, 2].semilogx(x, np.asarray(d['vtrms']) / 1e5, color=col, lw=0.8, ls='--')
        Lin = float(H['L_in'][-1])
        ax[1, 0].semilogx(x, np.asarray(d['Lr']) / Lin, color=col, lw=1.1, label=lab)
        ax[1, 1].semilogx(x, np.asarray(d['Lconv']) / Lin, color=col, lw=1.1, label=lab)
        ax[1, 2].semilogx(x, d['dlrho_rms'], color=col, lw=1.1, label=lab)
    ax[0, 0].loglog(icx, icrho, color='w', lw=0.9, ls=':', label='IC column')
    ax[0, 1].loglog(icx, icT, color='w', lw=0.9, ls=':', label='IC T_gas')
    ax[1, 1].semilogx(icx, icfm, color='w', lw=0.9, ls=':', label='IC scaffold F_MLT/F')
    ax[0, 0].set_ylabel('⟨ρ⟩ [g cm⁻³]'); ax[0, 0].set_title('density (shell mean)')
    ax[0, 1].set_ylabel('⟨T⟩ [K]'); ax[0, 1].set_title('gas temperature T_gas (solid), T_rad (dashed); shell means')
    ax[0, 2].set_ylabel('km/s'); ax[0, 2].set_title('ρ-weighted rms v_r (solid), rms |v_t| (dashed)')
    ax[1, 0].set_ylabel('L_rad(r) / L_in'); ax[1, 0].set_title('radiative luminosity 4πr²⟨F_r⟩'); ax[1, 0].set_ylim(-0.2, 1.6)
    ax[1, 1].set_ylabel('L_conv,res(r) / L_in'); ax[1, 1].set_title('resolved enthalpy flux (ideal-gas h)'); ax[1, 1].set_ylim(-0.3, 1.2)
    ax[1, 2].set_yscale('log'); ax[1, 2].set_ylabel('rms δ log ρ'); ax[1, 2].set_title('lateral density fluctuation per shell')
    for a in ax.flat:
        a.axvline(1, color='#a3adb8', lw=0.5, ls=':')
        a.set_xlabel('r / R_ph,IC'); a.grid(True, lw=0.5); a.legend(fontsize=7)
    fig.savefig(os.path.join(OUT, 'fig_profiles.png'), dpi=100)
    plt.close(fig)


# ---------------- slices (native cells) ----------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import csplot as cs  # noqa: E402
ZR = float(os.environ.get('ZOOM_R', '0.25'))
PB = (0.9, 1.2)                     # photosphere band, unwrapped


def pct(a, lo, hi):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return float(np.percentile(a, lo)), float(np.percentile(a, hi))


def ringq(d, key):
    """native ring quantities: lrho, dlrho (minus shell mean of log rho of the ring), vr, Tg, dT (Tg / shell mean - 1)"""
    lr = np.asarray(d[key + '_lrho'], float)
    Tg = np.asarray(d[key + '_Tg'], float)
    ns = lr.shape[1]
    return dict(lr=lr, dlr=lr - np.log10(np.asarray(d['rho'][:ns], float))[None, :], vr=np.asarray(d[key + '_vr'], float),
                Tg=Tg, dT=100 * (Tg / np.asarray(d['Tgas'][:ns], float)[None, :] - 1), fl=np.asarray(d[key + '_Tfl']))


def slicefig(d, fn):
    rf = np.asarray(d['rf'], float) / RPH0
    ns = int(d['ns'])
    rc = 0.5 * (rf[1:ns + 1] + rf[:ns])
    fig, ax = plt.subplots(3, 4, figsize=(18, 13.6), constrained_layout=True)
    for row, (key, zoom) in enumerate((('eqa', False), ('mea', False), ('eqa', True))):
        q = ringq(d, key)
        ang = np.asarray(d['eqa_phi' if key == 'eqa' else 'mea_psi'], float)
        i0, i1 = cs.rrange(rf, 0, ZR) if zoom else cs.rrange(rf, 0, 1.5)
        i1 = min(i1, ns)
        reg = (rc < ZR) if zoom else (rc > ZR) & (rc < 1.0)
        vv = max(1e-3, np.nanpercentile(abs(q['vr'][:, reg]), 99.5))
        dv = max(1e-5, np.nanpercentile(abs(q['dlr'][:, reg]), 99.5))
        Tl = pct(q['Tg'][:, reg & (rc > 0)] if zoom else q['Tg'][:, :i1][~q['fl'][:, :i1]], 0.5, 99.5) if zoom else (2000 * 1.003, max(2e3 * 1.01, np.nanmax(q['Tg'][:, :i1])))
        lrl = pct(q['lr'][:, reg], 0.5, 99.5) if zoom else (-17.2, float(np.nanmax(q['lr'])))
        pans = [(q['lr'], 'viridis', cs.norm_of('lin', lrl), 'log ρ [g cm⁻³]', None),
                (q['dlr'], 'RdBu_r', cs.norm_of('sym', dv), 'log ρ − ⟨log ρ⟩_shell', None),
                (q['Tg'], cs.cmap_floor('inferno'), cs.norm_of('log', Tl), 'T_gas [K] (grey: floor 2000 K)', q['fl']),
                (q['vr'], 'RdBu_r', cs.norm_of('sym', vv), 'v_r [km/s]', None)]
        for c, (qq, cm, nm, lab, fl) in enumerate(pans):
            a = ax[row, c]
            m = cs.ring_xy(a, ang, rf, cs.with_floor(qq, fl, nm), i0, i1, cm, nm)
            fig.colorbar(m, ax=a, label=lab, shrink=0.85, extend='both')
            a.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec='w', lw=0.4, ls=':'))
            if zoom:
                a.add_patch(plt.Circle((0, 0), 0.0494, fill=False, ec='w', lw=0.5, ls='--'))
                a.set_xlim(-ZR, ZR); a.set_ylim(-ZR, ZR)
            else:
                a.set_xlim(-1.5, 1.5); a.set_ylim(-1.5, 1.5)
            a.set_aspect('equal')
            a.set_xlabel('x / R_ph,IC'); a.set_ylabel('y / R_ph,IC' if key == 'eqa' else 'z / R_ph,IC')
        ax[row, 0].set_title(('equatorial plane z = 0' if key == 'eqa' else 'meridional plane y = 0 (x–z)') +
                             (' — zoom r < %.2f R_ph (dashed: MESA CZ base)' % ZR if zoom else ''), loc='left')
    fig.suptitle('t = %.1f d — native cells (512 angles × %d radial cells; dotted circle: R_ph of the IC, %s)'
                 % (float(d['t']) / DAY, ns, RPHS), fontsize=13)
    fig.savefig(os.path.join(OUT, fn), dpi=120)
    plt.close(fig)


def photfig(d, fn):
    rf = np.asarray(d['rf'], float) / RPH0
    ns = int(d['ns'])
    i0, i1 = cs.rrange(rf, *PB); i1 = min(i1, ns)
    fig, ax = plt.subplots(4, 1, figsize=(18, 15), constrained_layout=True, sharex=True)
    q = ringq(d, 'eqa'); ang = np.asarray(d['eqa_phi'], float)
    sl = slice(i0, i1)
    vv = max(1e-3, np.nanpercentile(abs(q['vr'][:, sl]), 99.5))
    dT = max(1e-3, np.nanpercentile(abs(q['dT'][:, sl][~q['fl'][:, sl]]), 99.5))
    Tl = (2000 * 1.003, max(2e3 * 1.01, float(np.nanmax(q['Tg'][:, sl]))))
    pans = [(q['lr'], 'viridis', cs.norm_of('lin', pct(q['lr'][:, sl], 0.5, 99.5)), 'log ρ [g cm⁻³]', None),
            (q['vr'], 'RdBu_r', cs.norm_of('sym', vv), 'v_r [km/s]', None),
            (q['Tg'], cs.cmap_floor('inferno'), cs.norm_of('log', Tl), 'T_gas [K] (grey: floor)', q['fl']),
            (q['dT'], cs.cmap_floor('PuOr_r'), cs.norm_of('sym', dT), 'δT_gas / ⟨T_gas⟩ [%]', q['fl'])]
    for a, (qq, cm, nm, lab, fl) in zip(ax, pans):
        m = cs.ring_unwrap(a, ang, rf, cs.with_floor(qq, fl, nm), i0, i1, cm, nm)
        fig.colorbar(m, ax=a, label=lab, pad=0.01, extend='both')
        a.set_ylabel('r / R_ph,IC'); a.set_ylim(PB)
        for e in (-135, -45, 45, 135):
            a.axvline(e, color='w', lw=0.4, ls=':')
    ax[-1].set_xlabel('φ [deg] (equatorial plane; dotted: panel edges)'); ax[-1].set_xlim(-135, 225)
    fig.suptitle('Photosphere band %.1f–%.1f R_ph,IC unwrapped, native cells (%d radial cells), t = %.1f d'
                 % (PB[0], PB[1], i1 - i0, float(d['t']) / DAY), fontsize=13)
    fig.savefig(os.path.join(OUT, fn), dpi=120)
    plt.close(fig)


def mapfig(d, fn):
    tags = [k[6:] for k in d.files if k.startswith('sh_vr_')]
    want = ['011', '030', '090', '098', '1085']
    tags = [t for t in want if t in tags]
    lab = lambda t: '%.3g' % (int(t) / (100 if len(t) == 3 else 1000))  # noqa: E731
    fig, axs = plt.subplots(4, 3, figsize=(18, 13.5), constrained_layout=True)
    pans = [('sh_vr_' + t, 'v_r at %s R_ph [km/s]' % lab(t)) for t in tags[:5]] + [('sh_Rph', 'R_ph / R_ph,IC')] + \
           [('sh_dT_' + t, 'δT_gas/⟨T_gas⟩ at %s R_ph [%%]' % lab(t)) for t in tags[:5]] + [('sh_Teff', 'T_eff [K]')]
    for a, (k, lb) in zip(axs.flat, pans):
        q = np.asarray(d[k], float)
        fl = None
        if k.startswith('sh_vr'):
            mu = cs.shell_mean(q)
            if abs(mu) > 2 * np.nanstd(q):              # shell moving as a whole: show v_r - <v_r>
                q = q - mu; lb = lb.replace('v_r', 'v_r − ⟨v_r⟩').replace('[km/s]', '[km/s]; ⟨v_r⟩ = %.3g km/s' % mu)
            v = max(1e-4, np.nanpercentile(abs(q), 99.5)); cm, nm = 'RdBu_r', cs.norm_of('sym', v)
        elif k.startswith('sh_dT'):
            q = 100 * q; fl = np.asarray(d['sh_Tfl_' + k[6:]])
            v = max(1e-4, np.nanpercentile(abs(q[~fl]), 99.5) if (~fl).any() else 1); cm, nm = cs.cmap_floor('PuOr_r'), cs.norm_of('sym', v)
        elif k == 'sh_Rph':
            cm, nm = 'cividis', cs.norm_of('lin', pct(q, 1, 99))
        else:
            cm, nm = 'inferno', cs.norm_of('lin', pct(q, 1, 99))
        m = cs.moll(a, q, cm, nm, W=900, cache=OUT, floor=fl)
        cb = fig.colorbar(m, ax=a, orientation='horizontal', shrink=0.7, pad=0.02, label=lb, extend='both')
        cb.ax.xaxis.set_major_locator(MaxNLocator(5))
    fig.suptitle('Mollweide maps on native cells (each 0.70° cell drawn; graticule every 45°, cube edges at lon ±45°, ±135°), '
                 't = %.1f d' % (float(d['t']) / DAY), fontsize=13)
    fig.savefig(os.path.join(OUT, fn), dpi=100)
    plt.close(fig)


for f in glob.glob(os.path.join(OUT, 'fig_slice_*')) + glob.glob(os.path.join(OUT, 'fig_map_*')) + glob.glob(os.path.join(OUT, 'fig_phot_*')):
    os.remove(f)
SL, MP, PH = [], [], []
R2 = [i for i, d in enumerate(R) if 'ver' in d.files and int(d['ver']) >= 2]
if R2:
    p = R2[-1]
    slicefig(R[p], 'fig_slice_%05d.png' % p); SL.append(('fig_slice_%05d.png' % p, float(tr[p])))
    photfig(R[p], 'fig_phot_%05d.png' % p); PH.append(('fig_phot_%05d.png' % p, float(tr[p])))
    mapfig(R[p], 'fig_map_%05d.png' % p); MP.append(('fig_map_%05d.png' % p, float(tr[p])))

# ---------------- status ----------------
_lg = {}
for dd in [RUN] + [x for x in os.environ.get('EXTRA_RUN', '').split(':') if x]:   # same log name in two dirs: RUN wins
    for f in glob.glob(os.path.join(dd, 'run*.log')):
        _lg.setdefault(os.path.basename(f), f)
LOGS = sorted(_lg.values(), key=os.path.getmtime)
log = ''.join(open(f, errors='replace').read() for f in LOGS)
LOGDIR = os.environ.get('LOGDIR', '')
lks = sorted(set(glob.glob(os.path.join(RUN, 'link.*.out')) + sum((glob.glob(os.path.join(x, 'link.*.out')) for x in LOGDIR.split(':') if x), [])), key=os.path.getmtime)
lk = ''.join(open(f, errors='replace').read() for f in lks)
el = re.findall(r'elapsed=([0-9.eE+-]+) cycle=([0-9]+) time=([0-9.eE+-]+)', log)
spc = np.nan
if len(el) > 3:
    (e0, c0, _), (e1, c1, _) = el[-min(len(el), 20)], el[-1]
    spc = (float(e1) - float(e0)) / max(1, int(c1) - int(c0))
cyc = int(el[-1][1]) if el else len(td)
try:
    sq = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', 'jinma'), '-n', os.environ.get('JOBNAME', 'rsgE21'), '-o', '%i|%T|%M|%l'],
                        capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
except Exception:
    sq = []
EV = {}
try:   # last row of the event counters (rsg.log; counts per 100 cycles, cell-stages)
    evl = [ln for ln in open(os.path.join(RUN, 'rsg.log')) if ln.strip()]
    hdr = [ln for ln in evl if ln.startswith('#  cycle')][-1][1:].split()
    rows = [ln.split() for ln in evl if not ln.startswith('#') and len(ln.split()) == len(hdr)]
    EV = dict(zip(hdr, rows[-1])) if rows else {}
except Exception:
    EV = {}
nonfin = [int(d['nonfinite']) for d in R]
T = dict(t_d=float(td[-1]), t_s=float(H['time'][-1]), cyc=cyc, spc=float(spc), dt=float(H['dt'][-1]),
         pic=float(np.mean(H['Picard'][1:])), picmax=float(np.max(H['Picard'])),
         LL=float(H['L_top'][-1] / H['L_in'][-1]), Lbot=float(H['L_bot'][-1] / H['L_in'][-1]),
         ke=float(H['KE_int'][-1]), dE=float((H['Etot'][-1] - E0) / abs(E0)), dEb=float(ct[-1] / abs(E0)),
         dM=float((H['M_tot'][-1] - M0) / M0),
         ndump=len(R), t_dump=float(tr[-1]) if len(tr) else 0.,
         rph=float(R[-1]['Rph_med']) / RPH0 if R else 0., teff=float(R[-1]['Teff_med']) if R else 0.,
         nonfinite=int(sum(nonfin)), vmax=float(R[-1]['vmax']) / 1e5 if R else 0.,
         fatal=len(re.findall(r'FATAL', log + lk)), div=len(re.findall(r'DIVERGED', log + lk)),
         nan=len(re.findall(r'\bnan\b|\bNaN\b', log + lk)), nc=len(re.findall(r'NON-CONVERGED=[1-9]', log + lk)),
         jobs=[s.split('|') for s in sq], slices=SL, maps=MP, phot=PH, ev=EV, nfloorT=(int(R[R2[-1]]['nfloorT']) if R2 else -1), ncell=(int(np.prod(R[R2[-1]]['sh_Rph'].shape)) * len(R[R2[-1]]['r']) if R2 else 0), jobname=os.environ.get('JOBNAME', 'rsgE21'))
json.dump(T, open(os.path.join(OUT, 'table.json'), 'w'), indent=1)
print('figs done', T)
