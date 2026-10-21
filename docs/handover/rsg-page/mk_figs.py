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

RUN = os.environ.get('RUN', '/raven/ptmp/jinma/rsg_g25_1019/prod/R')
ICF = os.environ.get('ICF', '/raven/ptmp/jinma/rsg_g25_1019/ic/pass3b/ic_rsg_g25.txt')
OLD = os.environ.get('OLDRUN', '/raven/ptmp/jinma/rsg_g25_1019/prod/R_oldseed')
RED, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
RS, DAY, LSUN = 6.957e10, 86400., 3.828e33
RPH0 = 5.80481e13
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
if os.path.exists(os.path.join(OLD, 'rsg.user.hst')):
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
a.set_ylabel('R_ph / R_ph,IC (834 R_sun)'); a.legend(loc='lower left'); a.set_title('photosphere (from the 1 d dumps)')
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
        ax[0, 1].loglog(x, d['Trad'], color=col, lw=1, label=lab)
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
    ax[0, 1].set_ylabel('⟨T_rad⟩ [K]'); ax[0, 1].set_title('radiation temperature (shell mean)')
    ax[0, 2].set_ylabel('km/s'); ax[0, 2].set_title('ρ-weighted rms v_r (solid), rms |v_t| (dashed)')
    ax[1, 0].set_ylabel('L_rad(r) / L_in'); ax[1, 0].set_title('radiative luminosity 4πr²⟨F_r⟩'); ax[1, 0].set_ylim(-0.2, 1.6)
    ax[1, 1].set_ylabel('L_conv,res(r) / L_in'); ax[1, 1].set_title('resolved enthalpy flux (ideal-gas h)'); ax[1, 1].set_ylim(-0.3, 1.2)
    ax[1, 2].set_yscale('log'); ax[1, 2].set_ylabel('rms δ log ρ'); ax[1, 2].set_title('lateral density fluctuation per shell')
    for a in ax.flat:
        a.axvline(1, color='#a3adb8', lw=0.5, ls=':')
        a.set_xlabel('r / R_ph,IC'); a.grid(True, lw=0.5); a.legend(fontsize=7)
    fig.savefig(os.path.join(OUT, 'fig_profiles.png'), dpi=85)
    plt.close(fig)


# ---------------- slices ----------------
def slicefig(d, fn):
    rf = np.asarray(d['rf'], float) / RPH0
    ns = int(d['ns'])
    ang = np.asarray(d['ang'])
    na = len(ang)
    af = np.linspace(0, 2 * np.pi, na + 1)
    AA, RR = np.meshgrid(af, rf[:ns + 1], indexing='ij')
    fig, ax = plt.subplots(3, 4, figsize=(17, 12.6), constrained_layout=True)
    rc = 0.5 * (rf[1:ns + 1] + rf[:ns])
    for row, key in enumerate(('eq', 'me', 'eq')):
        zoom = row == 2
        if key == 'eq':
            X, Y = RR * np.cos(AA), RR * np.sin(AA)
        else:
            X, Y = RR * np.sin(AA), RR * np.cos(AA)
        lr = np.asarray(d[key + '_lrho'], float); dr = np.asarray(d[key + '_drho'], float)
        T = np.asarray(d[key + '_T'], float); vr = np.asarray(d[key + '_vr'], float)
        reg = (rc < ZR) if zoom else (rc > ZR) & (rc < 1.0)
        vv = max(0.02, np.nanpercentile(abs(vr[:, reg]), 99.5))
        dv = max(0.0005, np.nanpercentile(abs(dr[:, reg]), 99.5))
        lrl = (np.nanpercentile(lr[:, reg], 0.5), np.nanpercentile(lr[:, reg], 99.5)) if zoom else (-17, -6)
        Tl = (np.nanpercentile(T[:, reg], 0.5), np.nanpercentile(T[:, reg], 99.5)) if zoom else (2, 1000)
        pans = [(lr, 'viridis', None, 'log ρ [g cm⁻³]', lrl),
                (dr, 'RdBu_r', TwoSlopeNorm(0, -dv, dv), 'log ρ − log⟨ρ⟩_shell', None),
                (T, 'inferno', LogNorm(*Tl), 'T_rad [kK]', None),
                (vr, 'RdBu_r', TwoSlopeNorm(0, -vv, vv), 'v_r [km/s]', None)]
        for c, (q, cm, nm, lab, lim) in enumerate(pans):
            a = ax[row, c]
            kw = dict(cmap=cm, shading='flat', rasterized=True)
            if nm is not None:
                kw['norm'] = nm
            else:
                kw['vmin'], kw['vmax'] = lim
            m = a.pcolormesh(X, Y, q, **kw)
            fig.colorbar(m, ax=a, label=lab, shrink=0.85)
            a.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec='w', lw=0.4, ls=':'))
            if zoom:
                a.add_patch(plt.Circle((0, 0), 0.0494, fill=False, ec='w', lw=0.5, ls='--'))
                a.set_xlim(-ZR, ZR); a.set_ylim(-ZR, ZR)
            a.set_aspect('equal')
            a.set_xlabel('x / R_ph,IC'); a.set_ylabel('y / R_ph,IC' if key == 'eq' else 'z / R_ph,IC')
        ax[row, 0].set_title(('equatorial plane z = 0' if key == 'eq' else 'meridional plane y = 0 (x–z)') +
                             (' — zoom r < %.2f R_ph (dashed: MESA CZ base)' % ZR if zoom else ''))
    fig.suptitle('t = %.1f d  (dotted circle: R_ph of the IC, 834 R_sun)' % (float(d['t']) / DAY), fontsize=13)
    fig.savefig(os.path.join(OUT, fn), dpi=60, pil_kwargs={'quality': 84})
    plt.close(fig)


def mapfig(d, fn):
    keys = [k for k in ('map_vr_007', 'map_vr_030', 'map_vr_090', 'map_vr_098') if k in d.files][:4]
    lat = np.radians(np.linspace(-90, 90, 181)); lon = np.radians(np.linspace(-180, 180, 361))
    fig, axs = plt.subplots(2, 3, figsize=(16, 6.8), subplot_kw=dict(projection='mollweide'), constrained_layout=True)
    pans = [(k, 'v_r at %.2f R_ph [km/s]' % (int(k[-3:]) / 100)) for k in keys] + [('map_Rph', 'R_ph / R_ph,IC'), ('map_Teff', 'T_eff [K]')]
    for a, (k, lab) in zip(axs.flat, pans):
        q = np.asarray(d[k], float)
        if k.startswith('map_vr'):
            v = max(0.02, np.nanpercentile(abs(q), 99.5)); kw = dict(cmap='RdBu_r', norm=TwoSlopeNorm(0, -v, v))
        elif k == 'map_Rph':
            lo, hi = np.nanpercentile(q, [1, 99]); kw = dict(cmap='cividis', vmin=lo, vmax=hi)
        else:
            lo, hi = np.nanpercentile(q, [1, 99]); kw = dict(cmap='inferno', vmin=lo, vmax=hi)
        m = a.pcolormesh(lon, lat, q, shading='flat', rasterized=True, **kw)
        a.set_xticklabels([]); a.set_yticklabels([]); a.grid(True, lw=0.3, color='#555')
        fig.colorbar(m, ax=a, orientation='horizontal', shrink=0.7, pad=0.03, label=lab)
    fig.suptitle('Mollweide maps, t = %.1f d' % (float(d['t']) / DAY), fontsize=13)
    fig.savefig(os.path.join(OUT, fn), dpi=60, pil_kwargs={'quality': 84})
    plt.close(fig)


for f in glob.glob(os.path.join(OUT, 'fig_slice_*.jpg')) + glob.glob(os.path.join(OUT, 'fig_map_*.jpg')):
    os.remove(f)
SL, MP = [], []
if len(R):
    picks = sorted(set([len(R) - 1, int(np.argmin(abs(tr - (tr[-1] - 10)))), int(np.argmin(abs(tr - (tr[-1] - 20))))]))
    for p in picks:
        fn = 'fig_slice_%05d.jpg' % p
        slicefig(R[p], fn)
        SL.append((fn, float(tr[p])))
    for p in picks[-2:]:
        fn = 'fig_map_%05d.jpg' % p
        mapfig(R[p], fn)
        MP.append((fn, float(tr[p])))

# ---------------- status ----------------
LOGS = sorted(glob.glob(os.path.join(RUN, 'run*.log')), key=os.path.getmtime)
log = ''.join(open(f, errors='replace').read() for f in LOGS)
lks = sorted(glob.glob(os.path.join(RUN, 'link.*.out')), key=os.path.getmtime)
lk = ''.join(open(f, errors='replace').read() for f in lks)
el = re.findall(r'elapsed=([0-9.eE+-]+) cycle=([0-9]+) time=([0-9.eE+-]+)', log)
spc = np.nan
if len(el) > 3:
    (e0, c0, _), (e1, c1, _) = el[-min(len(el), 20)], el[-1]
    spc = (float(e1) - float(e0)) / max(1, int(c1) - int(c0))
cyc = int(el[-1][1]) if el else len(td)
try:
    sq = subprocess.run(['squeue', '-h', '-u', os.environ.get('USER', 'jinma'), '-n', 'rsgP25', '-o', '%i|%T|%M|%l'],
                        capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
except Exception:
    sq = []
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
         jobs=[s.split('|') for s in sq], slices=SL, maps=MP)
json.dump(T, open(os.path.join(OUT, 'table.json'), 'w'), indent=1)
print('figs done', T)
