#!/usr/bin/env python3
"""HGPAGE 10-18: figures + run table for the He giant 14p5M production page (READ ONLY on the run).
usage: mk_figs.py <red_dir> <img_dir>    (reads hst/log of the Raven run and the red_*.npz reductions)"""
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
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

RUN = os.environ.get('RUN', '/raven/ptmp/jinma/hg145_prod_1018/prod/R')          # dir with hegiant.user.hst + run.log
ICF = os.environ.get('ICF', '/raven/ptmp/jinma/ic_star_runs/hegiant_14p5M_mk/ic_hegiant_14p5M_mk.txt')  # IC column (deck he_ic_file)
RED, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
RS, DAY, LSUN = 6.957e10, 86400., 3.828e33
RPH0 = 7.94653e12 / RS          # IC photosphere [Rsun]
OMEGA = np.pi / 3               # wedge solid angle: (pi/3) x (cos pi/3 - cos 2pi/3)
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
                     'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'axes.titlesize': 11,
                     'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
                     'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'})
TC = ['#9aa4ae', '#4cc3ff', '#7dff8a', '#ffb347', '#ff5c7a']   # time colours, IC first


def readhst(f):
    lab = None
    for ln in open(f):
        if ln.startswith('#  ['):
            lab = [x.split('=')[1] for x in ln[1:].split()]
            break
    d = np.atleast_2d(np.loadtxt(f, comments='#'))
    # restarts rewrite times: keep the last occurrence (scan backwards, keep strictly decreasing t)
    t = d[:, 0]
    keep = np.zeros(len(t), bool)
    tmin = np.inf
    for i in range(len(t) - 1, -1, -1):
        if t[i] < tmin:
            keep[i] = True
            tmin = t[i]
    d = d[keep]
    return {n: d[:, i] for i, n in enumerate(lab)}


H = readhst(os.path.join(RUN, 'hegiant.user.hst'))
td = H['time'] / DAY
reds = sorted(glob.glob(os.path.join(RED, 'red_*.npz')))
R = [np.load(f) for f in reds]
tr = np.array([float(d['t']) for d in R]) / DAY
rph_med = np.array([np.nanmedian(d['Rph']) for d in R]) / RS
rph_lo = np.array([np.nanpercentile(d['Rph'], 2) for d in R]) / RS
rph_hi = np.array([np.nanpercentile(d['Rph'], 98) for d in R]) / RS
teff_med = np.array([np.nanmedian(d['Teff']) for d in R])


def smooth(y, n):
    if len(y) < 2 * n:
        return y
    c = np.cumsum(np.insert(y, 0, 0.))
    s = (c[n:] - c[:-n]) / n
    return np.concatenate([np.full(n - 1, np.nan), s])


# ---------------- history ----------------
fig, ax = plt.subplots(2, 3, figsize=(15, 7.6), constrained_layout=True)
a = ax[0, 0]
a.plot(td, H['L_top'] / H['L_in'], color='#ffb347', lw=1.2, label='L_top / L_in')
a.plot(td, H['L_bot'] / H['L_in'], color='#4cc3ff', lw=1.0, label='L_bot / L_in')
a.plot(td, H['L_MLT'] / H['L_in'], color='#7dff8a', lw=1.0, label='L_MLT / L_in')
a.axhline(1, color='#a3adb8', lw=0.6, ls=':')
a.set_ylim(0.8, 1.1)
a.set_ylabel('luminosity / L_in'); a.legend(loc='lower right'); a.set_title('luminosities (wedge, L_in = imposed bottom)')
a = ax[0, 1]
a.semilogy(td, np.maximum(H['KE_int'], 1e30), color='#ff5c7a', lw=1.2, label='KE_int')
a.semilogy(td, np.maximum(H['KEr_int'], 1e30), color='#c77dff', lw=1.0, label='KE_r,int (radial)')
a.set_ylabel('erg (wedge)'); a.legend(); a.set_title('kinetic energy inside the star')
a = ax[0, 2]
if len(tr):
    a.fill_between(tr, rph_lo, rph_hi, color='#4cc3ff', alpha=0.25, lw=0, label='2-98 % over the wedge')
    a.plot(tr, rph_med, 'o-', ms=2.5, color='#4cc3ff', lw=1.2, label='median R_ph')
a.axhline(RPH0, color='#a3adb8', lw=0.6, ls=':', label='IC R_ph 114.2')
a.set_ylabel('R_ph [R_sun]'); a.legend(loc='lower left'); a.set_title('photosphere radius (from the 0.25 d dumps)')
a2 = a.twinx()
if len(tr):
    a2.plot(tr, teff_med / 1e3, color='#ffb347', lw=0.9, ls='--')
a2.set_ylabel('median T_eff [kK] (dashed)', color='#ffb347')
a = ax[1, 0]
n = max(1, min(200, len(td) // 50))
a.plot(td, H['Picard'], color='#3a424c', lw=0.4)
a.plot(td, smooth(H['Picard'], n), color='#e8ecf1', lw=1.2, label='running mean (%d cycles)' % n)
a.set_ylabel('Picard iterations per cycle'); a.legend(); a.set_title('implicit radiation solve effort')
a = ax[1, 1]
a.plot(td, H['dt'], color='#7dff8a', lw=1.0)
a.set_ylabel('dt [s]'); a.set_title('time step (cap 100 s)')
a = ax[1, 2]
a.semilogy(td, np.maximum(np.abs(H['Mdot_top']), 1e8), color='#66d9c6', lw=1.0, label='|Mdot_top| [g/s]')
a.semilogy(td, np.maximum(H['P_sponge'], 1e20), color='#c77dff', lw=1.0, label='P_sponge [erg/s]')
a.set_title('outer boundary: mass loss and sponge power'); a.legend()
for a in ax.flat:
    a.set_xlabel('t [d]'); a.grid(True, lw=0.5)
fig.savefig(os.path.join(OUT, 'fig_hist.png'), dpi=95)
plt.close(fig)

# ---------------- radial profiles ----------------
ic = np.loadtxt(ICF, comments='#')
icr, icrho, icT, icfm = ic[:, 0] / RS, ic[:, 1], ic[:, 5], ic[:, 6]
if len(R):
    nt = len(R)
    pick = sorted(set([0] + [int(round(x * (nt - 1))) for x in (0.1, 0.4, 1.0)]))
    pick = [p for p in pick if p < nt]
    if len(pick) > 1 and tr[pick[1]] < 0.9 and nt > 5:   # first shown time >= ~1 d when available
        k = int(np.argmin(abs(tr - 1.0)))
        pick = sorted(set([0, k] + pick[2:]))
    d0 = R[0]
    fig, ax = plt.subplots(2, 4, figsize=(16, 7.8), constrained_layout=True)
    for row, (xlo, xhi) in enumerate([(8, 1.5 * RPH0), (0.85 * RPH0, 1.12 * RPH0)]):
        for c, p in enumerate(pick):
            d = R[p]
            r = d['r'] / RS
            lab = 't = %.2f d' % tr[p]
            col = TC[c % len(TC)]
            Lin = np.interp(float(d['t']), H['time'], H['L_in'])
            Ftot = Lin / (OMEGA * d['r'] ** 2)
            if row == 0:
                ax[row, 0].semilogy(r, d['rho'], color=col, lw=1.1, label=lab)
                ax[row, 1].semilogy(r, d['Trad'], color=col, lw=1.1, label=lab)
            else:
                ax[row, 0].plot(r, d['rho'] / d0['rho'], color=col, lw=1.1, label=lab)
                ax[row, 1].plot(r, d['Trad'] / d0['Trad'], color=col, lw=1.1, label=lab)
            if p > 0:
                ax[row, 2].semilogy(r, np.maximum(d['vrrms'] / 1e5, 1e-5), color=col, lw=1.1, label=lab)
            ax[row, 3].semilogy(r, np.maximum(abs(d['Fenth'] / Ftot), 1e-7), color=col, lw=1.1, label=lab + ' (resolved)')
        if row == 0:
            ax[row, 0].semilogy(icr, icrho, color='#ffffff', lw=0.8, ls='--', label='IC column (MESA-anchored)')
            ax[row, 1].semilogy(icr, icT, color='#ffffff', lw=0.8, ls='--', label='IC column T_gas')
        ax[row, 3].semilogy(icr, np.maximum(icfm, 1e-7), color='#ffffff', lw=0.8, ls='--', label='IC scaffold F_MLT/F')
        for a in ax[row]:
            a.set_xlim(xlo, xhi); a.grid(True, lw=0.5); a.set_xlabel('r [R_sun]')
            a.axvline(RPH0, color='#a3adb8', lw=0.6, ls=':')
        ax[row, 3].set_ylim(1e-6, 2)
    ax[0, 0].set_ylabel('shell-mean ρ [g cm⁻³]'); ax[0, 1].set_ylabel('T_rad = (E/a)^¼ [K]')
    ax[1, 0].set_ylabel('ρ / ρ(t=0)'); ax[1, 1].set_ylabel('T_rad / T_rad(t=0)')
    ax[0, 2].set_ylabel('mass-weighted rms v_r [km/s]'); ax[1, 2].set_ylabel('rms v_r [km/s]')
    ax[0, 3].set_ylabel('|F_enth| / (L_in / Ω r²)'); ax[1, 3].set_ylabel('|F_enth| / F_tot')
    ax[1, 0].set_ylim(0.0, 2.5); ax[1, 1].set_ylim(0.7, 1.3)
    ax[0, 0].set_ylim(1e-15, 1e-3)
    ax[0, 0].set_title('density'); ax[0, 1].set_title('radiation temperature'); ax[0, 2].set_title('radial velocity rms')
    ax[0, 3].set_title('convective (enthalpy) flux fraction')
    ax[0, 0].legend(loc='lower left'); ax[0, 3].legend(loc='upper left')
    fig.savefig(os.path.join(OUT, 'fig_profiles.png'), dpi=90)
    plt.close(fig)

# ---------------- photosphere maps ----------------
if len(R):
    lst = [R[int(round(0.5 * (len(R) - 1)))], R[-1]] if len(R) > 2 else [R[-1]]
    fig, ax = plt.subplots(len(lst), 2, figsize=(12, 5.0 * len(lst)), constrained_layout=True, squeeze=False)
    for i, d in enumerate(lst):
        lat = 90 - np.degrees(d['th']); lon = np.degrees(d['ph'])
        Rp = d['Rph'] / RS
        med = np.nanmedian(Rp)
        dev = 100 * (Rp / med - 1)
        v = max(0.05, np.nanpercentile(abs(dev), 99.5))
        m = ax[i, 0].pcolormesh(lon, lat, dev.T, cmap='RdBu_r', norm=TwoSlopeNorm(0, -v, v), shading='nearest')
        fig.colorbar(m, ax=ax[i, 0], label='R_ph / median − 1 [%]')
        ax[i, 0].set_title('R_ph map, t = %.2f d (median %.2f R_sun)' % (float(d['t']) / DAY, med))
        dT = d['Teff'] - np.nanmedian(d['Teff'])
        vt = max(0.5, np.nanpercentile(abs(dT), 99.5))
        m = ax[i, 1].pcolormesh(lon, lat, dT.T, cmap='RdBu_r', norm=TwoSlopeNorm(0, -vt, vt), shading='nearest')
        fig.colorbar(m, ax=ax[i, 1], label='T_eff − median [K]')
        ax[i, 1].set_title('T_eff map, t = %.2f d (median %.0f K)' % (float(d['t']) / DAY, np.nanmedian(d['Teff'])))
        for a in ax[i]:
            a.set_aspect('equal'); a.set_xlabel('φ [deg]'); a.set_ylabel('latitude [deg]')
    fig.savefig(os.path.join(OUT, 'fig_rph.png'), dpi=85)
    plt.close(fig)


# ---------------- Cartesian slices ----------------
def edges(c):
    e = np.empty(len(c) + 1)
    e[1:-1] = 0.5 * (c[1:] + c[:-1]); e[0] = c[0] - (e[1] - c[0]); e[-1] = c[-1] + (c[-1] - e[-2])
    return e


def slicefig(d, d0, fn):
    ns = int(d['ns'])
    rf = d['rf'][:ns + 1] / RS
    the, phe = edges(d['th']), edges(d['ph'])
    fig, ax = plt.subplots(4, 4, figsize=(17, 17), constrained_layout=True)
    for row in range(4):
        eqp = row < 2
        A, RR = np.meshgrid(phe if eqp else the, rf, indexing='ij')
        if eqp:
            X, Y = RR * np.cos(A), RR * np.sin(A)
            key = 'eq'
        else:
            X, Y = RR * np.sin(A), RR * np.cos(A)
            key = 'me'
        lr = np.asarray(d[key + '_lrho'], float); lr0 = np.asarray(d0[key + '_lrho'], float)
        T = np.asarray(d[key + '_T'], float); vr = np.asarray(d[key + '_vr'], float)
        zoom = row % 2 == 1
        if zoom:
            vv = max(0.05, np.nanpercentile(abs(vr[:, rf[1:] > 0.85 * RPH0]), 99.5))
        else:
            vv = max(0.05, np.nanpercentile(abs(vr), 99.8))
        dl = lr - lr0
        dv = max(0.01, np.nanpercentile(abs(dl[:, rf[1:] < 1.12 * RPH0]), 99.5))
        pans = [(lr, 'viridis', None, 'log ρ [g cm⁻³]', (-15, -3.5) if not zoom else (-13, -6)),
                (dl, 'RdBu_r', TwoSlopeNorm(0, -dv, dv), 'log ρ/ρ(t=0)', None),
                (T, 'inferno', None, 'T_rad [kK]', (5, 200) if not zoom else (5, 30)),
                (vr, 'RdBu_r', TwoSlopeNorm(0, -vv, vv), 'v_r [km/s]', None)]
        for c, (q, cm, nm, lab, lim) in enumerate(pans):
            a = ax[row, c]
            kw = dict(cmap=cm, shading='flat', rasterized=True)
            if nm is not None:
                kw['norm'] = nm
            elif lim is not None:
                kw['vmin'], kw['vmax'] = lim
            if c == 2:
                from matplotlib.colors import LogNorm
                kw.pop('vmin'); kw.pop('vmax'); kw['norm'] = LogNorm(*lim)
            m = a.pcolormesh(X, Y, q, **kw)
            fig.colorbar(m, ax=a, label=lab, shrink=0.85)
            a.set_aspect('equal')
            if zoom:
                if eqp:   # box around the photosphere on the mid-phi ray
                    cx, cy = RPH0 * np.cos(np.pi / 6), RPH0 * np.sin(np.pi / 6)
                else:
                    cx, cy = RPH0, 0.0
                a.set_xlim(cx - 18, cx + 18); a.set_ylim(cy - 18, cy + 18)
            a.set_xlabel('x [R_sun]'); a.set_ylabel('y [R_sun]' if eqp else 'z [R_sun]')
        ax[row, 0].set_title(('equatorial plane θ = π/2' if eqp else 'meridional plane φ = 30°') + (' (zoom at R_ph)' if zoom else ''))
    fig.suptitle('t = %.2f d' % (float(d['t']) / DAY), fontsize=14)
    fig.savefig(os.path.join(OUT, fn), dpi=62, pil_kwargs={'quality': 85})
    plt.close(fig)


for f in glob.glob(os.path.join(OUT, 'fig_slice_*.jpg')):
    os.remove(f)
SL = []
if len(R) > 1:
    for p in sorted(set([int(round(0.5 * (len(R) - 1))), len(R) - 1])):
        fn = 'fig_slice_%05d.jpg' % p
        slicefig(R[p], R[0], fn)
        SL.append((fn, tr[p]))

# ---------------- run table / status ----------------
LOGS = sorted(glob.glob(os.path.join(RUN, os.environ.get('LOGGLOB', 'run*.log'))), key=os.path.getmtime)
log = ''.join(open(f, errors='replace').read() for f in LOGS)
el = re.findall(r'elapsed=([0-9.eE+-]+) cycle=([0-9]+) time=([0-9.eE+-]+)', log)
spc = np.nan
if len(el) > 3:
    (e0, c0, _), (e1, c1, _) = el[-min(len(el), 20)], el[-1]
    spc = (float(e1) - float(e0)) / max(1, int(c1) - int(c0))
try:
    sq = subprocess.run(['sacct', '-X', '-n', '-P', '--name=' + os.environ.get('JOBNAME', 'hgP145'), '-S', '2026-10-18', '-o', 'JobID,State,Elapsed'],
                        capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
except Exception:
    sq = []
T = dict(t_d=float(td[-1]), cyc=int(len(td)), spc=float(spc), pic=float(np.mean(H['Picard'][1:])) if len(td) > 1 else 0.,
         picmax=float(np.max(H['Picard'])), LL=float(H['L_top'][-1] / H['L_in'][-1]), dt=float(H['dt'][-1]),
         ndump=len(R), t_dump=float(tr[-1]) if len(tr) else 0., rph=float(rph_med[-1]) if len(tr) else 0.,
         teff=float(teff_med[-1]) if len(tr) else 0., ke=float(H['KE_int'][-1]),
         fatal=log.count('FATAL'), div=log.count('DIVERGED'), nc=len(re.findall(r'NON-CONVERGED=[1-9]', log)),
         jobs=[s.split('|') for s in sq], slices=SL)
json.dump(T, open(os.path.join(OUT, 'table.json'), 'w'), indent=1)
print('figs done', T['t_d'], T['ndump'])
