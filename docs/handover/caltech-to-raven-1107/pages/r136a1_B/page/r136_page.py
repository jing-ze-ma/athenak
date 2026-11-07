"""r136_page.py: figures, movie frames and status.json for the R136a1 phase-1 run page (r136a1-1104).

Reads the run directory (default run128/prod): r136.user.hst, run.log, bin/r136.hydro_w.NNNNN.bin + r136.m1.NNNNN.bin.
Every new dump pair is reduced once to cache/red_NNNNN.npz (shell profiles, Mdot_r, slices); figures are redrawn
from the cache. CPU only (run it inside a held allocation, see refresh.sh).

usage: python3 r136_page.py PAGEDIR [--full]   (--full: re-choose colour ranges and redraw every movie frame)

Definitions (shell means are solid-angle weighted over the wedge, Omega = pi/3 sr; wedge -> star x 4 pi / Omega):
  T          from (rho, e_gas) with the run's general EOS (the IC tool's dump of the code table, X 0.58 Y 0.4136)
  kappa_R    the run's Rosseland table (MESA a09 / TOPS), bilinear in (log T, log rho)
  tau        per column from the top inward, sum kappa_R rho dr; R_ph = r where the shell-mean tau crosses 1
             (as r_photo,3D on the BSG pages); T_eff = (L_top / (4 pi sigma R_ph^2))^(1/4), L_top from the history
  F0         comoving radial flux F0 = F_r - v_r (E + P_rr), F_r = lab-frame M1 flux m1_f1, P_rr = D_rr E from the m1_vet dump
             at the same time where one exists, else E/3 (11-05 USER: the lab flux carries the closed-wall odd-even v_r mode
             as a +-25-30 % zigzag through (4/3) E v; F0 is smooth).  Deep (G = (4/3) kappa E dt/c > 1e-3, r136kf-1106): from the
             E gradient, -(c / 3 chi_f) dE/dr (harmonic faces), the solver's own thick-limit F0: after a backward-Euler-path
             step F holds the velocity before the radiative kick while the dump's v holds it, and F - (E + P_rr) v is then
             F0 - G (F0 - F_wb), a grid-looking noise of up to 50 x at the wall (see r136_f0.py)
  Gamma_rad  shell mean of kappa_R F0 / (c g), g = G m(r) / r^2, m(r) = m(x1min) + mass of the dump below r (lab: faint)
  Gamma_L    kappa_R(shell) L / (4 pi c G m(r)), L = the star's luminosity
  P_rad/P    (E/3) / (E/3 + P_gas)
  L_rad      4 pi r^2 <F0> (comoving; lab 4 pi r^2 <F_r> faint);  L_conv = 4 pi r^2 [<(e_gas + P_gas)(v_r - v~)> + <(E + P_rr) v_r>],
             v~ = <rho v_r>/<rho>: the advected radiation enthalpy is booked with the convective luminosity (HG/RSG convention),
             so L_rad + L_conv = L_rad,lab + L_conv,gas exactly
  L_MLT      the frozen scaffold of the IC column (F_MLT/F x L) times the history's MLT ramp weight w_mlt
  Mdot_r     4 pi r^2 <rho v_r>  (> 0 outward), shown over M_env = 4 pi / Omega x wedge mass above x1min
"""
import glob
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pipe'))
import bin_convert_repro as bc  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r136_f0 import f0_solver  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm, SymLogNorm  # noqa: E402

W = os.path.dirname(os.path.abspath(__file__))
R = '/resnick/groups/carnegie_poc/jingze/r136a1_1104'
RUN = os.environ.get('R136_RUN', R + '/run128/prod')
CACHE = os.path.join(W, 'cache', os.path.basename(os.path.dirname(RUN)) + '_' + os.path.basename(RUN))
OPAC = R + '/opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase.txt'
EOSF = R + '/ic_cache/tables/eos_8f658626c61e/eos_dump.txt'
ICF = R + '/ic/ic_r136a1.txt'
RSUN, DAY, C, G, SIG = 6.957e10, 86400.0, 2.99792458e10, 6.674e-8, 5.670374e-5
LSTAR = 2.60446120e40           # MESA photosphere L (IC column, L const)
GM0 = float(os.environ.get('R136_GM0', '3.1056270511e28'))   # he_gm = G m(x1min) of the run's deck
ICF = os.environ.get('R136_IC', ICF)
X1MIN = 23.1                    # Rsun, set from the run's grid in main()
TLIM = float(os.environ.get('R136_TLIM', '1.728e6'))
LABEL = os.environ.get('R136_LABEL', 'R136a1 phase 1 (128$^2$)')
OMEGA = np.pi / 3 * (np.cos(np.pi / 3) - np.cos(2 * np.pi / 3))   # = pi/3 sr
FAC = 4 * np.pi / OMEGA
RMAX_PLOT = 32.0   # Rsun, slices/profiles outer limit: auto = max(32, 1.25 R_ph newest, 5-Rsun steps); R136_RMAX overrides


# ------------------------------------------------------------------ microphysics
class Table:
    """AthenaK stellar opacity table (log10 kappa on uniform log T x log rho)."""
    def __init__(self, fn):
        hdr = [ln for ln in open(fn) if ln.startswith('#')]
        g = [ln for ln in hdr if len(ln.split()) == 7 and ln.split()[1].isdigit()][0].split()[1:]
        self.nT, self.nD = int(g[0]), int(g[1])
        self.t0, self.dt, self.d0, self.dd = map(float, g[2:])
        self.v = np.loadtxt(fn, comments='#', dtype=np.float32).reshape(self.nT, self.nD)

    def __call__(self, rho, T):
        x = np.clip((np.log10(T) - self.t0) / self.dt, 0, self.nT - 1.001)
        y = np.clip((np.log10(rho) - self.d0) / self.dd, 0, self.nD - 1.001)
        i = x.astype(np.int32); j = y.astype(np.int32); fx = (x - i).astype(np.float32); fy = (y - j).astype(np.float32)
        v = self.v
        lk = ((1 - fx) * (1 - fy) * v[i, j] + fx * (1 - fy) * v[i + 1, j]
              + (1 - fx) * fy * v[i, j + 1] + fx * fy * v[i + 1, j + 1])
        return np.power(np.float32(10.0), lk)


class EOS:
    """the code's general-EOS table dump: log10(e/rho), log10(p/rho) on (log rho, log T)."""
    def __init__(self, fn):
        hdr = open(fn).read(2000).splitlines()
        nx, ny, x0, dx, y0, dy = hdr[2].lstrip('#').split()
        self.nx, self.ny = int(nx), int(ny)
        self.x0, self.dx, self.y0, self.dy = float(x0), float(dx), float(y0), float(dy)
        d = np.loadtxt(fn, comments='#', dtype=np.float64)
        self.le = d[:, 0].reshape(self.ny, self.nx).T.copy()      # [ix(rho), iy(T)]
        self.lp = d[:, 1].reshape(self.ny, self.nx).T.copy()
        self.lT = self.y0 + self.dy * np.arange(self.ny)

    def temp_pgas(self, rho, eint):
        """T and P_gas of cells (rho, e_gas [erg/cm3]); linear in log rho between table rows."""
        lr = np.log10(rho.astype(np.float64)).ravel()
        le = np.log10(np.maximum(eint.astype(np.float64), 1e-300)).ravel() - lr
        x = np.clip((lr - self.x0) / self.dx, 0, self.nx - 1.000001)
        i0 = x.astype(np.int64); f = x - i0
        T = np.empty_like(lr); P = np.empty_like(lr)
        for ii in np.unique(i0):
            s = np.where(i0 == ii)[0]
            lt0 = np.interp(le[s], self.le[ii], self.lT)
            lt1 = np.interp(le[s], self.le[min(ii + 1, self.nx - 1)], self.lT)
            lt = (1 - f[s]) * lt0 + f[s] * lt1
            T[s] = lt
            # p/rho at (rho, T): bilinear
            y = np.clip((lt - self.y0) / self.dy, 0, self.ny - 1.000001); j = y.astype(np.int64); g = y - j
            i1 = min(ii + 1, self.nx - 1)
            lp = ((1 - f[s]) * ((1 - g) * self.lp[ii, j] + g * self.lp[ii, j + 1])
                  + f[s] * ((1 - g) * self.lp[i1, j] + g * self.lp[i1, j + 1]))
            P[s] = lp
        T = (10 ** T).reshape(rho.shape).astype(np.float32)
        P = (10 ** P * 10 ** lr).reshape(rho.shape).astype(np.float32)
        return T, P


def assemble(d):
    geo = np.asarray(d['mb_geometry'])
    nx1, nx2, nx3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_mb'], d['nx2_mb'], d['nx3_mb']
    out = {}
    for v in d['var_names']:
        a = np.empty((nx3, nx2, nx1), np.float32)
        for m, blk in enumerate(d['mb_data'][v]):
            j0 = int(round((geo[m, 2] - d['x2min']) / (d['x2max'] - d['x2min']) * nx2))
            k0 = int(round((geo[m, 4] - d['x3min']) / (d['x3max'] - d['x3min']) * nx3))
            i0 = int(round((geo[m, 0] - d['x1min']) / (d['x1max'] - d['x1min']) * nx1)) if b1 < nx1 else 0
            a[k0:k0 + b3, j0:j0 + b2, i0:i0 + b1] = blk
        out[v] = a
    return out


def dump_time(path):
    with open(path, 'rb') as f:
        return float(re.search(r'time=\s*([-+0-9.eE]+)', f.read(400).decode('latin-1')).group(1))


_VMAP = None


def vet_file(t):
    global _VMAP
    if _VMAP is None:
        _VMAP = {}
        for f in glob.glob(RUN + '/bin/r136.m1_vet.*.bin'):
            try:
                _VMAP[round(dump_time(f))] = f
            except Exception:
                pass
    return _VMAP.get(round(t))


RED_VER = 3      # 2 (11-05): comoving F0 for Gamma / L_rad, radiation enthalpy in L_conv; 3 (11-06): solver-consistent F0 (r136_f0.py)


def dumps():
    out = []
    for h in sorted(glob.glob(RUN + '/bin/r136.hydro_w.*.bin')):
        n = int(h.split('.')[-2]); m = RUN + '/bin/r136.m1.%05d.bin' % n
        if os.path.exists(m) and time.time() - max(os.path.getmtime(h), os.path.getmtime(m)) > 180:
            out.append((n, h, m))
    return out


def x1faces():
    rows = [ln.split() for ln in open(RUN + '/r136.x1grid.txt') if ln.strip() and not ln.startswith('#')]
    i = np.array([int(x[0]) for x in rows]); x1f = np.array([float(x[1]) for x in rows])
    hdr = open(RUN + '/r136.x1grid.txt').readline()
    is_, ie = int(re.search(r'is = (\d+)', hdr).group(1)), int(re.search(r'ie = (\d+)', hdr).group(1))
    f = x1f[(i >= is_) & (i <= ie + 1)]
    return f


def ic_column():
    d = np.loadtxt(ICF, comments='#')
    return dict(r=d[:, 0], rho=d[:, 1], eint=d[:, 2], Fr=d[:, 3], E=d[:, 4], T=d[:, 5], fmlt=d[:, 6])


# ------------------------------------------------------------------ reduction
def reduce_dump(n, hf, mf, x1f, kap, eos):
    os.makedirs(CACHE, exist_ok=True)
    cf = os.path.join(CACHE, 'red_%05d.npz' % n)
    fv = vet_file(dump_time(hf))
    st = np.array([os.path.getsize(hf), os.path.getsize(mf), os.path.getsize(fv) if fv else 0], float)
    if os.path.exists(cf):
        c = np.load(cf)
        if 'ver' in c and int(c['ver']) >= RED_VER and np.array_equal(c['stamp'], st):
            return dict(c)
    dh = bc.read_binary(hf); t = dh['time']; w = assemble(dh); del dh
    dm = bc.read_binary(mf); assert abs(dm['time'] - t) < 1e-6 * max(1, t); m = assemble(dm); del dm
    nx3, nx2, nx1 = w['dens'].shape
    x2f = np.linspace(np.pi / 3, 2 * np.pi / 3, nx2 + 1); x3f = np.linspace(0, np.pi / 3, nx3 + 1)
    r = 0.5 * (x1f[1:] + x1f[:-1]); dr = np.diff(x1f); dvr = (x1f[1:] ** 3 - x1f[:-1] ** 3) / 3
    dom = np.diff(x3f)[:, None] * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :]
    om = dom.sum()
    sh = lambda a: np.einsum('kji,kj->i', a.astype(np.float64), dom) / om   # noqa: E731
    rho, vr = w['dens'], w['velx']
    T, Pg = eos.temp_pgas(rho, w['eint'])
    E, Flab = m['m1_e'], m['m1_f1']
    del m
    frr, vet = None, False
    if fv:
        try:
            dv = bc.read_binary(fv)
            if abs(float(dv['time']) - t) < 1e-6 * max(abs(t), 1.0):
                frr = assemble(dv)['m1_d11']; vet = True
            del dv
        except Exception as e:
            print('m1_vet unreadable', n, e, flush=True)
    hrad = E * ((1.0 + frr) if vet else np.float32(4. / 3.))     # E + P_rr
    del frr
    kr = kap(rho, T)
    Fr, _ = f0_solver(Flab - vr * hrad, E, kr, rho, x1f, r)        # comoving F0 (solver-consistent, r136_f0.py)
    dmass = np.einsum('kji,kj->i', rho.astype(np.float64), dom) * dvr          # wedge shell masses
    mcum = GM0 / G + FAC * (np.cumsum(dmass) - 0.5 * dmass)                    # m(r) at centres
    g = G * mcum / r ** 2
    res = dict(ver=RED_VER, stamp=st, time=t, n=n, r=r, vet=vet)
    res['rho'] = sh(rho); res['T'] = sh(T); res['Pg'] = sh(Pg); res['E'] = sh(E); res['Fr'] = sh(Fr)
    res['Fr_lab'] = sh(Flab)
    res['kap'] = sh(kr)
    res['gam_rad'] = sh(kr * Fr) / (C * g)
    res['gam_rad_lab'] = sh(kr * Flab) / (C * g)
    del Flab
    res['gam_L'] = res['kap'] * LSTAR / (4 * np.pi * C * G * mcum)
    res['prad_frac'] = sh((E / 3) / (E / 3 + Pg))
    vt = sh(rho * vr) / np.maximum(res['rho'], 1e-300)
    v2 = w['velx'] ** 2 + w['vely'] ** 2 + w['velz'] ** 2
    res['vr_rms'] = np.sqrt(np.maximum(sh(vr ** 2) - vt ** 2, 0))
    res['v_rms'] = np.sqrt(sh(v2))
    res['L_rad'] = 4 * np.pi * r ** 2 * res['Fr']
    res['L_rad_lab'] = 4 * np.pi * r ** 2 * res['Fr_lab']
    res['L_conv_gas'] = 4 * np.pi * r ** 2 * sh((w['eint'] + Pg) * (vr - vt[None, None, :].astype(np.float32)))
    res['L_radadv'] = 4 * np.pi * r ** 2 * (res['Fr_lab'] - res['Fr'])     # = <(E + P_rr) v_r> of the velocity in F; keeps L_rad + L_conv
    res['L_conv'] = res['L_conv_gas'] + res['L_radadv']
    del hrad
    res['mdot'] = 4 * np.pi * r ** 2 * sh(rho * vr)
    res['M_env'] = FAC * dmass.sum()
    res['KE'] = 0.5 * FAC * (np.einsum('kji,kj->i', (rho * v2).astype(np.float64), dom) * dvr).sum()
    dtau = kr * rho * dr[None, None, :].astype(np.float32)
    tau = np.flip(np.cumsum(np.flip(dtau.astype(np.float64), axis=2), axis=2), axis=2)
    res['tau'] = sh(tau)
    del tau, dtau
    # slices: equatorial = theta cell nearest pi/2 (j = nx2/2), meridional = phi cell nearest pi/6 (k = nx3/2)
    je, km = nx2 // 2, nx3 // 2
    gloc = kr * Fr / (C * g[None, None, :]).astype(np.float32)
    for nm, a in (('rho', rho), ('T', T), ('vr', vr), ('gam', gloc)):
        res['eq_' + nm] = a[:, je, :].copy()       # (k, i)
        res['me_' + nm] = a[km, :, :].copy()       # (j, i)
    res['x2f'], res['x3f'], res['x1f'] = x2f, x3f, x1f
    np.savez(cf, **res)
    print('reduced dump %05d t = %.3f d' % (n, t / DAY), flush=True)
    return res


def r_photo(r, tau):
    idx = np.where(tau >= 1.0)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[-1]
    if i + 1 >= len(r):
        return r[i]
    a, b = np.log(tau[i]), np.log(max(tau[i + 1], 1e-300))
    return r[i] + (r[i + 1] - r[i]) * (0 - a) / (b - a)


def read_hst():
    f = RUN + '/r136.user.hst'
    names = None
    for ln in open(f):
        if ln.startswith('#') and '[1]=' in ln:
            names = [x.split('=')[1] for x in re.findall(r'\[\d+\]=\S+', ln)]
    d = np.loadtxt(f, comments='#')
    # restarts rewind time and append: keep the last run's rows only (11-06: 5.0-5.33 d doubled)
    t_ = d[:, 0]; keep = t_ < np.append(np.minimum.accumulate(t_[::-1])[::-1][1:], np.inf)
    d = d[keep]
    return {k: d[:, i] for i, k in enumerate(names)}


# ------------------------------------------------------------------ figures
def fig_timeseries(h, red, out):
    t = h['time'] / DAY
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.2), constrained_layout=True)
    a = ax[0, 0]
    for k, lab in (('L_in', 'L$_{in}$ (inner BC)'), ('L_bot', 'L$_{bot}$'), ('L_top', 'L$_{top}$ (Marshak top)')):
        a.plot(t, FAC * h[k] / LSTAR, lw=1, label=lab)
    a.axhline(1, color='k', lw=0.6, ls=':')
    a.set_ylabel('L / L$_\\star$ (wedge $\\times$ 4$\\pi$/$\\Omega$)'); a.legend(fontsize=9); a.set_xlabel('t [d]')
    a = ax[0, 1]
    tr = np.array([x['time'] for x in red]) / DAY
    rph = np.array([r_photo(x['r'], x['tau']) for x in red]) / RSUN
    a.plot(tr, rph, 'o-', ms=3, color='C3'); a.set_ylabel('R$_{ph}$ [R$_\\odot$] (shell-mean $\\tau_R$ = 1)', color='C3')
    a.set_xlabel('t [d]')
    b = a.twinx()
    Lt = np.interp(tr * DAY, h['time'], FAC * h['L_top'])
    teff = (Lt / (4 * np.pi * SIG * (rph * RSUN) ** 2)) ** 0.25 / 1e3
    b.plot(tr, teff, 's--', ms=3, color='C0'); b.set_ylabel('T$_{eff}$ [kK] (L$_{top}$, R$_{ph}$)', color='C0')
    a = ax[1, 0]
    a.semilogy(t, FAC * np.maximum(h['KE_int'], 1e30), lw=1, label='KE')
    a.set_ylabel('kinetic energy [erg]'); a.set_xlabel('t [d]'); a.legend(fontsize=9)
    a = ax[1, 1]
    a.plot(t, FAC * (h['Etot'] - h['Etot'][0]), lw=1, color='C4')
    a.set_ylabel('E$_{tot}$(t) $-$ E$_{tot}$(0) [erg]'); a.set_xlabel('t [d]')
    fig.savefig(out, dpi=110); plt.close(fig)
    return rph, teff, tr


def pick_times(red):
    n = len(red)
    if n <= 4:
        return list(range(n))
    return sorted(set([0, n // 4, n // 2, (3 * n) // 4, n - 1]))


def fig_profiles(red, ic, h, out):
    sel = pick_times(red)
    cm = plt.cm.viridis(np.linspace(0, 0.9, len(sel)))
    r_ic = ic['r'] / RSUN
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.6), constrained_layout=True)
    a = ax.ravel()
    mic = (r_ic >= X1MIN) & (r_ic <= RMAX_PLOT)
    a[0].semilogy(r_ic[mic], ic['rho'][mic], 'k--', lw=1.2, label='IC column')
    a[1].semilogy(r_ic[mic], ic['T'][mic], 'k--', lw=1.2)
    a[3].plot(r_ic[mic], (ic['E'] / 3 / (ic['E'] / 3 + ic['Pg']))[mic], 'k--', lw=1.2)
    wm = np.interp([red[i]['time'] for i in sel], h['time'], h['w_mlt']) if 'w_mlt' in h else np.ones(len(sel))
    for c, i, wmi in zip(cm, sel, wm):
        x = red[i]; r = x['r'] / RSUN; s = r <= RMAX_PLOT; lab = 't = %.2f d' % (x['time'] / DAY)
        a[0].semilogy(r[s], x['rho'][s], color=c, lw=1, label=lab)
        a[1].semilogy(r[s], x['T'][s], color=c, lw=1)
        a[2].plot(r[s], x['gam_rad_lab'][s], color=c, lw=0.5, alpha=0.3)
        a[2].plot(r[s], x['gam_rad'][s], color=c, lw=1)
        a[2].plot(r[s], x['gam_L'][s], color=c, lw=1, ls=':')
        a[3].plot(r[s], x['prad_frac'][s], color=c, lw=1)
        a[4].semilogy(r[s], np.maximum(x['vr_rms'][s] / 1e5, 1e-4), color=c, lw=1)
        a[5].plot(r[s], x['L_rad_lab'][s] / LSTAR, color=c, lw=0.5, alpha=0.3)
        a[5].plot(r[s], x['L_rad'][s] / LSTAR, color=c, lw=1)
        a[5].plot(r[s], x['L_conv'][s] / LSTAR, color=c, lw=1, ls='--')
        fm = np.interp(x['r'][s], ic['r'], ic['fmlt'])
        a[5].plot(r[s], wmi * fm, color=c, lw=0.8, ls=':')
    a[2].axhline(1, color='k', lw=0.6, ls=':')
    a[0].set_ylabel('$\\rho$ [g cm$^{-3}$]'); a[0].legend(fontsize=8)
    a[1].set_ylabel('T [K]')
    a[2].set_ylabel('$\\Gamma$: solid $\\langle\\kappa_R F_0\\rangle$/(c g) comoving, faint lab, dotted $\\kappa_R$L/(4$\\pi$cGm)', fontsize=9)
    a[3].set_ylabel('P$_{rad}$ / (P$_{rad}$ + P$_{gas}$)'); a[3].set_ylim(0, 1)
    a[4].set_ylabel('v$_r$ rms [km/s]')
    a[5].set_ylabel('L$_{rad}$/L comoving (solid; faint lab), L$_{conv}$/L incl. (E+P$_{rr}$)v$_r$ (dashed), L$_{MLT}$ (dotted)', fontsize=8)
    a[2].set_ylim(0, max(2.0, min(6, np.nanmax([red[i]['gam_rad'][red[i]['r'] / RSUN <= RMAX_PLOT].max() for i in sel]) * 1.1)))
    for b in a:
        b.set_xlabel('r [R$_\\odot$]'); b.set_xlim(X1MIN, RMAX_PLOT)
    fig.savefig(out, dpi=105); plt.close(fig)


def fig_mdot(red, out):
    t = np.array([x['time'] for x in red]) / DAY
    r = red[0]['r'] / RSUN
    M = np.array([x['mdot'] / x['M_env'] * DAY for x in red])    # 1/d
    fig, ax = plt.subplots(figsize=(12, 5.2), constrained_layout=True)
    if len(t) >= 2:
        lim = np.nanpercentile(np.abs(M), 99.5) or 1e-6
        tf = np.concatenate([[t[0] - 0.5 * (t[1] - t[0])], 0.5 * (t[1:] + t[:-1]), [t[-1] + 0.5 * (t[-1] - t[-2])]])
        x1f = red[0]['x1f'] / RSUN
        pc = ax.pcolormesh(tf, x1f, M.T, cmap='RdBu_r', norm=SymLogNorm(linthresh=lim * 1e-3, vmin=-lim, vmax=lim), shading='flat')
        fig.colorbar(pc, ax=ax, label='$\\dot M_r$ / M$_{env}$ [d$^{-1}$]  (> 0 outward)')
    ax.set_ylim(X1MIN, max(RMAX_PLOT, 40)); ax.set_xlabel('t [d]'); ax.set_ylabel('r [R$_\\odot$]')
    ax.set_title('4$\\pi$r$^2\\langle\\rho v_r\\rangle$ over the envelope mass above x1min (M$_{env}$ = %.3g M$_\\odot$ at t = 0)'
                 % (red[0]['M_env'] / 1.989e33), fontsize=10)
    fig.savefig(out, dpi=110); plt.close(fig)


FIELDS = (('rho', '$\\rho$ [g cm$^{-3}$]', 'log', 'cividis'), ('T', 'T [K]', 'log', 'inferno'),
          ('vr', 'v$_r$ [km/s]', 'div', 'RdBu_r'), ('gam', '$\\Gamma_{loc}$ = $\\kappa_R$F$_0$/(c g), comoving', 'lin', 'magma'))


def ranges(x, full, fn):
    if os.path.exists(fn) and not full:
        return json.load(open(fn))
    rr = {}
    s = x['r'] / RSUN <= RMAX_PLOT
    for k, _, kind, _ in FIELDS:
        a = np.concatenate([x['eq_' + k][:, s].ravel(), x['me_' + k][:, s].ravel()])
        if k == 'vr':
            a = a / 1e5; v = float(np.nanpercentile(np.abs(a), 99.5)) or 1.0; rr[k] = [-v, v]
        elif kind == 'log':
            a = a[a > 0]; rr[k] = [float(np.nanpercentile(a, 0.5)), float(np.nanpercentile(a, 99.9))]
        else:
            rr[k] = [0.0, float(min(np.nanpercentile(a, 99.5), 8.0))]
    json.dump(rr, open(fn, 'w'))
    return rr


def frame(x, rr, out, which):
    r = x['x1f'] / RSUN; s = r <= RMAX_PLOT + 1e-9; rf = r[s]; ni = len(rf) - 1
    fig, ax = plt.subplots(1, 4, figsize=(17, 6.2) if which == 'eq' else (17, 6.6), constrained_layout=True)
    ang = x['x3f'] if which == 'eq' else x['x2f']
    A, Rr = np.meshgrid(ang, rf, indexing='ij')
    if which == 'eq':
        X, Y = Rr * np.cos(A), Rr * np.sin(A)
    else:
        X, Y = Rr * np.sin(A), Rr * np.cos(A)
    for a, (k, lab, kind, cmap) in zip(ax, FIELDS):
        d = x[which + '_' + k][:, :ni]
        if k == 'vr':
            d = d / 1e5; norm = TwoSlopeNorm(0, rr[k][0], rr[k][1])
        elif kind == 'log':
            norm = LogNorm(*rr[k])
        else:
            norm = Normalize(*rr[k])
        pc = a.pcolormesh(X, Y, d, cmap=cmap, norm=norm, shading='flat', rasterized=True)
        a.set_aspect('equal'); a.set_xlabel('[R$_\\odot$]')
        fig.colorbar(pc, ax=a, orientation='horizontal', label=lab, shrink=0.9, pad=0.06)
    ttl = ('equatorial plane ($\\theta$ = %.1f$^\\circ$), $\\phi$ 0-60$^\\circ$' % np.degrees(0.5 * (x['x2f'][len(x['x2f']) // 2 - 1] + x['x2f'][len(x['x2f']) // 2]))
           if which == 'eq' else 'meridional plane ($\\phi$ = %.1f$^\\circ$), $\\theta$ 60-120$^\\circ$' % np.degrees(0.5 * (x['x3f'][len(x['x3f']) // 2 - 1] + x['x3f'][len(x['x3f']) // 2])))
    fig.suptitle(LABEL + '  t = %.2f d   %s   r $\\leq$ %.0f R$_\\odot$' % (x['time'] / DAY, ttl, RMAX_PLOT), fontsize=12)
    fig.savefig(out, dpi=100); plt.close(fig)


# ------------------------------------------------------------------ status
def run_status(red, rph, teff, h):
    lg = open(RUN + '/run.log', 'rb').read().decode('latin-1', 'replace')
    el = re.findall(r'elapsed=\s*([0-9.eE+-]+)\s+cycle=(\d+)\s+time=([0-9.eE+-]+)\s+dt=([0-9.eE+-]+)', lg)
    nc = re.findall(r'NON-CONVERGED=\s*([0-9.eE+-]+)', lg)
    fatal = lg.count('FATAL')
    ll = open(os.path.dirname(RUN) + '/launch.log').read()
    m = re.findall(r'PROD start on (\S+) GPUs (\S+) \(job (\d+)\)', ll)
    node, gpus, job = m[-1] if m else ('?', '?', '?')
    ngpu = len(gpus.split(','))
    st = {}
    if el:
        e, cyc, tt, dt = el[-1]; e = float(e); cyc = int(cyc); tt = float(tt); dt = float(dt)
        k = max(0, len(el) - 20); e0, c0, t0 = float(el[k][0]), int(el[k][1]), float(el[k][2])
        spc = (e - e0) / max(cyc - c0, 1); rate = (tt - t0) / max(e - e0, 1e-9)   # model s per wall s (recent)
        tl = TLIM
        st['s-time'] = '%.2f d' % (tt / DAY); st['s-cyc'] = '%d' % cyc; st['s-dt'] = '%.2f s' % dt
        st['c-wall'] = '%.1f h' % (e / 3600); st['c-gpuh'] = '%.0f' % (e / 3600 * ngpu)
        st['c-spc'] = '%.2f s' % spc
        st['c-proj'] = '%.0f h' % (e / 3600 + max(tl - tt, 0) / max(rate, 1e-9) / 3600)
        st['c-perday'] = '%.2f h' % (DAY / max(rate, 1e-9) / 3600)
    st['s-hw'] = '%d H200 on %s (alloc %s, GPUs %s)' % (ngpu, node, job, gpus)
    st['s-health'] = 'FATAL %d, NON-CONVERGED (last summary) %s' % (fatal, nc[-1] if nc else 'n/a')
    if len(rph):
        st['s-rph'] = '%.2f R<sub>&#9737;</sub>' % rph[-1]; st['s-teff'] = '%.1f kK' % teff[-1]
    st['s-updated'] = time.strftime('%Y-%m-%d %H:%M %Z')
    st['s-ndump'] = '%d' % len(red)
    st['__v'] = str(int(time.time()))
    banner = (os.environ.get('R136_BANNERHEAD', 'Phase 1, first look at 128&times;128 lateral.') + ' Start-up transient: the IC is MESA&rsquo;s MLT++ structure '
              '(super-Eddington layers held by MESA&rsquo;s enforced gradient) re-solved with our EOS and opacities; the frozen '
              'MLT scaffold ramps off. <b>No line force: any outflow, &#7742; or v<sub>&infin;</sub> here is NOT physical.</b>')
    st['s-banner'] = banner
    return st


def main():
    global RMAX_PLOT
    page = sys.argv[1]; full = '--full' in sys.argv
    os.makedirs(page + '/img', exist_ok=True)
    global X1MIN
    x1f = x1faces(); X1MIN = x1f[0] / RSUN; kap = Table(OPAC); eos = EOS(EOSF); ic = ic_column(); h = read_hst()
    ic['Pg'] = eos.temp_pgas(ic['rho'], ic['eint'])[1].astype(np.float64)
    red = [reduce_dump(n, hf, mf, x1f, kap, eos) for n, hf, mf in dumps()]
    red.sort(key=lambda x: float(x['time']))
    rph, teff, tr = fig_timeseries(h, red, page + '/img/timeseries.png')
    rf = os.path.join(CACHE, 'ranges.json')
    old = json.load(open(rf)).get('rmax', 0) if os.path.exists(rf) else 0
    want = float(os.environ.get('R136_RMAX', 0)) or max(32.0, 5 * np.ceil(1.25 * np.nanmax(rph) / 5))
    if want > old or full:            # the envelope grew past the frame: new limit, re-chosen ranges, all frames redrawn
        full = True; RMAX_PLOT = want
    else:
        RMAX_PLOT = old
    rr = ranges(red[-1], full, rf)
    rr['rmax'] = RMAX_PLOT; json.dump(rr, open(rf, 'w'))
    fig_profiles(red, ic, h, page + '/img/profiles.png')
    fig_mdot(red, page + '/img/mdot.png')
    tlist = []
    for which in ('eq', 'me'):
        os.makedirs(page + '/movie/' + which, exist_ok=True)
    for i, x in enumerate(red):
        tlist.append(round(float(x['time']) / DAY, 3))
        for which in ('eq', 'me'):
            fn = page + '/movie/%s/%s_%05d.png' % (which, which, i)
            if full or not os.path.exists(fn):
                frame(x, rr, fn, which)
    json.dump({'n': len(red), 't': tlist}, open(page + '/movie/n.json', 'w'))
    st = run_status(red, rph, teff, h)
    st['c-rmax'] = '%.0f' % RMAX_PLOT
    json.dump(st, open(page + '/status.json', 'w'), indent=1)
    print('OK', st.get('s-time'), 'R_ph', st.get('s-rph'), 'dumps', len(red))


if __name__ == '__main__':
    main()
