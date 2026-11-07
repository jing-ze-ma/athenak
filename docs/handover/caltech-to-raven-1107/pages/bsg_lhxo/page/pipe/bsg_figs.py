"""REPRO copy (scripts_repro, 10-06): Figs 2, 3, 4, 6 of Ma, Bildsten & Jiang 2026
(arXiv 2609.28656) from the REPRODUCTION run $BSG_REPRO_RUN (default
/raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n, seeded, force_reference_work = split; 10-06)
(the paper's IC as is, 20-80 Rsun, nx1 267 stretched, Ma's padded TOPS tables, vet_gd;
READ ONLY).  Outputs, cache and movie ranges live under figs/scripts_repro/.

usage: python3 bsg_figs.py [--arm arm1|arm2] [--snap N] [--run DIR] [--rmax-plot R]

Per-dump reductions are cached in figs/cache/<arm>/red_NNNNN.npz (keyed by the dump's
size and
mtime), so re-running only reads new dumps. Definitions follow the paper's Appendix A:
  energies  = sums over the simulated wedge (NOT scaled to 4 pi)            (A1-A5)
  Mdot_r    = 4 pi/Omega sum rho v_r r^2 dOmega, / M_env, M_env = 4 pi/Omega * M_wedge
  rho, T    = solid-angle (dOmega) weighted shell means                     (A13-A14)
  tau_ijk   = int_r^rmax kappa_R rho dr (column-wise, cell-centre values),
  r_photo3D = r(tau_bar = 1), tau_bar = dOmega-weighted mean of tau_ijk     (A9-A11)
"""
import argparse
import glob
import json
import os
import re
import time
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert_repro as bc  # noqa: E402  (repro: '=' inside header comments)

import matplotlib  # noqa: E402
import matplotlib.ticker  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# handover copy (10-09): site paths from env, defaults = viper
BSG = os.environ.get('BSG_DATA_DIR', '/viper/ptmp2/jinma/bsg_1001')
FIGS = os.environ.get('BSG_WORK', os.path.dirname(os.path.abspath(__file__)))   # all writes here
CACHE = os.path.join(FIGS, 'cache')   # set per arm in main(): cache/<arm>
RSUN = 6.957e10
DAY = 86400.0
YR = 3.15576e7
MU_MU_OVER_K = 0.62 * 1.66053906660e-24 / 1.380649e-16   # T[K] = (P/rho) * this
GAMMA = 5.0 / 3.0
GM = 2.6542488e27
# Linhao Ma's tables (TOPS X 0.7 Z 0.008, padded with edge values), md5 450bc0c1 = run's
OPAC = os.environ.get('BSG_OPAC', os.path.join(
    BSG, 'ma2026_files/athenak_tables/rosseland_ma2026_x0.7_z0.008_padded.txt'))
L_STAR = 5.581461e38     # the IC file's L (implicit_flux_x1min = L/4 pi r_in^2)
KAPPA_EDD = 4 * np.pi * 2.99792458e10 * GM / L_STAR   # cm^2/g
# per-arm run dir, input name and 1-D IC (Fig 3 right panel)
ARMS = {
    'arm1': dict(run=os.path.join(BSG, 'step4/prod'), inp='bsg3d_prod.athinput',
                 ic=os.path.join(BSG, 'step3/ic_bsg/ic_bsg_mlt.txt')),
    'arm2': dict(run=os.path.join(BSG, 'arm2/prod_1002'), inp='bsg3d_arm2.athinput',
                 ic=os.path.join(BSG, 'arm2/ic3d/ic_bsg_arm2.txt')),
    'arm2s2': dict(run='/raven/ptmp/jinma/bsg_raven_1002/arm2r2b_4n',
                   inp='bsg3d_arm2.athinput',
                   ic=os.path.join(BSG, 'arm2/ic3d/ic_bsg_arm2.txt')),
    # reproduction run: the paper's IC as is (= repro_4n/bundle/ic_ma2026.txt, md5 21672b25)
    # (10-06) 'repro' = the run shown on the page: env BSG_REPRO_RUN (default = the seeded
    # split run repro_sf_4n); 'repro_full' = the unseeded force_reference_work = full run
    # repro_4n (stopped), kept for the KE comparison.  Same grid, IC and tables.
    'repro': dict(run=os.environ.get('BSG_REPRO_RUN',
                                     '/raven/ptmp/jinma/bsg_raven_1002/repro_sf_4n'),
                  inp=os.environ.get('BSG_INP', 'bsg3d_repro.athinput'),
                  ic=os.environ.get('BSG_IC', os.path.join(BSG, 'repro/ic3d/ic_ma2026.txt'))),
    'repro_full': dict(run=os.environ.get('BSG_REPRO_FULL_RUN',
                                          '/raven/ptmp/jinma/bsg_raven_1002/repro_4n'),
                       inp=os.environ.get('BSG_INP', 'bsg3d_repro.athinput'),
                       ic=os.environ.get('BSG_IC', os.path.join(BSG, 'repro/ic3d/ic_ma2026.txt'))),
}


def cache_dir(run):
    """per-run cache: cache/<basename of the run dir> (dump numbers repeat across runs)"""
    return os.path.join(FIGS, 'cache', os.path.basename(os.path.normpath(run)))
ARM = 'repro'
IC = ARMS[ARM]['ic']
INP = ARMS[ARM]['inp']
GRID_CANDIDATES = ['{run}/bsg3d.x1grid.txt']   # repro: the run's own grid dump only
STRETCH_KEYS = ['f_stretch_r_p_amp', 'f_stretch_r_p_xa', 'f_stretch_r_p_xb',
                'f_stretch_r_p_w', 'nx1', 'x1min', 'x1max']
# Fig 4 radius [Rsun]: the paper's 36 (its Fe peak 35 + 1), USER 11-06.  The run's own Fe
# kappa_R peak (shell mean) is measured per dump (fe_peak in numbers.json) and quoted in the caption.
R_SLICE = 36.0


# ---------------------------------------------------------------------------------------
# grid, opacity
def input_keys(path):
    out = {}
    blk = None
    for line in open(path):
        s = line.split('#')[0].strip()
        if s.startswith('<'):
            blk = s
        elif '=' in s and blk == '<mesh>':
            k, v = [x.strip() for x in s.split('=', 1)]
            out[k] = v
    return out


def radial_faces(run):
    """True x1 faces of the stretched grid from the code's grid dump; the dump's run must
    have the same stretch keys as this run (checked)."""
    want = input_keys(os.path.join(run, INP))
    for pat in GRID_CANDIDATES:
        f = pat.format(run=run)
        if not os.path.exists(f):
            continue
        d = os.path.dirname(f)
        inps = glob.glob(os.path.join(d, '*.athinput')) + glob.glob(
            os.path.join(os.path.dirname(d), '*.athinput'))
        ok = d == run
        for inp in inps:
            k = input_keys(inp)
            if all(float(k.get(x, 'nan')) == float(want[x]) for x in STRETCH_KEYS):
                ok = True
                break
        if not ok:
            continue
        hdr = open(f).readline()
        m = re.search(r'is = (\d+) ie = (\d+)', hdr)
        is_, ie = int(m.group(1)), int(m.group(2))
        a = np.loadtxt(f, usecols=(0, 1))
        x1f = a[is_:ie + 2, 1]
        assert abs(x1f[-1] / float(want['x1max']) - 1) < 1e-12
        assert len(x1f) == int(want['nx1']) + 1, (f, len(x1f))
        assert abs(x1f[0] / float(want['x1min']) - 1) < 1e-12
        return x1f, f
    sys.exit('no x1 grid dump with matching stretch keys found')


class Kappa:
    """Bilinear log10 kappa_R(log T, log rho) on the TOPS table (edge-filled, as the
    code)."""
    def __init__(self, path=OPAC):
        hdr = [ln for ln in open(path) if ln.startswith('#')]
        g = [ln for ln in hdr if re.match(r'#\s*\d', ln)][0]
        nT, nD, lT0, dlT, lD0, dlD = g[1:].split()
        self.nT, self.nD = int(nT), int(nD)
        self.lT0, self.dlT, self.lD0, self.dlD = map(float, (lT0, dlT, lD0, dlD))
        self.tab = np.loadtxt(path, comments='#').reshape(self.nT, self.nD)

    def __call__(self, rho, T):
        x = (np.log10(np.maximum(T, 1.0)) - self.lT0) / self.dlT
        y = (np.log10(np.maximum(rho, 1e-30)) - self.lD0) / self.dlD
        x = np.clip(x, 0, self.nT - 1 - 1e-9)
        y = np.clip(y, 0, self.nD - 1 - 1e-9)
        i = x.astype(np.int32)
        j = y.astype(np.int32)
        fx = (x - i).astype(np.float32)
        fy = (y - j).astype(np.float32)
        t = self.tab.astype(np.float32)
        v = ((1 - fx) * (1 - fy) * t[i, j] + fx * (1 - fy) * t[i + 1, j]
             + (1 - fx) * fy * t[i, j + 1] + fx * fy * t[i + 1, j + 1])
        return np.power(np.float32(10.0), v)


# ---------------------------------------------------------------------------------------
# dumps
def assemble(d):
    """dict var -> (nx3, nx2, nx1) float32 array from read_binary output."""
    geo = np.asarray(d['mb_geometry'])
    nx1, nx2, nx3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_mb'], d['nx2_mb'], d['nx3_mb']
    out = {}
    for v in d['var_names']:
        a = np.empty((nx3, nx2, nx1), np.float32)
        for m, blk in enumerate(d['mb_data'][v]):
            i0 = int(round((geo[m, 0] - d['x1min']) / (d['x1max'] - d['x1min']) * nx1))
            j0 = int(round((geo[m, 2] - d['x2min']) / (d['x2max'] - d['x2min']) * nx2))
            k0 = int(round((geo[m, 4] - d['x3min']) / (d['x3max'] - d['x3min']) * nx3))
            a[k0:k0 + b3, j0:j0 + b2, i0:i0 + b1] = blk
        out[v] = a
    return out


def dump_time(path):
    """time= of a bin dump from its text header (no data read)."""
    with open(path, 'rb') as f:
        hdr = f.read(400).decode('latin-1')
    return float(re.search(r'time=\s*([-+0-9.eE]+)', hdr).group(1))


def dump_list(run):
    """(n, hydro_w, m1, m1_vet or None).  m1_vet dumps (variable = m1_vet: the closure's
    f_K and D_11..D_33, added by a restart overlay) are numbered from their own first
    dump, so they are matched to the m1 dumps by time."""
    hs = sorted(glob.glob(os.path.join(run, 'bin/bsg3d.hydro_w.*.bin')))
    vets = {}
    for v in sorted(glob.glob(os.path.join(run, 'bin/bsg3d.m1_vet.*.bin'))):
        vets[round(dump_time(v), 3)] = v
    out = []
    for h in hs:
        n = int(h.split('.')[-2])
        m = os.path.join(run, 'bin/bsg3d.m1.%05d.bin' % n)
        # Caltech 10-09: live run: skip a dump still being written (modified < 3 min ago)
        if os.path.exists(m) and time.time() - max(os.path.getmtime(h), os.path.getmtime(m)) > 180:
            out.append((n, h, m, vets.get(round(dump_time(m), 3))))
    return out


def stamp(*files):
    return np.array([os.path.getsize(f) + 1e-9 * os.path.getmtime(f) for f in files])


def reduce_dump(n, hfile, mfile, x1f, kap, vfile=None):
    """Shell profiles, totals and the two slices of one dump (cached)."""
    os.makedirs(CACHE, exist_ok=True)
    cf = os.path.join(CACHE, 'red_%05d.npz' % n)
    st = stamp(hfile, mfile) if vfile is None else stamp(hfile, mfile, vfile)
    if os.path.exists(cf):
        c = np.load(cf)
        if (c['stamp'].shape == st.shape and np.allclose(c['stamp'], st, rtol=0, atol=0)
                and float(c.get('r_slice', 36.0)) == R_SLICE):
            return dict(c)
    dh = bc.read_binary(hfile)
    t = dh['time']
    w = assemble(dh)
    del dh
    nx3, nx2, nx1 = w['dens'].shape
    x2f = np.linspace(np.pi / 3, 2 * np.pi / 3, nx2 + 1)
    x3f = np.linspace(0.0, np.pi / 3, nx3 + 1)
    r = 0.5 * (x1f[1:] + x1f[:-1])
    dr = np.diff(x1f)
    dvr = (x1f[1:] ** 3 - x1f[:-1] ** 3) / 3.0
    dmu = np.cos(x2f[:-1]) - np.cos(x2f[1:])
    dphi = np.diff(x3f)
    dom = dphi[:, None] * dmu[None, :]              # (k, j)
    omega = dom.sum()
    rho = w['dens']
    pg = (GAMMA - 1.0) * w['eint']
    T = pg / rho * np.float32(MU_MU_OVER_K)
    v2 = w['velx'] ** 2 + w['vely'] ** 2 + w['velz'] ** 2
    dm = np.einsum('kji,kj->i', rho.astype(np.float64), dom) * dvr      # shell masses
    res = dict(stamp=st, time=t, n=n, omega=omega, r=r, r_slice=R_SLICE)
    res['rho_sh'] = np.einsum('kji,kj->i', rho.astype(np.float64), dom) / omega
    res['T_sh'] = np.einsum('kji,kj->i', T.astype(np.float64), dom) / omega
    res['mdot_r'] = (4 * np.pi / omega) * r ** 2 * np.einsum(
        'kji,kj->i', (rho * w['velx']).astype(np.float64), dom)
    res['M_wedge'] = dm.sum()
    res['E_k'] = 0.5 * (np.einsum('kji,kj->i', (rho * v2).astype(np.float64), dom)
                        * dvr).sum()
    res['U_gas'] = (np.einsum('kji,kj->i', w['eint'].astype(np.float64), dom) * dvr).sum()
    res['E_grav'] = -(GM * dm / r).sum()
    res['vrms_sh'] = np.sqrt(np.einsum('kji,kj->i', v2.astype(np.float64), dom) / omega)
    kr = kap(rho, T)
    # tau from the outer edge inward, cell-centre sampling (A9):
    # tau_i = sum_{i'>=i} k rho dr
    dtau = kr * rho * dr[None, None, :].astype(np.float32)
    tau = np.flip(np.cumsum(np.flip(dtau.astype(np.float64), axis=2), axis=2), axis=2)
    res['tau_mean'] = np.einsum('kji,kj->i', tau, dom) / omega
    res['kap_sh'] = np.einsum('kji,kj->i', kr.astype(np.float64), dom) / omega
    del tau, dtau
    # slices: phi = 0 plane (k = 0, first cell) and the shell nearest R_SLICE
    i36 = int(np.argmin(np.abs(r / RSUN - R_SLICE)))
    res['i36'] = i36
    for nm, a in (('rho', rho), ('T', T), ('vr', w['velx']), ('v2', v2), ('pg', pg),
                  ('kap', kr)):
        res['p0_' + nm] = a[0].copy()
        res['s36_' + nm] = a[:, :, i36].copy()
    del w, rho, T, v2, pg, kr
    dm1 = bc.read_binary(mfile)
    assert abs(dm1['time'] - t) < 1e-6 * max(1.0, t), (hfile, mfile)
    e = assemble(dm1)['m1_e']
    del dm1
    res['U_rad'] = (np.einsum('kji,kj->i', e.astype(np.float64), dom) * dvr).sum()
    res['p0_erad'] = e[0].copy()
    del e
    if vfile is not None:
        # the closure's Eddington tensor diagonal: P_rad,aa = D_aa E exactly
        dv = bc.read_binary(vfile)
        assert abs(dv['time'] - t) < 1e-6 * max(1.0, t), (mfile, vfile)
        v = assemble(dv)
        del dv
        for nm in ('m1_fk', 'm1_d11', 'm1_d22', 'm1_d33'):
            res['p0_' + nm[3:]] = v[nm][0].copy()
        res['fk_sh'] = np.einsum('kji,kj->i', v['m1_fk'].astype(np.float64), dom) / omega
        del v
    res['x2f'], res['x3f'], res['x1f'] = x2f, x3f, x1f
    np.savez(cf, **res)
    print('reduced dump %05d  t = %.3f d' % (n, t / DAY), flush=True)
    return res


def r_photo(r, taum):
    """r where the mean tau crosses 1 (searching inward from the top)."""
    idx = np.where(taum >= 1.0)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[-1]
    if i + 1 >= len(r):
        return r[i]
    lt0, lt1 = np.log(taum[i]), np.log(max(taum[i + 1], 1e-300))
    return r[i] + (r[i + 1] - r[i]) * (0 - lt0) / (lt1 - lt0)


def read_hst(path):
    hdr = open(path).readlines()[1]
    names = re.findall(r'\[\d+\]=(\S+)', hdr)
    # Caltech 10-09: the run is live; drop a partially written last row
    rows = [ln.split() for ln in open(path) if not ln.startswith('#')]
    rows = [r for r in rows if len(r) == len(names)]
    a = np.array(rows, dtype=float).reshape(-1, len(names))
    return {nm: a[:, i] for i, nm in enumerate(names)}


# ---------------------------------------------------------------------------------------
# figures: USER 11-06, every figure mimics Ma+2026 (panels, pixel layout, ranges, colormaps, labels).
# The axes rectangles below are the pixel boxes measured on the page's paper_figN.png, so ours and
# his have the same size at 100 dpi.  All limits are fixed (read off his figures); values beyond saturate.
matplotlib.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix',
                            'font.size': 15, 'axes.labelsize': 15, 'axes.titlesize': 16,
                            'xtick.labelsize': 15, 'ytick.labelsize': 15, 'legend.fontsize': 15,
                            'axes.linewidth': 1.2, 'legend.handlelength': 2.4})


def fs(n, title=None, leg=None, pad=3.5):
    """rc for one figure: font sizes [pt at 100 dpi] and tick-label pad, tuned against his raster"""
    return {'font.size': n, 'axes.labelsize': n, 'axes.titlesize': title or n, 'xtick.labelsize': n,
            'ytick.labelsize': n, 'legend.fontsize': leg or n, 'xtick.major.pad': pad}


def px_axes(fig, x0, x1, ytop, ybot, **kw):
    """axes from pixel box (origin top left), figure of fig.get_size_inches() x 100 px"""
    W, H = fig.get_size_inches() * 100
    return fig.add_axes([x0 / W, 1 - ybot / H, (x1 - x0) / W, (ybot - ytop) / H], **kw)


def time_stamp(fig, t_d, color='k'):
    """movie frames only: the dump's time, top left (his figures have no time label)"""
    fig.text(0.008, 0.975, 't = %.2f d' % t_d, ha='left', va='top', fontsize=16, color=color)


F2_T1 = 57.4          # his Fig 2 time axis
F2_RY = (20.0, 63.0)  # his Fig 2 radial axis [Rsun]
F2_EY = (-9.0, 0.97)  # his energy panel (read off: 1e44 erg)
F2_KY = (5.9e-5, 0.93)   # his kinetic-energy panel (read off)


def fig2(reds, hst, out, rmax_plot=None):
    with matplotlib.rc_context({**fs(17, leg=16, pad=8), 'axes.labelpad': 9}):
        return _fig2(reds, hst, out)


def _fig2(reds, hst, out):
    t = np.array([c['time'] for c in reds]) / DAY
    rph = np.array([r_photo(c['r'], c['tau_mean']) for c in reds]) / RSUN
    eg = np.array([c['E_grav'] for c in reds])
    ui = np.array([c['U_gas'] + c['U_rad'] for c in reds])
    ek = np.array([c['E_k'] for c in reds])
    menv = (4 * np.pi / reds[0]['omega']) * np.array([c['M_wedge'] for c in reds])
    tmax = max(t.max(), hst['time'].max() / DAY)
    fig = plt.figure(figsize=(14.44, 14.21), dpi=100)
    boxes = ((27.5, 243.0), (296.5, 513.5), (567.5, 783.5), (837.5, 1053.5), (1108.0, 1324.5))
    axs = [px_axes(fig, 126.5, 1271.5, a, b) for a, b in boxes]
    ax = axs[0]
    ax.plot(t, (eg - eg[0]) / 1e44, '-', color='green', lw=2.5,
            label=r'change of gravitational binding energy $\Delta E_{\rm grav}$')
    ax.plot(t, (ui - ui[0]) / 1e44, '-', color='orange', lw=2.5,
            label=r'change of internal energy $\Delta U_{\rm rad}+\Delta U_{\rm gas}$')
    ax.set_ylabel(r'Energy [$10^{44}$ erg]')
    ax.set_ylim(*F2_EY)
    ax.set_yticks([0, -2.5, -5, -7.5])
    ax.set_yticklabels(['0.0', r'$-2.5$', r'$-5.0$', r'$-7.5$'])
    ax.legend(loc='upper right', framealpha=0.8)
    ax = axs[1]
    ax.semilogy(t[t > 0], ek[t > 0] / 1e44, '-', color='blue', lw=2.5,
                label=r'gas kinetic energy $E_{\rm k}$')
    ax.set_ylabel(r'Energy [$10^{44}$ erg]')
    ax.set_ylim(*F2_KY)
    ax.set_yticks([1e-1, 1e-2, 1e-3, 1e-4])
    ax.legend(loc='lower right', framealpha=0.8)
    for ax in axs[:2]:
        ax.tick_params(axis='y', which='both', right=True, labelright=True)
    # colour maps: cell edges in time = midpoints between dumps
    if len(t) > 1:
        te = np.concatenate([[t[0] - 0.5 * (t[1] - t[0])], 0.5 * (t[1:] + t[:-1]),
                             [t[-1] + 0.5 * (t[-1] - t[-2])]])
    else:
        te = np.array([t[0] - 0.25, t[0] + 0.25])
    te = np.clip(te, 0, None)
    rf = reds[0]['x1f'] / RSUN
    mdot = np.array([c['mdot_r'] for c in reds]) / menv[:, None] * YR
    rho = np.log10(np.array([c['rho_sh'] for c in reds]))
    T = np.log10(np.array([c['T_sh'] for c in reds]))
    specs = [(axs[2], mdot, 'RdBu', -5, 5, [-4, -2, 0, 2, 4], (577, 775),
              r'$\dot{M}_r/M_{\rm envelope}$ [yr$^{-1}$]',
              r'radial mass flow $\dot{M}_r/M_{\rm envelope}$'),
             (axs[3], rho, 'Blues', -10.5, -3.5, [-10, -8, -6, -4], (847, 1045),
              r'$\log(\rho)$ [g cm$^{-3}$]', 'angular averaged gas density'),
             (axs[4], T, 'Blues', 3.8, 6.2, [4.0, 4.5, 5.0, 5.5, 6.0], (1117, 1315),
              r'$\log(T)$ [K]', 'angular averaged gas temperature')]
    for ax, z, cm, lo, hi, cbt, (cy0, cy1), cl, tl in specs:
        pc = ax.pcolormesh(te, rf, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat')
        cax = px_axes(fig, 1294.5, 1323.5, cy0 - 9.5, cy1 + 9.5)
        cb = fig.colorbar(pc, cax=cax, extend='both')
        cb.set_ticks(cbt)
        cb.set_label(cl)
        ax.plot(t, rph, '--', color='blue', lw=3, label=r'$r_{\rm photo,\,3D}$')
        ax.set_ylim(*F2_RY)
        ax.set_yticks([20, 40, 60])
        ax.set_ylabel(r'Radial Coordinate [$R_\odot$]')
        ax.text(0.02, 0.955, tl, transform=ax.transAxes, va='top')
        ax.text(43.7, 0.905, 'steady state', transform=ax.get_xaxis_transform(), color='purple',
                ha='center', va='center')
    axs[4].legend(loc='lower left', framealpha=0.8)
    for ax in axs:
        ax.set_xlim(0, F2_T1)
        ax.set_xticks([0, 10, 20, 30, 40, 50, 57.4])
        ax.set_xticklabels(['0', '10', '20', '30', '40', '50', '57.4'])
        ax.axvline(30.0, color='purple', ls='--', lw=2.5)
    axs[4].set_xlabel('Time since Simulation Starts [days]')
    fig.savefig(out, dpi=100)
    plt.close(fig)
    return dict(t_dumps_d=[round(float(x), 4) for x in t],
                t_hst_max_d=round(float(tmax), 4),
                r_photo_3D_Rsun=[round(float(x), 3) for x in rph],
                M_envelope_whole_star_Msun=float(menv[-1] / 1.98847e33),
                dE_grav_1e44=float((eg[-1] - eg[0]) / 1e44),
                dU_int_1e44=float((ui[-1] - ui[0]) / 1e44), E_k_1e44=float(ek[-1] / 1e44),
                dE_grav_min_1e44=float((eg - eg[0]).min() / 1e44),
                dE_grav_max_1e44=float((eg - eg[0]).max() / 1e44),
                dU_min_1e44=float((ui - ui[0]).min() / 1e44),
                dU_max_1e44=float((ui - ui[0]).max() / 1e44),
                E_k_max_1e44=float(ek.max() / 1e44), rph_max=float(np.nanmax(rph)),
                mdot_max=float(np.abs(mdot[:, rf[:-1] < 40]).max()))


# Fig 3: wedge geometry measured on paper_fig3.png (circle fit + the arcs at 25, 35, 45 Rsun):
# 7.2 px/Rsun, centre y 509 px, panel centres x 252, 718, 1184 px; wedge r = 21-52 Rsun.
F3_RMIN, F3_RMAX, F3_S, F3_CY = 21.0, 52.0, 7.2, 509.0
F3_CX = (252.0, 718.0, 1184.0)
F3_W, F3_H = 1444, 561
RAM_CMAP = LinearSegmentedColormap.from_list('ramp', ['white', 'blue', (0, 0, 77 / 255)])
F3_RAM = (0.0, 0.5)
F3_KAP = (0.005, 2.4)   # his colour bar: lower end ~0 (no tick there), upper end 2.4


def wedge_axes(fig, cx):
    """axes in Rsun whose origin is the wedge centre, at his pixel scale; his frame, arcs and labels"""
    x0, x1, yt, yb = cx - 232, cx + 232, 90, 420
    ax = px_axes(fig, x0, x1, yt, yb)
    ax.set_xlim(-232 / F3_S, 232 / F3_S)
    ax.set_ylim((F3_CY - yb) / F3_S, (F3_CY - yt) / F3_S)
    ax.axis('off')
    return ax


def wedge_lines(ax, cx):
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)
    for rr in (25, 35, 45):
        ax.plot(rr * np.cos(th), rr * np.sin(th), '-.', color='0.7', lw=1.2)
    for deg in (75, 90, 105):
        a = np.radians(deg)
        ax.plot([F3_RMIN * np.cos(a), F3_RMAX * np.cos(a)],
                [F3_RMIN * np.sin(a), F3_RMAX * np.sin(a)], '-.', color='0.7', lw=1.2)
    for rr in (F3_RMIN, F3_RMAX):
        ax.plot(rr * np.cos(th), rr * np.sin(th), 'k-', lw=1.5)
    for a in (np.pi / 3, 2 * np.pi / 3):
        ax.plot([F3_RMIN * np.cos(a), F3_RMAX * np.cos(a)],
                [F3_RMIN * np.sin(a), F3_RMAX * np.sin(a)], 'k-', lw=1.5)
    for deg in (60, 75, 90, 105, 120):
        a = np.radians(deg)
        ax.text(56.4 * np.cos(a), 56.4 * np.sin(a), '%d\u00b0' % deg, ha='center', va='center',
                fontsize=18, clip_on=False)
    # r labels at his pixel positions
    for rr, xp, yp in ((45, 453, 247), (35, 420, 309), (25, 384, 372)):
        ax.text((xp - 252) / F3_S, (F3_CY - yp) / F3_S, r'$%d\,R_\odot$' % rr, ha='center',
                va='center', fontsize=18, clip_on=False)


def prad_rr(red, est=None):
    """P_rad,rr on the phi = 0 plane and its provenance: D_11 E from an m1_vet dump; else
    f_K E from an offline estimate (prad/est_NNNNN.npz of the same dump, the Python port
    of the vet_col formal solution); else E/3."""
    e = red['p0_erad']
    if 'p0_d11' in red:
        return red['p0_d11'] * e, 'm1_vet'
    if est is not None and os.path.exists(est):
        z = np.load(est)
        if abs(float(z['time']) - float(red['time'])) < 1e-6 * float(red['time']):
            return z['fk_plane'].astype(np.float32) * e, 'offline vet_col port'
    return e / 3.0, 'E/3'


_IC_CACHE = {}


def ic_profile(kap):
    if IC not in _IC_CACHE:
        ic = np.loadtxt(IC)
        ric, rhoic, Tic, fmlt = ic[:, 0], ic[:, 1], ic[:, 5], ic[:, 6]
        _IC_CACHE[IC] = (ric, rhoic, Tic, fmlt, kap(rhoic, Tic))
    return _IC_CACHE[IC]


def fig3(red, rmax_plot, kap, out, est=None, stamp=False):
    with matplotlib.rc_context(fs(16.5, pad=4)):
        return _fig3(red, kap, out, est, stamp)


def _fig3(red, kap, out, est=None, stamp=False):
    rf = red['x1f'] / RSUN
    x2f = red['x2f']
    r_c = 0.5 * (rf[1:] + rf[:-1])
    i0 = int(np.searchsorted(rf, F3_RMIN) - 1)
    i0 = max(i0, 0)
    ni = int((rf <= F3_RMAX + 1e-9).sum())      # one past the last face below the cut
    fi = np.concatenate([[F3_RMIN], rf[i0 + 1:ni], [F3_RMAX]])
    Ra, TH = np.meshgrid(fi, x2f)
    X, Y = Ra * np.cos(TH), Ra * np.sin(TH)
    rho, pg, v2 = red['p0_rho'], red['p0_pg'], red['p0_v2']
    pram = 0.5 * rho * v2
    prr, psrc = prad_rr(red, est)
    ptot = pg + prr
    ratio_all = pram / ptot
    i1 = i0 + len(fi) - 1
    ratio = ratio_all[:, i0:i1]
    kp = red['p0_kap'][:, i0:i1]
    rph = r_photo(red['r'], red['tau_mean']) / RSUN
    ric, rhoic, Tic, fmlt, kic = ic_profile(kap)
    tau_ic = np.concatenate([np.cumsum((kic * rhoic * np.gradient(ric))[::-1])[::-1]])
    tau_ic_top = tau_ic - np.interp(red['x1f'][-1], ric, tau_ic)   # tau from OUR top
    j1 = np.where(tau_ic_top >= 1)[0][-1]
    rph_ic = ric[j1] / RSUN
    fig = plt.figure(figsize=(F3_W / 100, F3_H / 100), dpi=100)
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)

    def arc(ax, r, **kw):
        if F3_RMIN <= r <= F3_RMAX:
            ax.plot(r * np.cos(th), r * np.sin(th), **kw)
    axs = [wedge_axes(fig, cx) for cx in F3_CX]
    for ax, z, cm, (lo, hi) in ((axs[0], ratio, RAM_CMAP, F3_RAM), (axs[1], kp, 'GnBu', F3_KAP)):
        pc = ax.pcolormesh(X, Y, z, cmap=cm, vmin=lo, vmax=hi, shading='flat', rasterized=True)
        wedge_lines(ax, 0)
        arc(ax, rph, color='b', ls='--', lw=3)
        if ax is axs[0]:
            pc0 = pc
        else:
            pc1 = pc
    for cx, t1, t2 in ((F3_CX[0], r'3D $\phi=0$ plane', 'Ram Pressure vs Total Pressure'),
                       (F3_CX[1], r'3D $\phi=0$ plane', 'Rosseland Opacity'),
                       (F3_CX[2], '1D MESA profile', 'Rosseland Opacity')):
        fig.text(cx / F3_W, 1 - 38 / F3_H, t1, ha='center', va='center', fontsize=17)
        fig.text(cx / F3_W, 1 - 63 / F3_H, t2, ha='center', va='center', fontsize=17)
    cax = px_axes(fig, 126.5, 380.5, 430, 452)
    cb = fig.colorbar(pc0, cax=cax, orientation='horizontal', extend='max')
    cb.set_ticks([0, .1, .2, .3, .4, .5])
    cb.set_label(r'$\frac{1}{2}\rho v^2/(P_{\rm rad}+P_{\rm gas})$')
    cax = px_axes(fig, 592, 844, 430, 452)
    cb = fig.colorbar(pc1, cax=cax, orientation='horizontal', extend='both')
    cb.set_ticks([.5, 1, 1.5, 2])
    cb.set_label(r'$\kappa_{\rm Rosseland}$ [cm$^2$g$^{-1}$]')
    ax = axs[2]
    half = x2f >= np.pi / 2 - 1e-12
    Xl, Yl = X[half], Y[half]          # left half: theta >= 90 deg is x <= 0
    kic_f = np.interp(0.5 * (fi[1:] + fi[:-1]), ric / RSUN, kic)
    ax.pcolormesh(Xl, Yl, np.tile(kic_f, (half.sum() - 1, 1)), cmap='GnBu', vmin=F3_KAP[0],
                  vmax=F3_KAP[1], shading='flat', rasterized=True)
    wedge_lines(ax, 0)
    conv = fmlt > 1e-3     # the supplied profile has no MLT flux: stays empty (no grey zones)
    edges = np.flatnonzero(np.diff(conv.astype(int)))
    starts = list(ric[edges[np.diff(conv.astype(int))[edges] > 0] + 1] / RSUN)
    ends = list(ric[edges[np.diff(conv.astype(int))[edges] < 0]] / RSUN)
    if conv[0]:
        starts = [ric[0] / RSUN] + starts
    if conv[-1]:
        ends = ends + [ric[-1] / RSUN]
    thr = np.linspace(np.pi / 3, np.pi / 2, 100)
    arc(ax, rph_ic, color='k', ls='--', lw=3)
    arc(ax, rph, color='b', ls='--', lw=3)
    h = [plt.Line2D([], [], color='b', ls='--', lw=3, label=r'$r_{\rm photo,\,3D}$'),
         plt.Line2D([], [], color='k', ls='--', lw=3, label=r'$r_{\rm photo,\,MESA}$')]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(F3_CX[2] / F3_W, 1 - 462 / F3_H), fontsize=16,
               framealpha=0.8)
    # Fe kappa_R peak (max where 1.2e5 < T < 2.5e5 K; the He II bump near 50 Rsun, T ~4e4 K,
    # outgrows it after ~3 d) of the IC and of this dump's shell-mean kappa_R
    wi = (Tic > 1.2e5) & (Tic < 2.5e5)
    fe_ic = float(ric[wi][np.argmax(kic[wi])] / RSUN)
    fe_ic_k = float(kic[wi].max())
    wr = (red['T_sh'] > 1.2e5) & (red['T_sh'] < 2.5e5)
    fe_3d = float(red['r'][wr][np.argmax(red['kap_sh'][wr])] / RSUN)
    fe_3d_k = float(red['kap_sh'][wr].max())
    wh = (red['T_sh'] > 3e4) & (red['T_sh'] < 6e4)
    he_3d = float(red['r'][wh][np.argmax(red['kap_sh'][wh])] / RSUN) if wh.any() else np.nan
    he_3d_k = float(red['kap_sh'][wh].max()) if wh.any() else np.nan
    if stamp:
        time_stamp(fig, red['time'] / DAY)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    cz = [[round(float(a), 2), round(float(b), 2)] for a, b in zip(starts, ends)]
    return dict(t_d=round(float(red['time'] / DAY), 4),
                r_photo_3D_Rsun=round(float(rph), 3),
                prad=psrc,
                r_photo_IC_Rsun=round(float(rph_ic), 3), ic_convective_zones_Rsun=cz,
                kappa_edd=float(KAPPA_EDD),
                fe_peak_ic_Rsun=round(fe_ic, 3), fe_peak_ic_kappa=round(fe_ic_k, 4),
                fe_peak_3d_Rsun=round(fe_3d, 3), fe_peak_3d_kappa=round(fe_3d_k, 4),
                heii_bump_3d_Rsun=round(he_3d, 3), heii_bump_3d_kappa=round(he_3d_k, 4),
                kap_over_edd_max_3d=float(np.max(red['p0_kap'] / KAPPA_EDD)),
                ram_ratio_max_below_rphoto=float(np.nanmax(
                    ratio_all[:, red['r'] / RSUN < rph])) if rph == rph else None)


# BSG_RECHOOSE=1: redraw all movie frames (colour limits are fixed constants since 11-06)
RECHOOSE = os.environ.get('BSG_RECHOOSE', '0') == '1'

# Fig 4 (static figure and movie frames): his colour limits, dark plate as on the page (paper_fig4_dark)
F4_V = 60.0            # km/s
F4_RHO = (5.0, 8.0)    # 1e-8 g/cm^3
F4_T = (1.5, 1.7)      # 1e5 K
FIG4_BG = '#000000'
FIG4_RC = {'figure.facecolor': FIG4_BG, 'axes.facecolor': FIG4_BG, 'savefig.facecolor': FIG4_BG,
           'savefig.edgecolor': FIG4_BG, 'text.color': 'white', 'axes.labelcolor': 'white',
           'axes.edgecolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white'}
F4_W, F4_H = 1444, 669
F4_AX = ((118.0, 494.5), (570.5, 947.5), (1022.5, 1400.5))
F4_CB = ((132.0, 480.0), (587.0, 931.0), (1040.0, 1383.0))


def fig4_style():
    """context manager: the dark Fig 4 style (use around figure creation .. savefig)"""
    return matplotlib.rc_context(FIG4_RC)


def fig4(red, out, red_newest=None, stamp=False):
    with fig4_style(), matplotlib.rc_context(fs(16.5, title=20, pad=7)):
        return _fig4(red, out, stamp)


def theta_phi_axes(ax, ylabel):
    deg = ['0'] + ['%d\u00b0' % x for x in (15, 30, 45, 60)]
    ax.set_xlim(0, 60)
    ax.set_ylim(120, 60)
    ax.set_xticks([0, 15, 30, 45, 60])
    ax.set_xticklabels(deg)
    ax.set_yticks([60, 75, 90, 105, 120])
    ax.set_yticklabels(['%d\u00b0' % x for x in (60, 75, 90, 105, 120)])
    ax.set_xlabel(r'$\phi$ Coordinate')
    if ylabel:
        ax.set_ylabel(r'$\theta$ Coordinate', labelpad=9)


def _fig4(red, out, stamp=False):
    x2f, x3f = np.degrees(red['x2f']), np.degrees(red['x3f'])
    r36 = red['r'][int(red['i36'])] / RSUN
    vr = red['s36_vr'] / 1e5
    rho = red['s36_rho'] / 1e-8
    T = red['s36_T'] / 1e5
    fig = plt.figure(figsize=(F4_W / 100, F4_H / 100), dpi=100)
    info = {}
    specs = [(vr, 'seismic', -F4_V, F4_V, r'$v_r$ [km/s]', 'Radial Velocity',
              [-50, -25, 0, 25, 50], None),
             (rho, 'Blues', F4_RHO[0], F4_RHO[1], r'$\rho$ [$10^{-8}\,$g$\,$cm$^{-3}$]', 'Gas Density',
              [5, 6, 7, 8], None),
             (T, 'Blues_r', F4_T[0], F4_T[1], r'$T$ [$10^{5}\,$K]', 'Gas Temperature',
              [1.5, 1.55, 1.6, 1.65, 1.7], '%.2f')]
    for i, (z, cm, lo, hi, cl, tl, tk, fm) in enumerate(specs):
        ax = px_axes(fig, F4_AX[i][0], F4_AX[i][1], 59.5, 436.5)
        zmean_phys = float(z.mean())
        # z is (k = phi, j = theta): transpose for (x = phi, y = theta)
        pc = ax.pcolormesh(x3f, x2f, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat', rasterized=True)
        theta_phi_axes(ax, i == 0)
        ax.set_title(tl, pad=11)
        cax = px_axes(fig, F4_AX[i][0], F4_AX[i][1], 549, 567.5)
        cb = fig.colorbar(pc, cax=cax, orientation='horizontal', extend='both')
        cb.set_ticks(tk)
        if fm:
            cb.ax.set_xticklabels([fm % x for x in tk])
        cb.set_label(cl)
        info[tl] = dict(range=[lo, hi], range_source='paper', plotted='value',
                        min=float(z.min()), max=float(z.max()), mean_physical=zmean_phys,
                        frac_beyond=float(np.mean((z < lo) | (z > hi))))
    if stamp:
        time_stamp(fig, red['time'] / DAY, 'white')
    fig.savefig(out, dpi=100)
    plt.close(fig)
    return dict(t_d=round(float(red['time'] / DAY), 4), r_cell_Rsun=round(float(r36), 4),
                maps=info)


# ---------------------------------------------------------------------------------------
# emergent-flux ("intensity") map.  Ma+2026 show no intensity map; they use the radial
# radiative flux of the outermost layer only for L_bol (A6-A7).  In our M1 run the
# cell-centred F_r above the photosphere is NOT the escaping flux: in optically thin cells
# (c dt/dx >> 1) the implicit cell flux is not a transport flux (streaks, F_r < 0), and the
# mean-tau layer (tau_bar = 1) cuts the thin infalling channels where most columns are
# already thin (sk67 dump 57: median column tau 0.03 at that layer).  Since 10-06 the map
# is therefore F_r at each column's OWN tau_R = 1 point: the first cell outward (from the
# column's interior) whose column tau, integrated from the top (the same cell-centre tau
# as tau_mean / A9), is < 1.  See /viper/ptmp2/jinma/intmap_1006/RESULTS.md.
SIGMA_SB = 5.670374419e-5


def int_reduce(n, hfile, mfile, x1f, kap, red=None):
    """F_r (k, j) at each column's own tau_R = 1 cell, and that cell's r and tau; cached in
    cache/<arm>/int2_NNNNN.npz (keyed by both dumps' size and mtime)."""
    cf = os.path.join(CACHE, 'int2_%05d.npz' % n)
    st = stamp(hfile, mfile)
    if os.path.exists(cf):
        c = np.load(cf)
        if c['stamp'].shape == st.shape and np.allclose(c['stamp'], st, rtol=0, atol=0):
            return dict(c)
    dh = bc.read_binary(hfile)
    t = dh['time']
    w = assemble(dh)
    del dh
    rho, eint = w['dens'], w['eint']
    T = (GAMMA - 1.0) * eint / rho * np.float32(MU_MU_OVER_K)
    kr = kap(rho, T)
    nx3, nx2, nx1 = rho.shape
    dr = np.diff(x1f)
    tau = np.flip(np.cumsum(np.flip((kr * rho).astype(np.float64) * dr[None, None, :],
                                    axis=2), axis=2), axis=2)
    del w, rho, eint, T, kr
    thin = tau < 1.0
    iown = np.where(thin.any(axis=2), np.argmax(thin, axis=2), nx1 - 1)
    tau_own = np.take_along_axis(tau, iown[..., None], 2)[..., 0]
    del tau, thin
    d = bc.read_binary(mfile)
    assert abs(d['time'] - t) < 1e-6 * max(1.0, t), (hfile, mfile)
    a = assemble(d)
    del d
    fr = np.take_along_axis(a['m1_f1'], iown[..., None], 2)[..., 0].astype(np.float64)
    del a
    r = 0.5 * (x1f[1:] + x1f[:-1])
    res = dict(stamp=st, n=n, time=t, i_own=iown, r_own=r[iown], tau_own=tau_own,
               fr_own=fr, x2f=np.linspace(np.pi / 3, 2 * np.pi / 3, nx2 + 1),
               x3f=np.linspace(0.0, np.pi / 3, nx3 + 1))
    np.savez(cf, **res)
    return res


def wedge_mean(z, x2f):
    w = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :] * np.ones_like(z)
    return float((z * w).sum() / w.sum())


# not in the paper: fixed limits set by hand 11-06 (no his to copy), Blues_r as his T map
INT_RANGES = dict(int_lo=0.0, int_hi=2.0, teff_lo=11000.0, teff_hi=19000.0)


def int_ranges(ired_newest=None):
    return dict(INT_RANGES)


def fig_int(ired, rng, out, dpi=100):
    """Left: F_r/<F_r> at each column's own tau_R = 1 cell; right: T_eff = (F_r/sigma)^(1/4)
    [K].  Same panel style as his Fig 4.  Arrays are (k = phi, j = theta), plotted transposed."""
    with matplotlib.rc_context(fs(16.5, title=20, pad=7)):
        return _fig_int(ired, out, dpi)


def _fig_int(ired, out, dpi=100):
    x2f, x3f = np.degrees(ired['x2f']), np.degrees(ired['x3f'])
    fr = ired['fr_own'].astype(np.float64)
    fm = wedge_mean(fr, ired['x2f'])
    q = fr / fm
    te = (np.maximum(fr, 0.0) / SIGMA_SB) ** 0.25
    ro = ired['r_own'].astype(np.float64)
    lrat = 4 * np.pi * wedge_mean(ro ** 2 * fr, ired['x2f'])
    r1, r50, r99 = np.percentile(ro / RSUN, [1, 50, 99])
    fig = plt.figure(figsize=(9.92, 6.69), dpi=100)
    for i, (z, lo, hi, tl, cl, tk, x0) in enumerate((
            (q, INT_RANGES['int_lo'], INT_RANGES['int_hi'], 'Emergent Flux',
             r'$F_r/\langle F_r\rangle$', [0, 0.5, 1, 1.5, 2], 118.0),
            (te, INT_RANGES['teff_lo'], INT_RANGES['teff_hi'], 'Local Effective Temperature',
             r'$T_{\rm eff}=(F_r/\sigma_{\rm SB})^{1/4}$ [K]', [11000, 13000, 15000, 17000, 19000], 570.5))):
        ax = px_axes(fig, x0, x0 + 377, 59.5, 436.5)
        pc = ax.pcolormesh(x3f, x2f, z.T, cmap='Blues_r', vmin=lo, vmax=hi, shading='flat',
                           rasterized=True)
        theta_phi_axes(ax, i == 0)
        ax.set_title(tl, pad=11)
        cax = px_axes(fig, x0, x0 + 377, 549, 567.5)
        cb = fig.colorbar(pc, cax=cax, orientation='horizontal', extend='both')
        cb.set_ticks(tk)
        cb.set_label(cl)
    time_stamp(fig, float(ired['time']) / DAY)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return dict(t_d=round(float(ired['time'] / DAY), 4), r_ph_Rsun=round(float(r50), 3),
                r_own_p1_p99_Rsun=[round(float(r1), 3), round(float(r99), 3)],
                tau_own_median=float(np.median(ired['tau_own'])), F_mean=fm,
                L_4pir2F=lrat, L_over_Lstar=lrat / L_STAR,
                Teff_mean=(fm / SIGMA_SB) ** 0.25, q_min=float(q.min()),
                q_max=float(q.max()), q_sd=float(q.std()), fneg=float((fr < 0).mean()))


OMEGA_WEDGE = np.pi / 3.0     # solid angle of our wedge [sr]


def psd_ls(t, flux):
    """Lomb-Scargle PSD in ppm^2/muHz, normalised so int PSD df (0..Nyquist) =
    variance."""
    from astropy.timeseries import LombScargle
    y = (flux - 1.0) * 1e6
    T = t.max() - t.min()
    dtm = np.median(np.diff(t))
    fmin, fmax = 1.0 / T, 0.5 / dtm
    f = np.arange(fmin, fmax, 0.1 / T)                   # Hz, oversample 10
    p = LombScargle(t, y, normalization='psd').power(f)
    df = f[1] - f[0]
    p = p * np.var(y) / (p.sum() * df)                    # ppm^2/Hz, Parseval
    return f * 1e6, p * 1e-6                              # muHz, ppm^2/muHz


def hst_series(hst, key):
    t = hst['time']
    y = hst[key]
    m = t > 0
    t, y = t[m], y[m]
    _, ui = np.unique(t, return_index=True)               # restarts may repeat rows
    return t[ui], y[ui]


def coherent_mode(t, x, pmin_d=3.0, pmax_d=4.5):
    """Least-squares sine (free period in pmin..pmax d) to the detrended relative
    fluctuation x(t): the l = 0 radial mode (RADIAL_MODE_1003.md)."""
    best = None
    for P in np.linspace(pmin_d, pmax_d, 1501) * DAY:
        A = np.vstack([np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P)]).T
        c = np.linalg.lstsq(A, x, rcond=None)[0]
        r = np.sum((x - A @ c) ** 2)
        if best is None or r < best[0]:
            best = (r, P, c, A @ c)
    r, P, c, model = best
    return dict(P_d=P / DAY, f_muHz=1e6 / P, amp=float(np.hypot(*c)),
                var_frac=float(1.0 - r / np.sum(x ** 2)), model=model)


def lc_treatments(t, L):
    """Paper eq. 8 (linear detrend, fluctuation x sqrt(Omega/4pi)) and the variant
    with the coherent radial mode unscaled (coherent over the whole star) and only the
    residual (stochastic, ~independent patches) scaled."""
    s = np.sqrt(OMEGA_WEDGE / (4 * np.pi))
    x = L / np.polyval(np.polyfit(t, L, 1), t) - 1.0
    cm = coherent_mode(t, x)
    f_eq8 = 1.0 + s * x
    f_mode = 1.0 + cm['model'] + s * (x - cm['model'])
    return x, f_eq8, f_mode, cm


def logbin(f, p, nper=10):
    e = 10 ** (np.arange(np.floor(np.log10(f.min()) * nper),
                         np.ceil(np.log10(f.max()) * nper) + 1) / nper)
    i = np.digitize(f, e)
    ks = np.unique(i)
    return (np.array([f[i == k].mean() for k in ks]),
            np.array([p[i == k].mean() for k in ks]))


def psd_stats(f, p, fmode=None, T=None, fit=(1.25, 400.0)):
    """f [muHz], p [(dF/F)^2/muHz]. Peak, cumulative 50/80 %, red-noise slope and a
    Harvey-type fit P0/(1+(f/nu_char)^gamma) + C (SLF convention) on the log-binned
    PSD; the coherent-mode bins (|f - fmode| < 2/T) are left out of peak and fits."""
    from scipy.optimize import curve_fit
    sel = f >= 1.25
    c = np.cumsum(p[sel])
    out = dict(f50_muHz=float(f[sel][np.searchsorted(c / c[-1], 0.5)]),
               f80_muHz=float(f[sel][np.searchsorted(c / c[-1], 0.8)]),
               rms=float(np.sqrt(np.sum(p) * (f[1] - f[0]))))
    keep = np.ones_like(f, bool)
    if fmode is not None:
        keep = np.abs(f - fmode) > 2e6 / T
    fb, pb = logbin(f[keep & sel], p[keep & sel])
    out['peak_raw_muHz'] = float(f[sel][np.argmax(p[sel])])
    out['peak_binned_noMode_muHz'] = float(fb[np.argmax(pb)])
    out['peak_fP_noMode_muHz'] = float(fb[np.argmax(fb * pb)])
    m = (fb >= 10) & (fb <= 100)
    out['slope_10_100'] = float(np.polyfit(np.log10(fb[m]), np.log10(pb[m]), 1)[0])
    m = (fb >= fit[0]) & (fb <= fit[1])

    def harvey(lf, lp0, lnu, g, lc):
        return np.log10(10 ** lp0 / (1 + (10 ** lf / 10 ** lnu) ** g) + 10 ** lc)
    try:
        q = curve_fit(harvey, np.log10(fb[m]), np.log10(pb[m]),
                      p0=[np.log10(pb[m][:3].mean()), 1.0, 2.0, np.log10(pb[m].min())],
                      maxfev=20000)[0]
        out['harvey'] = dict(P0=float(10 ** q[0]), nu_char_muHz=float(10 ** q[1]),
                             gamma=float(q[2]), C=float(10 ** q[3]))
    except Exception as e:  # noqa: BLE001
        out['harvey'] = dict(error=str(e))
    return out


def ke_flatness(hst, t0_d=25.0, block_d=5.0):
    t, ke = hst_series(hst, 'KE_int')
    td = t / DAY
    if td.max() < t0_d + 2 * block_d:
        return None
    edges = np.arange(t0_d, td.max() + 1e-9, block_d)
    bm = np.array([ke[(td >= a) & (td < a + block_d)].mean() for a in edges[:-1]])
    w = td >= t0_d
    sl = np.polyfit(td[w], ke[w], 1)[0]
    return dict(t0_d=t0_d, t1_d=float(td.max()), block_d=block_d,
                block_means=[float(v) for v in bm],
                dev_max_frac=float(np.max(np.abs(bm / bm.mean() - 1))),
                slope_frac_per_10d=float(sl * 10 / ke[w].mean()),
                mean=float(ke[w].mean()))


F6_W, F6_H = 1446, 796
F6_FX = (1.18, 277.8)     # his PSD frequency axis [muHz] (read off; right end = Nyquist of 30 min)
F6_PY = (7e-9, 8e-4)      # his PSD axis
F6_TX = (30.0, 57.4)
F6_FY = (0.7, 1.2)
F6_BAND = (0.97, 1.03)    # his cyan band (read off: 0.972-1.030)


def fig6(hst, out, t0_d=30.0, t1_d=None):
    """Light curve + PSD over the analysis window t0_d..t1_d (default 30 d to the end
    of the hst, the paper's steady window), in his layout; the TESS panel is his raster."""
    with matplotlib.rc_context(fs(15, title=17, leg=14, pad=7)):
        return _fig6(hst, out, t0_d, t1_d)


def _fig6(hst, out, t0_d=30.0, t1_d=None):
    t, L = hst_series(hst, 'L_top')
    tmax_d = t.max() / DAY
    t1 = tmax_d if t1_d is None else t1_d
    w = (t >= t0_d * DAY) & (t <= t1 * DAY + 1.0)
    if w.sum() < 64:                                     # window not reached yet
        w = t <= t1 * DAY + 1.0
    t, L = t[w], L[w]
    s = np.sqrt(OMEGA_WEDGE / (4 * np.pi))
    x, flux, fluxm, cm = lc_treatments(t, L)
    fig = plt.figure(figsize=(F6_W / 100, F6_H / 100), dpi=100)
    ax = px_axes(fig, 110.5, 1351.5, 32.5, 272.5)
    td = t / DAY
    ax.plot(td, flux, 'b-', lw=2, label='3D Simulation, bolometric')
    ax.axhspan(F6_BAND[0], F6_BAND[1], color='c', alpha=0.3, lw=0,
               label='SK -67 133 range, TESS band')
    ax.set_xlim(*F6_TX)
    ax.set_ylim(*F6_FY)
    ax.set_xticks([30, 35, 40, 45, 50, 55])
    ax.set_yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    ax.set_xlabel('Time since Simulation Starts [days]')
    ax.set_ylabel('Normalized Flux')
    ax.legend(loc='lower left', framealpha=0.8)
    T = t.max() - t.min()
    info = dict(t_range_d=[round(float(td.min()), 4), round(float(td.max()), 4)],
                n_samples=int(len(td)), flux_min=float(flux.min()),
                flux_max=float(flux.max()), cadence_s=float(np.median(np.diff(t))),
                L_top_mean_wedge_erg_s=float(L.mean()),
                L_top_mean_whole_star_erg_s=float(L.mean() * 4 * np.pi / OMEGA_WEDGE),
                rms_wedge_detrended=float(x.std()),
                rms_eq8=float((flux - 1).std()), rms_mode=float((fluxm - 1).std()),
                mode=dict(P_d=cm['P_d'], f_muHz=cm['f_muHz'], amp_rel=cm['amp'],
                          var_frac=cm['var_frac']),
                ke_flat=ke_flatness(hst))
    ax2 = px_axes(fig, 110.5, 627.5, 384.5, 704.5)
    if len(t) >= 8:
        tw = ax2.twinx()
        f, p = psd_ls(t, flux)
        p = p * 1e-12      # ppm^2/muHz -> (dF/F)^2/muHz, the paper's actual unit despite its axis label
        ax2.loglog(f, p, 'b-', lw=2)
        sel = (f >= F6_FX[0]) & (f <= F6_FX[1])
        c = np.cumsum(p[sel])
        tw.plot(f[sel], c / c[-1], '--', color='navy', lw=3.5)
        info['psd_eq8'] = psd_stats(f, p, cm['f_muHz'], T)
        f2, p2 = psd_ls(t, fluxm)
        info['psd_mode'] = psd_stats(f2, p2 * 1e-12, cm['f_muHz'], T)
        tw.set_ylim(0, 1)
        tw.set_yticks([0, .2, .4, .6, .8, 1])
        tw.set_ylabel('Cumulative Integrated PSD')
        info['psd_f_range_muHz'] = [float(f.min()), float(f.max())]
        info['f80_muHz'] = info['psd_eq8']['f80_muHz']
    ax2.set_xlim(*F6_FX)
    ax2.set_ylim(*F6_PY)
    ax2.set_yticks([1e-8, 1e-7, 1e-6, 1e-5, 1e-4])
    ax2.set_xlabel(r'Frequency [$\mu$Hz]')
    ax2.set_ylabel(r'PSD [ppm$^2$ / $\mu$Hz]')
    ax2.set_title('3D Simulation', pad=9)
    # lower-right panel: his raster, at its own pixels (paper_fig6_src.png = page paper_fig6.png)
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'paper_fig6_src.png')
    a = plt.imread(src)
    fig.figimage(a[345:796, 700:1446], xo=700, yo=0, origin='upper', zorder=3)
    info['pre_window'] = dict(n=0, flux_min=None, flux_max=None)
    fig.savefig(out, dpi=100)
    plt.close(fig)
    return info


# ---------------------------------------------------------------------------------------
def main():
    global CACHE, IC, INP, ARM
    ap = argparse.ArgumentParser()
    ap.add_argument('--arm', default='repro', choices=sorted(ARMS))
    ap.add_argument('--run', default=None, help='override the arm\'s run dir')
    ap.add_argument('--snap', type=int, default=None, help='dump index for Figs 3, 4')
    ap.add_argument('--rmax-plot', type=float, default=F3_RMAX)
    ap.add_argument('--psd-t0', type=float, default=30.0,
                    help='Fig 6 analysis window start [d] (paper: 30)')
    ap.add_argument('--psd-t1', type=float, default=None,
                    help='Fig 6 analysis window end [d] (default: hst end)')
    ap.add_argument('--prad-est', action='store_true',
                    help='without an m1_vet dump, use prad/est_NNNNN.npz (offline f_K) '
                    'for '
                         'P_rad in Fig 3 instead of E/3')
    a = ap.parse_args()
    ARM = a.arm
    if a.run is None:
        a.run = ARMS[ARM]['run']
    IC, INP = ARMS[ARM]['ic'], ARMS[ARM]['inp']
    CACHE = cache_dir(a.run)
    os.makedirs(CACHE, exist_ok=True)
    x1f, gfile = radial_faces(a.run)
    kap = Kappa()
    dumps = dump_list(a.run)
    reds = [reduce_dump(n, h, m, x1f, kap, v) for n, h, m, v in dumps]
    hst = read_hst(os.path.join(a.run, 'bsg3d.user.hst'))
    if a.snap is None:   # his Figs 3, 4: the snapshot at t = 40 d (the newest dump until then)
        red = min(reds, key=lambda c: abs(float(c['time']) / DAY - 40.0))
    else:
        red = [c for c in reds if int(c['n']) == a.snap][0]
    i2 = fig2(reds, hst, os.path.join(FIGS, 'out', 'ours_fig2.png'), a.rmax_plot)
    est = (os.path.join(BSG, 'prad', 'est_%05d.npz' % int(red['n'])) if a.prad_est
           else None)
    i3 = fig3(red, a.rmax_plot, kap, os.path.join(FIGS, 'out', 'ours_fig3.png'), est)
    i4 = fig4(red, os.path.join(FIGS, 'out', 'ours_fig4.png'))
    i6 = fig6(hst, os.path.join(FIGS, 'out', 'ours_fig6.png'), a.psd_t0, a.psd_t1)
    nn, hh, mm, _ = [x for x in dumps if x[0] == int(red['n'])][0]
    ired = int_reduce(nn, hh, mm, x1f, kap, red)
    irng = int_ranges()
    iint = fig_int(ired, irng, os.path.join(FIGS, 'out', 'ours_int.png'))
    meta = dict(arm=ARM, run=a.run, input=INP, ic=IC, grid_file=gfile,
                r_bot_Rsun=round(float(x1f[0] / RSUN), 3),
                r_top_Rsun=round(float(x1f[-1] / RSUN), 3),
                grid=[len(x1f) - 1, int(red['s36_rho'].shape[1]),
                      int(red['s36_rho'].shape[0])],
                dumps_used=[int(c['n']) for c in reds],
                snapshot_fig3_4=int(red['n']), fig2=i2, fig3=i3, fig4=i4, fig6=i6,
                fig_int=iint)
    with open(os.path.join(CACHE, 'numbers.json'), 'w') as fo:
        json.dump(meta, fo, indent=1, default=float)
    print(json.dumps(meta, indent=1, default=float))


if __name__ == '__main__':
    main()
