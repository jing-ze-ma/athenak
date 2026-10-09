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
# Fig 4 radius [Rsun]: kept at 40.64 as on the arm-2 page (the paper uses 36 = its Fe
# peak 35 + 1).  The repro run's own Fe kappa_R peak (shell mean) is measured per dump
# (fe_peak in numbers.json) and quoted in the caption.
R_SLICE = 40.64


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
        if os.path.exists(m):
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
    a = np.loadtxt(path, ndmin=2)
    return {nm: a[:, i] for i, nm in enumerate(names)}


# ---------------------------------------------------------------------------------------
# figures
def mark_steady(ax, tmax_d):
    ax.axvline(30.0, color='purple', ls='--', lw=2)
    if tmax_d < 30:
        ax.axvspan(0, 30, color='0.85', alpha=0.25, lw=0)


def fig2(reds, hst, out, rmax_plot):
    t = np.array([c['time'] for c in reds]) / DAY
    rph = np.array([r_photo(c['r'], c['tau_mean']) for c in reds]) / RSUN
    eg = np.array([c['E_grav'] for c in reds])
    ui = np.array([c['U_gas'] + c['U_rad'] for c in reds])
    ek = np.array([c['E_k'] for c in reds])
    menv = (4 * np.pi / reds[0]['omega']) * np.array([c['M_wedge'] for c in reds])
    tmax = max(t.max(), hst['time'].max() / DAY)
    xr = (0, max(tmax, 1e-3))
    fig, axs = plt.subplots(5, 1, figsize=(9.75, 9.594), constrained_layout=True)
    ax = axs[0]
    ax.plot(t, (eg - eg[0]) / 1e44, 'g-o', ms=3, label=r'change of gravitational binding '
            r'energy $\Delta E_{\rm grav}$ (dumps)')
    ax.plot(t, (ui - ui[0]) / 1e44, '-o', color='orange', ms=3,
            label=r'change of internal energy $\Delta U_{\rm rad}+\Delta U_{\rm gas}$'
            r' (dumps)')
    th = hst['time'] / DAY
    uh = hst['E_rad'] + hst['e_gas']
    ax.plot(th, (uh - uh[0]) / 1e44, '-', color='orange', lw=0.8, alpha=0.6,
            label=r'$\Delta$(E_rad + e_gas), user.hst')
    ax.set_ylabel(r'Energy [$10^{44}$ erg]')
    ax.legend(fontsize=7, loc='lower left')
    ax = axs[1]
    ax.semilogy(t[t > 0], ek[t > 0] / 1e44, 'b-o', ms=3,
                label=r'gas kinetic energy $E_{\rm k}$'
                ' (dumps)')
    m = hst['KE_int'] > 0
    ax.semilogy(th[m], hst['KE_int'][m] / 1e44, 'b-', lw=0.8, alpha=0.6,
                label=r'KE_int (r < r_int), user.hst')
    ax.set_ylabel(r'Energy [$10^{44}$ erg]')
    ax.legend(fontsize=7, loc='lower right')
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
    specs = [(axs[2], mdot, 'RdBu', -5, 5, r'$\dot M_r/M_{\rm envelope}$ [yr$^{-1}$]',
              r'radial mass flow $\dot M_r/M_{\rm envelope}$'),
             (axs[3], rho, 'Blues', -11, -3.5, r'$\log(\rho)$ [g cm$^{-3}$]',
              'angular averaged gas density'),
             (axs[4], T, 'Blues', 3.8, 6.2, r'$\log(T)$ [K]',
              'angular averaged gas temperature')]
    for ax, z, cm, lo, hi, cl, tl in specs:
        pc = ax.pcolormesh(te, rf, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat')
        cb = fig.colorbar(pc, ax=ax, extend='both', pad=0.01)
        cb.set_label(cl)
        ax.plot(t, rph, 'b--', lw=2, label=r'$r_{\rm photo,\,3D}$')
        if rf[-1] < rmax_plot - 0.1:
            ax.axhline(rf[-1], color='k', ls=':', lw=1)
            ax.text(0.98, (rf[-1] - 20) / (rmax_plot - 20) + 0.02,
                    r'our top %.0f $R_\odot$' % rf[-1], transform=ax.transAxes,
                    ha='right', va='bottom', fontsize=8)
        ax.set_ylim(20, rmax_plot)
        ax.set_ylabel(r'Radial Coordinate [$R_\odot$]')
        ax.text(0.02, 0.92, tl, transform=ax.transAxes, va='top', fontsize=9)
    axs[4].legend(loc='lower left', fontsize=8)
    for ax in axs:
        ax.set_xlim(*xr)
        mark_steady(ax, tmax)
    axs[4].set_xlabel('Time since Simulation Starts [days]')
    fig.suptitle('Ours (AthenaK reproduction run, 1 %% T seed, split work, vet_gd), '
                 't = 0 - %.2f d: start-up '
                 'transient (t < 30 d)' % tmax, fontsize=10)
    fig.savefig(out, dpi=148.1)
    plt.close(fig)
    return dict(t_dumps_d=[round(float(x), 4) for x in t],
                t_hst_max_d=round(float(tmax), 4),
                r_photo_3D_Rsun=[round(float(x), 3) for x in rph],
                M_envelope_whole_star_Msun=float(menv[-1] / 1.98847e33),
                dE_grav_1e44=float((eg[-1] - eg[0]) / 1e44),
                dU_int_1e44=float((ui[-1] - ui[0]) / 1e44), E_k_1e44=float(ek[-1] / 1e44))


def wedge_axes(ax, rmax_plot, rlabels=(25, 35, 45)):
    ax.set_aspect('equal')
    ax.axis('off')
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)
    for rr in (20, rmax_plot):
        ax.plot(rr * np.cos(th), rr * np.sin(th), 'k-', lw=1.2)
    for a in (np.pi / 3, 2 * np.pi / 3):
        ax.plot([20 * np.cos(a), rmax_plot * np.cos(a)],
                [20 * np.sin(a), rmax_plot * np.sin(a)],
                'k-', lw=1.2)
    for rr in rlabels:
        ax.plot(rr * np.cos(th), rr * np.sin(th), '-.', color='0.6', lw=0.8)
        a = np.pi / 3
        ax.text(rr * np.cos(a) + 1.5, rr * np.sin(a) - 1.0, r'$%d\,R_\odot$' % rr,
                fontsize=10)
    for deg in (75, 90, 105):
        a = np.radians(deg)
        ax.plot([20 * np.cos(a), rmax_plot * np.cos(a)],
                [20 * np.sin(a), rmax_plot * np.sin(a)],
                '-.', color='0.6', lw=0.8)
    for deg in (60, 75, 90, 105, 120):
        a = np.radians(deg)
        ax.text(1.06 * rmax_plot * np.cos(a), 1.06 * rmax_plot * np.sin(a),
                r'$%d^\circ$' % deg,
                ha='center', va='center', fontsize=10)
    ax.set_xlim(-0.62 * rmax_plot, 0.62 * rmax_plot)
    ax.set_ylim(0.5 * 20 * 0.9 + 0, 1.12 * rmax_plot)


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


def fig3(red, rmax_plot, kap, out, est=None):
    rf = red['x1f'] / RSUN
    x2f = red['x2f']
    sel = rf <= rmax_plot + 1e-9
    ni = int(sel.sum()) - 1
    R, TH = np.meshgrid(rf[:ni + 1], x2f)
    X, Y = R * np.cos(TH), R * np.sin(TH)
    rho, pg, v2 = red['p0_rho'], red['p0_pg'], red['p0_v2']
    pram = 0.5 * rho * v2
    prr, psrc = prad_rr(red, est)
    ptot = pg + prr
    ratio = (pram / ptot)[:, :ni]
    kp = red['p0_kap'][:, :ni]
    rph = r_photo(red['r'], red['tau_mean']) / RSUN
    # IC profile
    ic = np.loadtxt(IC)
    ric, rhoic = ic[:, 0], ic[:, 1]
    Tic = ic[:, 5]
    fmlt = ic[:, 6]
    kic = kap(rhoic, Tic)
    tau_ic = np.concatenate([np.cumsum((kic * rhoic * np.gradient(ric))[::-1])[::-1]])
    tau_ic_top = tau_ic - np.interp(red['x1f'][-1], ric, tau_ic)   # tau from OUR top
    j1 = np.where(tau_ic_top >= 1)[0][-1]
    rph_ic = ric[j1] / RSUN
    fig, axs = plt.subplots(1, 3, figsize=(14.44, 5.61))
    blues = LinearSegmentedColormap.from_list('ramp', ['white', 'blue', 'navy'])
    th = np.linspace(np.pi / 3, 2 * np.pi / 3, 200)
    ax = axs[0]
    pc = ax.pcolormesh(X, Y, ratio, cmap=blues, vmin=0, vmax=0.5, shading='flat',
                       rasterized=True)
    wedge_axes(ax, rmax_plot)
    ax.plot(rph * np.cos(th), rph * np.sin(th), 'b--', lw=2)
    ax.set_title('Ours 3D $\\phi=0$ plane (t = %.2f d)\nRam Pressure vs Total Pressure'
                 % (red['time'] / DAY), fontsize=11)
    cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend='max', shrink=0.75,
                      pad=0.02)
    plab = {'m1_vet': r'$P_{\rm rad}=D_{rr}E$ (vet_gd)',
            'offline vet_col port': r'$P_{\rm rad}=f_K E$ (offline vet_col)',
            'E/3': r'$P_{\rm rad}=E/3$'}[psrc]
    cb.set_label(r'$\frac{1}{2}\rho v^2/(P_{\rm rad}+P_{\rm gas})$,  ' + plab)
    ax = axs[1]
    pc = ax.pcolormesh(X, Y, kp, cmap='GnBu', vmin=0.05, vmax=2.5, shading='flat',
                       rasterized=True)
    wedge_axes(ax, rmax_plot)
    ax.plot(rph * np.cos(th), rph * np.sin(th), 'b--', lw=2)
    ax.set_title('Ours 3D $\\phi=0$ plane (t = %.2f d)\nRosseland Opacity'
                 % (red['time'] / DAY), fontsize=11)
    cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend='both', shrink=0.75,
                      pad=0.02)
    cb.set_label(r'$\kappa_{\rm Rosseland}$ [cm$^2$g$^{-1}$]  (Ma+2026 TOPS tables)')
    ax = axs[2]
    kic_f = np.interp(0.5 * (rf[1:ni + 1] + rf[:ni]), ric / RSUN, kic)
    half = x2f >= np.pi / 2 - 1e-12
    Xl, Yl = X[half], Y[half]          # left half: theta >= 90 deg is x <= 0
    pc = ax.pcolormesh(Xl, Yl, np.tile(kic_f, (half.sum() - 1, 1)), cmap='GnBu',
                       vmin=0.05,
                       vmax=2.5, shading='flat', rasterized=True)
    wedge_axes(ax, rmax_plot)
    # convective shells (F_MLT/F > 0) on the right half
    conv = fmlt > 1e-3
    edges = np.flatnonzero(np.diff(conv.astype(int)))
    starts = list(ric[edges[np.diff(conv.astype(int))[edges] > 0] + 1] / RSUN)
    ends = list(ric[edges[np.diff(conv.astype(int))[edges] < 0]] / RSUN)
    if conv[0]:
        starts = [ric[0] / RSUN] + starts
    if conv[-1]:
        ends = ends + [ric[-1] / RSUN]
    thr = np.linspace(np.pi / 3, np.pi / 2, 100)
    for a, b in zip(starts, ends):
        if b < 20 or a > rmax_plot:
            continue
        a, b = max(a, 20), min(b, rmax_plot)
        xs = np.concatenate([a * np.cos(thr), b * np.cos(thr[::-1])])
        ys = np.concatenate([a * np.sin(thr), b * np.sin(thr[::-1])])
        ax.fill(xs, ys, color='0.6', lw=0)
        if b - a > 2:
            ax.text(0.5 * (a + b) * np.cos(np.radians(70)),
                    0.5 * (a + b) * np.sin(np.radians(70)),
                    'convection', rotation=-20, ha='center', va='center', fontsize=10,
                    color='0.25')
    ax.plot([0, 0], [20, rmax_plot], '-.', color='0.6', lw=0.8)
    ax.plot(rph_ic * np.cos(th), rph_ic * np.sin(th), 'k--', lw=2,
            label=r'$r_{\rm photo,\,IC}$ = %.1f $R_\odot$' % rph_ic)
    ax.plot(rph * np.cos(th), rph * np.sin(th), 'b--', lw=2,
            label=r'$r_{\rm photo,\,3D}$ = %.1f $R_\odot$' % rph)
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
    ax.plot(fe_ic * np.cos(thr), fe_ic * np.sin(thr), '-', color='darkred', lw=1.5,
            label=r'IC Fe $\kappa$ peak %.1f $R_\odot$ ($\kappa/\kappa_{\rm Edd}$ = %.2f)'
            % (fe_ic, fe_ic_k / KAPPA_EDD))
    ax.set_title('The paper\'s 1D IC as run (%s, radiative only)\n'
                 'Rosseland Opacity ($F_{\\rm MLT}$ = 0: no convective shell marked)'
                 % os.path.basename(IC), fontsize=11)
    if rf[-1] < rmax_plot - 0.1:
        for k, a in enumerate(axs):
            a.plot(rf[-1] * np.cos(th), rf[-1] * np.sin(th), 'k:', lw=1.2,
                   label=(r'our top %.0f $R_\odot$' % rf[-1]) if k == 2 else None)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.06), fontsize=10)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.88, bottom=0.04, wspace=0.05)
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
                    ratio[:, red['r'][:ni] / RSUN < rph])) if rph == rph else None)


RANGES_JSON = os.path.join(FIGS, 'movie', 'ranges.json')
# repro (no seed): the Fig 4 shell is uniform to ~1e-6 at first; floor of the fixed relative
# rho / T colour range [%] so round-off is not painted as structure
REL_MIN_PCT = 0.1
# BSG_RECHOOSE=1: re-choose all fixed colour ranges from the newest dump in this process
# (overwrites ranges.json; update_page_repro.sh does this on every refresh)
RECHOOSE = os.environ.get('BSG_RECHOOSE', '0') == '1'
_DONE = set()


def shell_mean(red, key):
    """solid-angle-weighted mean of a shell slice s36_<key>"""
    x2f = red['x2f']
    w = (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :] * np.ones_like(red['s36_' + key])
    return float((red['s36_' + key] * w).sum() / w.sum())


def abs_ranges(red_newest):
    """Fig 4 absolute colour ranges for rho and T, fixed for all frames (movie and static
    figure): X_ref (1 -+ a), X_ref = shell mean at the Fig 4 shell in the newest dump,
    a = the fixed relative range in ranges.json (rho, T in %; 1-99 percentile of the
    newest dump when absent).  Stored in ranges.json (rho_ref [g/cm^3], T_ref [K],
    rho_abs, T_abs); delete those keys (or the file) to re-choose."""
    d = {}
    if os.path.exists(RANGES_JSON) and not (RECHOOSE and 'abs' not in _DONE):
        d = json.load(open(RANGES_JSON))
    if 'rho_abs' in d and 'T_abs' in d:
        return d
    _DONE.add('abs')
    for k in ('rho', 'T'):
        ref = shell_mean(red_newest, k)
        if k not in d:
            z = (red_newest['s36_' + k] / ref - 1) * 100
            lo, hi = np.percentile(z, [1, 99])
            m = max(float(max(abs(lo), abs(hi))), REL_MIN_PCT)
            d[k] = [-m, m]
        a = d[k][1] / 100.0
        d[k + '_ref'] = ref
        d[k + '_abs'] = [ref * (1 - a), ref * (1 + a)]
    d['abs_from_dump'] = int(red_newest['n'])
    d['abs_t_d'] = float(red_newest['time'] / DAY)
    os.makedirs(os.path.dirname(RANGES_JSON), exist_ok=True)
    json.dump(d, open(RANGES_JSON, 'w'), indent=1)
    return d


def fig4(red, out, red_newest=None):
    x2f, x3f = np.degrees(red['x2f']), np.degrees(red['x3f'])
    r36 = red['r'][int(red['i36'])] / RSUN
    vr = red['s36_vr'] / 1e5
    rho = red['s36_rho'] / 1e-8
    T = red['s36_T'] / 1e5
    fig, axs = plt.subplots(1, 3, figsize=(14.44, 6.69), constrained_layout=True)
    info = {}

    def rng(a, lo, hi):
        """paper range if our data overlap it, else our 2-98 percentiles"""
        p2, p98 = np.percentile(a, [2, 98])
        if p98 < lo or p2 > hi or (p98 - p2) < 0.05 * (hi - lo):
            return float(p2), float(p98), 'auto'
        return lo, hi, 'paper'
    vlim = 60.0
    vmax_ours = float(np.percentile(np.abs(vr), 99.5))
    if vmax_ours < 0.05 * vlim:
        vlim = max(vmax_ours, 1e-6)
        vmode = 'auto'
    else:
        vmode = 'paper'
    # rho, T: absolute values, fixed range = shell mean (newest dump) x (1 -+ a), the
    # same ranges as the movie frames (ranges.json)
    ar = abs_ranges(red if red_newest is None else red_newest)
    rlo, rhi = ar['rho_abs'][0] / 1e-8, ar['rho_abs'][1] / 1e-8
    tlo, thi = ar['T_abs'][0] / 1e5, ar['T_abs'][1] / 1e5
    specs = [(vr, 'seismic', -vlim, vlim, vmode, r'$v_r$ [km/s]', 'Radial Velocity',
              'both'),
             (rho, 'Blues', rlo, rhi, 'fixed',
              r'$\rho$ [$10^{-8}$ g cm$^{-3}$]  (range $\bar\rho\,(1\pm%.3g\%%)$, '
              r'$\bar\rho$ = %.3f)' % (ar['rho'][1], ar['rho_ref'] / 1e-8),
              'Gas Density', 'both'),
             (T, 'Blues_r', tlo, thi, 'fixed',
              r'$T$ [$10^5$ K]  (range $\bar T\,(1\pm%.3g\%%)$, $\bar T$ = %.4f)'
              % (ar['T'][1], ar['T_ref'] / 1e5), 'Gas Temperature', 'both')]
    for ax, (z, cm, lo, hi, mode, cl, tl, ext) in zip(axs, specs):
        zmean_phys = float(z.mean())
        # z is (k = phi, j = theta): transpose for (x = phi, y = theta)
        pc = ax.pcolormesh(x3f, x2f, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat',
                           rasterized=True)
        ax.set_xlim(0, 60)
        ax.set_ylim(120, 60)
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_xticklabels([r'$%d^\circ$' % x for x in (0, 15, 30, 45, 60)])
        ax.set_yticks([60, 75, 90, 105, 120])
        ax.set_yticklabels([r'$%d^\circ$' % x for x in (60, 75, 90, 105, 120)])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\phi$ Coordinate')
        ax.set_title(tl + ('' if mode == 'paper' else '  (own colour range)'))
        cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend=ext, shrink=0.9)
        cb.set_label(cl)
        cb.formatter.set_useOffset(False)
        cb.locator = matplotlib.ticker.MaxNLocator(5)
        cb.update_ticks()
        info[tl] = dict(range=[round(float(lo), 6), round(float(hi), 6)],
                        range_source=mode, plotted='value',
                        min=float(z.min()), max=float(z.max()), mean_physical=zmean_phys)
    axs[0].set_ylabel(r'$\theta$ Coordinate')
    fig.suptitle(r'Ours: $\theta-\phi$ plane at r = %.2f $R_\odot$ (nearest cell to %.2f), '
                 r't = %.2f d' % (r36, R_SLICE, red['time'] / DAY))
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


def int_ranges(ired_newest):
    """fixed ranges of the intensity map, stored in ranges.json: int_lo/int_hi (F/<F>) and
    teff_lo/teff_hi [K], 1-99 percentiles of the newest dump's own-tau=1 map."""
    d = json.load(open(RANGES_JSON)) if os.path.exists(RANGES_JSON) else {}
    if ('int_lo' in d and 'teff_lo' in d and d.get('int_def') == 'own_tau1'
            and not (RECHOOSE and 'int' not in _DONE)):
        return d
    _DONE.add('int')
    fr = ired_newest['fr_own'].astype(np.float64)
    q = fr / wedge_mean(fr, ired_newest['x2f'])
    lo, hi = np.percentile(q, [1, 99])
    d['int_lo'], d['int_hi'] = float(max(lo, 0.0)), float(hi)
    te = (np.maximum(fr, 0.0) / SIGMA_SB) ** 0.25
    lo, hi = np.percentile(te, [1, 99])
    d['teff_lo'], d['teff_hi'] = float(lo), float(hi)
    d['int_def'] = 'own_tau1'
    d['int_from_dump'] = int(ired_newest['n'])
    d['int_t_d'] = float(ired_newest['time'] / DAY)
    json.dump(d, open(RANGES_JSON, 'w'), indent=1)
    return d


def fig_int(ired, rng, out, dpi=100):
    """Left: F_r/<F_r> at each column's own tau_R = 1 cell; right: T_eff = (F_r/sigma)^(1/4)
    [K].  Arrays are (k = phi, j = theta): plotted transposed (theta down, phi across)."""
    x2f, x3f = np.degrees(ired['x2f']), np.degrees(ired['x3f'])
    fr = ired['fr_own'].astype(np.float64)
    fm = wedge_mean(fr, ired['x2f'])
    q = fr / fm
    te = (np.maximum(fr, 0.0) / SIGMA_SB) ** 0.25
    ro = ired['r_own'].astype(np.float64)
    lrat = 4 * np.pi * wedge_mean(ro ** 2 * fr, ired['x2f'])
    r1, r50, r99 = np.percentile(ro / RSUN, [1, 50, 99])
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 6.6), constrained_layout=True)
    for ax, z, (lo, hi), tl, cl in (
            (axs[0], q, (rng['int_lo'], rng['int_hi']), 'Emergent flux (intensity)',
             r'$F_r/\langle F_r\rangle$'),
            (axs[1], te, (rng['teff_lo'], rng['teff_hi']), 'Local effective temperature',
             r'$T_{\rm eff}=(F_r/\sigma_{\rm SB})^{1/4}$ [K]')):
        pc = ax.pcolormesh(x3f, x2f, z.T, cmap='Blues_r', vmin=lo, vmax=hi, shading='flat',
                           rasterized=True)
        ax.set_xlim(0, 60)
        ax.set_ylim(120, 60)
        ax.set_xticks([0, 15, 30, 45, 60])
        ax.set_yticks([60, 75, 90, 105, 120])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_title(tl, fontsize=12)
        cb = fig.colorbar(pc, ax=ax, orientation='horizontal', extend='both', shrink=0.9)
        cb.set_label(cl)
    axs[0].set_ylabel(r'$\theta$ [deg]')
    fig.suptitle(r'$F_r$ at own $\tau_R$=1 (r %.1f-%.1f $R_\odot$)  t = %6.2f d' % (
        r1, r99, ired['time'] / DAY) + '\n' + r'$\langle T_{\rm eff}\rangle$ = %.0f K,'
        r'  $4\pi\langle r^2F_r\rangle/L$ = %.3f' % ((fm / SIGMA_SB) ** 0.25,
                                                     lrat / L_STAR),
        fontsize=15, family='monospace')
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


def fig6(hst, out, t0_d=30.0, t1_d=None):
    """Light curve + PSD over the analysis window t0_d..t1_d (default 30 d to the end
    of the hst, the paper's steady window)."""
    t, L = hst_series(hst, 'L_top')
    tmax_d = t.max() / DAY
    t1 = tmax_d if t1_d is None else t1_d
    w = (t >= t0_d * DAY) & (t <= t1 * DAY + 1.0)
    if w.sum() < 64:                                     # window not reached yet
        w = t <= t1 * DAY + 1.0
    full = t <= t1 * DAY + 1.0                           # t = 0 .. window end
    tF, LF = t[full], L[full]
    t, L = t[w], L[w]
    s = np.sqrt(OMEGA_WEDGE / (4 * np.pi))
    x, flux, fluxm, cm = lc_treatments(t, L)
    fig = plt.figure(figsize=(14.46, 7.96))
    ax = fig.add_axes([0.06, 0.60, 0.90, 0.34])
    td = t / DAY
    # before the window: NOT detrended; 1 + sqrt(Omega/4pi) (L_top / <L_top>_window - 1)
    pre = tF < t.min()
    fpre = 1.0 + s * (LF[pre] / L.mean() - 1.0)
    if pre.any():
        ax.plot(tF[pre] / DAY, fpre, '-', color='#8fb3e8', lw=1.0,
                label='Ours before the window: 1 + $\\sqrt{\\Omega/4\\pi}$'
                '($L_{\\rm top}/\\langle L_{\\rm top}\\rangle_{\\rm window}$ - 1), '
                'not detrended')
    ax.plot(td, flux, 'b-', lw=1.5,
            label='Ours, bolometric, paper eq. 8 (detrended over the window, '
            '$\\sqrt{\\Omega/4\\pi}$ = %.3f)' % s)
    ax.plot(td, fluxm, '-', color='darkorange', lw=1.0, alpha=0.8,
            label='Ours, radial mode (P = %.2f d) unscaled + residual scaled' % cm['P_d'])
    ax.axhspan(0.97, 1.03, color='c', alpha=0.3, lw=0,
               label='SK -67 133 range, TESS band '
               '(read off the paper, $\\pm$3 %)')
    xr = (0, max(30.0, tmax_d))
    ax.axvspan(td.min(), td.max(), color='#ffe9a8', alpha=0.6, lw=0, zorder=0,
               label='analysis window %.1f - %.1f d' % (td.min(), td.max()))
    ax.set_xlim(0, max(tF.max() / DAY * 1.02, 0.1))
    fall = np.concatenate([flux, fpre])                  # auto-expand for the start-up
    ax.set_ylim(min(0.7, fall.min() - 0.02), max(1.2, fall.max() + 0.02))
    ax.set_xlabel('Time since Simulation Starts [days]')
    ax.set_ylabel('Normalized Flux')
    ax.legend(loc='lower left', fontsize=8)
    T = t.max() - t.min()
    info = dict(t_range_d=[round(float(td.min()), 4), round(float(td.max()), 4)],
                n_samples=int(len(td)), flux_min=float(flux.min()),
                flux_max=float(flux.max()), cadence_s=float(np.median(np.diff(t))),
                L_top_mean_wedge_erg_s=float(L.mean()),
                L_top_mean_whole_star_erg_s=float(L.mean() * 4 * np.pi / OMEGA_WEDGE),
                xr=xr, rms_wedge_detrended=float(x.std()),
                rms_eq8=float((flux - 1).std()), rms_mode=float((fluxm - 1).std()),
                mode=dict(P_d=cm['P_d'], f_muHz=cm['f_muHz'], amp_rel=cm['amp'],
                          var_frac=cm['var_frac']),
                ke_flat=ke_flatness(hst))
    ax2 = fig.add_axes([0.06, 0.08, 0.40, 0.40])
    if len(t) >= 8:
        tw = ax2.twinx()
        for fl, col, cc, key in ((flux, 'b', 'navy', 'psd_eq8'),
                                 (fluxm, 'darkorange', 'saddlebrown', 'psd_mode')):
            f, p = psd_ls(t, fl)
            p = p * 1e-12      # ppm^2/muHz -> (dF/F)^2/muHz, the paper's actual unit
            ax2.loglog(f, p, '-', color=col, lw=1.2)
            sel = (f >= 1.25) & (f <= 278)
            c = np.cumsum(p[sel])
            tw.plot(f[sel], c / c[-1], '--', color=cc, lw=2.5)
            info[key] = psd_stats(f, p, cm['f_muHz'], T)
        ax2.axvline(f.max(), color='0.5', ls=':', lw=1)
        tw.set_ylim(0, 1)
        tw.set_ylabel('Cumulative Integrated PSD')
        info['psd_f_range_muHz'] = [float(f.min()), float(f.max())]
        info['f80_muHz'] = info['psd_eq8']['f80_muHz']
    ax2.set_xlim(1.0, 300)
    ax2.set_ylim(7e-9, 8e-4)                               # the paper's panel range
    ax2.set_xlabel(r'Frequency [$\mu$Hz]')
    ax2.set_ylabel(r'PSD [$(\Delta F/F)^2$ / $\mu$Hz]')
    ax2.set_title('Ours, t = %.2f - %.2f d (blue eq. 8, orange mode unscaled)'
                  % (td.min(), td.max()), fontsize=10)
    ax3 = fig.add_axes([0.52, 0.05, 0.46, 0.47])
    ax3.axis('off')
    crop = os.path.join(FIGS, 'tess_s96_crop.png')   # lower-right panel of paper Fig 6
    if os.path.exists(crop):
        ax3.imshow(plt.imread(crop), interpolation='lanczos')
        ax3.text(0.5, -0.01, 'TESS sector 96, from Ma+2026 Fig. 6 (ppm$^2$/$\\mu$Hz)',
                 transform=ax3.transAxes, ha='center', va='top', fontsize=10)
    else:
        ax3.text(0.5, 0.5, 'TESS, SK -67 133:\n'
                 'sector 96 PSD peaks at ~2-4 uHz (Ma+2026)\n'
                 'SLF fit nu_peak = 2.58 uHz, nu_char = 3.14 uHz (Ma+2024)',
                 ha='center', va='center', fontsize=12)
    info['pre_window'] = dict(n=int(pre.sum()), flux_min=float(fpre.min()) if pre.any()
                              else None, flux_max=float(fpre.max()) if pre.any() else None)
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
    ap.add_argument('--rmax-plot', type=float, default=65.0)
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
    if a.snap is None:
        red = reds[-1]
    else:
        red = [c for c in reds if int(c['n']) == a.snap][0]
    i2 = fig2(reds, hst, os.path.join(FIGS, 'out', 'ours_fig2.png'), a.rmax_plot)
    est = (os.path.join(BSG, 'prad', 'est_%05d.npz' % int(red['n'])) if a.prad_est
           else None)
    i3 = fig3(red, a.rmax_plot, kap, os.path.join(FIGS, 'out', 'ours_fig3.png'), est)
    i4 = fig4(red, os.path.join(FIGS, 'out', 'ours_fig4.png'), reds[-1])
    i6 = fig6(hst, os.path.join(FIGS, 'out', 'ours_fig6.png'), a.psd_t0, a.psd_t1)
    nn, hh, mm, _ = [x for x in dumps if x[0] == int(red['n'])][0]
    ired = int_reduce(nn, hh, mm, x1f, kap, red)
    nl = dumps[-1]
    irng = int_ranges(int_reduce(nl[0], nl[1], nl[2], x1f, kap, reds[-1]))
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
