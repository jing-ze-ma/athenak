#!/usr/bin/env python3
"""bsg_fr2d_1015 fr2d.py RUN [RUN ...]: 2-D look at the radial radiation flux F_r near the photosphere.

Per run (one dump): F_r = m1_f1 (LAB-frame radial flux; rad_m1 is mixed-frame, docs/dev/rad_m1_design.md:35),
tau_R per column integrated inward from the outer boundary with the run's own Rosseland table (bilinear in
log T, log rho as radm1::M1TableOpacities, same he_opac_logd_min extension as HsReadOpacityTable), T from the run's
EOS (ideal mu 0.62 or the general-EOS table dump).  r_ph = exact tau_R = 2/3 crossing per column (cell-constant
kappa rho, as the code's tau), r10 likewise.  Outputs PNGs + <run>.json in out/."""
import sys, os, re, json, math
import numpy as np
PAGE = '/raven/ptmp/jinma/bsgrec_1011/ref/docs/handover/bsg-hrdet-1009/page'
sys.path.insert(0, PAGE)
import bin_convert_repro as bc  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy.ndimage import uniform_filter  # noqa: E402

W = '/raven/ptmp/jinma/bsg_fr2d_1015'
D = W + '/data'
OUT = W + '/out'
RSUN = 6.957e10
DAY = 86400.0
SIG = 5.670374419e-5
CL = 2.99792458e10
MU_K = 0.62 * 1.66053906660e-24 / 1.380649e-16

RUNS = {
    'sk67': dict(
        label='Viper BSG sk67 (Sk-67 133; TOPS, general EOS; d6164a6a) dump 116',
        h=D + '/sk67/bsg3d.hydro_w.00116.bin', m=D + '/sk67/bsg3d.m1.00116.bin',
        grid=dict(nx1=205, x1min=1.391400e12, x1max=3.854178e12, p=[(-7.725141242938e-01, 2.208282271879e-01,
                  9.263535133940e-01, 8.227593520313e-03)]),
        eos=D + '/sk67/eos_dump.txt',
        opac='/raven/u/jinma/bsg_raven_1002/bundle/rosseland_tops_x0.7_z0.008.txt', ldmin=-21.0, L=5.558683e38),
    'repro': dict(
        label='Viper BSG repro (Ma+2026 IC; Ma tables, ideal gas; d6164a6a) dump 99',
        h=D + '/repro/bsg3d.hydro_w.00099.bin', m=D + '/repro/bsg3d.m1.00099.bin',
        gridfile=D + '/repro/bsg3d.x1grid.txt',
        grid=dict(nx1=267, x1min=1.391400e12, x1max=5.565600e12, p=[(-6.808500000000e-01, 1.723088670756e-01,
                  7.035565721502e-01, 5.507821105970e-03)]),
        eos=None, opac='/raven/ptmp/jinma/bsg_truerepro_1008/files/rosseland_ma2026_x0.7_z0.008_padded.txt',
        ldmin=None, L=5.581461e38),
    'P': dict(
        label='Raven BSG prod P (Ma IC; FIN3 9972a56d) dump 42',
        h='/raven/ptmp/jinma/bsgprod_1012/prod/P/bin/bsg3d.hydro_w.00042.bin',
        m='/raven/ptmp/jinma/bsgprod_1012/prod/P/bin/bsg3d.m1.00042.bin',
        vet='/raven/ptmp/jinma/bsgprod_1012/prod/P/bin/bsg3d.m1_vet.00042.bin',
        face='/raven/ptmp/jinma/bsgprod_1012/prod/P/bin/bsg3d.m1_face.00042.bin',
        gridfile='/raven/ptmp/jinma/bsgprod_1012/prod/P/bsg3d.x1grid.txt',
        eos=None, opac='/raven/ptmp/jinma/bsg_truerepro_1008/files/rosseland_ma2026_x0.7_z0.008_padded.txt',
        ldmin=None, L=5.581461e38),
    'agcE1': dict(
        label='Viper AG Car E1 (apudev; general EOS X 0.36; d3bf42d7) dump 26',
        h=D + '/agcE1/agcar3d.hydro_w.00026.bin', m=D + '/agcE1/agcar3d.m1.00026.bin',
        gridfile=D + '/agcE1/agcar3d.x1grid.txt',
        eos='/raven/ptmp/jinma/rampA_1011/ana/eos_dump.txt',
        opac=D + '/agcE1/rosseland_ext2_gs98_x0.36_z0.02.txt', ldmin=-21.0,
        L=8.142577965749e12 * 4 * np.pi * 5.7743100000e12 ** 2),
    'agcL64': dict(
        label='Raven AG Car E1 reduced (pulse_1014 L268_64: 268x64x64, prod ramp) dump 8',
        h='/raven/ptmp/jinma/pulse_1014/runs/L268_64/bin/agcar3d.hydro_w.00079.bin',
        m='/raven/ptmp/jinma/pulse_1014/runs/L268_64/bin/agcar3d.m1.00008.bin',
        gridfile='/raven/ptmp/jinma/pulse_1014/runs/L268_64/agcar3d.x1grid.txt',
        eos='/raven/ptmp/jinma/rampA_1011/ana/eos_dump.txt',
        opac=D + '/agcE1/rosseland_ext2_gs98_x0.36_z0.02.txt', ldmin=-21.0,
        L=8.142577965749e12 * 4 * np.pi * 5.7743100000e12 ** 2),
    'agcS': dict(
        label='Viper AG Car E1 reduced post-ramp (rol2u_drift S_B: 268x32x32, w=0) dump 30',
        h='/raven/ptmp/jinma/rol2u_drift_1015/vruns/S_B/bin/agcar3d.hydro_w.00300.bin',
        m='/raven/ptmp/jinma/rol2u_drift_1015/vruns/S_B/bin/agcar3d.m1.00030.bin',
        gridfile='/raven/ptmp/jinma/rol2u_drift_1015/vruns/S_B/agcar3d.x1grid.txt',
        eos='/raven/ptmp/jinma/rampA_1011/ana/eos_dump.txt',
        opac=D + '/agcE1/rosseland_ext2_gs98_x0.36_z0.02.txt', ldmin=-21.0,
        L=8.142577965749e12 * 4 * np.pi * 5.7743100000e12 ** 2),
}

for _k in ('agcL64', 'agcS', 'agcE1'):
    RUNS[_k + '_idealT'] = dict(RUNS[_k], eos=None, label=RUNS[_k]['label'] + ' [WHAT-IF: ideal-gas T, mu 0.62]')


# ---------------------------------------------------------------------------------------------- grid
def logcosh(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2 * ax)) - math.log(2.0)


def plat_G(xa, xb, w, xi):
    return 0.5 * w * (logcosh((xi - xa) / w) - logcosh(-xa / w) - logcosh((xi - xb) / w) + logcosh(-xb / w))


def stretch_faces(nx1, x1min, x1max, c=(), p=()):
    """grid_stretch.hpp StretchRPoly on the uniform faces"""
    xi = np.arange(nx1 + 1) / nx1
    u = xi.copy()
    xik = xi.copy()
    for ck in c:
        u += ck * xik * (1 - xi)
        xik = xik * xi
    for a, xa, xb, w in p:
        u += a * (plat_G(xa, xb, w, xi) - xi * plat_G(xa, xb, w, 1.0))
    return x1min + (x1max - x1min) * u


def faces_from_file(f):
    m = re.search(r'is = (\d+) ie = (\d+)', open(f).readline())
    is_, ie = int(m.group(1)), int(m.group(2))
    a = np.loadtxt(f, usecols=(0, 1))
    return a[is_:ie + 2, 1]


def centroid(rl, rr):
    q = rl / rr
    return 0.25 * (q * q + 1) / ((q * q + q + 1) / 3) * (rr + rl)


# ---------------------------------------------------------------------------------------------- opacity / EOS
class Kappa:
    """radm1::M1TableOpacities Rosseland part: bilinear in (log T, log rho) of log10 kappa, clamped at the
    grid edges; optional he_opac_logd_min extension (HsReadOpacityTable: slope of the first rho interval
    clamped to [0, 1])."""
    def __init__(self, path, ldmin=None):
        g = None
        for ln in open(path):
            if not ln.startswith('#'):
                break
            mm = re.match(r'#\s*(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', ln)
            if mm and g is None:
                g = mm.groups()
        nT, nD = int(g[0]), int(g[1])
        lT0, dlT, lD0, dlD = map(float, g[2:])
        tab = np.loadtxt(path, comments='#').reshape(nT, nD)
        if ldmin is not None and ldmin < lD0 - 1e-9 * dlD:
            nadd = int(math.ceil((lD0 - ldmin) / dlD - 1e-9))
            sl = np.clip((tab[:, 1] - tab[:, 0]) / dlD, 0.0, 1.0)
            ext = tab[:, :1] - sl[:, None] * (nadd - np.arange(nadd))[None, :] * dlD
            tab = np.concatenate([ext, tab], axis=1)
            lD0 -= nadd * dlD
            nD += nadd
        self.nT, self.nD, self.lT0, self.dlT, self.lD0, self.dlD, self.tab = nT, nD, lT0, dlT, lD0, dlD, tab

    def __call__(self, rho, T):
        x = np.clip((np.log10(T) - self.lT0) / self.dlT, 0, self.nT - 1 - 1e-9)
        y = np.clip((np.log10(rho) - self.lD0) / self.dlD, 0, self.nD - 1 - 1e-9)
        i = x.astype(np.int32)
        j = y.astype(np.int32)
        fx = x - i
        fy = y - j
        t = self.tab
        v = ((1 - fx) * ((1 - fy) * t[i, j] + fy * t[i, j + 1]) + fx * ((1 - fy) * t[i + 1, j] + fy * t[i + 1, j + 1]))
        return 10.0 ** v


class GEOS:
    """general-EOS table dump (x = log rho, y = log T; cols log10 e/rho, log10 p/rho): invert e -> T"""
    def __init__(self, path):
        hdr = [ln for ln in open(path) if ln.startswith('#')][:6]
        nx, ny, xmin, dx, ymin, dy = hdr[2][1:].split()
        self.nx, self.ny = int(nx), int(ny)
        self.xmin, self.dx, self.ymin, self.dy = map(float, (xmin, dx, ymin, dy))
        a = np.loadtxt(path, comments='#').reshape(self.ny, self.nx, 2)
        self.le, self.lp = a[..., 0], a[..., 1]

    def T_P(self, rho, e):
        x = np.clip((np.log10(rho) - self.xmin) / self.dx, 0, self.nx - 1 - 1e-9)
        ix = x.astype(np.int32)
        wx = x - ix
        le = np.log10(np.maximum(e / rho, 1e-30))
        f = lambda j: (1 - wx) * self.le[j, ix] + wx * self.le[j, ix + 1]  # noqa: E731
        lo = np.zeros(rho.shape, np.int32)
        hi = np.full(rho.shape, self.ny - 2, np.int32)
        for _ in range(12):
            mid = (lo + hi + 1) // 2
            c = f(mid) <= le
            lo = np.where(c, mid, lo)
            hi = np.where(c, hi, mid - 1)
        j = lo
        f0, f1 = f(j), f(j + 1)
        fy = np.clip((le - f0) / (f1 - f0), 0.0, 1.0)
        self.nlow = int((le < f(np.zeros_like(j))).sum())
        lT = self.ymin + (j + fy) * self.dy
        lp = (1 - fy) * ((1 - wx) * self.lp[j, ix] + wx * self.lp[j, ix + 1]) + \
            fy * ((1 - wx) * self.lp[j + 1, ix] + wx * self.lp[j + 1, ix + 1])
        return 10.0 ** lT, rho * 10.0 ** lp


def assemble(d):
    geo = np.asarray(d['mb_geometry'])
    nx1, nx2, nx3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_mb'], d['nx2_mb'], d['nx3_mb']
    out = {}
    for v in d['var_names']:
        a = np.empty((nx3, nx2, nx1), np.float64)
        for m, blk in enumerate(d['mb_data'][v]):
            i0 = int(round((geo[m, 0] - d['x1min']) / (d['x1max'] - d['x1min']) * nx1))
            j0 = int(round((geo[m, 2] - d['x2min']) / (d['x2max'] - d['x2min']) * nx2))
            k0 = int(round((geo[m, 4] - d['x3min']) / (d['x3max'] - d['x3min']) * nx3))
            a[k0:k0 + b3, j0:j0 + b2, i0:i0 + b1] = blk
        out[v] = a
    return out, d


# ---------------------------------------------------------------------------------------------- helpers
def surf_r(tauf, x1f, kr, t0):
    """exact r where the inward-integrated tau (faces tauf[...,i] at x1f[i]; cell-constant kr = kappa rho)
    first reaches t0 from the top; nan if never (column optically thin)"""
    nx1 = kr.shape[-1]
    ok = tauf[..., 0] >= t0
    # last face index with tau >= t0 (tau decreases outward)
    ge = tauf[..., :nx1] >= t0
    i = nx1 - 1 - np.argmax(ge[..., ::-1], axis=-1)  # cell i: tauf[i] >= t0 > tauf[i+1]
    ta = np.take_along_axis(tauf, (i + 1)[..., None], -1)[..., 0]
    k = np.take_along_axis(kr, i[..., None], -1)[..., 0]
    r = x1f[i + 1] - (t0 - ta) / k
    return np.where(ok, r, np.nan), i


def at_r(Q, x1v, r):
    """linear interpolation of cell-centred Q (k, j, i) in r to r (k, j) per column"""
    ic = np.clip(np.searchsorted(x1v, np.nan_to_num(r, nan=x1v[0])) - 1, 0, len(x1v) - 2)
    w = np.clip((np.nan_to_num(r, nan=x1v[0]) - x1v[ic]) / (x1v[ic + 1] - x1v[ic]), 0.0, 1.0)
    q0 = np.take_along_axis(Q, ic[..., None], -1)[..., 0]
    q1 = np.take_along_axis(Q, (ic + 1)[..., None], -1)[..., 0]
    return np.where(np.isfinite(r), (1 - w) * q0 + w * q1, np.nan)


def wmean(z, dom):
    m = np.isfinite(z)
    return float((z * dom)[m].sum() / dom[m].sum())


def wpct(z, dom, ps):
    m = np.isfinite(z)
    zz, ww = z[m].ravel(), dom[m].ravel()
    o = np.argsort(zz)
    c = np.cumsum(ww[o])
    c /= c[-1]
    return [float(np.interp(p / 100.0, c, zz[o])) for p in ps]


def box3(z):
    """one 3x3 box pass, periodic in both angles (the wedge is periodic in theta and phi)"""
    return uniform_filter(z, size=3, mode='wrap')


# ---------------------------------------------------------------------------------------------- main per run
def run(name):
    R = RUNS[name]
    os.makedirs(OUT, exist_ok=True)
    g = R.get('grid')
    if 'gridfile' in R:
        x1f = faces_from_file(R['gridfile'])
        if g:
            xs = stretch_faces(g['nx1'], g['x1min'], g['x1max'], g.get('c', ()), g['p'])
            print('%s: stretch formula vs grid file max |dr|/dr = %.2e' % (name, np.max(np.abs(xs - x1f) / np.diff(x1f).min())))
    else:
        x1f = stretch_faces(g['nx1'], g['x1min'], g['x1max'], g.get('c', ()), g['p'])
    x1v = centroid(x1f[:-1], x1f[1:])
    dr = np.diff(x1f)
    w, dh = assemble(bc.read_binary(R['h']))
    t = dh['time']
    rho, eint = w['dens'], w['eint']
    vr = w['velx']
    v2 = vr ** 2 + w['vely'] ** 2 + w['velz'] ** 2
    del w
    nx3, nx2, nx1 = rho.shape
    assert nx1 == len(x1v), (nx1, len(x1v))
    x2f = np.linspace(dh['x2min'], dh['x2max'], nx2 + 1)
    x3f = np.linspace(dh['x3min'], dh['x3max'], nx3 + 1)
    dom = np.diff(x3f)[:, None] * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :]
    if R['eos'] is None:
        pg = (5.0 / 3.0 - 1.0) * eint
        T = pg / rho * MU_K
        nlow = 0
    else:
        ge = GEOS(R['eos'])
        T, pg = ge.T_P(rho, eint)
        nlow = ge.nlow
    kap = Kappa(R['opac'], R['ldmin'])
    kr = kap(rho, T) * rho
    dtau = kr * dr[None, None, :]
    tauf = np.zeros((nx3, nx2, nx1 + 1))
    tauf[..., :nx1] = np.flip(np.cumsum(np.flip(dtau, -1), -1), -1)
    tauc = tauf[..., 1:] + 0.5 * dtau
    rph, iph = surf_r(tauf, x1f, kr, 2.0 / 3.0)
    r1, _ = surf_r(tauf, x1f, kr, 1.0)
    r10, _ = surf_r(tauf, x1f, kr, 10.0)
    m1, dm = assemble(bc.read_binary(R['m']))
    assert abs(dm['time'] - t) <= 1e-6 * max(1.0, t), (dm['time'], t)
    E, Fr = m1['m1_e'], m1['m1_f1']
    del m1
    L4 = R['L'] / (4 * np.pi)
    D11 = None
    if R.get('vet'):
        vv, _ = assemble(bc.read_binary(R['vet']))
        D11 = vv['m1_d11']
        del vv
    Prr = E * (D11 if D11 is not None else 1.0 / 3.0)
    F0r = Fr - vr * (E + Prr)   # comoving radial flux, O(v/c): F0_r = F_r - v_r E - v_j P_jr (P_rr part only)
    q = Fr * x1v[None, None, :] ** 2 / L4
    # ---------------- photosphere numbers
    Fph = at_r(Fr, x1v, rph)
    F0ph = at_r(F0r, x1v, rph)
    qph = Fph * rph ** 2 / L4
    q0ph = F0ph * rph ** 2 / L4
    Tph = at_r(T, x1v, rph)
    qbox = box3(np.nan_to_num(qph, nan=np.nanmedian(qph)))
    # own-tau=1 cell (the page's int_reduce definition): innermost cell with tau(centre) < 1
    thin = tauc < 1.0
    iown = np.where(thin.any(-1), np.argmax(thin, -1), nx1 - 1)
    Fown = np.take_along_axis(Fr, iown[..., None], -1)[..., 0]
    rown = x1v[iown]
    qown = Fown * rown ** 2 / L4
    # top-of-grid flux scaled to r_ph (escaping-flux definition)
    qtop = Fr[..., -1] * x1v[-1] ** 2 / L4
    # tau-window mean of r^2 F_r over tau in [0.1, 2/3]  (radial smoothing along the column)
    win = (tauc >= 0.1) & (tauc <= 2.0 / 3.0)
    qwin = np.where(win.any(-1), (q * win * dr).sum(-1) / np.maximum((win * dr).sum(-1), 1e-300), np.nan)
    res = dict(run=name, label=R['label'], t_s=t, t_d=t / DAY, L=R['L'], nx=[nx1, nx2, nx3],
               eos_cells_below_table=nlow,
               rph_Rsun_mean=wmean(rph, dom) / RSUN, rph_Rsun_p2_50_98=[x / RSUN for x in wpct(rph, dom, [2, 50, 98])],
               r10_Rsun_p2_50_98=[x / RSUN for x in wpct(r10, dom, [2, 50, 98])],
               rtop_Rsun=x1f[-1] / RSUN, thin_columns=int(np.isnan(rph).sum()),
               frac_Fph_le0=float(((Fph <= 0) * dom)[np.isfinite(Fph)].sum() / dom[np.isfinite(Fph)].sum()),
               frac_Fph_le0_cols=float((Fph <= 0).mean()),
               q_ph_p2_50_98=wpct(qph, dom, [2, 50, 98]), q_ph_mean=wmean(qph, dom),
               q_ph_sd=float(np.nanstd(qph)), q_ph_min=float(np.nanmin(qph)), q_ph_max=float(np.nanmax(qph)),
               q0_ph_comoving_p2_50_98=wpct(q0ph, dom, [2, 50, 98]), frac_F0ph_le0=float((F0ph <= 0).mean()),
               q_box3_p2_50_98=wpct(qbox, dom, [2, 50, 98]), frac_box3_le0=float((qbox <= 0).mean()),
               q_own_p2_50_98=wpct(qown, dom, [2, 50, 98]), frac_own_le0=float((qown <= 0).mean()),
               q_top_p2_50_98=wpct(qtop, dom, [2, 50, 98]), frac_top_le0=float((qtop <= 0).mean()),
               q_win_p2_50_98=wpct(qwin, dom, [2, 50, 98]), frac_win_le0=float((qwin <= 0).mean()),
               Tph_K_p2_50_98=wpct(Tph, dom, [2, 50, 98]),
               Teff_local_p2_50_98=wpct((np.maximum(Fph, 0) / SIG) ** 0.25, dom, [2, 50, 98]),
               Teff_mean_flux=float((wmean(qph, dom) * L4 / (wmean(rph, dom) ** 2) / SIG) ** 0.25),
               L_ph_over_L=wmean(qph, dom), L_top_over_L=wmean(qtop, dom))
    # ---------------- radial structure: sign changes, flux decomposition (shell means, area weighted)
    om = dom.sum()
    shm = lambda a: np.einsum('kji,kj->i', a, dom) / om  # noqa: E731
    neg = shm((Fr <= 0).astype(np.float64))
    h = eint + pg
    vm, hm = shm(vr), shm(h)
    r2L = x1v ** 2 / L4
    Fh = (shm(h * vr) - hm * vm) * r2L        # enthalpy fluctuation flux = F_conv
    Fht = shm(h * vr) * r2L
    Fk = shm(0.5 * rho * v2 * vr) * r2L
    Frad = shm(Fr) * r2L
    Fadv = shm(vr * (E + Prr)) * r2L
    Frad0 = shm(F0r) * r2L
    taush = shm(tauc)
    idx = np.where(neg > 0)[0]
    idx1 = np.where(neg > 0.01)[0]
    res['sign_change_r_Rsun_any'] = [float(x1v[idx[0]] / RSUN), float(x1v[idx[-1]] / RSUN)] if len(idx) else None
    res['sign_change_r_Rsun_gt1pct'] = [float(x1v[idx1[0]] / RSUN), float(x1v[idx1[-1]] / RSUN)] if len(idx1) else None
    res['neg_frac_max'] = [float(neg.max()), float(x1v[np.argmax(neg)] / RSUN)]

    def intervals(mask):
        out, i = [], 0
        while i < len(mask):
            if mask[i]:
                j = i
                while j + 1 < len(mask) and mask[j + 1]:
                    j += 1
                out.append([round(float(x1v[i] / RSUN), 2), round(float(x1v[j] / RSUN), 2), round(float(neg[i:j + 1].max()), 4)])
                i = j + 1
            else:
                i += 1
        return out
    res['neg_intervals_gt1pct'] = intervals(neg > 0.01)
    res['neg_intervals_gt0'] = intervals(neg > 0)
    res['neg_profile_Rsun'] = [[round(float(x1v[i] / RSUN), 2), round(float(neg[i]), 4), float(taush[i])]
                               for i in range(0, nx1, max(1, nx1 // 60))]
    # in local-tau coordinates
    tb = np.array([1e-3, 1e-2, 0.1, 0.3, 2 / 3., 1, 3, 10, 30, 100, 1e3, 1e4, 1e5, 1e7])
    ng = []
    for a, b in zip(tb[:-1], tb[1:]):
        mk = (tauc >= a) & (tauc < b)
        ng.append([float(a), float(b), int(mk.sum()), float((Fr[mk] <= 0).mean()) if mk.any() else None])
    res['neg_frac_by_local_tau'] = ng
    # flux decomposition on tau surfaces (per column at r(tau0), fluctuation over the surface)
    surf = {}
    for t0 in (0.1, 2.0 / 3.0, 1.0, 3.0, 10.0, 100.0):
        rs, _ = surf_r(tauf, x1f, kr, t0)
        hs, vs, fs = at_r(h, x1v, rs), at_r(vr, x1v, rs), at_r(Fr, x1v, rs)
        ks = at_r(0.5 * rho * v2 * vr, x1v, rs)
        es = at_r(vr * (E + Prr), x1v, rs)
        sc = rs ** 2 / L4
        mh, mv = wmean(hs, dom), wmean(vs, dom)
        surf['%.3g' % t0] = dict(r_Rsun=wmean(rs, dom) / RSUN, Frad=wmean(fs * sc, dom),
                                 Fconv=wmean((hs - mh) * (vs - mv) * sc, dom), Fenth_tot=wmean(hs * vs * sc, dom),
                                 Fkin=wmean(ks * sc, dom), Fradadv=wmean(es * sc, dom),
                                 vr_rms_kms=float(np.sqrt(wmean(vs ** 2, dom)) / 1e5),
                                 Fconv_over_Frad=wmean((hs - mh) * (vs - mv) * sc, dom) / wmean(fs * sc, dom))
    res['surf'] = surf
    iph_m = int(np.searchsorted(x1v, res['rph_Rsun_mean'] * RSUN))
    res['shell_at_rph'] = dict(r_Rsun=float(x1v[iph_m] / RSUN), Frad=float(Frad[iph_m]), Fconv=float(Fh[iph_m]),
                               Fkin=float(Fk[iph_m]), Fradadv=float(Fadv[iph_m]), neg_frac=float(neg[iph_m]))
    res['shell_Lrad_over_L_minmax_above_r10'] = None
    np.savez_compressed(os.path.join(OUT, name + '_prof.npz'), r=x1v, neg=neg, Fh=Fh, Fht=Fht, Fk=Fk, Frad=Frad, Fadv=Fadv,
                        Frad0=Frad0, taush=taush, qph=qph, rph=rph, r10=r10, x2f=x2f, x3f=x3f)
    json.dump(res, open(os.path.join(OUT, name + '.json'), 'w'), indent=1)
    print(json.dumps({k: res[k] for k in res if k not in ('surf', 'neg_frac_by_local_tau')}, indent=1))
    print(json.dumps(res['surf'], indent=1))
    # ================================================================ figures
    lab = '%s   t = %.2f d' % (R['label'], t / DAY)
    th = np.degrees(x2f)
    ph = np.degrees(x3f)
    rr = x1f / RSUN
    k0 = nx3 // 2
    vmax = 2.0
    # ---- Fig 1: (r, theta) slice at one phi
    fig, axs = plt.subplots(2, 1, figsize=(13, 10.5), constrained_layout=True)
    rlo = max(x1f[0], np.nanmin(r10[k0]) - 0.25 * (x1f[-1] - np.nanmin(r10[k0]))) / RSUN
    for ax, (a, b), ttl in ((axs[0], (rr[0], rr[-1]), 'full radial range'),
                            (axs[1], (rlo, rr[-1]), 'zoom: photosphere region')):
        pc = ax.pcolormesh(rr, th, q[k0], cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='flat', rasterized=True)
        tc = 0.5 * (th[1:] + th[:-1])
        ax.plot(rph[k0] / RSUN, tc, 'k-', lw=1.3, label=r'$\tau_R=2/3$')
        ax.plot(r10[k0] / RSUN, tc, '-', color='lime', lw=1.3, label=r'$\tau_R=10$')
        ax.set_xlim(a, b)
        ax.set_ylim(th[-1], th[0])
        ax.set_xlabel(r'$r\ [R_\odot]$')
        ax.set_ylabel(r'$\theta$ [deg]')
        ax.set_title(ttl)
        ax.legend(loc='lower left', fontsize=10, framealpha=0.8)
        cb = fig.colorbar(pc, ax=ax, extend='both', pad=0.01)
        cb.set_label(r'$F_r r^2/(L_*/4\pi)$  (lab frame, linear)')
    fig.suptitle(lab + r'   slice $\phi$ = %.1f deg' % (0.5 * (ph[k0] + ph[k0 + 1])), fontsize=13)
    fig.savefig(os.path.join(OUT, name + '_rtheta.png'), dpi=110)
    plt.close(fig)
    # ---- Fig 2: (theta, phi) at r_ph
    teff = (np.maximum(Fph, 0) / SIG) ** 0.25
    teffb = (np.maximum(qbox, 0) * L4 / rph ** 2 / SIG) ** 0.25
    tlo, thi = np.nanpercentile(np.where(Fph > 0, teff, np.nan), [1, 99])
    tlo = min(tlo, np.nanpercentile(Tph, 1))
    thi = max(thi, np.nanpercentile(Tph, 99))
    panels = [(qph, 'RdBu_r', (-vmax, vmax), r'$F_r r^2/(L/4\pi)$ at $\tau_R=2/3$ (no smoothing)'),
              (qbox, 'RdBu_r', (-vmax, vmax), r'same, one 3x3 box pass'),
              (rph / RSUN, 'viridis', tuple(np.nanpercentile(rph / RSUN, [1, 99])), r'$r_{\rm ph}\ [R_\odot]$ ($\tau_R=2/3$)'),
              (teff, 'inferno', (tlo, thi), r'$T_{\rm eff}=(\max(F_r,0)/\sigma)^{1/4}$ at $r_{\rm ph}$ [K]'),
              (teffb, 'inferno', (tlo, thi), r'$T_{\rm eff}$ from the 3x3 box pass [K]'),
              (Tph, 'inferno', (tlo, thi), r'$T_{\rm gas}(\tau_R=2/3)$ [K]')]
    fig, axs = plt.subplots(2, 3, figsize=(16, 11.5), constrained_layout=True)
    for ax, (z, cm, (lo, hi), ttl) in zip(axs.ravel(), panels):
        pc = ax.pcolormesh(ph, th, z.T, cmap=cm, vmin=lo, vmax=hi, shading='flat', rasterized=True)
        if cm == 'RdBu_r':
            ax.contour(0.5 * (ph[1:] + ph[:-1]), 0.5 * (th[1:] + th[:-1]), z.T, levels=[0.0], colors='k', linewidths=0.6)
        ax.set_xlim(ph[0], ph[-1])
        ax.set_ylim(th[-1], th[0])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_ylabel(r'$\theta$ [deg]')
        ax.set_title(ttl, fontsize=11)
        fig.colorbar(pc, ax=ax, orientation='horizontal', extend='both', shrink=0.9, pad=0.02)
    fig.suptitle(lab + '\n' + r'$F_r\leq 0$ at $r_{\rm ph}$: %.2f%% of columns;  p2/50/98 of $F_r r^2/(L/4\pi)$ = %.3f / %.3f / %.3f;'
                 r'  $\langle r_{\rm ph}\rangle$ = %.2f $R_\odot$ (black: $F_r=0$)' % (
                     100 * res['frac_Fph_le0_cols'], *res['q_ph_p2_50_98'], res['rph_Rsun_mean']), fontsize=12)
    fig.savefig(os.path.join(OUT, name + '_thetaphi.png'), dpi=110)
    plt.close(fig)
    # ---- Fig 3: radial profiles
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)
    ax = axs[0]
    ax.plot(x1v / RSUN, Frad, 'k', lw=2, label=r'$\langle F_r\rangle$ (lab)')
    ax.plot(x1v / RSUN, Frad0, 'k:', lw=1.5, label=r'$\langle F_{0,r}\rangle$ (comoving, $P_{rr}$ term)')
    ax.plot(x1v / RSUN, Fadv, 'b', label=r'$\langle v_r(E+P_{rr})\rangle$ (rad. advection)')
    ax.plot(x1v / RSUN, Fh, 'r', lw=2, label=r"$F_{\rm conv}=\langle h'v_r'\rangle$, $h=e+P_g$")
    ax.plot(x1v / RSUN, Fk, 'g', label=r'$\langle \frac{1}{2}\rho v^2 v_r\rangle$')
    ax.plot(x1v / RSUN, Frad + Fh + Fk, color='0.5', ls='--', label='rad + conv + kin')
    ax.axhline(0, color='0.7', lw=0.8)
    ax.axhline(1, color='0.7', lw=0.8, ls=':')
    for xx, cc in ((res['rph_Rsun_mean'], 'k'), (res['r10_Rsun_p2_50_98'][1], 'lime')):
        ax.axvline(xx, color=cc, lw=1)
    ax.set_xlabel(r'$r\ [R_\odot]$')
    ax.set_ylabel(r'shell flux $\times r^2/(L/4\pi)$')
    ax.set_title(r'shell means (black line: $\langle r_{\rm ph}\rangle$, green: median $r(\tau_R=10)$)', fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[1]
    ax.semilogy(x1v / RSUN, np.maximum(neg, 1e-6), 'k', lw=2, label=r'shell fraction of cells with $F_r\leq 0$')
    ax.semilogy(x1v / RSUN, taush, 'b', lw=1, label=r'shell mean $\tau_R$ (cell centre)')
    ax.axhline(2 / 3., color='b', ls=':', lw=0.8)
    ax.axvline(res['rph_Rsun_mean'], color='k', lw=1)
    ax.axvline(res['r10_Rsun_p2_50_98'][1], color='lime', lw=1)
    ax.set_ylim(1e-6, 1e6)
    ax.set_xlabel(r'$r\ [R_\odot]$')
    ax.legend(fontsize=9)
    fig.suptitle(lab, fontsize=12)
    fig.savefig(os.path.join(OUT, name + '_prof.png'), dpi=110)
    plt.close(fig)
    return res


if __name__ == '__main__':
    for n in sys.argv[1:]:
        run(n)
