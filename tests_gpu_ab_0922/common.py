"""Shared readers for the 2026-09-22 apudev A/B smoke analysis (read-only)."""
import os
import glob
import numpy as np

B = '/viper/u2/jinma/ATHENAK/bench/'
OUT = '/viper/u2/jinma/ATHENAK/athenak/tests_gpu_ab_0922/'
EOSTAB = B + 'cs_ens/analysis/eos_table.txt'
PROT = 3.0501e5          # 2 pi / Omega, Omega = 2.06e-5
CPOLY = [-0.068392, -2.191487, 2.464818, -1.366698]
CFL = 0.3
DFLOOR = 5.0e-14
PFLOOR = 1.0e-3

# comparison 1: restart arms (rot 283)
ARMS1 = [('prodbin', B + 'ck_sph_ab/prodbin'),
         ('sph',     B + 'ck_sph_ab/sph'),
         ('sphbeam', B + 'ck_sph_ab/sphbeam'),
         ('noang',   B + 'cs_noang/noang'),
         ('offfix',  B + 'ck_sph_ab/offfix')]
# comparison 2: cold-start reconstruction arms
ARMS2 = [('plm_ctl', B + 'cs_recon_ab/plm_ctl'),
         ('ppmx',    B + 'cs_recon_ab/ppmx'),
         ('wenoz',   B + 'cs_recon_ab/wenoz')]
PROD3 = B + 'cs_mhd_prod3'

COL = {'time': 0, 'dt': 1, 'mass': 2, 'm1': 3, 'm2': 4, 'm3': 5, 'E': 6,
       'ke1': 7, 'ke2': 8, 'ke3': 9, 'me1': 10, 'me2': 11, 'me3': 12}


def hst(d):
    a = np.loadtxt(os.path.join(d, 'dhj.mhd.hst'))
    return np.atleast_2d(a)


def evlog(d):
    """event-counter log -> (names, cycles, counts) ; rows are PER INTERVAL."""
    fn = os.path.join(d, 'dhj.log')
    with open(fn) as f:
        lines = f.readlines()
    names = lines[1].replace('#', '').split()
    rows = [[float(v) for v in ln.split()] for ln in lines[2:]
            if ln.strip() and not ln.lstrip().startswith('#')]
    a = np.array(rows)
    return names, a[:, 0], a


def cyclelines(d):
    """stdout elapsed/cycle/time/dt lines -> array (n,4)."""
    fn = glob.glob(os.path.join(d, 'log.out.*'))[0]
    rows = []
    with open(fn) as f:
        for ln in f:
            if ln.startswith('elapsed='):
                p = dict(kv.split('=') for kv in ln.split())
                rows.append([float(p['elapsed']), float(p['cycle']),
                             float(p['time']), float(p['dt'])])
    return np.array(rows), fn


def desum(d):
    """rt_desum lines -> (cycle, src, de, rel) arrays, all lines."""
    fn = glob.glob(os.path.join(d, 'log.out.*'))[0]
    out = []
    with open(fn) as f:
        for ln in f:
            if 'rt_desum' in ln:
                t = ln.split()
                c = float(t[2].split('=')[1])
                src = float(t[5])
                de = float(t[8])
                rel = float(t[11])
                out.append([c, src, de, rel])
    return np.array(out) if out else np.zeros((0, 4))


# ---------------- cubed-sphere geometry, exactly as coordinates.cpp ----------------
def stretch_r(x, r0, r1):
    xi = (x - r0) / (r1 - r0)
    u = xi.copy()
    xik = xi.copy()
    for k in range(1, 5):
        u = u + CPOLY[k - 1] * xik * (1 - xi)
        xik = xik * xi
    return r0 + (r1 - r0) * u


def radial_centroid(rl, rr):
    return 0.75 * (rr**4 - rl**4) / (rr**3 - rl**3)


def cs_metric(raw):
    """Reproduce Coordinates dx1/dx2/dx3, sin_cell, r_c for a cubed-sphere dump.

    Returns dict with rc (ni,), rf (ni+1,), dx1 (ni,), dx2 (nmb,nk,nj,ni),
    dx3 (same), sin_cell (nmb,nk,nj), lat/lon (nmb,nk,nj).
    """
    g = np.asarray(raw['mb_geometry'])
    nmb = g.shape[0]
    ni, nj, nk = raw['nx1_mb'], raw['nx2_mb'], raw['nx3_mb']
    r0, r1 = raw['x1min'], raw['x1max']
    rf = stretch_r(g[0, 0] + (g[0, 1] - g[0, 0]) * np.arange(ni + 1) / ni, r0, r1)
    rl, rr = rf[:-1], rf[1:]
    rc = radial_centroid(rl, rr)
    dx1 = rr - rl
    sin_cell = np.zeros((nmb, nk, nj))
    dth_xi = np.zeros((nmb, nk, nj))
    dth_eta = np.zeros((nmb, nk, nj))
    for m in range(nmb):
        xic = np.pi / 4 * (g[m, 2] + (g[m, 3] - g[m, 2]) * (np.arange(nj) + 0.5) / nj)
        xif = np.pi / 4 * (g[m, 2] + (g[m, 3] - g[m, 2]) * np.arange(nj + 1) / nj)
        etac = np.pi / 4 * (g[m, 4] + (g[m, 5] - g[m, 4]) * (np.arange(nk) + 0.5) / nk)
        etaf = np.pi / 4 * (g[m, 4] + (g[m, 5] - g[m, 4]) * np.arange(nk + 1) / nk)
        X, Y = np.meshgrid(np.tan(xic), np.tan(etac))         # (nk,nj)
        XL, YL = np.meshgrid(np.tan(xif[:-1]), np.tan(etaf[:-1]))
        XR, YR = np.meshgrid(np.tan(xif[1:]), np.tan(etaf[1:]))
        C = np.sqrt(1 + X * X)
        D = np.sqrt(1 + Y * Y)
        sin_cell[m] = np.sqrt(1 + X * X + Y * Y) / (C * D)
        dth_xi[m] = np.arccos((1 + XL * XR + Y * Y)
                              / np.sqrt(1 + XL * XL + Y * Y)
                              / np.sqrt(1 + XR * XR + Y * Y))
        dth_eta[m] = np.arccos((1 + X * X + YL * YR)
                               / np.sqrt(1 + X * X + YL * YL)
                               / np.sqrt(1 + X * X + YR * YR))
    dx2 = rc[None, None, None, :] * dth_xi[..., None]
    dx3 = rc[None, None, None, :] * dth_eta[..., None]
    return dict(rc=rc, rf=rf, dx1=dx1, dx2=dx2, dx3=dx3, sin_cell=sin_cell,
                dth_xi=dth_xi, dth_eta=dth_eta, nmb=nmb, ni=ni, nj=nj, nk=nk)


# ---------------- EOS: T, p, Gamma_1 from the run's own table ----------------
class EOS2:
    """dhjcs.EOS plus Gamma_1, from the tabulated log10(e/rho), log10(p/rho)."""

    def __init__(self, fn=EOSTAB):
        with open(fn) as f:
            lines = [next(f) for _ in range(6)]
        nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2][1:].split()]
        self.nx, self.ny = int(nx), int(ny)
        self.xmin, self.dx, self.ymin, self.dy = xmin, dx, ymin, dy
        d = np.loadtxt(fn)
        self.le = d[:, 0].reshape(self.ny, self.nx)
        self.lp = d[:, 1].reshape(self.ny, self.nx)
        self.lt = self.ymin + self.dy * np.arange(self.ny)
        self.ld = self.xmin + self.dx * np.arange(self.nx)

    def invert(self, rho, eint):
        """rho, eint (erg/cm^3) -> T [K], p [erg/cm^3], Gamma_1 [-]."""
        rho = np.asarray(rho, float)
        eint = np.asarray(eint, float)
        sh = rho.shape
        rho = rho.ravel()
        eint = eint.ravel()
        x = (np.log10(rho) - self.xmin) / self.dx
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        fx = np.clip(x - ix, 0, 1)
        le = self.le[:, ix] * (1 - fx) + self.le[:, ix + 1] * fx
        lp = self.lp[:, ix] * (1 - fx) + self.lp[:, ix + 1] * fx
        target = np.log10(eint / rho)
        fin = np.isfinite(le)
        lef = np.where(fin, le, -np.inf)
        first = np.argmax(fin, axis=0)
        cnt = (lef < target[None, :]).sum(axis=0)
        k = np.clip(cnt, first + 1, self.ny - 1)
        ar = np.arange(len(k))
        e0 = le[k - 1, ar]
        e1 = le[k, ar]
        f = np.clip((target - e0) / (e1 - e0), 0, 1)
        lT = self.lt[k - 1] + f * self.dy
        lP = lp[k - 1, ar] * (1 - f) + lp[k, ar] * f
        T = 10**lT
        p = 10**lP * rho
        # Gamma_1 = (dlnP/dlnrho)_T + (dlnP/dlnT)_rho * (dlnT/dlnrho)_s
        # with (dlnT/dlnrho)_s from de = (p/rho^2) drho at constant entropy.
        iy = np.clip(k - 1, 1, self.ny - 2)
        jx = np.clip(ix, 1, self.nx - 2)
        LP = self.lp + self.ld[None, :]            # log10 p
        LU = self.le                               # log10 (e/rho) = log10 u
        dlnp_dlnr = (LP[iy, jx + 1] - LP[iy, jx - 1]) / (2 * self.dx)
        dlnp_dlnT = (LP[iy + 1, jx] - LP[iy - 1, jx]) / (2 * self.dy)
        dlnu_dlnr = (LU[iy, jx + 1] - LU[iy, jx - 1]) / (2 * self.dx)
        dlnu_dlnT = (LU[iy + 1, jx] - LU[iy - 1, jx]) / (2 * self.dy)
        u = eint / rho                             # specific internal energy
        with np.errstate(invalid='ignore', divide='ignore'):
            # (du/dlnrho)_s = p/rho  ->  (dlnT/dlnrho)_s
            dlnT_dlnr_s = (p / rho - u * dlnu_dlnr) / (u * dlnu_dlnT)
            g1 = dlnp_dlnr + dlnp_dlnT * dlnT_dlnr_s
        g1 = np.where(np.isfinite(g1), g1, 5. / 3.)
        g1 = np.clip(g1, 1.01, 2.0)
        return T.reshape(sh), p.reshape(sh), g1.reshape(sh)


def fast_speed(rho, p, g1, bn, bt1, bt2):
    """Newtonian MHD fast speed with bn the NORMAL field component."""
    a2 = g1 * p / rho
    bsq = bn * bn + bt1 * bt1 + bt2 * bt2
    ca2 = bsq / rho
    qsq = a2 + ca2
    tmp = np.sqrt(np.maximum(qsq * qsq - 4.0 * a2 * bn * bn / rho, 0.0))
    return np.sqrt(0.5 * (qsq + tmp))
