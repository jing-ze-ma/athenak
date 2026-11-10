"""csio.py: AthenaK cubed-sphere dump I/O for fe_hires (self-contained; no /viper paths).

Verbatim copies of the pieces fe_hires used on viper:
  EOS        dhjcs.EOS (docs/handover/scripts/dhjcs.py), invert() chunked as hr_obs._inv_chunked
  panel_map  dhjcs.panel_map / FRAMES
  geometry   w121_obs_1003/obs.py geometry()
  radii      hr_obs.radii(): true faces/centres from a version-1.2 .bin (mb_x1f/mb_x1v); for
             version-1.1 dumps the radial grid is rebuilt from the run deck (8-coefficient
             polynomial stretch only)
bin_convert.py is the athenak vis/python reader (v1.1 + v1.2), vendored next to this file.
"""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bin_convert as bc  # noqa: E402


def read_dump(fn):
    return bc.read_binary(fn)


from numba import njit, prange  # noqa: E402


@njit(cache=True, parallel=True)
def _kscan(LE, LP, ix, fx, target):
    """per cell: the EOS._invert_np bracket k (count of finite le < target, clipped to
    [first finite + 1, ny - 1]) and le, lp at k - 1 and k; le/lp interpolated in rho exactly as
    there (a*(1 - f) + b*f, no FMA)."""
    ny = LE.shape[0]
    n = ix.size
    K = np.empty(n, np.int64)
    E0 = np.empty(n)
    E1 = np.empty(n)
    P0 = np.empty(n)
    P1 = np.empty(n)
    for i in prange(n):
        a = ix[i]
        f = fx[i]
        g = 1 - f
        t = target[i]
        first = -1
        cnt = 0
        for j in range(ny):
            le = LE[j, a]*g + LE[j, a + 1]*f
            if np.isfinite(le):
                if first < 0:
                    first = j
                if le < t:
                    cnt += 1
            elif -np.inf < t:
                cnt += 1
        if first < 0:
            first = 0
        k = min(max(cnt, first + 1), ny - 1)
        K[i] = k
        E0[i] = LE[k - 1, a]*g + LE[k - 1, a + 1]*f
        E1[i] = LE[k, a]*g + LE[k, a + 1]*f
        P0[i] = LP[k - 1, a]*g + LP[k - 1, a + 1]*f
        P1[i] = LP[k, a]*g + LP[k, a + 1]*f
    return K, E0, E1, P0, P1


class EOS:
    """AthenaK general-EOS table dump (eos_table*.txt): rho, eint -> T [K], p [erg/cm^3]."""
    def __init__(self, fn):
        with open(fn) as f:
            lines = [next(f) for _ in range(6)]
        nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2][1:].split()]
        self.nx, self.ny = int(nx), int(ny)
        self.xmin, self.dx, self.ymin, self.dy = xmin, dx, ymin, dy
        d = np.loadtxt(fn)
        self.le = d[:, 0].reshape(self.ny, self.nx)
        self.lp = d[:, 1].reshape(self.ny, self.nx)
        self.lt = self.ymin + self.dy*np.arange(self.ny)
        self.ld = self.xmin + self.dx*np.arange(self.nx)

    def _invert(self, rho, eint):
        """numba scan (fehires_fast): same arithmetic as _invert_np, without the (ny, N) arrays."""
        rho = np.asarray(rho, float)
        eint = np.asarray(eint, float)
        sh = rho.shape
        rho = rho.ravel()
        eint = eint.ravel()
        x = (np.log10(rho) - self.xmin)/self.dx
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        fx = np.clip(x - ix, 0, 1)
        target = np.log10(eint/rho)
        k, e0, e1, p0, p1 = _kscan(self.le, self.lp, ix, fx, target)
        f = np.clip((target - e0)/(e1 - e0), 0, 1)
        lT = self.lt[k - 1] + f*self.dy
        lP = p0*(1 - f) + p1*f
        T = 10**lT
        p = 10**lP*rho
        return T.reshape(sh), p.reshape(sh)

    def _invert_np(self, rho, eint):
        rho = np.asarray(rho, float)
        eint = np.asarray(eint, float)
        sh = rho.shape
        rho = rho.ravel()
        eint = eint.ravel()
        x = (np.log10(rho) - self.xmin)/self.dx
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        fx = np.clip(x - ix, 0, 1)
        le = self.le[:, ix]*(1 - fx) + self.le[:, ix + 1]*fx      # (ny, N)
        lp = self.lp[:, ix]*(1 - fx) + self.lp[:, ix + 1]*fx
        target = np.log10(eint/rho)
        fin = np.isfinite(le)
        lef = np.where(fin, le, -np.inf)
        first = np.argmax(fin, axis=0)
        cnt = (lef < target[None, :]).sum(axis=0)
        k = np.clip(cnt, first + 1, self.ny - 1)
        e0 = le[k - 1, np.arange(len(k))]
        e1 = le[k, np.arange(len(k))]
        f = np.clip((target - e0)/(e1 - e0), 0, 1)
        lT = self.lt[k - 1] + f*self.dy
        lP = lp[k - 1, np.arange(len(k))]*(1 - f) + lp[k, np.arange(len(k))]*f
        T = 10**lT
        p = 10**lP*rho
        return T.reshape(sh), p.reshape(sh)

    def invert(self, rho, eint, ch=200000):
        rho = np.asarray(rho, float)
        eint = np.asarray(eint, float)
        r, e = rho.ravel(), eint.ravel()
        T, p = np.empty_like(r), np.empty_like(r)
        for s in range(0, r.size, ch):
            T[s:s+ch], p[s:s+ch] = self._invert(r[s:s+ch], e[s:s+ch])
        return T.reshape(rho.shape), p.reshape(rho.shape)


FRAMES = {0: ([0, 1, 0], [0, 0, 1], [1, 0, 0]), 1: ([-1, 0, 0], [0, 0, 1], [0, 1, 0]),
          2: ([0, -1, 0], [0, 0, 1], [-1, 0, 0]), 3: ([0, 1, 0], [-1, 0, 0], [0, 0, 1]),
          4: ([1, 0, 0], [0, 0, 1], [0, -1, 0]), 5: ([0, 1, 0], [1, 0, 0], [0, 0, -1])}


def panel_map(p, xi, eta):
    a, b, n = [np.array(v, float) for v in FRAMES[p]]
    x = np.tan(xi)
    y = np.tan(eta)
    dl = np.sqrt(1 + x*x + y*y)
    return (a[None, None, :]*x[..., None] + b[None, None, :]*y[..., None]
            + n[None, None, :])/dl[..., None]


def geometry(raw):
    """unit radial vector X and unit tangent vectors along xi/eta per column; solid angle."""
    g = np.asarray(raw['mb_geometry'])
    nmb, nj, nk = g.shape[0], raw['nx2_mb'], raw['nx3_mb']
    per = nmb//6
    X = np.zeros((nmb, nk, nj, 3))
    EXI = np.zeros_like(X)
    EET = np.zeros_like(X)
    w = np.zeros((nmb, nk, nj))
    for m in range(nmb):
        p = m//per
        xi = np.pi/4*(g[m, 2] + (g[m, 3] - g[m, 2])*(np.arange(nj) + 0.5)/nj)
        eta = np.pi/4*(g[m, 4] + (g[m, 5] - g[m, 4])*(np.arange(nk) + 0.5)/nk)
        XI, ETA = np.meshgrid(xi, eta)
        X[m] = panel_map(p, XI, ETA)
        d = 1e-6
        e = panel_map(p, XI + d, ETA) - panel_map(p, XI - d, ETA)
        EXI[m] = e/np.linalg.norm(e, axis=-1)[..., None]
        e = panel_map(p, XI, ETA + d) - panel_map(p, XI, ETA - d)
        EET[m] = e/np.linalg.norm(e, axis=-1)[..., None]
        x, y = np.tan(XI), np.tan(ETA)
        w[m] = (1 + x*x)*(1 + y*y)/(1 + x*x + y*y)**1.5
    dOm = 4*np.pi*w/w.sum()
    return X, EXI, EET, dOm


def _deck(fn):
    kv = {}
    for ln in open(fn):
        s = ln.split('#')[0].split('=')
        if len(s) == 2 and s[0].strip() not in kv:
            kv[s[0].strip()] = s[1].strip()
    return kv


def radii(raw, deck=None):
    """radial faces rf (n+1) and centres rc (n) [cm]."""
    if 'mb_x1f' in raw:
        rf = np.asarray(raw['mb_x1f'], float)[0]
        rc = np.asarray(raw['mb_x1v'], float)[0]
        return rf, rc
    if deck is None:
        raise SystemExit('version-1.1 dump without mb_x1f: pass --deck <run.athinput>')
    kv = _deck(deck)
    if any(k.startswith('f_stretch_r_b') for k in kv):
        raise SystemExit('deck uses the f_stretch_r_b* bump terms; only version-1.2 dumps '
                         '(with mb_x1f) are supported for it')
    n = int(kv['nx1'])
    r0, r1 = float(kv['x1min']), float(kv['x1max'])
    cs = [float(kv['f_stretch_r_c%d' % k]) for k in range(1, 9)]

    def stretch(xi):
        u = xi.copy()
        xk = xi.copy()
        for ck in cs:
            u += ck*xk*(1 - xi)
            xk *= xi
        return r0 + (r1 - r0)*u
    return stretch(np.arange(n + 1)/n), stretch((np.arange(n) + 0.5)/n)
