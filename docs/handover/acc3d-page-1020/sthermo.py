"""Entropy of the RY Per accretor gas (X 0.70, Y 0.28; read-only analysis).
Two independent gas entropies, both in k_B per baryon (s m_u / k_B):
  s_saha : analytic ideal H/He mixture with Saha ionization, Sackur-Tetrode per species (H0 g=2, p, He0, He+ g=2, He++, e g=2).
  s_tab  : thermodynamic integration of the run's own general-EOS table (log e/rho, log p/rho on log rho x log T):
           at T0 = table T_min, Helmholtz f(rho) = f_ref + int (p/rho) dln rho; s = (e - f)/T0; then s(rho,T) = s(rho,T0) + int de/T.
           The single constant (f_ref) is fixed by s_tab = s_saha at (log rho -8, T_min).  This checks the Saha model against the table.
Bundle copy of page/entropy/sthermo.py (10-20): table path from common.EOSTBL, cache in OUT/data.
Radiation: s_rad = 4 a T^3 / (3 rho) (LTE, T = T_gas) and from the M1 E_r: s_rad = (4/3) a^(1/4) E_r^(3/4) / rho."""
import numpy as np
k = 1.380649e-16; h = 6.62607015e-27; me = 9.1093837e-28; mu = 1.66053907e-24; mH = 1.6735575e-24; mp = 1.67262192e-24
mHe = 6.6464731e-24; eV = 1.602176634e-12; a_rad = 7.565733e-15
X, Y = 0.70, 0.28
import os as _os
from common import EOSTBL as TBL, OUT as _OUT


def nQ(m, T):
    return (2 * np.pi * m * k * T / h ** 2) ** 1.5


def s_saha(rho, T, full=False):
    """rho [g/cc], T [K] -> s_gas [k_B/baryon] (and H ionized fraction)."""
    rho = np.asarray(rho, float); T = np.asarray(T, float)
    nH = X * rho / mH; nHe = Y * rho / mHe; q = nQ(me, T); kT = k * T
    AH = q * np.exp(-13.598434 * eV / kT); A1 = 4 * q * np.exp(-24.587389 * eV / kT); A2 = q * np.exp(-54.417765 * eV / kT)
    lo = np.full(np.broadcast(rho, T).shape, -60.0); hi = np.log(nH + 2 * nHe) + 1e-9
    for _ in range(70):  # bisection on ln n_e (charge neutrality); f is monotone decreasing in n_e
        mid = 0.5 * (lo + hi); ne = np.exp(mid); r1 = A1 / ne; r2 = A2 / ne
        f = nH * AH / (AH + ne) + nHe * (r1 + 2 * r1 * r2) / (1 + r1 + r1 * r2) - ne
        pos = f > 0; lo = np.where(pos, mid, lo); hi = np.where(pos, hi, mid)
    ne = np.exp(0.5 * (lo + hi)); xH = AH / (AH + ne); r1 = A1 / ne; r2 = A2 / ne; z = 1 + r1 + r1 * r2
    S = np.zeros_like(lo)
    for n, g, m in ((nH * (1 - xH), 2, mH), (nH * xH, 1, mp), (nHe / z, 1, mHe), (nHe * r1 / z, 2, mHe),
                    (nHe * r1 * r2 / z, 1, mHe), (ne, 2, me)):
        n = np.maximum(n, 1e-300); S = S + n * (2.5 + np.log(g * nQ(m, T) / n))
    s = S * mu / rho
    return (s, xH) if full else s


def s_rad_T(rho, T):
    return 4 * a_rad * T ** 3 / (3 * rho) * mu / k


def s_rad_E(rho, Er):
    return 4.0 / 3.0 * a_rad ** 0.25 * np.maximum(Er, 0) ** 0.75 / rho * mu / k


class STable:
    def __init__(self, fn=TBL, cache=_os.path.join(_OUT, 'data', 'stable_cache.npz')):
        import os
        if cache and os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(fn):
            c = np.load(cache); self.S = c['S']; self.const = float(c['const'])
            self.nx, self.ny, self.x0, self.dx, self.y0, self.dy = [c['g'][i] for i in range(6)]
            self.nx, self.ny = int(self.nx), int(self.ny); return
        hdr = [ln for ln in open(fn) if ln.startswith('#')][:6]
        nx, ny, x0, dx, y0, dy = [float(v) for v in hdr[2][1:].split()]
        self.nx, self.ny, self.x0, self.dx, self.y0, self.dy = int(nx), int(ny), x0, dx, y0, dy
        d = np.loadtxt(fn).reshape(self.ny, self.nx, 2)
        le, lp = d[..., 0], d[..., 1]                     # [T, rho]
        lrho = x0 + dx * np.arange(self.nx); lT = y0 + dy * np.arange(self.ny)
        e = 10 ** le; T = 10 ** lT[:, None]
        # s along T at fixed rho: trapezoid of de/T
        ds = np.diff(e, axis=0) * 0.5 * (1 / T[1:] + 1 / T[:-1])
        sT = np.vstack([np.zeros((1, self.nx)), np.cumsum(ds, axis=0)])
        # Helmholtz at T0 along rho: df = (p/rho) dln rho
        prho = 10 ** lp[0]; dl = dx * np.log(10)
        f0 = np.concatenate([[0.0], np.cumsum(0.5 * (prho[1:] + prho[:-1]) * dl)])
        s0 = (e[0] - f0) / T[0]
        S = (s0[None, :] + sT) * mu / k                     # k_B/baryon, up to a constant
        iref = int(round((-8.0 - x0) / dx))
        self.const = float(s_saha(10 ** lrho[iref], 10 ** lT[0]) - S[0, iref])
        self.S = S + self.const; self.lrho = lrho; self.lT = lT
        if cache: os.makedirs(os.path.dirname(cache), exist_ok=True); np.savez(cache, S=self.S, const=self.const, g=np.array([self.nx, self.ny, x0, dx, y0, dy]))

    def __call__(self, rho, T):
        fx = np.clip((np.log10(rho) - self.x0) / self.dx, 0, self.nx - 1.000001)
        fy = np.clip((np.log10(T) - self.y0) / self.dy, 0, self.ny - 1.000001)
        i = fx.astype(np.int64); j = fy.astype(np.int64); wx = fx - i; wy = fy - j; S = self.S
        return ((1 - wy) * ((1 - wx) * S[j, i] + wx * S[j, i + 1]) + wy * ((1 - wx) * S[j + 1, i] + wx * S[j + 1, i + 1]))


if __name__ == '__main__':
    st = STable()
    print('const', st.const)
    for lr in (-14, -12, -10, -9, -8, -7, -6, -4):
        row = []
        for T in (5e3, 1e4, 2e4, 3e4, 5e4, 1e5, 3e5, 1e6, 1e7):
            row.append('%7.2f/%7.2f' % (st(10.0 ** lr, T), s_saha(10.0 ** lr, T)))
        print('logrho %3d  tab/saha at T 5e3..1e7:' % lr, ' '.join(row))
