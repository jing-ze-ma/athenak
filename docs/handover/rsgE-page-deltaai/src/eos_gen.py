#!/usr/bin/env python3
"""rsgE page: numpy port of AthenaK's general EOS (src/eos/eos_composition.hpp EvaluateAt + eos_table.cpp grid), for
T_gas(rho, e_int) of the dumps.  Only the configuration of this run is supported: eos_h2 = false, eos_ionization = true
(H, He, He+ Saha), eos_metal_ionization = false, eos_radiation = false (e_int in the dump is the GAS energy per volume).
Energy zero point as in the code with H2 off: e = 3/2 n_tot k T + n_HII chi_H + n_HeII chi_He1 + n_HeIII (chi_He1 + chi_He2).
The table is the code's node grid (log10 rho, log10 T; dlogd = eos_dlog, dlogt = 0.2 eos_dlog) of log10(e/rho); the
inversion is bilinear in (log rho, log10 e/rho) per cell (the code uses bicubic Hermite + Newton; at dlogT = 0.004 the
difference is far below the plot accuracy).  Below the table's lowest T the code clamps; T is returned clipped to the
table range and the caller marks floor cells."""
import numpy as np

KB = 1.380649e-16; HP = 6.62607015e-27; MU = 1.66053906660e-24; ME = 9.1093837015e-28; EV = 1.602176634e-12
A_H, A_HE = 1.008, 4.002602
CHI_H, CHI_HE1, CHI_HE2 = 13.598434 * EV, 24.587389 * EV, 54.417765 * EV


def evaluate(rho, t, xh=0.653848, yhe=0.326395, a_metal=16.0):
    """vectorised EvaluateAt (H2 off, metal ion off): returns e_spec, p_spec, mu, xe"""
    rho = np.asarray(rho, float); t = np.asarray(t, float)
    rho, t = np.broadcast_arrays(rho, t)
    nh_tot = xh * rho / (A_H * MU); nhe_tot = yhe * rho / (A_HE * MU)
    nz = (1 - xh - yhe) * rho / (a_metal * MU) if a_metal > 0 else 0 * rho
    sf = (2 * np.pi * ME * KB * t / HP ** 2) ** 1.5
    bh, b1, b2 = -CHI_H / (KB * t), -CHI_HE1 / (KB * t), -CHI_HE2 / (KB * t)
    kh = np.where(bh > -700, sf * np.exp(np.maximum(bh, -700)), 0.)
    khe1 = np.where(b1 > -700, 4 * sf * np.exp(np.maximum(b1, -700)), 0.)
    khe2 = np.where(b2 > -700, sf * np.exp(np.maximum(b2, -700)), 0.)
    nmax = nh_tot + 2 * nhe_tot

    def species(nel):
        yion = np.where(kh > 0, kh / nel, 0.)
        nh = nh_tot / (1 + yion); nhii = yion * nh
        r1 = np.where(khe1 > 0, khe1 / nel, 0.); r2 = np.where(khe2 > 0, khe2 / nel, 0.)
        nhe = nhe_tot / (1 + r1 + r1 * r2)
        return nh, nhii, nhe, r1 * nhe, r1 * r2 * nhe

    lo = np.log(nmax) - 40 * np.log(10); hi = np.log(nmax)
    for _ in range(80):                                 # same bisection as the code (80 halvings or 1e-13)
        mid = 0.5 * (lo + hi)
        ntry = np.exp(mid)
        nh, nhii, nhe, nheii, nheiii = species(ntry)
        pos = (nhii + nheii + 2 * nheiii - ntry) > 0
        lo = np.where(pos, mid, lo); hi = np.where(pos, hi, mid)
        if np.all(hi - lo < 1e-13):
            break
    anyion = (kh > 0) | (khe1 > 0)
    nel = np.where(anyion, np.exp(0.5 * (lo + hi)), 1e-40 * nmax)
    nh, nhii, nhe, nheii, nheiii = species(nel)
    nel = nhii + nheii + 2 * nheiii
    ntot = nh + nhii + nhe + nheii + nheiii + nz + nel
    pres = ntot * KB * t
    eint = 1.5 * pres + nhii * CHI_H + nheii * CHI_HE1 + nheiii * (CHI_HE1 + CHI_HE2)
    return eint / rho, pres / rho, rho / (ntot * MU), nel / ntot


class Table:
    def __init__(self, logd=(-21.0, -1.0), logt=(3.3, 6.8), dlog=0.02, **comp):
        dlogd, dlogt = dlog, 0.2 * dlog
        nx = int(round((logd[1] - logd[0]) / dlogd)) + 1; ny = int(round((logt[1] - logt[0]) / dlogt)) + 1
        self.nx, self.ny = nx, ny
        self.x0, self.dx = logd[0], (logd[1] - logd[0]) / (nx - 1)
        self.y0, self.dy = logt[0], (logt[1] - logt[0]) / (ny - 1)
        X, Y = np.meshgrid(self.x0 + np.arange(nx) * self.dx, self.y0 + np.arange(ny) * self.dy, indexing='xy')  # [j(T), i(rho)]
        es, ps, mu, xe = evaluate(10 ** X, 10 ** Y, **comp)
        self.fe = np.log10(es)                            # [ny, nx]
        self.fp = np.log10(ps)

    def temperature(self, rho, eint):
        """T [K] from rho [g/cm3] and gas internal energy per volume [erg/cm3]; any shape. Clipped to the table range.
        Returns (T, below) with below = True where e/rho is under the table's lowest-T row (floor / clamp)."""
        shp = np.shape(rho)
        x = np.clip((np.log10(np.ravel(rho)) - self.x0) / self.dx, 0, self.nx - 1 - 1e-9)
        i = x.astype(np.int64); fx = x - i
        ft = np.log10(np.maximum(np.ravel(eint), 1e-300) / np.ravel(rho))
        row = lambda j: (1 - fx) * self.fe[j, i] + fx * self.fe[j, i + 1]   # noqa: E731
        jlo = np.zeros(i.shape, np.int64); jhi = np.full(i.shape, self.ny - 1, np.int64)
        while True:                                       # bisection on the (monotone) rows
            act = jhi - jlo > 1
            if not act.any():
                break
            jm = (jlo + jhi) // 2
            up = row(jm) > ft
            jhi = np.where(act & up, jm, jhi); jlo = np.where(act & ~up, jm, jlo)
        flo, fhi = row(jlo), row(jhi)
        fy = np.clip((ft - flo) / np.where(fhi > flo, fhi - flo, 1.0), 0, 1)
        logT = self.y0 + (jlo + fy) * self.dy
        below = ft < row(np.zeros_like(i))
        return (10 ** logT).reshape(shp).astype(np.float32), below.reshape(shp)


if __name__ == '__main__':
    import sys
    import time
    t0 = time.time(); tb = Table()
    print('table %d x %d built in %.1f s' % (tb.nx, tb.ny, time.time() - t0))
    es, ps, mu, xe = evaluate(10 ** -11.0, 10 ** 5.05)
    print('grid centre mu %.6f xe %.6f   (code run log: mu 0.648009, n_e/n_tot 0.526021)' % (mu, xe))
    if len(sys.argv) > 1:                                 # IC column check: r rho eint F E | T F_MLT/F
        ic = np.loadtxt(sys.argv[1], comments='#')
        T, below = tb.temperature(ic[:, 1], ic[:, 2])
        rel = T / ic[:, 5] - 1
        ok = ic[:, 5] > 2100
        print('IC column: %d rows, max |T_inv/T_IC - 1| = %.2e (T > 2100 K), median %.2e; floored rows %d' %
              (len(ic), np.abs(rel[ok]).max(), np.median(np.abs(rel[ok])), (~ok).sum()))
