"""Vectorised python port of src/eos/eos_composition.hpp -> the free-electron
fraction x_e = n_e/n_tot that ResistivityEOS() reads.

Built on the SAME (log10 rho, log10 T) node grid the run's EOS table uses, then
bilinearly interpolated (the code uses a bicubic Hermite patch on the same nodes;
on this grid, dlogT = 0.01, the difference is far below the decade-scale features
we care about here).  Read-only; nothing here writes into the run directories.
"""
import numpy as np

kboltz = 1.380649e-16
hplanck = 6.62607015e-27
m_u = 1.66053906660e-24
m_el = 9.1093837015e-28
ev = 1.602176634e-12
a_hyd = 1.008
a_hel = 4.002602
chi_h = 13.598434 * ev
chi_he1 = 24.587389 * ev
chi_he2 = 54.417765 * ev
chi_d = 4.478007 * ev
theta_rot = 85.4
theta_vib = 6332.0

# abun, chi, gfac, tc_a, tc_b, tc_c, tc_d           (Na K Ca Al Mg Fe)
DONOR = np.array([
    [2.0e-6, 5.139, 1.0,     10.045, 0.72,   1.08,  0.0],
    [1.2e-7, 4.341, 1.0,     12.479, 0.879,  0.879, 0.0],
    [2.2e-6, 6.113, 4.0,      4.990, 0.2394, 0.595, 1.398e-3],
    [3.0e-6, 5.986, 1.0 / 3,  5.014, 0.2179, 0.580, 2.264e-3],
    [3.8e-5, 7.646, 4.0,      5.89,  0.37,   0.73,  0.0],
    [3.2e-5, 7.902, 2.4,      5.44,  0.48,   0.48,  0.0]])
DNAME = ['Na', 'K', 'Ca', 'Al', 'Mg', 'Fe']

XH, YHE, A_METAL, MH = 0.7381, 0.2485, 16.0, 0.0


def _rotpart(t):
    z = np.empty_like(t)
    lo = t < 10.0 * theta_rot
    if lo.any():
        tl = t[lo]
        zs = np.zeros_like(tl)
        for jj in range(200):
            gns = 1.0 if jj % 2 == 0 else 3.0
            ej = jj * (jj + 1.0) * theta_rot
            zs += gns * (2.0 * jj + 1.0) * np.exp(-ej / tl)
        z[lo] = zs
    hi = ~lo
    if hi.any():
        u = theta_rot / t[hi]
        z[hi] = 2.0 * (t[hi] / theta_rot) * (1 + u / 3 + u * u / 15
                                             + 4 * u ** 3 / 315)
    return z


def _vibpart(t):
    u = theta_vib / t
    ex = np.where(u > 300.0, 0.0, np.exp(-np.minimum(u, 300.0)))
    return 1.0 / (1.0 - ex)


def _kdis(t):
    zrot = _rotpart(t)
    zvib = _vibpart(t)
    mh = a_hyd * m_u
    x = np.pi * mh * kboltz * t / hplanck**2
    expo = -chi_d / (kboltz * t)
    return np.where(expo < -700.0, 0.0,
                    x**1.5 * (16.0 / (zrot * zvib)) * np.exp(np.maximum(expo, -700.0)))


def _species(nel, nh_tot, nhe_tot, kdis, kh, khe1, khe2, kmet, nmet_base, fcond):
    """returns nh2, nh, nhii, nhe, nheii, nheiii, nmion (all same shape as nel)."""
    yion = np.where(kh > 0, kh / nel, 0.0)
    b = 1.0 + yion
    disc = b * b + 8.0 * nh_tot / np.where(kdis > 0, kdis, 1.0)
    nh = np.where(kdis > 0, 2.0 * nh_tot / (b + np.sqrt(disc)), 0.0)
    nh2 = np.where(kdis > 0, 0.5 * (nh_tot - nh * b), 0.5 * nh_tot)
    nhii = yion * nh
    r1 = np.where(khe1 > 0, khe1 / nel, 0.0)
    r2 = np.where(khe2 > 0, khe2 / nel, 0.0)
    nhe = nhe_tot / (1.0 + r1 + r1 * r2)
    nheii = r1 * nhe
    nheiii = r1 * r2 * nhe
    nmion = np.zeros_like(nel)
    for s in range(6):
        rm = np.where(kmet[s] > 0, kmet[s] / nel, 0.0)
        nmion += DONOR[s, 0] * nmet_base * fcond[s] * rm / (1.0 + rm)
    return nh2, nh, nhii, nhe, nheii, nheiii, nmion


def _cond_factors(t, p_cgs):
    """gas-phase fraction per donor; eos_metal_condensation = true path."""
    lgp = np.log10(np.maximum(p_cgs, 1e-30) / 1.0e6)     # bar
    out = []
    for s in range(6):
        a, bq, c, dq = DONOR[s, 3], DONOR[s, 4], DONOR[s, 5], DONOR[s, 6]
        denom = a - bq * lgp + dq * lgp * lgp - c * MH
        tc = np.where(denom > 0, 1.0e4 / np.where(denom > 0, denom, 1.0), 0.0)
        out.append(0.5 * (1.0 + np.tanh((t - tc) / (0.05 * np.where(tc > 0, tc, 1.0)))))
    return out


def _eval_at(rho, t, fcond):
    nh_tot = XH * rho / (a_hyd * m_u)
    nhe_tot = YHE * rho / (a_hel * m_u)
    z = 1.0 - XH - YHE
    nz = z * rho / (A_METAL * m_u)
    kdis = _kdis(t)
    sf = (2.0 * np.pi * m_el * kboltz * t / hplanck**2)**1.5

    def sah(chi, g):
        b = -chi / (kboltz * t)
        return np.where(b > -700.0, g * sf * np.exp(np.maximum(b, -700.0)), 0.0)
    kh = sah(chi_h, 1.0)
    khe1 = sah(chi_he1, 4.0)
    khe2 = sah(chi_he2, 1.0)
    kmet = [sah(DONOR[s, 1] * ev, DONOR[s, 2]) for s in range(6)]
    nmet_base = nh_tot * 10.0**MH

    nmax = nh_tot + 2.0 * nhe_tot
    lo = np.log(nmax) - 40.0 * np.log(10.0)
    hi = np.log(nmax).copy()
    lo = lo.copy()
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        ntry = np.exp(mid)
        _, _, nhii, _, nheii, nheiii, nmion = _species(
            ntry, nh_tot, nhe_tot, kdis, kh, khe1, khe2, kmet, nmet_base, fcond)
        res = nhii + nheii + 2.0 * nheiii + nmion - ntry
        lo = np.where(res > 0.0, mid, lo)
        hi = np.where(res > 0.0, hi, mid)
    nel = np.exp(0.5 * (lo + hi))
    nh2, nh, nhii, nhe, nheii, nheiii, nmion = _species(
        nel, nh_tot, nhe_tot, kdis, kh, khe1, khe2, kmet, nmet_base, fcond)
    nel = nhii + nheii + 2.0 * nheiii + nmion
    ntot = nh2 + nh + nhii + nhe + nheii + nheiii + nz + nel
    p = ntot * kboltz * t
    return nel / ntot, p


def evaluate(rho, t):
    """x_e(rho, T) with the two-pass condensation of EOSCompositionModel::Evaluate."""
    one = [np.ones_like(rho) for _ in range(6)]
    _, p0 = _eval_at(rho, t, one)
    fc = _cond_factors(t, p0)
    xe, _ = _eval_at(rho, t, fc)
    return xe


class XeTable:
    """log10 x_e on the run's own EOS node grid, bilinearly interpolated."""

    def __init__(self, xmin=-14.0, dx=0.05, nx=281, ymin=1.5, dy=0.01, ny=451):
        self.xmin, self.dx, self.nx = xmin, dx, nx
        self.ymin, self.dy, self.ny = ymin, dy, ny
        ld = xmin + dx * np.arange(nx)
        lt = ymin + dy * np.arange(ny)
        LD, LT = np.meshgrid(ld, lt)               # (ny, nx)
        xe = evaluate(10.0**LD, 10.0**LT)
        self.lxe = np.log10(np.maximum(xe, 1e-300))

    def __call__(self, rho, t):
        sh = np.asarray(rho).shape
        x = (np.log10(np.asarray(rho).ravel()) - self.xmin) / self.dx
        y = (np.log10(np.asarray(t).ravel()) - self.ymin) / self.dy
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        iy = np.clip(np.floor(y).astype(int), 0, self.ny - 2)
        fx = x - ix
        fy = y - iy                                 # linear continuation off-table
        f = (self.lxe[iy, ix] * (1 - fx) * (1 - fy)
             + self.lxe[iy, ix + 1] * fx * (1 - fy)
             + self.lxe[iy + 1, ix] * (1 - fx) * fy
             + self.lxe[iy + 1, ix + 1] * fx * fy)
        return (10.0**f).reshape(sh)


def eta_eos(xe, T, max_eta):
    """ResistivityEOS(), exactly."""
    xemin = 230.0 * np.sqrt(T) / max_eta
    xu = np.maximum(xe, xemin)
    eta = 230.0 * np.sqrt(T) / xu + 5.2e11 * 20.0 / (T * np.sqrt(T))
    return np.minimum(eta, max_eta)
