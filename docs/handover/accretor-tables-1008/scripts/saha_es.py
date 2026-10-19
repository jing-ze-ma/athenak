#!/usr/bin/env python3
"""LTE (Saha) electron-scattering opacity kappa_es = sigma_T n_e / rho [cm^2/g] for the TOPS
mixture (number fractions read from tables/mixture_x0.36_z0.02.txt).

Species: H (13.598 eV), He I/II (24.587, 54.418 eV), and the 19 metals of the TOPS mixture
with their FIRST ionisation only (they matter for n_e only where H is neutral, T < ~5000 K;
higher metal stages add <= Z/2 ~ 1 % of the electrons and are absorbed into the edge match
of the tabulated opacity).  Partition-function ratios: ground-state statistical weights
(2 U+/U0 below).  n_e is found by bisection in log n_e (charge conservation)."""
import numpy as np

KB_EV = 8.617333262e-5
MU = 1.66053906660e-24
SIGT = 6.6524587321e-25
SAHA_C = 2.4146830e15          # (2 pi m_e k / h^2)^{3/2} [cm^-3 K^-3/2]

AMASS = dict(h=1.00794, he=4.002602, c=12.0107, n=14.0067, o=15.9994, ne=20.1797,
             na=22.98977, mg=24.305, al=26.98154, si=28.0855, p=30.97376, s=32.065,
             cl=35.453, ar=39.948, k=39.0983, ca=40.078, ti=47.867, cr=51.9961,
             mn=54.93805, fe=55.845, co=58.9332, ni=58.6934)
# first ionisation potential [eV], 2*g(+)/g(0) (ground terms)
ION1 = dict(c=(11.260, 2*6/9), n=(14.534, 2*9/4), o=(13.618, 2*4/9), ne=(21.565, 2*6/1),
            na=(5.139, 2*1/2), mg=(7.646, 2*2/1), al=(5.986, 2*1/6), si=(8.152, 2*6/9),
            p=(10.487, 2*9/4), s=(10.360, 2*4/9), cl=(12.968, 2*9/4), ar=(15.760, 2*6/1),
            k=(4.341, 2*1/2), ca=(6.113, 2*2/1), ti=(6.828, 2*28/21), cr=(6.767, 2*6/7),
            mn=(7.434, 2*7/6), fe=(7.902, 2*30/25), co=(7.881, 2*21/28), ni=(7.640, 2*10/21))


def read_mixture(fn):
    v = open(fn).read().split()
    f = {v[i+1].lower(): float(v[i]) for i in range(0, len(v), 2)}
    s = sum(f.values())
    return {k: x/s for k, x in f.items()}


class Saha:
    def __init__(self, mixfile):
        self.f = read_mixture(mixfile)
        self.mbar = sum(self.f[k]*AMASS[k] for k in self.f)*MU   # g per nucleus
        self.X = self.f['h']*AMASS['h']*MU/self.mbar

    def ne(self, T, rho):
        """n_e [cm^-3] for scalar T [K] and scalar or array rho [g/cm^3] (vectorised)."""
        rho = np.atleast_1d(np.asarray(rho, float))
        n = rho/self.mbar
        kT = KB_EV*T
        g = SAHA_C*T**1.5
        stages = [('h', [(13.598, 1.0)]), ('he', [(24.587, 4.0), (54.418, 1.0)])]
        stages += [(k, [ION1[k]]) for k in self.f if k in ION1]
        phis = [(self.f[k], [g*w*np.exp(-chi/kT) for chi, w in st]) for k, st in stages]

        def charge(ne):
            q = np.zeros_like(ne)
            for fk, ph in phis:
                r = [np.ones_like(ne)]
                for p in ph:
                    r.append(r[-1]*p/ne)
                s = sum(r)
                q += fk*n*sum(j*r[j] for j in range(len(r)))/s
            return q - ne
        lo = np.full_like(n, -40.0)
        hi = np.log10(n*2.0)
        for _ in range(70):
            mid = 0.5*(lo + hi)
            pos = charge(10**mid) > 0
            lo = np.where(pos, mid, lo)
            hi = np.where(pos, hi, mid)
        return 10**(0.5*(lo + hi))

    def kes(self, T, rho):
        return SIGT*self.ne(T, rho)/rho


if __name__ == '__main__':
    s = Saha('/viper/ptmp2/jinma/lbv_1008/agcar/tables/mixture_x0.36_z0.02.txt')
    print('X from mixture %.4f, mbar/m_u %.4f, full-ionisation 0.2(1+X)-ish:' % (
        s.X, s.mbar/MU), s.kes(1e7, 1e-10))
    for T in (3000, 4000, 5000, 5800, 7000, 1e4, 2e4, 5e4):
        print(T, ['%.3e' % s.kes(T, r) for r in (1e-20, 1e-18, 1e-16, 1e-14, 1e-12)])
