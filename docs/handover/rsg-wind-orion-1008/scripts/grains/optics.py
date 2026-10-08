"""Mie optics of amorphous Mg2SiO4 grains (Jaeger et al. 2003 sol-gel n,k via the
Kitzmann & Heng 2018 LX-MIE compilation) + absorption variants.
Writes out/optics.npz: a, lam, Qabs[v,a,lam], Qsca, g, Qpr, Planck means.
usage: venv/bin/python optics.py"""
import numpy as np
import miepython as mp

D = '/orion/ptmp/jinma/rsg_wind_1008/grains/'
VARIANTS = ['lowk', 'pure', 'Fe3e-4', 'Fe1e-3']   # see nk_variant
h, c, kB = 6.62607e-27, 2.99792458e10, 1.380649e-16


def load_nk():
    d = np.loadtxt(D + 'data/Mg2SiO4_amorph_sol-gel.dat', comments='#')
    lam, i = np.unique(d[:, 0], return_index=True)
    return lam, d[i, 1], d[i, 2]


def nk_variant(var, lam, n, k):
    """lowk : k x 0.1 below 5 um (assumed: sol-gel visible/NIR k is often read as an upper
              limit; crystalline forsterite is cleaner) ;
       pure : Jaeger et al. 2003 as measured ;
       FeXX : k = max(k, XX) below 8 um, a grey Fe-inclusion absorption floor (assumed;
              k = 1e-3 ~ Mg_{2-2x}Fe_{2x}SiO4 with x ~ 0.01 by linear k scaling of the
              Dorschner et al. 1995 MgFeSiO4 glass, k(1 um) = 0.069 at x = 0.5)."""
    k = k.copy()
    if var == 'lowk':
        k[lam < 5] *= 0.1
    elif var.startswith('Fe'):
        kf = float(var[2:])
        k[lam < 8] = np.maximum(k[lam < 8], kf)
    return n, k


def planck_lam(lam_um, T):
    l = lam_um * 1e-4
    return 2 * h * c**2 / l**5 / np.expm1(np.minimum(h * c / (l * kB * T), 700.))


def mie_check():
    """Bohren & Huffman (1983) App. A test: m = 1.55, x = 5.213 -> Qsca = Qext = 3.1054;
    small-x limit Qabs = 4 x Im[(m^2-1)/(m^2+2)]."""
    qe, qs, qb, g = mp.efficiencies_mx(1.55, 5.213)
    m = 1.6 - 1e-3j
    x = 1e-3
    qe2, qs2, _, _ = mp.efficiencies_mx(m, x)
    ray = 4 * x * np.imag((m**2 - 1) / (m**2 + 2)) * -1     # miepython: m = n - ik
    return qe, qs, (qe2 - qs2), abs(ray)


if __name__ == '__main__':
    print('Mie check (BH83: 3.1054):', mie_check())
    lam0, n0, k0 = load_nk()
    lam = np.geomspace(0.2, 900., 400)
    a = np.geomspace(1e-7, 5e-4, 60)              # cm, 1 nm .. 5 um
    nv = len(VARIANTS)
    Qa = np.zeros((nv, len(a), len(lam)))
    Qs = np.zeros_like(Qa)
    gg = np.zeros_like(Qa)
    kv = np.zeros((nv, len(lam)))
    for iv, v in enumerate(VARIANTS):
        n, k = nk_variant(v, lam0, n0, k0)
        nn = np.interp(np.log(lam), np.log(lam0), n)
        kk = np.exp(np.interp(np.log(lam), np.log(lam0), np.log(k)))
        kv[iv] = kk
        m = nn - 1j * kk
        for ia, aa in enumerate(a):
            x = 2 * np.pi * aa / (lam * 1e-4)
            qe, qs, qb, g = mp.efficiencies_mx(m, x)
            Qa[iv, ia] = qe - qs
            Qs[iv, ia] = qs
            gg[iv, ia] = g
    Qpr = Qa + (1 - gg) * Qs
    T = np.geomspace(100., 5000., 200)
    B = np.array([planck_lam(lam, t) for t in T])            # [T, lam]
    wl = np.gradient(lam)
    Bn = (B * wl).sum(1)
    QaT = np.einsum('val,tl->vat', Qa, B * wl) / Bn          # Planck-mean Qabs(a, T)
    QprT = np.einsum('val,tl->vat', Qpr, B * wl) / Bn
    np.savez(D + 'out/optics.npz', a=a, lam=lam, Qabs=Qa, Qsca=Qs, g=gg, Qpr=Qpr,
             k=kv, T=T, QaT=QaT, QprT=QprT, variants=np.array(VARIANTS))
    print('done', Qa.shape)
