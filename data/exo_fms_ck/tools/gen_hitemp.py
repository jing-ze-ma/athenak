"""Build the high-temperature extensions of the Exo-FMS 11-band correlated-k data.

  python3 gen_hitemp.py [--met 1x|10x] [--v1]

  CE_tables/FastChem_ck_<met>_int_hiT2.txt : 100..10100 K (25 K), rows <= 6100 K verbatim
  ck/Premixed_<met>_g8_11_hiT2.txt         : 100..6100 K verbatim + 6125..7075 K (25 K)
                                             + 7100..10100 K (200 K)
Inputs ck/Premixed_<met>_g8_11.txt and CE_tables/FastChem_ck_<met>_int.txt under CK_DATA,
outputs under CK_OUT (default CK_DATA). FastChem abundances from FASTCHEM_INPUT (for 10x:
Asplund 2009 with the metals +1 dex). --v1 rebuilds the first version (_hiT, value-
continuous only) byte-identically.

Composition above 6100 K: FastChem 3.1.3 (pyfastchem, gas phase), times the
per-(p, column) ratio upstream/FastChem at 6100 K.

Line opacity above 6100 K (per p, band b, all g alike):
  k_b(T) = k(6100) * [X_b(T)/X_b(6100)]^alpha_b  +  kH_b(T)
X_b = FastChem VMR/mu summed over the neutral line carriers of the band (molecules
for b <= 7, i.e. >= 0.85 um; neutral metals + TiO/VO/FeH for b >= 8). alpha_b in
[0, 1] is fitted on the table's own 5100..6100 K points (log-log least squares through
6100 K). kH is the neutral-hydrogen bound-free (n = 1..40, hydrogenic, g = 1) + H ff
(H+ e-, g = 1) continuum at the band's representative wavelength (the one the code uses
for H-), with stimulated emission; it is absent from the upstream data and from
ck_continuum.

C1 join (v2): between 6100 and 7100 K every k (per p, b, g) and every composition column
is a log-space blend  log q = (1-w) log q_a + w log q_b  of the upstream table continued
at its own last-interval slope (q_a = q(6100) (T/6100)^s, s = the 5900..6100 K secant for
k, 6075..6100 K for the composition) and the physics above (q_b), with the quintic
smootherstep w(x), x = ln(T/6100)/ln(7100/6100), so value, first and second derivative
of log q in log T are continuous at both ends. The 25 K nodes over the window keep the
piecewise-linear lookup's slope steps small. v1 instead added kH(T) - kH(6100) (clipped
at 0) and joined the composition by the ratio only: value-continuous, slope jump ~650 in
d log k/d log T at low p."""
from fc_run import run as fc_vmr
from carriers import group_vmr
from ck_lib import read_ktable, read_ce, DATA
import pyfastchem
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

T_CE_NEW = np.arange(6125., 10100.1, 25.)
T_K_NEW_V1 = np.arange(6300., 10100.1, 200.)
T_K_NEW_V2 = np.concatenate([np.arange(6125., 7099., 25.),
                             np.arange(7100., 10100.1, 200.)])
TW0, TW1 = 6100., 7100.          # C1 blend window (v2)
KB = 1.380649e-16
H = 6.62607015e-27
C = 2.99792458e10


def kH_neutral(lam_um, T, P, vmrH, vmre, vmrHp, rho):
    """H bf + ff opacity [cm^2/g] at lam (um), hydrogenic with Gaunt factor 1."""
    nu = C/(lam_um*1e-4)
    x = H*nu/(KB*T)
    stim = -np.expm1(-x)
    ntot = P*1e6/(KB*T)
    nH = vmrH*ntot
    chi = 13.598*1.602176634e-12
    bf = 0.0
    for n in range(1, 41):
        if nu >= 3.28984e15/n**2:
            pop = n**2*np.exp(-chi*(1-1/n**2)/(KB*T)) / \
                2.0   # N_n/N_H, U_H ~ 2 (ground g=2)
            pop *= 2.0                                      # g_n = 2 n^2
            bf += 2.815e29/(n**5*nu**3)*pop
    ff = 3.692e8/np.sqrt(T)/nu**3*vmre*ntot*vmrHp*ntot      # per H atom-independent
    return (bf*nH*stim + ff*stim)/rho


def smootherstep(T):
    """Quintic 0 -> 1 in x = ln(T/TW0)/ln(TW1/TW0), C2 at both ends."""
    x = np.clip(np.log(np.asarray(T, float)/TW0)/np.log(TW1/TW0), 0.0, 1.0)
    return x**3*(10.0 - 15.0*x + 6.0*x*x)


def main(met='1x', v1=False):
    OUT = os.environ.get('CK_OUT', DATA).rstrip('/') + '/'
    suf = '_hiT' if v1 else '_hiT2'
    out_ce = OUT+f'CE_tables/FastChem_ck_{met}_int{suf}.txt'
    out_k = OUT+f'ck/Premixed_{met}_g8_11{suf}.txt'
    in_ce = DATA+f'CE_tables/FastChem_ck_{met}_int.txt'
    kt = read_ktable(DATA+f'ck/Premixed_{met}_g8_11.txt')
    cT, cP, CE = read_ce(in_ce)
    raw_ce = open(in_ce).read().split()
    assert cT[-1] == 6100. and kt['T'][-1] == 6100.
    T_K_NEW = T_K_NEW_V1 if v1 else T_K_NEW_V2
    # ---- composition -------------------------------------------------------------
    fc, flag = fc_vmr(np.concatenate([[6100.], T_CE_NEW]), cP)
    assert np.count_nonzero(flag) == 0, 'FastChem did not converge everywhere'
    ratio = CE[-1]/fc[0]                            # (nP, 6)
    CEn = fc[1:]*ratio[None]
    if not v1:
        # C1: log blend of the upstream columns continued at their last slope
        sC = np.log(CE[-1]/CE[-2])/np.log(cT[-1]/cT[-2])            # (nP, 6)
        la = np.log(CE[-1])[None] + sC[None]*np.log(T_CE_NEW/6100.)[:, None, None]
        w = smootherstep(T_CE_NEW)[:, None, None]
        CEn = np.exp((1.0-w)*la + w*np.log(CEn))
    # H+ is not a table column; get it for the H ff term from the same FastChem call
    nT0 = len(cT)
    nP = len(cP)
    with open(out_ce, 'w') as f:
        f.write(f"{nT0+len(T_CE_NEW)} {nP} {(nT0+len(T_CE_NEW))*nP} 5\nH2 He H e- H-\n")
        f.write(' '.join(raw_ce[9:9+nT0]+[f"{t:.1f}" for t in T_CE_NEW])+'\n')
        f.write(' '.join(raw_ce[9+nT0:9+nT0+nP])+'\n')
        p = 9+nT0+nP
        for i in range(nT0*nP):
            f.write(' '.join(raw_ce[p+6*i:p+6*i+6])+'\n')
        for it in range(len(T_CE_NEW)):
            for ip in range(nP):
                f.write(' '.join(f"{v:.10e}" for v in CEn[it, ip])+'\n')
    # ---- k-table -------------------------------------------------------------------
    T = kt['T']
    P = kt['P']
    K = kt['K']
    gw = kt['gw']
    wl = kt['wl']
    nb = len(wl)-1
    lam = 2/(1/wl[:-1]+1/wl[1:])
    fit = [5100., 5300., 5500., 5700., 5900.]
    ifit = [list(T).index(t) for t in fit]
    Knew = np.empty((len(T_K_NEW), len(P), nb, K.shape[-1]))
    alpha = np.empty((len(P), nb))
    from carriers import D as FCD
    fcH = pyfastchem.FastChem(FCD+'element_abundances/asplund_2009.dat',
                              FCD+'logK/logK.dat', 0)
    iH, ie, iHp = (fcH.getGasSpeciesIndex(s) for s in ('H', 'e-', 'H1+'))

    def kH_all(Tv, Pv):
        inp = pyfastchem.FastChemInput()
        out = pyfastchem.FastChemOutput()
        inp.temperature = [Tv]
        inp.pressure = [Pv]
        fcH.calcDensities(inp, out)
        n = np.array(out.number_densities)[0]
        ng = Pv*1e6/(KB*Tv)
        mu = out.mean_molecular_weight[0]
        rho = Pv*1e6*mu*1.66053907e-24/(KB*Tv)
        return np.array([kH_neutral(lm, Tv, Pv, n[iH]/ng, n[ie]/ng, n[iHp]/ng, rho)
                         for lm in lam])
    # upstream last-interval slope d log k / d log T, per (p, b, g)
    lK0 = np.log(np.maximum(K[-1], 1e-99))
    sK = (lK0 - np.log(np.maximum(K[-2], 1e-99)))/np.log(T[-1]/T[-2])
    wK = smootherstep(T_K_NEW)
    for ip, Pv in enumerate(P):
        X = group_vmr(np.concatenate([[6100.], fit, T_K_NEW]), Pv, nb)   # (1+5+nnew, nb)
        km = (K[:, ip]*gw).sum(-1)                                      # (nT, nb)
        for b in range(nb):
            x = np.log(X[1:6, b]/X[0, b])
            y = np.log(np.maximum(km[ifit, b], 1e-300)/max(km[-1, b], 1e-300))
            a = (x*y).sum()/max((x*x).sum(), 1e-30)
            alpha[ip, b] = min(max(a, 0.0), 1.0)
        S = (X[6:]/X[0])**alpha[ip]                                      # (nnew, nb)
        kH0 = kH_all(6100., Pv)
        for it, Tv in enumerate(T_K_NEW):
            if v1:
                dH = np.maximum(kH_all(Tv, Pv)-kH0, 0.0)
                Knew[it, ip] = K[-1, ip]*S[it][:, None] + dH[:, None]
            else:
                kb = K[-1, ip]*S[it][:, None] + kH_all(Tv, Pv)[:, None]
                la = lK0[ip] + sK[ip]*np.log(Tv/6100.)
                Knew[it, ip] = np.exp((1.0-wK[it])*la + wK[it]*np.log(kb))
    Knew = np.maximum(Knew, 1e-99)
    if not v1:
        nfix = np.count_nonzero(np.diff(Knew, axis=-1) < 0)
        Knew = np.maximum.accumulate(Knew, axis=-1)
        print('g-monotonicity fixes in the blend window:', nfix)
    assert np.all(np.diff(Knew, axis=-1) >= 0), 'k-distribution not monotone in g'
    raw = kt['raw']
    nT0 = len(T)
    ng = K.shape[-1]
    hdr_end = 4+nT0+len(P)+2*(nb+1)+2*ng
    tag = (' | hiT: T > 6100 K extended by gen_hitemp.py '
           '(carrier-abundance taper + H bf/ff), see HITEMP.md\n') if v1 else \
        (' | hiT2: T > 6100 K extended by gen_hitemp.py (carrier-abundance taper'
         ' + H bf/ff, C1 log blend 6100-7100 K), see HITEMP.md\n')
    with open(out_k, 'w') as f:
        f.write(kt['title'].rstrip('\n') + tag)
        f.write(f"{nT0+len(T_K_NEW)} {len(P)} {nb} {ng}\n")
        f.write(' '.join(raw[4:4+nT0]+[f"{t:.1f}" for t in T_K_NEW])+'\n')
        q = 4+nT0
        for n in (len(P), nb+1, nb+1):
            f.write(' '.join(raw[q:q+n])+'\n')
            q += n
        f.write(' \n')
        for n in (ng, ng):
            f.write(' '.join(raw[q:q+n])+'\n')
            q += n
        f.write(' \n')
        assert q == hdr_end
        body = raw[hdr_end:]
        for r in range(nT0*len(P)*nb):
            f.write(' '.join(body[ng*r:ng*r+ng])+'\n')
        for it in range(len(T_K_NEW)):
            for ip in range(len(P)):
                for b in range(nb):
                    f.write(' '.join(f"{v:.10e}" for v in Knew[it, ip, b])+'\n')
    np.savez(os.environ.get('CK_META', 'gen_meta.npz'), alpha=alpha, P=P, ratio=ratio)
    print('wrote', out_k, out_ce)
    print('alpha at p = 1, 10, 100, 230(~213), 1000 bar:')
    for Pv in (1., 10., 100., 213.8, 1000.):
        print(f"  {Pv:7g}", np.array2string(alpha[np.argmin(abs(P-Pv))], precision=2))
    print('CE ratio upstream/FastChem at 6100 K, min/max per column:',
          np.array2string(ratio.min(0), precision=3),
          np.array2string(ratio.max(0), precision=3))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--met', default='1x')
    ap.add_argument('--v1', action='store_true')
    a = ap.parse_args()
    main(a.met, a.v1)
