"""Python mirror of src/utils/correlated_k.hpp: k-table, FastChem composition, CIA,
Rayleigh, H- (John 1988), band Planck fractions (with the blue/red tail fold), and the
Rosseland/Planck means exactly as ck_build_rosseland_table builds them."""
import numpy as np
import os
DATA = os.environ.get('CK_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                              '..')).rstrip('/') + '/'
KB = 1.380649e-16
SIG = 5.670374419e-5


def read_ktable(fn=DATA+'ck/Premixed_1x_g8_11.txt'):
    with open(fn) as f:
        title = f.readline()
        tok = f.read().split()
    nT, nP, nb, ng = map(int, tok[:4])
    p = 4
    T = np.array(tok[p:p+nT], float)
    p += nT
    P = np.array(tok[p:p+nP], float)
    p += nP
    wl = np.array(tok[p:p+nb+1], float)
    p += nb+1
    wn = np.array(tok[p:p+nb+1], float)
    p += nb+1
    gx = np.array(tok[p:p+ng], float)
    p += ng
    gw = np.array(tok[p:p+ng], float)
    p += ng
    K = np.array(tok[p:p+nT*nP*nb*ng], float).reshape(nT, nP, nb, ng)
    assert p+nT*nP*nb*ng == len(tok)
    return dict(title=title, T=T, P=P, wl=wl, wn=wn, gx=gx, gw=gw, K=K, raw=tok)


def read_ce(fn=DATA+'CE_tables/FastChem_ck_1x_int.txt'):
    tok = open(fn).read().split()
    nT, nP, nrec, nsp = map(int, tok[:4])
    p = 4+nsp
    T = np.array(tok[p:p+nT], float)
    p += nT
    P = np.array(tok[p:p+nP], float)
    p += nP
    C = np.array(tok[p:p+nT*nP*6], float).reshape(nT, nP, 6)
    return T, P, C


def read_cia():
    out = []
    for s in ('H2-H2', 'H2-He', 'H2-H', 'He-H'):
        tok = open(DATA+f'cia/{s}_reform_11.txt').read().split()
        nT, nb = int(tok[0]), int(tok[1])
        p = 2
        T = np.array(tok[p:p+nT], float)
        p += nT + nb
        k = np.array(tok[p:p+nT*nb], float).reshape(nT, nb)
        out.append((T, k))
    return out


def read_ray():
    out = []
    for s in ('H2', 'He', 'H', 'e-'):
        tok = open(DATA+f'ray/Ray_{s}_11.txt').read().split()
        out.append(np.array(tok[1:12], float))
    return np.array(out)


def planck_below(lamT):
    lamT = np.atleast_1d(np.asarray(lamT, float))
    xi = 1.4387769e4/lamT
    n = np.arange(1, 501)[:, None]
    nx = n*xi
    s = np.where(nx < 700, np.exp(-np.minimum(nx, 700)) /
                 n*(xi**3+3*xi**2/n+6*xi/n**2+6/n**3), 0)
    return 15/np.pi**4*s.sum(0)


def planck_frac(T, wl):
    T = np.atleast_1d(T)
    F = np.array([planck_below(w*T) for w in wl])   # (nb+1, nT)
    f = F[:-1]-F[1:]
    f[0] += 1-F[0]
    f[-1] += F[-1]
    return f.T   # (nT, nb)


def interp_lin(x, xg, y):   # clamped linear, like ck_tp_index
    i = np.clip(np.searchsorted(xg, x)-1, 0, len(xg)-2)
    fr = np.clip((x-xg[i])/(xg[i+1]-xg[i]), 0, 1)
    return i, fr


def comp(ce, T, P):
    cT, cP, C = ce
    i, a = interp_lin(np.log10(T), np.log10(cT), None)
    j, b = interp_lin(np.log10(P), np.log10(cP), None)
    return (1-a)*((1-b)*C[i, j]+b*C[i, j+1]) + a*((1-b)*C[i+1, j]+b*C[i+1, j+1])


def continuum(ce, cia, ray, wl, T, P):
    v = comp(ce, T, P)
    mu = v[0]
    ntot = P*1e6/(KB*T)
    rho = P*1e6*mu*1.6726e-24/(KB*T)
    kc = np.zeros(len(wl)-1)
    for (Tg, k), (a, b) in zip(cia, ((1, 1), (1, 2), (1, 3), (2, 3))):
        if T <= Tg[0]:
            it, ft = 0, 0.
        elif T >= Tg[-1]:
            it, ft = len(Tg)-2, 1.
        else:
            it = np.searchsorted(Tg, T)-1
            ft = (T-Tg[it])/(Tg[it+1]-Tg[it])
        kc += ((1-ft)*k[it]+ft*k[it+1])*v[a]*ntot*v[b]*ntot/rho
    for s in range(4):
        kc += ray[s]*v[s+1]*ntot/rho
    kc += hminus(wl, T, v[5]*ntot, v[4]*ntot*v[3]*ntot*KB*T)/rho
    return kc, rho, v


def hminus_xsec(lam, T):
    """John (1988) H- bf [cm^2 per H-] and ff [cm^4/dyn per (Pe nH)] at lam [um]."""
    lam0 = 1.6419
    Cbf = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]
    # rows A..F, columns n = 1..6 (set 1: 0.1823-0.3645 um, set 2: longward)
    A1 = [[518.1021, 472.2636, -482.2089, 115.5291, 0, 0],
          [-734.8666, 1443.4137, -737.1616, 169.6374, 0, 0],
          [1021.1775, -1977.3395, 1096.8827, -245.6490, 0, 0],
          [-479.0721, 922.3575, -521.1341, 114.2430, 0, 0],
          [93.1373, -178.9275, 101.7963, -21.9972, 0, 0],
          [-6.4285, 12.3600, -7.0571, 1.5097, 0, 0]]
    A2 = [[0, 2483.3460, -3449.8890, 2200.0400, -696.2710, 88.2830],
          [0, 285.8270, -1158.3820, 2427.7190, -1841.4000, 444.5170],
          [0, -2054.2910, 8746.5230, -13651.1050, 8642.9700, -1863.8640],
          [0, 2827.7760, -11485.6320, 16755.5240, -10051.5300, 2095.2880],
          [0, -1341.5370, 5303.6090, -7510.4940, 4400.0670, -901.7880],
          [0, 208.9520, -812.9390, 1132.7380, -655.0200, 132.9850]]
    xbf = 0.
    if lam < lam0:
        dk = 1/lam-1/lam0
        xbf = 1e-18*lam**3*dk**1.5*sum(Cbf[n]*dk**(n/2) for n in range(6))
    sff = 0.
    if lam > 0.1823:
        A = A2 if lam >= 0.3645 else A1
        th = 5040./T
        sff = sum(th**((n+2)/2)*(lam**2*A[0][n] + A[1][n] + A[2][n]/lam + A[3][n]/lam**2
                                 + A[4][n]/lam**3 + A[5][n]/lam**4) for n in range(6))
    return xbf, 1e-29*sff


def hminus(wl, T, nHm, Pe_nH):
    out = np.zeros(len(wl)-1)
    for b in range(len(wl)-1):
        lam = 2/(1/wl[b]+1/wl[b+1])
        xbf, xff = hminus_xsec(lam, T)
        out[b] = xbf*nHm + xff*Pe_nH
    return out


def means(kt, ce, cia, ray, T, P, lineK=None):
    """kappa_R (code definition), kappa_P, per-band g-mean k_line and continuum at
    (T, P).
    lineK: (nb, ng) line k at this point; default = bilinear log interp of kt
    (clamped)."""
    wl = kt['wl']
    gw = kt['gw']
    if lineK is None:
        i, a = interp_lin(np.log10(T), np.log10(kt['T']), None)
        j, b = interp_lin(np.log10(P), np.log10(kt['P']), None)
        L = np.log10(np.maximum(kt['K'], 1e-99))
        lk = (1-a)*((1-b)*L[i, j]+b*L[i, j+1]) + a*((1-b)*L[i+1, j]+b*L[i+1, j+1])
        lineK = 10**lk
    kc, rho, v = continuum(ce, cia, ray, wl, T, P)
    fp = planck_frac(np.array([1.01*T, 0.99*T]), wl)
    wb = SIG*(1.01*T)**4*fp[0] - SIG*(0.99*T)**4*fp[1]
    ktot = lineK + kc[:, None]
    inv = (gw/ktot).sum(1)
    ok = wb > 0
    kR = wb[ok].sum()/(wb[ok]*inv[ok]).sum()
    fP = planck_frac(T, wl)[0]
    kP = (fP*(gw*ktot).sum(1)).sum()
    return dict(kR=kR, kP=kP, kline=(gw*lineK).sum(1), kc=kc, rho=rho, v=v)
