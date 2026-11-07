#!/usr/bin/env python3
"""RSG E page (11-01): Lomb-Scargle power spectrum of the light curve, relaxed phase only. Called from mk_figs.py on every
refresh with arrays it already holds (hst time / L_top / L_in, dump times / median R_ph, KE_int); numpy only (no astropy in
the venv): generalised LS with a floating mean (Zechmeister & Kuerster 2009), 'standard' normalisation, Baluev (2008) FAP as
in astropy's LombScargle.false_alarm_probability(method='baluev').  ~0.1 s for a 1 d binned series of a few thousand days."""
import numpy as np
from scipy.special import gammaln
from scipy.optimize import brentq

DAY = 86400.
BIN = 1.0       # d, the light curve is binned to daily means (hst rows every ~1000 s): anti-alias, Nyquist 2 d << 20 d
WRUN = 100      # d, boxcar of the running means used for the relaxed-phase test
KSIG = 3.0      # running mean must stay within KSIG x (spread of the running mean in the reference half)
GRID = 50.      # d, t_start is rounded up to this grid
PMIN = 20.      # d, shortest period plotted / searched
OVERS = 10      # frequency oversampling per 1/T


def daily(t_d, y, t0=0.):
    """daily means of y on [t0, ...): returns bin centres, means (bins without hst rows are interpolated)"""
    b = np.floor((t_d - t0) / BIN).astype(int)
    m = b >= 0
    b, yy = b[m], y[m]
    n = b.max() + 1
    c = np.bincount(b, minlength=n)
    s = np.bincount(b, yy, minlength=n)
    tc = t0 + (np.arange(n) + 0.5) * BIN
    ok = c > 0
    ym = np.interp(tc, tc[ok], s[ok] / c[ok])
    if c[-1] < 0.5 * np.median(c[ok]):   # drop a partly filled last day
        tc, ym = tc[:-1], ym[:-1]
    return tc, ym


def running(tc, y, w=WRUN):
    n = int(round(w / BIN))
    cs = np.cumsum(np.insert(y, 0, 0.))
    return tc[:len(y) - n + 1] + 0.5 * (w - BIN), (cs[n:] - cs[:-n]) / n


def relaxed_start(tc, series):
    """earliest t after which the WRUN-d running mean of EVERY series stays within KSIG sigma of its reference level;
    reference level / sigma = median / std of the running mean over the second half of the run (t > t_end/2); only
    excursions before t_end/2 count (the reference half defines 'stationary').  Returns t_start (rounded up to GRID) and
    per-series details."""
    th = tc[-1] / 2.
    det = {}
    tl = 0.
    for nm, y in series.items():
        tr, rm = running(tc, y)
        late = tr > th
        ref, sd = np.median(rm[late]), np.std(rm[late])
        out = (np.abs(rm - ref) > KSIG * sd) & (tr <= th)
        last = tr[out].max() + 0.5 * WRUN if out.any() else 0.
        det[nm] = dict(ref=float(ref), sd=float(sd), last=float(last))
        tl = max(tl, last)
    return float(np.ceil(tl / GRID) * GRID), det


def trend_t(t, y, nb=100.):
    """t-statistic of a linear trend fitted to the nb-day block means (red-noise-safe-ish significance)"""
    k = int((t[-1] - t[0]) // nb)
    if k < 4:
        return np.nan
    bc = t[0] + (np.arange(k) + 0.5) * nb
    bm = np.array([y[(t >= t[0] + i * nb) & (t < t[0] + (i + 1) * nb)].mean() for i in range(k)])
    p = np.polyfit(bc, bm, 1)
    sd = np.std(bm - np.polyval(p, bc), ddof=2)
    return float(p[0] / (sd / np.sqrt(((bc - bc.mean()) ** 2).sum())))


def gls(t, y, f):
    """generalised LS (floating mean), standard normalisation P = 1 - chi2(f)/chi2_const, unit weights"""
    w = np.full(len(t), 1. / len(t))
    Y = (w * y).sum()
    YY = (w * y * y).sum() - Y * Y
    P = np.empty(len(f))
    for i0 in range(0, len(f), 512):   # chunks: N x 512 doubles
        ff = f[i0:i0 + 512]
        a = 2 * np.pi * np.outer(ff, t)
        c, s = np.cos(a), np.sin(a)
        C, S = c @ w, s @ w
        YC = (c * y) @ w - Y * C
        YS = (s * y) @ w - Y * S
        CC = (c * c) @ w - C * C
        SS = (s * s) @ w - S * S
        CS = (c * s) @ w - C * S
        D = CC * SS - CS * CS
        P[i0:i0 + 512] = (SS * YC * YC + CC * YS * YS - 2 * CS * YC * YS) / (YY * D)
    return P


def _gam(n):
    return np.sqrt(2. / n) * np.exp(gammaln(n / 2.) - gammaln((n - 1) / 2.))


def fap_baluev(z, t, fmax):
    """astropy _statistics.fap_baluev, standard normalisation, dH = 1, dK = 3"""
    N = len(t)
    NH, NK = N - 1, N - 3
    z = np.clip(z, 0., 1. - 1e-15)
    fs = 1 - (1 - z) ** (0.5 * NK)          # 1 - FAP_single
    W = fmax * np.sqrt(4 * np.pi * np.var(t))
    tau = W * _gam(NH) * (1 - z) ** (0.5 * (NK - 1)) * np.sqrt(0.5 * NH * z)
    return 1 - fs * np.exp(-tau)


def fal_baluev(p, t, fmax):
    return brentq(lambda z: fap_baluev(z, t, fmax) - p, 1e-12, 1 - 1e-12)


def spectrum(tc, y, detrend=True):
    yy = y - y.mean()
    tt = tc - tc[0]
    lin = np.polyfit(tt, yy, 1)
    if detrend:
        yy = yy - np.polyval(lin, tt)
    T = tt[-1] - tt[0] + BIN
    fmin, fmax = 2. / T, 1. / PMIN
    f = np.arange(fmin, fmax, 1. / (OVERS * T))
    return f, gls(tt, yy, f), T, lin, fmax


def peaks(f, P, k=3):
    """top-k local maxima of P, separated by at least one resolution element"""
    i = np.where((P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]))[0] + 1
    i = i[np.argsort(P[i])[::-1]]
    sel = []
    for j in i:   # one resolution element 1/T = OVERS grid steps
        if all(abs(j - s) > OVERS for s in sel):
            sel.append(j)
        if len(sel) == k:
            break
    return sel


def window(f, tc, f0):
    """spectral window |sum exp(2 pi i (f - f0) t)|^2 / N^2 centred on f0"""
    a = 2 * np.pi * np.outer(f - f0, tc - tc[0])
    return (np.cos(a).sum(1) ** 2 + np.sin(a).sum(1) ** 2) / len(tc) ** 2
