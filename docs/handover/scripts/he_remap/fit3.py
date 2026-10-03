# flake8: noqa
"""remap128/grid/fit3.py: fit2.py with the target from the pp2 70 ks profiles.\nhe_mltpp_1002 (copy of hepresn_ext_1001/grid/fit.py): caps also from the MLT++ IC and the 3-D
periodic textbook wedge at 17.6/23.5/29.4/35 ks (collapse3d/hp3d.npz).
Fit the code's StretchRPoly (8 poly coefs + 2 sech^2 bumps + tanh plateau) to a target
dr(r) for x1min 0.40 R, x1max 3 R, with dr <= Hp/5 (IC + beon128 first/last envelopes)."""
import numpy as np
from scipy.optimize import least_squares
import os
R = 2.3717e11
GROW = float(os.environ.get('GROW', '1.03'))
R0, R1 = 0.40 * R, 3.0 * R
FAC = float(os.environ.get('FAC', '0.85'))
DR1 = float(os.environ.get('DR1', '0.0018'))   # plateau dr/R, 0.40-0.55 and 0.85-RFINE
DR2 = float(os.environ.get('DR2', '0.0015'))   # photosphere/porous band 0.55-0.85 R
RFINE = float(os.environ.get('RFINE', '1.0'))  # fine to here, geometric growth beyond
# cap(r) = min over the 12 pp2 dumps 64-70 ks of the shell-mean TOTAL-pressure H_p / 5,
# running min over +-0.01 R (remap128/grid/prof70.npz)
q = np.load('/viper/ptmp2/jinma/he_mltpp_1002/remap128/grid/prof70.npz')
rr = np.linspace(R0, R1, 20001)
hmin = q['Ht'].min(axis=0)
cap = np.interp(rr, q['rv'], hmin / 5.0)
w = int(0.01 * R / (rr[1] - rr[0]))
capm = np.array([cap[max(0, i - w):i + w + 1].min() for i in range(rr.size)])
capm = np.minimum(capm, 0.012 * rr)
x = rr / R
tgt = np.minimum(FAC * capm, np.where((x > 0.55) & (x < 0.85), DR2 * R, DR1 * R))
grow = DR1 * R * np.exp(np.log(GROW) * np.maximum(0, (rr - RFINE * R)) / (DR1 * R * 1.6))
tgt = np.where(x > RFINE, np.minimum(grow, 0.012 * rr), tgt)
tgt = np.minimum(tgt, 0.97 * capm)
N_t = np.trapz(1.0 / tgt, rr)
print('target cells %.1f' % N_t)
NSTR = 8


def lch(z):
    a = np.abs(z)
    return a + np.log1p(np.exp(-2 * a)) - np.log(2.0)


def G(xa, xb, w, xi):
    return 0.5 * w * (lch((xi - xa) / w) - lch(-xa / w) -
                      lch((xi - xb) / w) + lch(-xb / w))


def ufun(p, xi):
    u = xi.copy()
    xik = xi.copy()
    for k in range(8):
        u = u + p[k] * xik * (1 - xi)
        xik = xik * xi
    for b in range(2):
        a, xb, w = p[8 + 3 * b:11 + 3 * b]
        if a != 0:
            u = u + a * w * (np.tanh((xi - xb) / w) - (1 - xi) *
                             np.tanh(-xb / w) - xi * np.tanh((1 - xb) / w))
    a, xa, xb, w = p[14:18]
    if a != 0:
        u = u + a * (G(xa, xb, w, xi) - xi * G(xa, xb, w, 1.0))
    return u


def faces(p, N):
    xi = np.linspace(0, 1, N + 1)
    return R0 + (R1 - R0) * ufun(p, xi)


if __name__ == '__main__':
    import sys
    import os
    N = int(sys.argv[1]) if len(sys.argv) > 1 else int(round(N_t))
    # target faces from the target density
    cum = np.concatenate(
        [[0], np.cumsum(0.5 * (1 / tgt[1:] + 1 / tgt[:-1]) * np.diff(rr))])
    xi_t = np.linspace(0, 1, N + 1)
    rf_t = np.interp(xi_t * cum[-1], cum, rr)

    def res(p):
        rf = faces(p, N)
        dr = np.diff(rf)
        if np.any(dr <= 0):
            return np.full(2 * N + N + 1 + N - 1, 1e3)
        rcn = 0.5 * (rf[1:] + rf[:-1])
        cp = np.interp(rcn, rr, capm)
        tg = np.interp(rcn, rr, tgt)
        r1 = np.log(dr / tg)
        r2 = float(os.environ.get('PW', '20')) * np.maximum(0, np.log(dr / (0.985 * cp)))
        r3 = 0.3 * (rf - rf_t) / (0.003 * R)
        lr = np.log(dr[1:] / dr[:-1])
        r4 = float(os.environ.get('SW', '50')) * np.maximum(0, np.abs(lr) -
                                                            np.log(float(os.environ.get('RMAX', '1.05'))))
        return np.concatenate([r1, r2, r3, r4])

    def xp(r): return (r - R0) / (R1 - R0)
    # initial: polynomial fit to the target u, plateau over 0.53..1.05 R
    p0 = np.zeros(18)
    p0[14:18] = [-0.6, -0.05, np.interp(RFINE * R, rf_t, xi_t), 0.01]
    p0[8:14] = [-0.05, np.interp(0.55 * R, rf_t, xi_t),
                0.03, -0.05, np.interp(0.85 * R, rf_t, xi_t), 0.02]
    if os.environ.get('P0'):
        p0 = np.load(os.environ['P0'])
    best = None
    for it in range(6):
        lb = np.full(18, -np.inf)
        ub = np.full(18, np.inf)
        for b in range(2):
            lb[8 + 3 * b], ub[8 + 3 * b] = -3.0, 3.0
            lb[9 + 3 * b], ub[9 + 3 * b] = 0.0, 1.0
            lb[10 + 3 * b], ub[10 + 3 * b] = 0.005, 0.5
        lb[14], ub[14] = -0.99, 0.0
        lb[15], ub[15] = -0.2, 1.0
        lb[16], ub[16] = 0.0, 1.0
        lb[17], ub[17] = 0.003, 0.2
        p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)
        sol = least_squares(
            res,
            p0,
            method='trf',
            max_nfev=4000,
            x_scale='jac',
            bounds=(
                lb,
                ub))
        if best is None or sol.cost < best.cost:
            best = sol
        p0 = best.x * (1 + 0.02 * np.random.default_rng(it).standard_normal(18))
    p = best.x
    rf = faces(p, N)
    dr = np.diff(rf)
    rcn = 0.5 * (rf[1:] + rf[:-1])
    cp = np.interp(rcn, rr, capm)
    print('N %d cost %.4g  max dr/cap %.4f at r/R %.3f' % (N, best.cost, (dr / cp).max(),
          rcn[np.argmax(dr / cp)] / R))
    np.save('p_%d.npy' % N, p)
    for xx in (
            0.40,
            0.42,
            0.45,
            0.48,
            0.5,
            0.505,
            0.51,
            0.52,
            0.55,
            0.58,
            0.6,
            0.65,
            0.7,
            0.75,
            0.8,
            0.85,
            0.9,
            1.0,
            1.05,
            1.1,
            1.2,
            1.5,
            2.0,
            2.5,
            2.99):
        i = np.argmin(abs(rcn - xx * R))
        print('  r/R %.2f dr/R %.5f dr/r %.5f  cap/R %.5f tgt/R %.5f' % (
            xx, dr[i] / R, dr[i] / rcn[i], cp[i] / R, np.interp(rcn[i], rr, tgt) / R))
    print('p =', ' '.join('%.12e' % v for v in p))
