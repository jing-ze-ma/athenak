"""Von Neumann analysis of the multi-D implicit M1 step (rad_m1_implicit.cpp, bicgstab,
offdiag operator/lagged) about a uniform state E0 = 1, F = f0 x^ (c = chat = 1), frozen
coefficients, rho kappa_t = K, no absorption/emission, v = 0.

Unknowns per Fourier mode (alpha = kx dx, beta = ky dy): e (cell E), p (x1-face F1),
q (x2-face F2), each at its own position.  The step (lagged closure = one step lag,
implicit_closure_lag = step; = one Picard pass under `pass`):
  p = th [ p^n - 2i nu sa dP11 - i (nu/r) ca sin(beta) dP12 ]
  q = th [ wm q^n - 2i (nu/r) sb dP22 - i nu cb sin(alpha) dP12 ]
  e = e^n - 2i nu sa p - 2i (nu/r) sb q
  dP11 = chi e + chi' df,  dP22 = a e - chi'/2 df,  dP12 = (b/f0) cb q_c
  df = ca p_c - f0 e      (p_c, q_c: the cell F from the face means ca p, cb q)
th = 1/(1 + nu tau), nu = c dt/dx, tau = K dx, r = dy/dx, sa = sin(alpha/2), ca = cos(alpha/2).
closure terms (df, q_c) from the OLD state (lag) or the NEW state (implicit = the fixed
point of converged Picard, i.e. nonlinear backward Euler).
"""
import numpy as np


def chi_m1(f):   # Levermore / Minerbo-type M1 closure used by the code (M1Chi)
    f2 = f * f
    return (3 + 4 * f2) / (5 + 2 * np.sqrt(4 - 3 * f2))


def chi_ker(f):
    return (1 + 2 * f * f) / 3


def dchi(ch, f, h=1e-6):
    return (ch(f + h) - ch(f - h)) / (2 * h)


def gmat(al, be, f0, tau, nu, r=1.0, wm=1.0, closure='m1', mode='lag', od=True,
         dchi_on=True, ups=0.0):
    """ups: extra 'streaming opacity' on the x2 faces (cure c3 variant), in units 1/dx."""
    ch = {'m1': chi_m1, 'ker': chi_ker, 'edd': lambda f: 0 * f + 1 / 3}[closure]
    chi = ch(f0)
    cp = dchi(ch, f0) if dchi_on else 0.0
    a, b = (1 - chi) / 2, (3 * chi - 1) / 2
    sa, ca, sb, cb = np.sin(al / 2), np.cos(al / 2), np.sin(be / 2), np.cos(be / 2)
    th = 1 / (1 + nu * tau)
    thq = 1 / (1 + nu * tau + nu * ups)
    od = 1.0 if od else 0.0
    # closure perturbations as linear functionals of a state (e,p,q)
    dfv = np.array([-f0, ca, 0.0], complex)            # df
    q12 = np.array([0.0, 0.0, od * (b / f0) * cb], complex)   # dP12
    # rows: residual R(X_new; Y_old) = A X - B Y = 0
    A = np.zeros((3, 3), complex)
    B = np.zeros((3, 3), complex)
    C = B if mode == 'lag' else A          # where closure terms go (sign handled below)
    sgn = -1.0 if mode == 'lag' else 1.0   # moved to the rhs: B gets -coef
    # p row: p + th*2i nu sa (chi e + cp df) + th*i(nu/r) ca sin(be) dP12 = th p^n
    A[0, 0] += th * 2j * nu * sa * chi
    A[0, 1] += 1.0
    B[0, 1] += th
    C[0] += sgn * (th * 2j * nu * sa * cp * dfv + th * 1j * (nu / r) * ca * np.sin(be) * q12)
    # q row
    A[1, 0] += thq * 2j * (nu / r) * sb * a
    A[1, 2] += 1.0
    B[1, 2] += thq * wm
    C[1] += sgn * (thq * 2j * (nu / r) * sb * (-cp / 2) * dfv
                   + thq * 1j * nu * cb * np.sin(al) * q12)
    # e row
    A[2, 0] += 1.0
    A[2, 1] += 2j * nu * sa
    A[2, 2] += 2j * (nu / r) * sb
    B[2, 0] += 1.0
    return np.linalg.solve(A, B)


def gmax(f0, tau, nu, n=48, **kw):
    best, arg = 0.0, None
    for al in np.linspace(0, np.pi, n + 1):
        for be in np.linspace(0, np.pi, n + 1):
            if al == 0 and be == 0:
                continue
            g = np.max(np.abs(np.linalg.eigvals(gmat(al, be, f0, tau, nu, **kw))))
            if g > best:
                best, arg = g, (al, be)
    return best, arg


def pde_rate(kx, ky, f0, K, closure='m1'):
    """max Re(lambda) of the linearised continuous M1 system (c = 1) with opacity K"""
    ch = {'m1': chi_m1, 'ker': chi_ker}[closure]
    chi, cp = ch(f0), dchi(ch, f0)
    a, b = (1 - chi) / 2, (3 * chi - 1) / 2
    # U = (E, F1, F2); dE/dt = -i(kx F1 + ky F2); dF/dt = -i k.dP - K F
    # dP11 = chi E + cp df, dP22 = a E - cp/2 df, dP12 = (b/f0) F2, df = F1 - f0 E
    M = np.zeros((3, 3), complex)
    M[0, 1], M[0, 2] = -1j * kx, -1j * ky
    dP11 = np.array([chi - cp * f0, cp, 0])
    dP22 = np.array([a + cp / 2 * f0, -cp / 2, 0])
    dP12 = np.array([0, 0, b / f0])
    M[1] = -1j * (kx * dP11 + ky * dP12)
    M[2] = -1j * (kx * dP12 + ky * dP22)
    M[1, 1] -= K
    M[2, 2] -= K
    return np.max(np.linalg.eigvals(M).real)


if __name__ == '__main__':
    import sys
    for f0 in (0.2, 0.5, 0.8, 0.95):
        for tau in (0.05, 0.125, 0.5, 1.0, 2.0):
            row = []
            for mode in ('lag', 'imp'):
                g, arg = gmax(f0, tau, 1e6, mode=mode)
                row.append('%s %.3f @(%.2f,%.2f)' % (mode, g, arg[0], arg[1]))
            print('f0=%.2f tau=%.3f  ' % (f0, tau) + '  '.join(row))
