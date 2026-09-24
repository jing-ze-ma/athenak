"""Step-to-step relaxation of the lagged closure: the closure perturbations USED in step
n+1 are c_{n+1} = (1-w) c_n + w c(state_n).  Augmented map on (e, p, q, df_used, q_used)
(q_used = the dP12 functional).  Prints max |eig| over modes for w = 1/(1+G^2) with
G = cg/tau, and the plain lagged gain, nu = 1e4, r = 1."""
import numpy as np
from fourier import chi_m1, dchi


def aug(al, be, f0, tau, nu, w, r=1.0):
    chi, cp = chi_m1(f0), dchi(chi_m1, f0)
    a, b = (1 - chi) / 2, (3 * chi - 1) / 2
    sa, ca, sb, cb = np.sin(al / 2), np.cos(al / 2), np.sin(be / 2), np.cos(be / 2)
    th = 1 / (1 + nu * tau)
    # new (e,p,q) from old (e,p,q) and USED closure (d, s): d = df, s = dP12
    A = np.zeros((3, 3), complex)
    Bs = np.zeros((3, 3), complex)
    Cc = np.zeros((3, 2), complex)
    A[0, 0] = th * 2j * nu * sa * chi
    A[0, 1] = 1
    Bs[0, 1] = th
    Cc[0, 0] = -th * 2j * nu * sa * cp
    Cc[0, 1] = -th * 1j * (nu / r) * ca * np.sin(be)
    A[1, 0] = th * 2j * (nu / r) * sb * a
    A[1, 2] = 1
    Bs[1, 2] = th
    Cc[1, 0] = th * 2j * (nu / r) * sb * cp / 2
    Cc[1, 1] = -th * 1j * nu * cb * np.sin(al)
    A[2, 0] = 1
    A[2, 1] = 2j * nu * sa
    A[2, 2] = 2j * (nu / r) * sb
    Bs[2, 0] = 1
    Ai = np.linalg.inv(A)
    # closure from a state X: d = [-f0, ca, 0].X, s = [0, 0, (b/f0) cb].X
    L = np.array([[-f0, ca, 0], [0, 0, (b / f0) * cb]], complex)
    # step: used' = (1-w) used + w L X ; X' = Ai (Bs X + Cc used')
    M = np.zeros((5, 5), complex)
    U = (1 - w) * np.eye(2)
    M[3:, :3] = w * L
    M[3:, 3:] = U
    M[:3, :3] = Ai @ (Bs + Cc @ (w * L))
    M[:3, 3:] = Ai @ Cc @ U
    return M


def gm(f0, tau, nu, w, n=24):
    return max(np.max(np.abs(np.linalg.eigvals(aug(al, be, f0, tau, nu, w))))
               for al in np.linspace(0, np.pi, n + 1) for be in np.linspace(0, np.pi, n + 1)
               if (al, be) != (0, 0))


if __name__ == '__main__':
    for f0 in (0.3, 0.5):
        for tau in (0.01, 0.03, 0.0625, 0.125, 0.25, 0.5, 1.0):
            g1 = gm(f0, tau, 1e4, 1.0)
            row = []
            for cg in (0.5, 1.0):
                G = cg / tau
                w = min(1.0, 1.0 / (1.0 + G * G) * (1 + 1.0 / (cg * cg)))
                row.append('cg=%.1f w=%.2e g=%.5f' % (cg, w, gm(f0, tau, 1e4, w)))
            print('f0=%.1f tau=%.4f lag g=%.3f | ' % (f0, tau, g1) + ' | '.join(row))
