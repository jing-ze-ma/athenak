"""weight used in the code: G = C max(chi', b/f)/tau + f chi'/min(chi, 1-chi),
w = min(1, 2/(1+G^2)); max |eig| of the relaxed map (m1 closure), and the plain lag."""
import numpy as np
from fourier import chi_m1, dchi
from relax import gm


def weight(f0, tau, C):
    chi, cp = chi_m1(f0), dchi(chi_m1, f0)
    b = (3 * chi - 1) / 2
    G = C * max(cp, b / f0) / tau + f0 * cp / min(chi, 1 - chi)
    return min(1.0, 2.0 / (1.0 + G * G))


if __name__ == '__main__':
    for C in (1.0, 1.5):
        for f0 in (0.1, 0.3, 0.5, 0.65, 0.8, 0.95):
            row = []
            for tau in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0):
                w = weight(f0, tau, C)
                row.append('%.4f(%.0e)' % (gm(f0, tau, 1e4, w, n=16), w))
            print('C=%.1f f0=%.2f: ' % (C, f0) + ' '.join(row))
