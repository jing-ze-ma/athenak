"""the implemented weight: w = min(1, 2/(1+(C/tau)^2)); max |eig| of the relaxed map"""
from relax import gm
for C in (1.0, 2.0):
    for f0 in (0.3, 0.5, 0.65, 0.8, 0.95):
        row = []
        for tau in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0):
            w = min(1.0, 2.0 / (1.0 + (C / tau) ** 2))
            row.append('%.4f' % gm(f0, tau, 1e4, w, n=16))
        print('C=%.0f f0=%.2f  tau=.003..3: ' % (C, f0) + ' '.join(row))
