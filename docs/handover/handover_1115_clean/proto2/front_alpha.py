"""front_alpha.py: dissipation weight alpha at the steepest face (the front) after 1, 4 steps
for rayflux, A3 and S1 (beam into vacuum, pulse in chi dx 0.1 scattering, slab chi dx 10)."""
import numpy as np
import scheme as S
S.SLIN = True
n = 100
x = (np.arange(n) + 0.5)/n
for name, opt in (('rayflux', dict(cand='rayflux', recon=None, npic=1, Dguard=True)),
                  ('A3', dict(cand='A3', recon=None, npic=1, Dguard=True)),
                  ('S1', dict(cand='S1', recon=None, npic=1, Dguard=True))):
    for case in ('beam vac', 'pulse 0.1', 'slab 10'):
        row = []
        for cfl in (0.5, 4.0, 100.0, 1000.0):
            if case == 'beam vac':
                pb = S.Problem(n, bc='inflow'); E = np.zeros(n); F = np.zeros(n)
            elif case == 'slab 10':
                pb = S.Problem(n, ks=1000.0*(x > 0.5), bc='inflow'); E = np.zeros(n); F = np.zeros(n)
            else:
                pb = S.Problem(n, ks=10.0); m = (x > 0.2) & (x < 0.4)
                E = 1.0*m; F = 1.0*m
            U = np.r_[E, F]
            al = []
            for k in range(4):
                U, _ = S.step(pb, U, cfl/n, opt)
                e = np.r_[U[:n], U[n - 1]] if pb.bc != 'periodic' else U[:n]
                d = np.abs(np.diff(np.r_[e[-1] if pb.bc == 'periodic' else pb.ghost[0], e]))
                j = np.argmax(d[:len(S.RAYX['alpha'])])
                if k in (0, 3):
                    al.append(S.RAYX['alpha'][j])
            row.append('cfl %-5g %.2f/%.2f' % (cfl, al[0], al[1]))
        print('%-8s %-10s alpha at front (step1/step4): %s' % (name, case, '  '.join(row)))
