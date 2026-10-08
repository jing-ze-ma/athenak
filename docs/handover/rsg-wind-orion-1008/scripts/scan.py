"""(1) Warm/cool stable roots vs the thermalisation RATIO R = eps_el/eps_IR (eps_IR fixed,
continuum eps = 1 kept) -- the roots depend on eps only through R and the continuum share.
(2) t_th / t_dyn map over q x A_IR, (a) scattering electronic bands, isobaric.
Run: TAB=lowP|old python scan.py  -> scan_<TAB>.txt"""
import os
import pickle
import sys
import numpy as np
sys.dont_write_bytecode = True
import nlte as N   # noqa: E402  (TAB env selects the table)
L = N.L
RS = (0., 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.2, 0.3, 0.5, 1., 3., 10., 100.)
EIR = (1., 1e-2, 1e-4)
out = []
d = N.tables()
for s in N.STARS:
    st = L.star(s)
    for x in (2., 3., 4., 5.):
        W = N.dilution(x)
        for iv, P in enumerate(N.PBAR):
            kl, kc = d['P_kl'][iv], d['P_kc'][iv]
            for e in EIR:
                row = []
                for R in RS:
                    eps = np.full((len(N.TG), L.NB), e)
                    eps[:, N.IEL] = min(R*e, 1.)
                    H = N.heat(kl, kc, eps, W, st['Teff'])[0]
                    stv = sorted([r[0] for r in N.roots(H) if r[1]], reverse=True)
                    row.append('/'.join(f'{t:.0f}' for t in stv) or 'none')
                out.append(f'{s} x={x:.0f} P={P:.0e} eps_IR={e:.0e} | ' +
                           ' | '.join(f'R={R:g}: {r}' for R, r in zip(RS, row)))
# t_th map
res = pickle.load(open(N.OUT + f'res_{N.TAB}.pkl', 'rb'))
R0 = {'golden16': 669.875*L.RSUN, 'm20lgl5.5': 1311.634*L.RSUN}
VC = {'golden16': 7.46e5, 'm20lgl5.5': 8.298e5}
out.append('\nt_th/min(t_flow(30 km/s), t_sh) at the cool stable root, (a), A_el=1e6; '
           'columns q=1e-12,1e-11,1e-10 for A_IR=10 | A_IR=100')
for r0 in res:
    if r0['key'] != 'P' or r0['scen'] != 'lte':
        continue
    s, x, P = r0['star'], r0['x'], r0['val']
    td = min(x*R0[s]/3e6, R0[s]/VC[s])
    cells = []
    for air in (10., 100.):
        for q in N.QS:
            name = f'a q{q:.0e} AIR{air:.0f} Ael1e+06'
            r = [z for z in res if (z['star'], z['x'], z['key'], z['val'], z['scen'])
                 == (s, x, 'P', P, name)][0]
            stv = [z for z in r['roots'] if z['stable']]
            cells.append(f"{min(stv, key=lambda z: z['T'])['tth']/td:.1e}" if stv else '-')
    lte = [z for z in r0['roots'] if z['stable']]
    out.append(f"{s} x={x:.0f} P={P:.0e} t_dyn={td:.1e} s | LTE max t_th/t_dyn "
               f"{max(z['tth'] for z in lte)/td:.1e} | " + ' '.join(cells[:3]) + ' | '
               + ' '.join(cells[3:]))
open(N.OUT + f'scan_{N.TAB}.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
