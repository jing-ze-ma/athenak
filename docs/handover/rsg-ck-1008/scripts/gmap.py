"""(a) Thin-limit Gamma map: Gamma_F (B_nu(Teff) flux spectrum), Gamma_P, Gamma_R (ck),
Gamma_R (OPLIB+AESOPUS table of the red_giant runs).  usage: python gmap.py"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import rsglib as L  # noqa: E402

Tg = np.arange(1500., 6001., 250.)
lrg = np.arange(-17., -8.99, 0.5)
out = {}
for mode in ('clamp', 'extrap'):
    kF = {s: np.zeros((len(Tg), len(lrg))) for s in ('golden16', 'betelgeuse')}
    kP = np.zeros((len(Tg), len(lrg)))
    kR = np.zeros_like(kP)
    kA = np.zeros_like(kP)
    Pb = np.zeros_like(kP)
    for i, T in enumerate(Tg):
        for j, lr in enumerate(lrg):
            st = L.state(T, 10**lr, mode)
            kR[i, j], kP[i, j] = L.means(T, st['K'])
            for s in kF:
                kF[s][i, j] = L.kappa_thin_flux(st['K'], L.star(s)['Teff'])
            kA[i, j] = L.kappa_aesopus(T, 10**lr)[0]
            Pb[i, j] = st['P']*1e-6
    out[mode] = dict(kF=kF, kP=kP, kR=kR, kA=kA, Pb=Pb)
sv = {}
for m in out:
    for k, v in out[m].items():
        if k == 'kF':
            for s_ in v:
                sv[f'{m}_kF_{s_}'] = v[s_]
        else:
            sv[f'{m}_{k}'] = v
np.savez('gmap.npz', Tg=Tg, lrg=lrg, **sv)

kEg = L.star('golden16')['kE']
kEb = L.star('betelgeuse')['kE']
with open('gmap_table.md', 'w') as f:
    for m in ('clamp', 'extrap'):
        o = out[m]
        f.write(f'\n### Gamma_F thin limit, golden16 (Teff {L.star("golden16")["Teff"]:.0f} K'
                f' spectrum, kE {kEg:.3f}), mode {m}\n\n')
        sel = [0, 2, 4, 6, 8, 10, 12, 14, 16]
        f.write('| T [K] \\ log rho | ' + ' | '.join(f'{lrg[j]:.0f}' for j in sel) + ' |\n')
        f.write('|---' * (len(sel)+1) + '|\n')
        for i in range(0, len(Tg), 2):
            f.write(f'| {Tg[i]:.0f} | ' + ' | '.join(
                f'{o["kF"]["golden16"][i, j]/kEg:.2g}' for j in sel) + ' |\n')
    o = out['clamp']
    for nm, arr, kE in (('Gamma_P (gas-T Planck, ck)', o['kP'], kEg),
                        ('Gamma_R (ck Rosseland)', o['kR'], kEg),
                        ('Gamma_R (OPLIB+AESOPUS table of the red_giant runs)', o['kA'], kEg),
                        ('Gamma_F thin, Betelgeuse set (Teff 3600 K, kE %.3f)' % kEb,
                         o['kF']['betelgeuse'], kEb)):
        f.write(f'\n### {nm}, golden16 kE unless stated, mode clamp\n\n')
        sel = [0, 4, 8, 12, 16]
        f.write('| T [K] \\ log rho | ' + ' | '.join(f'{lrg[j]:.0f}' for j in sel) + ' |\n')
        f.write('|---' * (len(sel)+1) + '|\n')
        for i in range(0, len(Tg), 2):
            f.write(f'| {Tg[i]:.0f} | ' + ' | '.join(f'{arr[i, j]/kE:.2g}' for j in sel)
                    + ' |\n')
    f.write('\n### P [bar] of the grid points (clamp mode; table edge 1e-8 bar)\n\n')
    sel = [0, 4, 8, 12, 16]
    f.write('| T [K] \\ log rho | ' + ' | '.join(f'{lrg[j]:.0f}' for j in sel) + ' |\n')
    f.write('|---' * (len(sel)+1) + '|\n')
    for i in (0, 4, 8, 12, 16, 18):
        f.write(f'| {Tg[i]:.0f} | ' + ' | '.join(f'{o["Pb"][i, j]:.1e}' for j in sel)
                + ' |\n')

fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
ext = [lrg[0]-0.25, lrg[-1]+0.25, Tg[0]-125, Tg[-1]+125]
panels = [('Gamma_F thin, golden16, clamp', out['clamp']['kF']['golden16']/kEg),
          ('Gamma_F thin, golden16, extrap (P<1e-8 bar)',
           out['extrap']['kF']['golden16']/kEg),
          ('Gamma_F thin, Betelgeuse set, clamp', out['clamp']['kF']['betelgeuse']/kEb),
          ('Gamma_P (gas T), ck, clamp', out['clamp']['kP']/kEg),
          ('Gamma_R, ck, clamp', out['clamp']['kR']/kEg),
          ('Gamma_R, OPLIB+AESOPUS (runs)', out['clamp']['kA']/kEg)]
for a, (t, z) in zip(ax.flat, panels):
    im = a.imshow(np.log10(z), origin='lower', extent=ext, aspect='auto', cmap='viridis',
                  vmin=-6, vmax=1)
    cs = a.contour(lrg, Tg, np.log10(z), levels=[-2, -1, np.log10(0.3), 0], colors='w',
                   linewidths=0.8)
    a.clabel(cs, fmt={-2: '0.01', -1: '0.1', np.log10(0.3): '0.3', 0: '1'}, fontsize=8)
    a.set_title(t, fontsize=10)
    a.set_xlabel('log rho [g/cm3]')
    a.set_ylabel('T [K]')
fig.colorbar(im, ax=ax, label='log10 Gamma (kappa / kappa_Edd)', shrink=0.8)
fig.savefig('gmap.png', dpi=110)
print(open('gmap_table.md').read())
sys.exit(0)
