"""Cold-fraction maps: 6 cubed-sphere panels side by side, cube vertices marked."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = '/viper/u2/jinma/ATHENAK/bench/cs_ens/analysis/coldmap/'
z = np.load(OUT + 'coldfrac.npz')
NC = 32


def fig(arm, rots, fname, title):
    n = len(rots)
    f, ax = plt.subplots(n, 6, figsize=(16.5, 2.9*n + 0.9), squeeze=False)
    for r, nn in enumerate(rots):
        F = z['%s_%d' % (arm, nn)]
        for p in range(6):
            a = ax[r, p]
            m = a.pcolormesh(np.arange(NC), np.arange(NC), F[p].T, cmap='viridis',
                             vmin=0, vmax=0.6, shading='nearest')
            ov = F[p].T > 0.30
            a.contour(np.arange(NC), np.arange(NC), ov.astype(float), levels=[0.5],
                      colors='r', linewidths=0.8)
            for jc in (0, NC-1):
                for kc in (0, NC-1):
                    a.plot(jc, kc, marker='*', color='w', ms=11, mec='k', mew=0.6)
            for e in (2.5, NC-3.5):
                a.axhline(e, color='w', ls=':', lw=0.6)
                a.axvline(e, color='w', ls=':', lw=0.6)
            a.set_xticks([0, 15, 31])
            a.set_yticks([0, 15, 31])
            a.set_aspect('equal')
            if r == 0:
                a.set_title('panel %d' % p, fontsize=10)
            if p == 0:
                a.set_ylabel('rot %d\n K (eta index)' % nn, fontsize=9)
            a.set_xlabel('J (xi index)', fontsize=8)
    f.colorbar(m, ax=ax, shrink=0.85,
               label='fraction of the 128 radial cells with T < 1500 K')
    f.suptitle(title + '   (white stars = cube vertices; dotted = 2-cell edge band; '
               'red = cold columns >30%)',
               fontsize=11)
    f.savefig(OUT + fname, dpi=125, bbox_inches='tight')
    plt.close(f)


fig('ctl2', [2, 3, 4, 5], 'fig_coldmap_ctl2.png',
    'ctl2/s01  EXPLICIT RT source (dies rot 5.71): cold-column map')
fig('si', [5], 'fig_coldmap_si_rot5.png',
    'si/s01  SEMI-IMPLICIT RT source (survives): cold-column map at rot 5')
fig('si', [1, 2, 3, 4, 5], 'fig_coldmap_si_all.png',
    'si/s01  SEMI-IMPLICIT RT source: cold-column map, rot 1-5')
print('ok')
