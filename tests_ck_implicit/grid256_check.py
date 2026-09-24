"""The nx1 = 256, 8-coefficient radial grid as the code builds it (LeftEdgeX +
StretchRPoly, same arithmetic order), its monotonicity and end radii, and its cells per
scale height against the design (nx1time_0924/design256.py numbers,
radres_0924/meas.npz, read only).
Run: python3 gridcheck.py > gridcheck.txt   (writes gridcheck.png)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

R0, R1 = 9.44e9, 2.0556e10          # prod4 mesh/x1min, x1max
C8 = [0.418017, 4.818585, -84.625150, 397.754052, -955.667782, 1270.732651,
      -887.750373, 254.765233]
TODAY = [-0.068392, -2.191487, 2.464818, -1.366698]
OUT = '/viper/ptmp2/jinma/grid256_0925/'


def code_faces(c, n):
    """x1f exactly as Coordinates builds it (cell_locations.hpp, grid_stretch.hpp)."""
    c = list(c) + [0.0] * (8 - len(c))
    f = np.empty(n + 1)
    for i in range(n + 1):
        x = float(i) / float(n)
        r = (x * R1 - x * R0) - (0.5 * R1 - 0.5 * R0) + (0.5 * R0 + 0.5 * R1)
        xi = (r - R0) / (R1 - R0)
        u, xik = xi, xi
        for k in range(8):
            u += c[k] * xik * (1.0 - xi)
            xik *= xi
        f[i] = R0 + (R1 - R0) * u
    return f


def centroids(f):
    q = f[:-1] / f[1:]
    return 0.25 * (q * q + 1.0) / ((1.0 / 3.0) * (q * q + q + 1.0)) * (f[1:] + f[:-1])


Z = np.load('/viper/ptmp2/jinma/radres_0924/meas.npz', allow_pickle=True)
rc0 = Z['rc']
TAGS = ['s0', 's1', 's2', 's3', 's4']
P = np.concatenate([Z[t + '_p'].astype(float) for t in TAGS])
H = np.concatenate([Z[t + '_Hhse'].astype(float) for t in TAGS])
CZ = np.concatenate([Z[t + '_cosz'].astype(float) for t in TAGS])
SUB = slice(0, None, 7)
REG = [('top 1e-9..1e-6', 1e-9, 1e-6), ('upper 1e-6..1e-3', 1e-6, 1e-3),
       ('phot/jet-base 1e-3..1e-1', 1e-3, 1e-1), ('jet 1e-1..10', 1e-1, 10),
       ('deep 10..300', 10, 300)]


def cph(rc, dr):
    nH = np.exp(np.array([np.interp(rc, rc0, x) for x in np.log(H[SUB])]))
    pb = np.exp(np.array([np.interp(rc, rc0, x) for x in np.log(P[SUB])])) / 1e6
    return nH / dr[None, :], pb


def table(name, f):
    rc, dr = centroids(f), np.diff(f)
    x, pb = cph(rc, dr)
    cz = CZ[SUB]
    s = f'{name:14s}'
    for nm, lo, hi in REG:
        for side, m in (('d', cz > 0.5), ('n', cz < -0.5)):
            sel = (pb >= lo) & (pb < hi) & m[:, None]
            s += f' {np.percentile(x[sel], 10):5.1f}/{np.median(x[sel]):5.1f}'
            s += ' |' if side == 'd' else '  '
    print(s)
    return rc, x, pb


f256 = code_faces(C8, 256)
f128 = code_faces(TODAY, 128)
dr = np.diff(f256)
print(f'# 8-coefficient grid, nx1 = 256: x1f[0] = {f256[0]!r} (x1min {R0!r}), '
      f'x1f[256] = {f256[-1]!r} (x1max {R1!r})'.replace('np.float64', ''))
print(f'# monotonic: {bool(np.all(dr > 0))}; dr min {dr.min():.4e} (i = {dr.argmin()}), '
      f'max {dr.max():.4e} (i = {dr.argmax()}); max |dr ratio - 1| = '
      f'{np.max(np.abs(dr[1:] / dr[:-1] - 1)):.4f}')
print('# i  x1f[i] [cm]  dr[i] = x1f[i+1]-x1f[i] [cm]')
for i in range(257):
    print(f'{i:4d} {f256[i]:.9e} ' + (f'{dr[i]:.6e}' if i < 256 else ''))
print('\n# cells per H, p10/median over columns, day | night; design (design256.txt) '
      'A 8c @256 = 16.4/17.2 | 2.6/5.3, 14.9/15.5 | 4.7/7.8, 7.5/15.0 | 4.8/5.6, '
      '4.7/5.2 | 4.6/5.1, 4.6/4.9 | 4.6/4.9')
print('# grid          ' + '  '.join(r[0] for r in REG))
res = {'today 128': table('today 128', f128), 'code 256 8c': table('code 256 8c', f256)}

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for (nm, (rc, x, pb)), col in zip(res.items(), ['0.5', 'C0']):
    cz = CZ[SUB]
    for side, m, ls in (('day', cz > 0.5, '-'), ('night', cz < -0.5, '--')):
        lp = np.median(np.log10(pb[m]), axis=0)
        ax[0].plot(np.percentile(x[m], 10, axis=0), lp, ls, color=col,
                   label=f'{nm} {side} p10')
    ff = f256 if '256' in nm else f128
    ax[1].plot(rc / 1e9, np.diff(ff) / 1e7, color=col, label=nm)
ax[0].axvline(5, color='C3', lw=0.8)
ax[0].set_xscale('log')
ax[0].invert_yaxis()
ax[0].set_xlabel('cells per scale height (p10 over columns)')
ax[0].set_ylabel('log10 p [bar]')
ax[0].legend(fontsize=7)
ax[1].set_xlabel('r [1e9 cm]')
ax[1].set_ylabel('dr [1e7 cm]')
ax[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT + 'gridcheck.png', dpi=110)
