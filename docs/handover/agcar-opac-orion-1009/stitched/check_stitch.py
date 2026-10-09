#!/usr/bin/env python3
"""Checks + plots for the stitched Planck table.
usage: check_stitch.py <outdir> <ext2_table>   (outdir holds the table and stitch_fields.npz)"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402

od, e2 = sys.argv[1], sys.argv[2]
z = np.load(od + '/stitch_fields.npz')
lT, lD, E, N, dl, Df, jf, wT, O, i0 = (z[k] for k in
                                      ('lT', 'lD', 'E', 'N', 'delta', 'Dfull', 'jf', 'wT',
                                       'O', 'i0'))
i0 = int(i0)
nTo = O.shape[0]
new = [s for s in open(od + '/planck_ext2orion_gs98_x0.36_z0.02.txt')]
old = [s for s in open(e2)]
print('lines: new %d ext2 %d' % (len(new), len(old)))
hn = [k for k, s in enumerate(new) if s.startswith('#')]
ho = [k for k, s in enumerate(old) if s.startswith('#')]
print('comment line positions identical:', hn == ho, '; differing header lines:',
      [k + 1 for k in ho if new[k] != old[k]])
dn = [s for s in new if not s.startswith('#')]
do = [s for s in old if not s.startswith('#')]
chg = (dl != 0).ravel()
ndiff = sum(a != b for a, b in zip(dn, do))
ndiff_unch = sum(a != b for a, b, c in zip(dn, do, chg) if not c)
print('data lines differing: %d; cells with delta != 0: %d; differing among unchanged: %d'
      % (ndiff, chg.sum(), ndiff_unch))
# value of the written text vs computed
Nw = np.array([float(s) for s in dn]).reshape(N.shape)
print('max |written - computed| %.1e' % abs(Nw - N).max())
# regions
rT = [('3.450-3.475 taper', 3.449, 3.476), ('3.500-4.000 Ferguson', 3.499, 4.001),
      ('4.025-4.175 blend', 4.024, 4.176), ('4.200-4.250 TOPS', 4.199, 4.251),
      ('4.275-4.300 taper', 4.274, 4.301)]
for name, a, b in rT:
    r = (lT >= a) & (lT <= b)
    d = dl[r]
    print('%-22s cells %4d  max|delta| %.3f  mean delta %+.3f  (rho_f %s)' % (
        name, (d != 0).sum(), abs(d).max(), d[d != 0].mean(),
        ' '.join('%.2f' % lD[j] for j in jf[r])))
for lr in (-21, -19, -17):
    j = int(round((lr + 21)/0.05))
    print('  delta at log rho %d: min %+.3f max %+.3f' % (lr, dl[i0:i0+nTo, j].min(),
                                                          dl[i0:i0+nTo, j].max()))
# (b) continuity across rho_f: jump between node jf-1 and jf
for name, A in (('ext2', E), ('new', N)):
    jr = [abs(A[i, jf[i]] - A[i, jf[i]-1]) for i in range(i0, i0+nTo)]
    k = int(np.argmax(jr))
    print('%s: max |kP(rho_f) - kP(rho_f - 1 node)| = %.4f at log T %.3f' %
          (name, max(jr), lT[i0+k]))
# jump in log T between adjacent T rows, in the rows' changed region and across the T edges
for name, A in (('ext2', E), ('new', N)):
    sl = slice(i0-2, i0+nTo+2)
    dT = abs(np.diff(A[sl, :], axis=0))
    jmax = jf[i0:i0+nTo].max()
    print('%s: max |T-row jump| over log T %.3f..%.3f, log rho < %.2f: %.4f' % (
        name, lT[i0-2], lT[i0+nTo+1], lD[jmax], dT[:, :jmax].max()))
    for lab, ii in (('low edge 3.400-3.500', range(i0-2, i0+2)),
                    ('high edge 4.250-4.350', range(i0+nTo-3, i0+nTo+1))):
        m = max(abs(A[i+1, :jmax] - A[i, :jmax]).max() for i in ii)
        print('   %s: max T-row jump %.4f' % (lab, m))
# full-correction (untapered) jump that the taper spreads
print('untapered delta at the edge rows: max|Dfull| 3.450 %.3f, 4.300 %.3f' % (
    abs(Df[i0]).max(), abs(Df[i0+nTo-1]).max()))
# plots
fig, ax = plt.subplots(figsize=(7.5, 5))
sl = slice(i0-4, i0+nTo+4)
TT, DD = np.meshgrid(lT[sl], lD, indexing='ij')
v = abs(dl).max()
pc = ax.pcolormesh(TT, DD, dl[sl], cmap='RdBu_r', vmin=-v, vmax=v, shading='nearest')
ax.plot(lT[i0:i0+nTo], lD[jf[i0:i0+nTo]], 'k.-', ms=3, lw=0.8, label=r'$\rho_f(T)$ (ext2 data floor)')
ax.set_xlabel(r'log T [K]')
ax.set_ylabel(r'log $\rho$ [g cm$^{-3}$]')
ax.set_ylim(-21, -12)
ax.legend(loc='upper left')
fig.colorbar(pc, label=r'log $\kappa_P^{new}$ - log $\kappa_P^{ext2}$ [dex]')
fig.tight_layout()
fig.savefig(od + '/stitch_delta_map.png', dpi=130)
fig, axs = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
for a, T in zip(axs.ravel(), (4000, 5000, 6000, 8000, 10000, 15000)):
    t = np.log10(T)
    i = int(np.argmin(abs(lT - t)))
    io = i - i0
    a.plot(lD, E[i], 'C0-', label='ext2')
    if 0 <= io < nTo:
        a.plot(lD[:O.shape[1]], O[io], 'C2--', label='orion LTE')
    a.plot(lD, N[i], 'C3:', lw=2, label='stitched')
    a.axvline(lD[jf[i]], color='k', lw=0.6)
    a.set_title('T = %d K (log T %.3f node)' % (T, lT[i]))
    a.set_xlim(-21, -8)
    a.set_ylabel(r'log $\kappa_P$ [cm$^2$ g$^{-1}$]')
for a in axs[1]:
    a.set_xlabel(r'log $\rho$ [g cm$^{-3}$]')
axs[0, 0].legend()
fig.tight_layout()
fig.savefig(od + '/stitch_cuts.png', dpi=130)
