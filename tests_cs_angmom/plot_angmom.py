import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

D = np.load('/viper/u2/jinma/ATHENAK/athenak/tests_cs_angmom/angmom.npz')
NAMES = [('cs_mhd_prod3', 'cs  6x32x32', 'C3'),
         ('sp_mhd_prod3', 'sp  64x128 (polar HLLE)', 'C0'),
         ('sp_mhd_nopole', 'sp  64x128 (polar HLLD)', 'C2')]

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
tab = []
for name, lab, c in NAMES:
    a = D[name]
    rot, M, Lrel, Lfrm = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    L = Lrel + Lfrm
    L0 = L[0]
    ax[0, 0].plot(rot, (L-L0)/L0, c, label=lab)
    ax[0, 1].plot(rot, (M-M[0])/M[0], c, label=lab)
    ax[1, 0].plot(rot, (Lfrm-Lfrm[0])/L0, c, label=lab+' frame')
    ax[1, 0].plot(rot, (Lrel-Lrel[0])/L0, c, ls='--', label=lab+' relative')
    ax[1, 1].plot(rot, Lrel/L0, c, label=lab)
    # secular drift: linear fit over rot > 20 (after the Rayleigh-drag spin-up)
    sel = rot > 20
    p = np.polyfit(rot[sel], (L[sel]-L0)/L0, 1)
    pm = np.polyfit(rot[sel], (M[sel]-M[0])/M[0], 1)
    spec = (L/M)/(L0/M[0]) - 1.0
    ps = np.polyfit(rot[sel], spec[sel], 1)
    ax[0, 0].plot(rot, spec, c, ls=':', lw=1)
    tab.append((name, lab, rot[-1], L0, (L[-1]-L0)/L0, p[0], pm[0],
                (Lfrm[-1]-Lfrm[0])/L0, (Lrel[-1]-Lrel[0])/L0, spec[-1], ps[0]))

for a_, t in zip(ax.ravel(), ['(L_z - L_z0)/L_z0, inertial',
                              'mass drift (M-M0)/M0',
                              'split: frame (solid) vs relative (dashed), /L_z0',
                              'L_rel / L_z0']):
    a_.set_title(t, fontsize=10)
    a_.set_xlabel('rotations')
    a_.axhline(0, color='k', lw=0.5)
    a_.grid(alpha=0.3)
ax[0, 0].legend(fontsize=8)
ax[1, 0].legend(fontsize=6)
fig.tight_layout()
fig.savefig('/viper/u2/jinma/ATHENAK/athenak/tests_cs_angmom/angmom.png', dpi=130)

print('run                    rot_end   L_z0[cgs]     dL/L0(end)   drift/rot   '
      'dM/M0/rot   dLfrm/L0   dLrel/L0    d(L/M)end   d(L/M)/rot')
for r in tab:
    print('%-22s %6.1f  %.6e  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e'
          % (r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]))

print()
print('L_z at ~10 times (relative to L_z0), and mass-weighted <u> [cm/s] at 3 shells')
for name, lab, c in NAMES:
    a = D[name]
    idx = np.unique(np.linspace(0, len(a)-1, 10).astype(int))
    print('==', lab)
    print('   rot     (L-L0)/L0    (M-M0)/M0     Lrel/L0      u(r=9.91e9) '
          'u(1.30e10)  u(1.84e10)')
    L = a[:, 2]+a[:, 3]
    for i in idx:
        print('  %6.1f  %+.4e  %+.4e  %+.4e   %9.2f  %9.2f  %9.2f'
              % (a[i, 0], (L[i]-L[0])/L[0], (a[i, 1]-a[0, 1])/a[0, 1],
                 a[i, 2]/L[0], a[i, 4], a[i, 5], a[i, 6]))
