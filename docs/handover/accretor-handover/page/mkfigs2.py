#!/usr/bin/env python3
"""Stage-2 figures for the RY Per accretor page: the hot-star (c_ph 85.3 km/s) impact
run gpu2/t6_s1 (face-on density zooms) and its mass / angular-momentum budget."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
sys.path.insert(0, '/viper/ptmp2/jinma/accretor_1006/page')
import mkfigs as mf  # noqa: E402  (ballistic orbit, load, facemap)

A = '/viper/ptmp2/jinma/accretor_1006'
OUT = A + '/page'
RUN = A + '/gpu2/t6_s1'
RACC, P_ORB = 4.06, 0.8519

# face-on zooms at four times
fig, axs = plt.subplots(1, 4, figsize=(17, 4.6), constrained_layout=True)
for ax, n in zip(axs, (2, 5, 9, 15)):
    m = mf.facemap(ax, f'{RUN}/bin/ryper.hydro_w.{n:05d}.bin', 'hot star (c_ph 85 km/s)',
                   -7, 0, rmax=8)
axs[0].set_ylabel('y [R$_\\odot$]')
fig.colorbar(m, ax=axs, shrink=0.85, label='log$_{10}$ $\\rho/\\rho_{max}$ (θ-mean)')
fig.savefig(OUT + '/fig_hot_maps.png', dpi=110)

# budget from the user history (cumulative columns)
h = np.loadtxt(f'{RUN}/ryper.user.hst')
t, menv, jenv, mstr, jstr, mR, jR = h[:, 0], h[:, 4], h[:, 5], h[:, 6], h[:, 14], h[:, 12], h[:, 13]
GM = 2.93e5*RACC
jk = np.sqrt(GM*RACC)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
ax = axs[0]
ax.plot(t/P_ORB, mstr, color='#2a78d6', lw=2, label='stream inflow (cumulative)')
ax.plot(t/P_ORB, mR, color='#eb6834', lw=2, label='mass through R$_{acc}$ (cumulative)')
ax.set_xlabel('time [orbits]')
ax.set_ylabel('mass [code units]')
ax.legend(frameon=False)
ax = axs[1]
ok = mstr > 1e-3*max(mstr.max(), 1e-30)
ax.plot(t[ok]/P_ORB, jstr[ok]/mstr[ok]/jk, color='#2a78d6', lw=2, label='j brought in by the stream')
okR = np.abs(mR) > 1e-2*max(np.abs(mR).max(), 1e-30)
ax.plot(t[okR]/P_ORB, jR[okR]/mR[okR]/jk, color='#eb6834', lw=2, label='j of gas through R$_{acc}$')
ax.axhline(0, color='0.6', lw=0.8)
ax.set_xlabel('time [orbits]')
ax.set_ylabel('specific AM / j$_{Kep}$(R$_{acc}$)')
ax.legend(frameon=False)
fig.savefig(OUT + '/fig_hot_budget.png', dpi=120)
print('ok', t[-1], mstr[-1], mR[-1])
