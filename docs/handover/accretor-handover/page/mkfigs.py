#!/usr/bin/env python3
"""Figures for the RY Per accretor results page (stage 1).
Face-on density maps of the GPU stream tests with the ballistic orbit, and the
angular-momentum budget of the CPU spin 1.0 / 7.2 arms.  Code units: Rsun, km/s."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc  # noqa: E402

A = '/viper/ptmp2/jinma/accretor_1006'
OUT = A + '/page'
RACC = 4.06
P_ORB = 0.8519      # orbital period in code time units (IMPL2.md)

# ballistic stream (same integrator as tests/ana.py)
Ma, Md, a_rs = 6.24, 1.69, 30.3
mu = Md/(Ma + Md)
xa, xd = -mu, 1 - mu
G, Msun, Rsun = 6.674e-8, 1.989e33, 6.957e10
om = np.sqrt(G*(Ma + Md)*Msun/(a_rs*Rsun)**3)
aom = a_rs*Rsun*om/1e5
eps = np.sqrt(1.381e-16*6250/(1.27*1.6726e-24))/1e5/aom


def gx(x):
    return -(1 - mu)*(x - xa)/abs(x - xa)**3 - mu*(x - xd)/abs(x - xd)**3 + x


xL1 = brentq(gx, xa + 1e-3, xd - 1e-3)


def rhs(t, s):
    x, y, vx, vy = s
    r1 = np.hypot(x - xa, y)
    r2 = np.hypot(x - xd, y)
    return [vx, vy, -(1 - mu)*(x - xa)/r1**3 - mu*(x - xd)/r2**3 + x + 2*vy,
            -(1 - mu)*y/r1**3 - mu*y/r2**3 + y - 2*vx]


sol = solve_ivp(rhs, [0, 3], [xL1 - 1e-4, 0, -eps, 0], rtol=1e-11, atol=1e-13,
                max_step=1e-3)
bx, by = (sol.y[0] - xa)*a_rs, sol.y[1]*a_rs


def load(fn):
    fd = bc.read_binary(fn)
    bl = []
    for b in range(fd['n_mbs']):
        bl.append((fd['mb_x3v'][b][0], fd['mb_x1v'][b], fd['mb_x3v'][b],
                   fd['mb_data']['dens'][b]))
    bl.sort(key=lambda t: t[0])
    r = bl[0][1]
    ph = np.concatenate([b[2] for b in bl])
    d = np.concatenate([b[3] for b in bl], axis=0).mean(axis=1)   # (phi, r)
    return fd.get('time', np.nan), r, ph, d


def facemap(ax, fn, title, vmin, vmax, rmax=None):
    t, r, ph, d = load(fn)
    rf = np.concatenate([[r[0] - 0.5*(r[1] - r[0])], 0.5*(r[1:] + r[:-1]),
                         [r[-1] + 0.5*(r[-1] - r[-2])]])
    dph = ph[1] - ph[0]
    pf = np.concatenate([ph - 0.5*dph, [ph[-1] + 0.5*dph]])
    R, P = np.meshgrid(rf, pf)
    X, Y = R*np.cos(P), R*np.sin(P)
    ld = np.log10(np.maximum(d/np.nanmax(d), 1e-12))
    m = ax.pcolormesh(X, Y, ld, cmap='magma', vmin=vmin, vmax=vmax, shading='flat',
                      rasterized=True)
    th = np.linspace(0, 2*np.pi, 400)
    rin = rf[0]
    # inside the inner boundary is not simulated: plain grey disk with the boundary drawn
    ax.fill(rin*np.cos(th), rin*np.sin(th), color='#cfd3d8', zorder=3)
    ax.plot(rin*np.cos(th), rin*np.sin(th), color='black', lw=1.0, zorder=4)
    lab = 'inner boundary\n(not simulated)' if rin > 3.5 else 'inner\nboundary'
    ax.text(0, 0, lab, ha='center', va='center', fontsize=8, color='0.25', zorder=5)
    if abs(rin - RACC) > 1e-3:
        ax.plot(RACC*np.cos(th), RACC*np.sin(th), color='#2a78d6', lw=1.2, zorder=4)
    rb = np.hypot(bx, by)
    bxm = np.where((rb <= rf[-1]) & (rb >= rin), bx, np.nan)
    bym = np.where((rb <= rf[-1]) & (rb >= rin), by, np.nan)
    ax.plot(bxm, bym, color='white', lw=0.9, ls='--', zorder=4)
    lim = rmax or rf[-1]
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_title(f'{title}\nt = {t/P_ORB:.2f} orbit', fontsize=10)
    ax.set_xlabel('x [R$_\\odot$] (donor to +x)')
    return m


# Fig 1: cold-stream test (r_in = R_acc), four times
fns = [f'{A}/gpu/a2_rin406/bin/ryper.hydro_w.{n:05d}.bin' for n in (1, 3, 6, 10)]
fig, axs = plt.subplots(1, 4, figsize=(17, 4.6), constrained_layout=True)
for ax, fn in zip(axs, fns):
    m = facemap(ax, fn, 'cold stream, absorbing star', -6, 0)
axs[0].set_ylabel('y [R$_\\odot$]')
fig.colorbar(m, ax=axs, shrink=0.85, label='log$_{10}$ $\\rho/\\rho_{max}$ (θ-mean)')
fig.savefig(OUT + '/fig_stream.png', dpi=110)

# Fig 2: zoom on the impact, last dump, both inner radii
fig, axs = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
facemap(axs[0], f'{A}/gpu/a2_rin406/bin/ryper.hydro_w.00010.bin',
        'r$_{in}$ = R$_{acc}$ (absorbing surface)', -6, 0, rmax=7)
m = facemap(axs[1], f'{A}/gpu/a2_rin2/bin/ryper.hydro_w.00010.bin',
            'r$_{in}$ = 2 R$_\\odot$ (no star: free stream)', -6, 0, rmax=7)
axs[0].set_ylabel('y [R$_\\odot$]')
fig.colorbar(m, ax=axs, shrink=0.85, label='log$_{10}$ $\\rho/\\rho_{max}$')
fig.savefig(OUT + '/fig_impact.png', dpi=110)

# Fig 3: specific AM of the accreted gas, spin 1.0 vs 7.2 (CPU low res)
GM = 2.93e5*RACC          # (km/s)^2 Rsun, GM/R = 2.93e5 (IMPL2.md)
jk = np.sqrt(GM*RACC)
fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.0), constrained_layout=True)
for lab, d, c in (('spin 1.0 × sync', 'c_spin1', '#2a78d6'), ('spin 7.2 × sync', 'c_spin72', '#eb6834')):
    h = np.loadtxt(f'{A}/tests/{d}/ryper.user.hst')
    t, macc, jacc = h[:, 0], h[:, 6], h[:, 7]
    ok = macc > 1e-3*macc.max()
    ax.plot(t[ok]/P_ORB, jacc[ok]/macc[ok]/jk, color=c, lw=2, label=lab)
ax.axhline(1.0647, color='0.5', ls=':', lw=1)
ax.text(0.02, 1.0647, ' ballistic value 1.065', va='bottom', fontsize=9, color='0.35',
        transform=ax.get_yaxis_transform())
ax.set_xlabel('time [orbits]')
ax.set_ylabel('accreted j / j$_{Kep}$(R$_{acc}$)')
ax.legend(frameon=False)
fig.savefig(OUT + '/fig_jacc.png', dpi=120)
print('ok')
