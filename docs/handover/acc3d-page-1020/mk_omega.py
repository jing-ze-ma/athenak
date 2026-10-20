#!/usr/bin/env python3
"""ACC3D page 10-20: equatorial angular velocity profiles (READ ONLY on the run). usage: mk_omega.py <imgdir>
Omega(r) = sum(dm R v_phi) / sum(dm R^2) over |lat| < 5 deg and all phi (mass weighted), + Omega_orb (inertial frame), in
units of Omega_orb (P 3.81 d; 13.28 code). TRUE faces (common.load: mb_x1f/mb_x2f). With tracers also the s0-weighted
Omega of the stream gas. Per-dump profiles cached in $OUT/data/omega/."""
import os
import sys
import json
import numpy as np
import matplotlib.ticker
import common as C

IMG = sys.argv[1]
CD = os.path.join(C.OUT, 'data', 'omega'); os.makedirs(CD, exist_ok=True)
GM = 3.053e6   # G * 16 Msun [Rsun (km/s)^2]


def prof(f):
    key = os.path.join(CD, '%s_%05d.npz' % (os.path.basename(os.path.dirname(os.path.dirname(f))), C.dnum(f)))
    if os.path.exists(key) and os.path.getmtime(key) > os.path.getmtime(f):
        z = np.load(key); return {k: z[k] for k in z.files}
    h = C.load(f, ['dens', 'velz', 's_00'])
    r, th = h['x1f'], h['x2f']
    rc = h['x1v']; tc = 0.5 * (th[1:] + th[:-1])
    eq = np.abs(tc - np.pi / 2) < np.radians(5.0)
    dm = (h['dens'] * C.dvol(h))[:, eq, :]
    R = (rc[None, None, :] * np.sin(tc)[None, :, None])[:, eq, :]
    vp = h['velz'][:, eq, :]
    J = (dm * R * vp).sum(axis=(0, 1)); I = (dm * R ** 2).sum(axis=(0, 1))
    o = dict(t=h['Time'], r=rc, om=(J / I + C.OM) / C.OM)
    if 's_00' in h:
        w = dm * np.clip(h['s_00'][:, eq, :], 0, 1)
        Js = (w * R * vp).sum(axis=(0, 1)); Is = (w * R ** 2).sum(axis=(0, 1))
        o['om_s0'] = np.where(Is > 0, Js / np.maximum(Is, 1e-300) + C.OM, np.nan) / C.OM
        o['m_s0'] = w.sum(axis=(0, 1))
    np.savez(key, **o)
    return o


plt = C.mpl()
L = C.dumps('hydro_w')
stt = [x for x in L if x[2] == 'static']
sm = [x for x in L if x[2] == 'stream']
pick = ([stt[-1]] if stt else []) + ([sm[i] for i in np.unique(np.linspace(0, len(sm) - 1, 8).round().astype(int))[1:]] if len(sm) > 1 else sm)
fig, ax = plt.subplots(figsize=(8.6, 5.8))
cm = plt.cm.viridis(np.linspace(0.05, 0.95, max(len(pick), 1)))
out = []
for c, (t, f, ph) in zip(cm, pick):
    p = prof(f)
    lab = 'static phase, t %.2f d' % (t * C.DAYS) if ph == 'static' else 'stream +%.2f d' % ((t - C.TS) * C.DAYS)
    m = p['r'] <= C.DL1
    ax.plot(p['r'][m], p['om'][m], color=c, lw=1.6, label=lab)
    out.append(dict(t=float(t), file=os.path.basename(f)))
    if 'om_s0' in p and f == pick[-1][1]:
        ok = m & (p['m_s0'] > 1e-6 * np.nanmax(p['m_s0']))
        ax.plot(p['r'][ok], p['om_s0'][ok], color=C.C5, lw=1.2, ls='--', label='latest, s0-weighted (stream gas only)')
rr = np.linspace(C.RIN, C.DL1, 300)
ax.axhline(1.0, color='#e8ecf1', ls='--', lw=1.0, label=r'corotation ($\Omega_\mathrm{orb}$, P 3.81 d)')
ax.plot(rr, np.sqrt(GM / rr ** 3) / C.OM, color=C.C1, ls=':', lw=1.4, label=r'Keplerian $\Omega_K(r)$ (16 M$_\odot$)')
ax.axvline(C.RACC, color='0.6', ls='-.', lw=1, label=r'accretor surface (9.0 R$_\odot$)')
ax.axvline(C.RMEAS, color='0.7', ls=':', lw=1, label=r'$r_\mathrm{meas}$ 9.57')
ax.set_xlim(7.5, C.DL1); ax.set_ylim(0.8, 3.0); ax.set_xscale('log'); ax.set_xticks([7.5, 8, 8.5, 9, 10, 12, 14, 16])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel(r'r [R$_\odot$] (to $d_{L1}$)'); ax.set_ylabel(r'$\Omega / \Omega_\mathrm{orb}$ (inertial)')
ax.set_title(f'Equatorial angular velocity, |lat| < 5°, mass weighted ({C.RUN_LABEL}, {C.SITE})', fontsize=10)
ax.grid(alpha=0.3); fig.legend(*ax.get_legend_handles_labels(), fontsize=7.5, loc='lower center', ncol=4, frameon=False)
fig.tight_layout(rect=(0, 0.19, 1, 1)); fig.savefig(os.path.join(IMG, 'fig_omega_eq.png'), dpi=110)
json.dump(out, open(os.path.join(IMG, 'omega.json'), 'w'))
print('omega: %d profiles' % len(pick))
