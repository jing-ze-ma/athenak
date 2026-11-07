#!/usr/bin/env python3
"""cmp_mdot.py [outpng]: Ma, Bildsten & Jiang 2026 radial mass flow and 3-D photosphere for both bsgli_1104 runs.
Paper definitions (2609.28656): Sec. 3.1 + Fig. 2 third panel, eq. (A16): Mdot_r,i = (4 pi/Omega) sum_jk rho v_r r_i^2 dOmega_jk
(signed net shell flux, each snapshot, no time averaging), plotted as Mdot_r / M_envelope [yr^-1] in (t, r);
text: "in the first 10 days ... below r ~ 40 Rsun ... often exceeding Mdot_r/M_envelope ~ 4 yr^-1".
r_photo,3D: eqs. (A9)-(A11): tau_jk(r) = int_r^rmax kappa_R rho dr along each (theta_j, phi_k), tau_bar = (1/Omega) sum tau_jk dOmega,
r_photo,3D = r(tau_bar = 1) (interpolated). M_envelope = (4 pi/Omega) x domain mass (each snapshot).
Per-dump reductions cached in cmp_mdot_cache/<run>_<n>.npz (incremental).  Run inside an allocation (srun --overlap -c 4)."""
import sys, glob, os, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/bsgli_1104/page_none/pipe'); import bin_convert_repro as bc
W = '/resnick/groups/carnegie_poc/jingze/bsgli_1104/'; C = W + 'cmp_mdot_cache/'; os.makedirs(C, exist_ok=True)
RUNS = [('none', W + 'prod'), ('athena', W + 'prod_athena')]
Rs, day, YR = 6.957e10, 86400.0, 3.15576e7
TAB = '/resnick/groups/carnegie_poc/jingze/bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt'
hdr = [l for l in open(TAB) if l.startswith('#')]; nT, nD, lT0, dlT, lD0, dlD = [float(x) for x in hdr[4].split()[1:7]]
KT = np.array([float(l) for l in open(TAB) if not l.startswith('#')]).reshape(int(nT), int(nD))
TU = 0.62*1.66053907e-24/1.380649e-16
def kap(rho, T):
    x = np.clip((np.log10(T) - lT0)/dlT, 0, KT.shape[0] - 1.001); y = np.clip((np.log10(rho) - lD0)/dlD, 0, KT.shape[1] - 1.001)
    i, j = x.astype(int), y.astype(int); fx, fy = x - i, y - j
    return 10.0**((1-fx)*(1-fy)*KT[i, j] + fx*(1-fy)*KT[i+1, j] + (1-fx)*fy*KT[i, j+1] + fx*fy*KT[i+1, j+1])
def reduce(f, out):
    h = bc.read_binary(f); nb = len(h['mb_x1v'])
    rho = np.array(h['mb_data']['dens']); vr = np.array(h['mb_data']['velx']); ei = np.array(h['mb_data']['eint'])
    x1f = np.array(h['mb_x1f'][0]); x1v = np.array(h['mb_x1v'][0]); x2f = np.array(h['mb_x2f']); x3f = np.array(h['mb_x3f'])
    dOm = np.array([np.outer(x3f[b][1:] - x3f[b][:-1], np.cos(x2f[b][:-1]) - np.cos(x2f[b][1:])) for b in range(nb)])
    Om = dOm.sum()
    mdot = (4*np.pi/Om)*x1v**2*np.einsum('bkji,bkj->i', rho*vr, dOm)
    mwed = np.einsum('bkji,bkj->i', rho, dOm)*np.diff(x1f**3)/3.0
    menv = (4*np.pi/Om)*mwed.sum()
    T = (2.0/3.0)*ei/rho*TU
    dtau = kap(rho, T)*rho*np.diff(x1f)[None, None, None, :]
    tau_f = np.concatenate([np.cumsum(dtau[..., ::-1], axis=-1)[..., ::-1], np.zeros(dtau.shape[:-1] + (1,))], axis=-1)   # at faces
    tbar = np.einsum('bkji,bkj->i', tau_f, dOm)/Om
    i = np.where(tbar >= 1.0)[0][-1]   # last face with tau_bar >= 1
    lt0, lt1 = np.log(tbar[i]), np.log(max(tbar[i+1], 1e-300)); w = lt0/(lt0 - lt1)
    rph = x1f[i] + w*(x1f[i+1] - x1f[i])
    np.savez(out, t=h['time'], r=x1v, mdot=mdot, menv=menv, rph=rph)
def series(name, d):
    rows = []
    for f in sorted(glob.glob(d + '/bin/bsg3d.hydro_w.*.bin')):
        n = f.split('.')[-2]; out = C + '%s_%s.npz' % (name, n)
        if not os.path.exists(out): reduce(f, out)
        rows.append(dict(np.load(out)))
    return rows
S = {n: series(n, d) for n, d in RUNS}
print('Mdot_r/M_env [1/yr] (paper eq. A16), signed; max|.| over r < 40 Rsun; r_photo,3D (A9-A11)')
print('  t[d] | none: max|Mdot|/M (r at max)  Mdot/M at 22 30 38 Rsun   r_ph | athena: same')
def line(c):
    r = c['r']/Rs; q = c['mdot']/c['menv']*YR; m = r < 40; k = np.argmax(np.abs(q[m]))
    at = [q[np.argmin(np.abs(r - x))] for x in (22, 30, 38)]
    return '%6.2f (%5.1f)  %+6.2f %+6.2f %+6.2f  %6.2f' % (abs(q[m][k]), r[m][k], *at, c['rph']/Rs)
for i in range(max(len(v) for v in S.values())):
    a = S['none'][i] if i < len(S['none']) else None; b = S['athena'][i] if i < len(S['athena']) else None
    t = (a or b)['t']/day
    print('%6.2f | %s | %s' % (t, line(a) if a else '-', line(b) if b else '-'))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
for j, (n, _) in enumerate(RUNS):
    rr = S[n]
    if not rr: continue
    t = np.array([c['t'] for c in rr])/day; r = rr[0]['r']/Rs; q = np.array([c['mdot']/c['menv']*YR for c in rr])
    pc = ax[0, j].pcolormesh(t, r, q.T, cmap='RdBu', vmin=-5, vmax=5, shading='nearest')
    ax[0, j].plot(t, [c['rph']/Rs for c in rr], 'k--', lw=1.5, label=r'$r_{\rm photo,3D}$')
    ax[0, j].set(title=r'%s: $\dot M_r/M_{\rm env}$ [yr$^{-1}$] (eq. A16)' % n, xlabel='t [d]', ylabel=r'r [R$_\odot$]', ylim=(20, 80))
    ax[0, j].legend(loc='upper right'); fig.colorbar(pc, ax=ax[0, j])
    mm = [np.abs(c['mdot']/c['menv']*YR)[c['r']/Rs < 40].max() for c in rr]
    ax[1, 0].plot(t, mm, 'o-', label=n); ax[1, 1].plot(t, [c['rph']/Rs for c in rr], 'o-', label=n)
ax[1, 0].axhline(4, color='k', ls=':', label='paper: often > 4 yr$^{-1}$ (first 10 d, r < 40)')
ax[1, 0].set(xlabel='t [d]', ylabel=r'max$_{r<40}|\dot M_r|/M_{\rm env}$ [yr$^{-1}$]', yscale='log'); ax[1, 0].legend()
ax[1, 1].axhline(43, color='k', ls=':', label=r'paper steady $r_{\rm photo,3D}\approx 43$ R$_\odot$')
ax[1, 1].set(xlabel='t [d]', ylabel=r'$r_{\rm photo,3D}$ [R$_\odot$]'); ax[1, 1].legend()
out = sys.argv[1] if len(sys.argv) > 1 else W + 'cmp_mdot.png'; fig.savefig(out, dpi=110); print('wrote', out)
