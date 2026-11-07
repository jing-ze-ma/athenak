#!/usr/bin/env python3
"""cmp_arms.py [outpng]: the Athena++-equivalent arm (doface_1105/prod_vsnmlat_ns: doppler_vel = vstar_n + lateral merge x2 below 35 Rsun, no sink) against
the stopped vstar_n-only arm (doface_1105/prod_doface, ~18 d) and the Linhao-exact production (bsgli_1104/prod, 'none', 57.4 d).
Panels: r_photo,3D; max_{r<40}|Mdot_r|/M_env; Mdot_r/M_env at 22 Rsun (near the 20.2 Rsun base); L_top and L_in (wall)
from user.hst (1-d running mean over the 1000 s rows); solid-angle-mean T of the first 3 wall cells.
Definitions as cmp_mdot.py (Ma+2026 eqs. A9-A11, A16).  Per-dump reductions cached in cache/arms/<arm>_<n>.npz
(incremental); the 'none' run is reduced only up to the newest doface dump + 2 d.  CPU only (run inside an allocation)."""
import sys, glob, os, re, numpy as np
S = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, S); import bin_convert_repro as bc
C = S + '/cache/arms/'; os.makedirs(C, exist_ok=True)
J = '/resnick/groups/carnegie_poc/jingze/'
RUNS = [('cfl035', os.environ.get('BSG_REPRO_RUN', J + 'doface_1105/prod_vsnmlat_cfl35'), 'C3', '-'),
        ('vsnmlat', J + 'doface_1105/prod_vsnmlat_ns', 'C2', '-.'),
        ('vstar_n', J + 'doface_1105/prod_doface', 'C1', ':'),
        ('none', J + 'bsgli_1104/prod', 'C0', '--'),
        ('vstarfix', J + 'doface_1105/prod_vstarfix', 'C4', '-'),
        ('thr1e7', J + 'doface_1105/prod_vsnmlat_thr1e7', 'C5', ':')]
LAB = {'cfl035': 'this arm: his deck values, wrong v, CFL 0.3, GM 19.78', 'vstarfix': 'kick fix, GM 20, stopped 15.4 d', 'thr1e7': 'kicked, GM 20, stopped 35 d', 'vsnmlat': 'CFL 0.3, same (stopped 23.7 d)', 'vstar_n': 'vstar_n only (stopped ~18 d)', 'none': 'Linhao-exact production (none)'}
Rs, day, YR = 6.957e10, 86400.0, 3.15576e7
LW = 5.581461e38 / 12.0        # L of the profile x Omega/4pi (wedge pi/3 sr), the hst luminosity unit
NW = 3                         # wall cells shown
TAB = os.environ.get('BSG_OPAC', J + 'bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt')
hdr = [l for l in open(TAB) if l.startswith('#')]; nT, nD, lT0, dlT, lD0, dlD = [float(x) for x in hdr[4].split()[1:7]]
KT = np.array([float(l) for l in open(TAB) if not l.startswith('#')]).reshape(int(nT), int(nD))
TU = 0.62*1.66053907e-24/1.380649e-16


def kap(rho, T):
    x = np.clip((np.log10(T) - lT0)/dlT, 0, KT.shape[0] - 1.001); y = np.clip((np.log10(rho) - lD0)/dlD, 0, KT.shape[1] - 1.001)
    i, j = x.astype(int), y.astype(int); fx, fy = x - i, y - j
    return 10.0**((1-fx)*(1-fy)*KT[i, j] + fx*(1-fy)*KT[i+1, j] + (1-fx)*fy*KT[i, j+1] + fx*fy*KT[i+1, j+1])


def reduce(f, out):
    h = bc.read_binary(f); nb = len(h['mb_x1v'])
    rho = np.array(h['mb_data']['dens'], dtype=np.float64); vr = np.array(h['mb_data']['velx']); ei = np.array(h['mb_data']['eint'])
    x1f = np.array(h['mb_x1f'][0]); x1v = np.array(h['mb_x1v'][0]); x2f = np.array(h['mb_x2f']); x3f = np.array(h['mb_x3f'])
    dOm = np.array([np.outer(x3f[b][1:] - x3f[b][:-1], np.cos(x2f[b][:-1]) - np.cos(x2f[b][1:])) for b in range(nb)])
    Om = dOm.sum()
    mdot = (4*np.pi/Om)*x1v**2*np.einsum('bkji,bkj->i', rho*vr, dOm)
    menv = (4*np.pi/Om)*(np.einsum('bkji,bkj->i', rho, dOm)*np.diff(x1f**3)/3.0).sum()
    T = (2.0/3.0)*ei/rho*TU
    dtau = kap(rho, T)*rho*np.diff(x1f)[None, None, None, :]
    tau_f = np.concatenate([np.cumsum(dtau[..., ::-1], axis=-1)[..., ::-1], np.zeros(dtau.shape[:-1] + (1,))], axis=-1)
    tbar = np.einsum('bkji,bkj->i', tau_f, dOm)/Om
    i = np.where(tbar >= 1.0)[0][-1]
    lt0, lt1 = np.log(tbar[i]), np.log(max(tbar[i+1], 1e-300)); w = lt0/(lt0 - lt1)
    rph = x1f[i] + w*(x1f[i+1] - x1f[i])
    Tw = np.einsum('bkji,bkj->i', T[..., :NW], dOm)/Om
    rw = np.einsum('bkji,bkj->i', rho[..., :NW], dOm)/Om
    Tw_min, Tw_max = T[..., 0].min(), T[..., 0].max()
    np.savez(out, t=h['time'], r=x1v, mdot=mdot, menv=menv, rph=rph, Tw=Tw, rhow=rw, Tw_min=Tw_min, Tw_max=Tw_max)


def dump_time(f):
    return float(re.search(r'time=\s*([-+0-9.eE]+)', open(f, 'rb').read(400).decode('latin-1')).group(1))


def series(name, d, tmax=None):
    rows = []
    for f in sorted(glob.glob(d + '/bin/bsg3d.hydro_w.*.bin')):
        n = f.split('.')[-2]; out = C + '%s_%s.npz' % (name, n)
        if not os.path.exists(out):
            if tmax is not None and dump_time(f) > tmax:
                continue
            reduce(f, out)
        rows.append(dict(np.load(out)))
    if tmax is not None:
        rows = [c for c in rows if c['t'] <= tmax]
    return rows


def hst(d):
    """user.hst rows, restart overlaps resolved (the later segment wins), sorted in time."""
    fs = sorted(glob.glob(d + '/bsg3d.user.hst'))
    if not fs:
        return None
    a = np.loadtxt(fs[0], comments='#', ndmin=2)
    keep, tmin = [], np.inf
    for k in range(len(a) - 1, -1, -1):
        if a[k, 0] < tmin:
            keep.append(k); tmin = a[k, 0]
    return a[keep[::-1]]


def runmean(t, x, w=day):
    c = np.concatenate([[0.0], np.cumsum(x)]); out = np.empty_like(x)
    lo = np.searchsorted(t, t - w/2); hi = np.searchsorted(t, t + w/2, side='right')
    return (c[hi] - c[lo])/np.maximum(hi - lo, 1)


A = series(*RUNS[0][:2])
tA = max([c['t'] for c in A] + [0.0])
H = {RUNS[0][0]: hst(RUNS[0][1])}
tH = H[RUNS[0][0]][-1, 0] if H[RUNS[0][0]] is not None else tA
tcut = max(max(tA, tH) + 2*day, 24.0*day)   # comparison runs: the whole stopped CFL 0.3 arm (23.7 d)
SR = {RUNS[0][0]: A}
for nm, dd, _, _ in RUNS[1:]:
    SR[nm] = series(nm, dd, tmax=tcut)
    hb = hst(dd); H[nm] = hb[hb[:, 0] <= tcut] if hb is not None else None
B = SR['none']; Bv = SR['vstar_n']

i22 = None
print('vsnmlat arm vs vstar_n-only vs Linhao-exact none: per dump  [Mdot/M_env in 1/yr; T in K]')
print('  t[d] |  vsnmlat: r_ph  max|Md|/M  Md/M@22  T1        T2        T3       | vstar_n: same | none: same')


def cols(c):
    r = c['r']/Rs; q = c['mdot']/c['menv']*YR; m = r < 40
    return '%6.2f  %7.3f  %+7.3f  %.3e %.3e %.3e' % (c['rph']/Rs, np.abs(q[m]).max(), q[np.argmin(np.abs(r - 22))], *c['Tw'][:3])


bt = {round(float(c['t'])/day*2): c for c in B}; bv = {round(float(c['t'])/day*2): c for c in Bv}
for c in A:
    k = round(float(c['t'])/day*2)
    print('%6.2f | %s | %s | %s' % (c['t']/day, cols(c), cols(bv[k]) if k in bv else '-', cols(bt[k]) if k in bt else '-'))

import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize=(17, 9.2), constrained_layout=True)
ax = ax.ravel()
xmax = max(tcut/day, 1.0)
YL = {3: [], 4: []}
for name, d, col, ls in RUNS:
    rr = SR[name]; lab = LAB[name]
    if rr:
        t = np.array([c['t'] for c in rr])/day
        ax[0].plot(t, [c['rph']/Rs for c in rr], ls, color=col, marker='o', ms=3, label=lab)
        mm = np.array([np.abs(c['mdot']/c['menv']*YR)[c['r']/Rs < 40].max() for c in rr]); pos = mm > 0   # t = 0: exactly 0 (log axis)
        ax[1].plot(t[pos], mm[pos], ls, color=col, marker='o', ms=3, label=lab)
        ax[2].plot(t, [(c['mdot']/c['menv']*YR)[np.argmin(np.abs(c['r']/Rs - 22))] for c in rr], ls, color=col, marker='o', ms=3, label=lab)
        Tw = np.array([c['Tw'] for c in rr])
        for k, a in zip(range(NW), (1.0, 0.6, 0.35)):
            ax[5].plot(t, Tw[:, k], ls, color=col, alpha=a, lw=1.6, label='%s, cell %d' % ('lhxo' if name == 'cfl035' else name, k + 1))
    h = H.get(name)
    if h is not None and len(h) > 1:
        th = h[:, 0]
        for j, (ia, lw_) in enumerate(((5, 'L_top'), (6, 'L_in (wall)'))):
            y = h[:, ia]/LW
            a_ = ax[3] if j == 0 else ax[4]
            a_.plot(th/day, y, color=col, lw=0.4, alpha=0.35)
            YL[3 if j == 0 else 4].append(y)
            a_.plot(th/day, runmean(th, y), ls, color=col, lw=1.8, label=lab + ' (1-d mean)')
for a_, ys in ((ax[3], YL[3]), (ax[4], YL[4])):   # y-range from the 1-99 % of the raw rows (start-up spikes of a few rows off-scale)
    if ys:
        y = np.concatenate(ys); lo, hi = np.percentile(y, [1, 99]); pad = 0.1*(hi - lo + 1e-9)
        a_.set_ylim(min(lo - pad, 0.0 if a_ is ax[4] else lo - pad), max(hi + pad, 1.1))
ax[0].axhline(43, color='k', ls=':', label=r'paper steady $\approx$ 43 R$_\odot$ (after 30 d)')
ax[0].set(ylabel=r'$r_{\rm photo,3D}$ [R$_\odot$]', title=r'3D photosphere $r(\bar\tau=1)$ (eqs. A9-A11)')
ax[1].axhline(4, color='k', ls=':', label=r'paper: often > 4 yr$^{-1}$ (first 10 d, r < 40)')
ax[1].set(ylabel=r'max$_{r<40}|\dot M_r|/M_{\rm env}$ [yr$^{-1}$]', yscale='log', title=r'Radial mass flow, max over $r < 40$ R$_\odot$ (eq. A16)')
ax[2].axhline(0, color='k', lw=0.6)
ax[2].set(ylabel=r'$\dot M_r/M_{\rm env}$ at 22 R$_\odot$ [yr$^{-1}$]', title=r'Near-base mass flow (r = 22 R$_\odot$; base 20.2), >0 outward')
ax[3].axhline(1, color='k', ls=':', lw=0.8)
ax[3].set(ylabel=r'$L_{\rm top}/L$', title='Emergent luminosity (top face, user.hst)')
ax[4].axhline(1, color='k', ls=':', lw=0.8, label='L of the profile')
ax[4].axhline(0, color='k', lw=0.6)
ax[4].set(ylabel=r'$L_{\rm in}/L$  (>0 into the envelope)', title='Wall luminosity L$_{\\rm in}$ (inner face, user.hst)')
ax[5].set(ylabel='solid-angle mean T [K]', title='First 3 cells above the wall: T (collapse shows here first)')
for a in ax:
    a.set_xlabel('t [d]'); a.set_xlim(0, xmax); a.grid(alpha=0.25)
for a in (ax[0], ax[1], ax[3], ax[4]):
    a.legend(fontsize=8, loc='best')
ax[5].legend(fontsize=6, ncol=4, loc='best')
fig.suptitle('This arm (his deck values, wrong v, CFL 0.3, GM 19.78) vs earlier GM 20 arms vs Linhao-exact none; this arm at t = %.2f d' % (max(tA, tH)/day), fontsize=12)
out = sys.argv[1] if len(sys.argv) > 1 else S + '/out/ours_arms.png'
fig.savefig(out, dpi=110); print('wrote', out)
