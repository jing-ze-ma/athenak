#!/usr/bin/env python3
"""H2H 1D: Athena++ (bsg_star.cpp) vs AthenaK (rad_m1, Linhao IC+BC) on the same 1-D setup.
usage: h2h.py  -> h2h/ts_*.json, figures h2h/*.png"""
import sys, glob, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')
import bin_convert as bc

W = '/raven/ptmp/jinma/athpp_audit_1015/h2h'
RS = 6.957e10; DAY = 86400.0; L0 = 5.581461e38; C = 2.99792458e10; A = 7.565733250033928e-15
KB = 1.380649e-16; MU = 1.66053906660e-24; MUW = 0.62; GM = 2.6542488e27
OMEGA = np.pi/3.0   # AthenaK wedge solid angle: (cos pi/3 - cos 2pi/3) * pi/3

# Rosseland table (same file both codes use)
def _tab(fn):
    for ln in open(fn):
        if ln.startswith('# ') and len(ln.split()) == 7 and ln.split()[1].isdigit():
            nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in ln.split()[1:]]
    return np.loadtxt(fn, comments='#').reshape(int(nT), int(nD)), lT0, dlT, lD0, dlD
TV, LT0, DLT, LD0, DLD = _tab('/raven/ptmp/jinma/bsg_raven_1002/ma2026_tables/rosseland_ma2026_x0.7_z0.008_padded.txt')
def kap(T, rho):
    x = (np.log10(T)-LT0)/DLT; y = (np.log10(np.maximum(rho, 1e-30))-LD0)/DLD
    i = np.clip(x.astype(int), 0, TV.shape[0]-2); j = np.clip(y.astype(int), 0, TV.shape[1]-2)
    fx = np.clip(x-i, 0, 1); fy = np.clip(y-j, 0, 1)
    return 10**((1-fx)*(1-fy)*TV[i, j]+fx*(1-fy)*TV[i+1, j]+(1-fx)*fy*TV[i, j+1]+fx*fy*TV[i+1, j+1])


def prof_athpp(f):
    """Athena++ tab (prim+rad) -> cgs profile dict"""
    t = float(open(f).readline().split('time=')[1].split()[0])
    x = np.loadtxt(f, comments='#')
    tu = 1.899769e4; Tu = 1e5; ru = 1e-9; cs = RS/tu
    r = x[:, 1]*RS; rho = x[:, 2]*ru; p = x[:, 3]*ru*cs**2; v = x[:, 4]*cs
    E = x[:, 7]*A*Tu**4; F = x[:, 8]*C*A*Tu**4
    return dict(t=t*tu, r=r, rho=rho, T=p*MUW*MU/(rho*KB), eint=1.5*p, v=v, E=E, F=F)


def prof_athk(fh, fm):
    out = {}
    for f in (fh, fm):
        d = bc.read_binary(f)
        for k in d['var_names']:
            out[k] = d['mb_data'][k][0].mean(axis=(0, 1))   # one block, lateral mean
        g = d['mb_geometry'][0]; t = d['time']
    n = out['dens'].size; dr = (g[1]-g[0])/n
    rf = g[0] + dr*np.arange(n+1); r = 0.75*(rf[1:]**4-rf[:-1]**4)/(rf[1:]**3-rf[:-1]**3)
    rho = out['dens']; eint = out['eint']
    return dict(t=t, r=r, rho=rho, T=eint*(2./3.)*MUW*MU/(rho*KB), eint=eint, v=out['velx'], E=out['m1_e'], F=out['m1_f1'])


def series(profs):
    ts = {k: [] for k in ('t', 'T0', 'T1', 'T2', 'T3', 'rho0', 'rho1', 'rho2', 'rho3', 'F0', 'F1', 'F2', 'F3',
                          'Lwall', 'Ltop', 'rsurf', 'U', 'Erad', 'Egrav', 'Ekin', 'Etot', 'Lout_adv')}
    for p in profs:
        r = p['r']; dr = np.gradient(r); dV = 4*np.pi*r**2*dr   # full sphere
        ts['t'].append(p['t']/DAY)
        for i in range(4):
            ts['T%d' % i].append(p['T'][i]); ts['rho%d' % i].append(p['rho'][i])
            ts['F%d' % i].append(4*np.pi*r[i]**2*p['F'][i]/L0)
        ts['Lwall'].append(4*np.pi*r[0]**2*p['F'][0]/L0); ts['Ltop'].append(4*np.pi*r[-1]**2*p['F'][-1]/L0)
        s = np.where(p['rho'] > 1e-10)[0]; ts['rsurf'].append(r[s[-1]]/RS if len(s) else np.nan)
        U = np.sum(p['eint']*dV); Er = np.sum(p['E']*dV); Eg = np.sum(-GM*p['rho']/r*dV); Ek = np.sum(0.5*p['rho']*p['v']**2*dV)
        for k, v in (('U', U), ('Erad', Er), ('Egrav', Eg), ('Ekin', Ek), ('Etot', U+Er+Eg+Ek)): ts[k].append(v)
        vt = max(p['v'][-1], 0.0)
        ts['Lout_adv'].append(4*np.pi*r[-1]**2*vt*(p['eint'][-1]*5/3 + 0.5*p['rho'][-1]*vt**2 - GM*p['rho'][-1]/r[-1]
                                                      + 4/3*p['E'][-1])/L0)
    return {k: np.array(v) for k, v in ts.items()}


def gamma(p):
    return kap(p['T'], p['rho'])*p['F']/(C*GM/p['r']**2)


KNAMES = ('K1', 'K2', 'K3', 'K0')
import os as _os
if _os.environ.get('EMUL'):
    KNAMES = ('K1', 'E0', 'R1', 'R1t1', 'R2', 'R1R2', 'R3', 'R1R3', 'R4', 'R5a', 'R5b', 'R5c')
runs = {}
ANAMES = ('A_tc4', 'A_nf4', 'A_af1')
if _os.environ.get('BCS'):
    ANAMES = ('A_tc4', 'A_af1', 'Ba0', 'Ba1', 'Bb0', 'Bb1', 'Bc0', 'Bc1', 'Bd0', 'Bd1', 'Be0', 'Be1')
    KNAMES = ('K1',)
for name in ANAMES:
    od = 'outb' if _os.path.isdir(W + '/runs/%s/outb' % name) else 'out'
    fs = sorted(glob.glob(W + '/runs/%s/%s/bsg.block0.out2.*.tab' % (name, od)))
    if fs:
        runs[name] = [prof_athpp(f) for f in fs]
for name in KNAMES:
    fh = sorted(glob.glob(W + '/runs/%s/bin/bsg3d.hydro_w.*.bin' % name))
    fm = sorted(glob.glob(W + '/runs/%s/bin/bsg3d.m1.*.bin' % name))
    n = min(len(fh), len(fm))
    if n:
        runs[name] = [prof_athk(fh[i], fm[i]) for i in range(n)]

# drop non-finite (crashed) dumps
for k in list(runs):
    runs[k] = [p for p in runs[k] if np.all(np.isfinite(p['T'])) and np.all(np.isfinite(p['F'])) and p['rho'][:20].min() > 1e-15]
TS = {k: series(v) for k, v in runs.items()}
LAB = {'Be0': 'Be closed wall + reflecting hydro (af 0)', 'Be1': 'Be closed wall + reflecting hydro (af 1)', 'Bd0': 'Bd closed fixed-flux wall (af 0)', 'Bd1': 'Bd closed fixed-flux wall (af 1)', 'Ba0': 'Ba our flux BC (angle_flag 0)', 'Ba1': 'Ba our flux BC (angle_flag 1)', 'Bb0': 'Bb diffusion bath (af 0)', 'Bb1': 'Bb diffusion bath (af 1)', 'Bc0': 'Bc = Bb, no clip (af 0)', 'Bc1': 'Bc = Bb, no clip (af 1)','E0': 'AthenaK emul build, keys off', 'R1': 'R1 Athena++ DO wall face (taucell 4)', 'R1t1': 'R1 taucell 1', 'R2': 'R2 per-angle RHS clip', 'R1R2': 'R1+R2', 'R3': 'R3 no curvature (angle_flag 0)', 'R1R3': 'R1+R3', 'R4': 'R4 HLL diffusion eps', 'R5a': 'R5a central flux', 'R5b': 'R5b GM x1.0078', 'R5c': 'R5c loose solve','A_af1': 'Athena++ angle_flag=1 (spherical angular flux)', 'A_tc4': 'Athena++ (his code, taucell 4)', 'A_nf4': 'Athena++ no RHS clip', 'K1': 'AthenaK prod keys',
       'K2': 'AthenaK @ Athena++ dt', 'K3': 'AthenaK no pos floors', 'K0': 'AthenaK, our flux BC (not his)'}
COL = {'Be0': '#8c564b', 'Be1': '#c49c94', 'Bd0': '#ff7f0e', 'Bd1': '#ffbb78', 'Ba0': '#1f77b4', 'Ba1': '#17becf', 'Bb0': '#2ca02c', 'Bb1': '#98df8a', 'Bc0': '#9467bd', 'Bc1': '#c5b0d5','E0': '#1f5fa8', 'R1': '#d62728', 'R1t1': '#ff7f0e', 'R2': '#9467bd', 'R1R2': '#8c564b', 'R3': '#e377c2', 'R1R3': '#bcbd22', 'R4': '#17becf', 'R5a': '#2ca02c', 'R5b': '#7f7f7f', 'R5c': '#393b79','A_af1': '#d81b60', 'A_tc4': '#c0392b', 'A_nf4': '#e67e22', 'K1': '#1f5fa8', 'K2': '#2ca7a0', 'K3': '#7b4fb0', 'K0': '#7f7f7f'}

# AthenaK own budget from user.hst (wedge -> full sphere)
HST = {}
for k in KNAMES:
    f = W + '/runs/%s/bsg3d.user.hst' % k
    if os.path.exists(f):
        hdr = [l for l in open(f) if l.startswith('#')][-1]
        cols = {c.split('=')[1]: int(c.split('=')[0].strip('[]'))-1 for c in hdr[1:].split()}
        a = np.loadtxt(f, ndmin=2); a = a[np.all(np.isfinite(a), axis=1)]
        sc = 4*np.pi/OMEGA
        HST[k] = dict(t=a[:, cols['time']]/DAY, Lin=a[:, cols['L_in']]*sc/L0, Ltop=a[:, cols['L_top']]*sc/L0,
                      Etot=a[:, cols['Etot']]*sc)

# ---------- figures ----------
fig, ax = plt.subplots(3, 3, figsize=(16, 13))
for k, ts in TS.items():
    kw = dict(color=COL[k], label=LAB[k], lw=1.4)
    ax[0, 0].semilogy(ts['t'], ts['T1'], **kw); ax[0, 0].semilogy(ts['t'], ts['T0'], color=COL[k], lw=0.7, ls=':')
    ax[0, 1].semilogy(ts['t'], ts['rho1'], **kw); ax[0, 1].semilogy(ts['t'], ts['rho0'], color=COL[k], lw=0.7, ls=':')
    ax[0, 2].plot(ts['t'], ts['F0'], **kw)
    ax[1, 0].plot(ts['t'], (ts['Etot']-ts['Etot'][0])/1e45, **kw)
    ax[1, 1].plot(ts['t'], (ts['U']-ts['U'][0])/1e45, **kw); ax[1, 1].plot(ts['t'], (ts['Egrav']-ts['Egrav'][0])/1e45, color=COL[k], ls='--', lw=1)
    ax[1, 2].plot(ts['t'], ts['rsurf'], **kw)
    ax[2, 0].plot(ts['t'], ts['Ltop'], **kw)
for k, h in HST.items():
    ax[2, 1].plot(h['t'], h['Lin'], color=COL[k], lw=1, label=LAB[k] + ' L_in (hst)')
    cum = np.concatenate([[0], np.cumsum(0.5*(h['Lin'][1:]+h['Lin'][:-1]-h['Ltop'][1:]-h['Ltop'][:-1])*np.diff(h['t']))])*L0*DAY
    ax[2, 2].plot(h['t'], (h['Etot']-h['Etot'][0])/1e45, color=COL[k], lw=1.4, label=LAB[k] + ' dEtot')
    ax[2, 2].plot(h['t'], cum/1e45, color=COL[k], lw=1, ls='--', label='  int(L_in-L_top)')
for k in ANAMES:
    if k in TS:
        ts = TS[k]
        ax[2, 1].plot(ts['t'], ts['F0'], color=COL[k], lw=1, label=LAB[k] + ' cell-0 4pi r^2 F')
        lin = ts['F0']; lt = ts['Ltop'] + ts['Lout_adv']
        cum = np.concatenate([[0], np.cumsum(0.5*(lin[1:]+lin[:-1]-lt[1:]-lt[:-1])*np.diff(ts['t']))])*L0*DAY
        ax[2, 2].plot(ts['t'], (ts['Etot']-ts['Etot'][0])/1e45, color=COL[k], lw=1.4, label=LAB[k] + ' dEtot')
        ax[2, 2].plot(ts['t'], cum/1e45, color=COL[k], lw=1, ls='--', label='  int(L_wall-L_top) (cell)')
tt = ['cell 1 T [K] (dotted: cell 0)', 'cell 1 rho [g/cc] (dotted: cell 0)', 'wall cell-0 4 pi r^2 F_r / L',
      'dE_tot (full sphere) [1e45 erg]', 'dU (solid), dE_grav (dashed) [1e45 erg]', 'r(rho = 1e-10 g/cc) [Rsun]',
      'L_top / L (last cell)', 'bath luminosity L_in / L', 'own budget: dE_tot vs int(L_in - L_top) [1e45 erg]']
for a, t in zip(ax.flat, tt):
    a.set_title(t, fontsize=10); a.set_xlabel('t [d]'); a.grid(alpha=0.3)
ax[0, 2].set_ylim(-500, 100); ax[2, 1].set_ylim(-500, 100); ax[2, 0].set_ylim(-1, 10)
ax[0, 0].legend(fontsize=7); ax[2, 2].legend(fontsize=6); ax[2, 1].legend(fontsize=6)
fig.suptitle('Linhao BSG, identical 1-D setup (20-80 Rsun, 384 cells, his IC + inner heat-bath BC): Athena++ vs AthenaK')
fig.tight_layout(); fig.savefig(W + ('/bcswap_timeseries.png' if _os.environ.get('BCS') else '/emul_timeseries.png' if _os.environ.get('EMUL') else '/h2h_timeseries.png'), dpi=110)

# profiles
fig, ax = plt.subplots(2, 2, figsize=(14, 10))
for k, pr in runs.items():
    ts = np.array([p['t'] for p in pr])/DAY
    for j, td in enumerate((0, 1, 5, 20)):
        if ts[-1] < td - 0.05: continue
        p = pr[int(np.argmin(abs(ts-td)))]
        kw = dict(color=COL[k], lw=1.2, alpha=[1, .8, .6, .45][j], ls=['-', '--', '-.', ':'][j],
                  label='%s %g d' % (k, td))
        ax[0, 0].semilogy(p['r']/RS, p['T'], **kw); ax[0, 1].semilogy(p['r']/RS, p['rho'], **kw)
        ax[1, 0].plot(p['r']/RS, 4*np.pi*p['r']**2*p['F']/L0, **kw); ax[1, 1].plot(p['r']/RS, gamma(p), **kw)
for a, t in zip(ax.flat, ('T [K]', 'rho [g/cc]', 'L_r / L', 'Gamma = kappa_R F / (c g)')):
    a.set_title(t); a.set_xlabel('r [Rsun]'); a.grid(alpha=0.3); a.set_xlim(20, 70)
ax[1, 0].set_ylim(-5, 8); ax[1, 1].set_ylim(-1, 3); ax[0, 0].legend(fontsize=6, ncol=2)
fig.tight_layout(); fig.savefig(W + ('/bcswap_profiles.png' if _os.environ.get('BCS') else '/emul_profiles.png' if _os.environ.get('EMUL') else '/h2h_profiles.png'), dpi=110)

USERK = ('A_af1', 'K1', 'R1', 'R1t1', 'R3', 'R1R3', 'R4', 'R5b', 'R5c', 'R2')
# user figure
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for k in (('A_tc4', 'A_af1', 'Ba1', 'Bb1', 'Bc1', 'Bd1', 'Be0', 'Be1', 'K1') if _os.environ.get('BCS') else USERK if _os.environ.get('EMUL') else ('A_tc4', 'A_af1', 'K1', 'K0')):
    if k not in TS: continue
    ts = TS[k]; kw = dict(color=COL[k], label=LAB[k], lw=1.8 if k in ('A_tc4', 'K1') else 1.0)
    ax[0].semilogy(ts['t'], ts['T1'], **kw)
    ax[1].plot(ts['t'], (ts['Etot']-ts['Etot'][0])/1e45, **kw)
    ax[2].plot(ts['t'], ts['rsurf'], **kw)
ax[0].set_title('Temperature of the 2nd cell above the inner boundary'); ax[0].set_ylabel('T [K]')
ax[1].set_title('Total energy change (gas+rad+kin+grav, full sphere)'); ax[1].set_ylabel('1e45 erg')
ax[2].set_title('Radius where rho = 1e-10 g/cc'); ax[2].set_ylabel('Rsun')
for a in ax: a.set_xlabel('time [days]'); a.grid(alpha=0.3)
ax[0].legend(fontsize=8)
fig.suptitle('Athena++ (his bsg_star.cpp, 1-D) with only the inner boundary swapped; AthenaK K1 for reference' if _os.environ.get('BCS') else 'Athena++ ingredients added to AthenaK one at a time (same 1-D Linhao setup); Athena++ (his pgen) for reference' if _os.environ.get('EMUL') else 'Same 1-D BSG setup (Linhao IC and inner boundary) in two codes')
fig.tight_layout(); fig.savefig(W + ('/bcswap_user.png' if _os.environ.get('BCS') else '/emul_user.png' if _os.environ.get('EMUL') else '/h2h_user.png'), dpi=130)

# Fe bump (38-44 Rsun) comparison
fe = {}
for k, pr in runs.items():
    rows = []
    for p in pr:
        m = (p['r']/RS > 38) & (p['r']/RS < 44)
        rows.append((p['t']/DAY, float(np.mean(p['T'][m])), float(np.mean(p['rho'][m])), float(np.max(gamma(p)[m])),
                     float(np.max(np.abs(p['v'][m]))/1e5)))
    fe[k] = rows
json.dump({k: {kk: vv.tolist() for kk, vv in ts.items()} for k, ts in TS.items()} | {'Fe': fe,
          'hst': {k: {kk: vv.tolist() for kk, vv in h.items()} for k, h in HST.items()}},
          open(W + ('/bcswap.json' if _os.environ.get('BCS') else '/emul.json' if _os.environ.get('EMUL') else '/h2h.json'), 'w'))
for k, ts in TS.items():
    i = np.searchsorted(ts['t'], [1, 5, 10, 20, 30])
    print('%-6s t_end %.2f d' % (k, ts['t'][-1]))
    for j in i:
        if j < len(ts['t']):
            print('   t %5.2f T1 %.3e rho1 %.3e Fwall/L %+8.2f Ltop/L %6.2f rsurf %6.2f dEtot %+.3e dU %+.3e dEg %+.3e'
                  % (ts['t'][j], ts['T1'][j], ts['rho1'][j], ts['F0'][j], ts['Ltop'][j], ts['rsurf'][j],
                     ts['Etot'][j]-ts['Etot'][0], ts['U'][j]-ts['U'][0], ts['Egrav'][j]-ts['Egrav'][0]))
    if k in fe:
        for row in fe[k][::max(1, len(fe[k])//6)]:
            print('   Fe(38-44) t %5.2f <T> %.3e <rho> %.3e maxGamma %.3f max|v| %.1f km/s' % row)
for k, h in HST.items():
    cum = np.trapz(h['Lin']-h['Ltop'], h['t'])*L0*DAY
    print('%s own budget to %.2f d: dEtot %.4e  int(Lin-Ltop) %.4e  mean L_in/L %.2f' % (k, h['t'][-1], h['Etot'][-1]-h['Etot'][0], cum, np.mean(h['Lin'])))
