#!/usr/bin/env python3
"""athpp_fix_1016: summary + figure for the Athena++ implicit-radiation convergence fix.
usage: fixana.py  -> fix_summary.txt, fix_user.png, fix.json"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

F = '/raven/ptmp/jinma/athpp_fix_1016/runs'
A = '/raven/ptmp/jinma/athpp_audit_1015/h2h/runs'
DAYC = 4.547922                     # code time per day
EU = 1e-9*(6.957e10)**3*(6.957e10/1.899769e4)**2
SR = 4*np.pi/((np.cos(np.pi/3)-np.cos(2*np.pi/3))*(np.pi/3))
E45 = EU*SR/1e45                    # code energy (wedge) -> 1e45 erg full sphere
TU = 0.62*1.66053906660e-24*(6.957e10/1.899769e4)**2/1.380649e-16  # p/rho code -> K


def ledger(dirs):
    rows = []
    for d in dirs:
        f = d + '/ledger.txt'
        if os.path.exists(f):
            a = np.loadtxt(f); a = a[np.all(np.isfinite(a), axis=1)]
            rows.append(a)
    if not rows: return None
    a = rows[0]
    for b in rows[1:]:      # restart continuation: offset to the first file's baseline
        b = b[b[:, 0] > a[-1, 0]]
        if len(b) == 0: continue
        j = np.argmin(np.abs(a[:, 0] - b[0, 0]))
        a = np.vstack([a, b])
    b0 = 1 + 2*12
    # stop at the crash (dt collapse / NaN): first row with |E_tot change| > 1e3 x 1e45 erg
    bad = np.where(np.abs((a[:, b0] - a[0, b0])*E45) > 1e3)[0]
    if len(bad): a = a[:bad[0]]
    t = a[:, 0]/DAYC
    E = (a[:, b0] - a[0, b0])*E45
    T = (a[:, b0+1:b0+12] - a[0, b0+1:b0+12])*E45
    return dict(t=t, dE=E, nontel=T[:, 3]-T[:, 4], bnd=T[:, 4], reset=T[:, 6])


def ledger_concat(d1, d2):
    """d2 restarted from d1's end: ledger in d2 starts from zero at the restart time"""
    a = ledger([d1]); b = ledger([d2])
    if b is None: return a
    m = a['t'] < b['t'][0] + 1e-9
    a = {k: v[m] for k, v in a.items()}
    # E_tot measured in d2 is relative to d2's first row: shift by a's value at that time
    off = a['dE'][-1]; offn = a['nontel'][-1]; offb = a['bnd'][-1]; offr = a['reset'][-1]
    return dict(t=np.r_[a['t'], b['t']], dE=np.r_[a['dE'], b['dE'] + off],
                nontel=np.r_[a['nontel'], b['nontel'] + offn], bnd=np.r_[a['bnd'], b['bnd'] + offb],
                reset=np.r_[a['reset'], b['reset'] + offr])


def tabs(d):
    fs = sorted(glob.glob(d + '/bsg.block0.out2.*.tab'))
    t, T1, rs = [], [], []
    for f in fs:
        try:
            tt = float(open(f).readline().split('time=')[1].split()[0])
            x = np.loadtxt(f, comments='#')
        except Exception:
            continue
        if x.ndim != 2 or not np.all(np.isfinite(x[:, 2:5])): continue
        t.append(tt/DAYC); T1.append(x[1, 3]/x[1, 2]*TU)
        s = np.where(x[:, 2] > 0.1)[0]; rs.append(x[s[-1], 1] if len(s) else np.nan)
    o = np.argsort(t)
    return dict(t=np.array(t)[o], T1=np.array(T1)[o], rsurf=np.array(rs)[o])


def conv(d):
    f = d + '/conv.txt'
    if not os.path.exists(f): return None
    a = np.loadtxt(f, comments='#', ndmin=2)
    a = a[np.isfinite(a[:, 5]) & (a[:, 3] > 1e-8)]
    return dict(t=a[:, 2]/DAYC, niter=a[:, 4], cpu=a[:, 7], cell=a[:, 8] if a.shape[1] > 8 else None,
                dt=a[:, 3])


RUNS = {  # name: (ledger dirs, tab dir, conv dir, label)
    'Bd1 default': dict(led=A + '/L_Bd1/out', tab=A + '/L_Bd1/out', conv=F + '/j2/out'),
    'Bd1 converged ref': dict(led=(A + '/C_Bd1/out', F + '/R_Bd1/outC'), tab=(A + '/C_Bd1/out', F + '/R_Bd1/outC'),
                              conv=F + '/R_Bd1/outC'),
    'Bd1 fix': dict(led=F + '/F_Bd1/out', tab=F + '/F_Bd1/out', conv=F + '/F_Bd1/out'),
    'tc4 default': dict(led=A + '/L_tc4/out', tab=A + '/L_tc4/out', conv=None),
    'tc4 fix': dict(led=F + '/F_tc4/out', tab=F + '/F_tc4/out', conv=F + '/F_tc4/out'),
    'af1 default': dict(led=F + '/L_af1/out', tab=F + '/L_af1/out', conv=F + '/L_af1/out'),
    'af1 fix': dict(led=F + '/F_af1/out', tab=F + '/F_af1/out', conv=F + '/F_af1/out'),
}
EXTRA = [a for a in sys.argv[1:]]


def load(k):
    r = RUNS[k]
    if isinstance(r['led'], tuple):
        L = ledger_concat(*r['led'])
        tb = [tabs(d) for d in r['tab']]
        Tb = {kk: np.r_[tb[0][kk], tb[1][kk][tb[1]['t'] > (tb[0]['t'][-1] if len(tb[0]['t']) else -1)]]
              for kk in tb[0]}
    else:
        L = ledger([r['led']]); Tb = tabs(r['tab'])
    C = conv(r['conv']) if r['conv'] else None
    return L, Tb, C


out = open('/raven/ptmp/jinma/athpp_fix_1016/fix_summary.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

D = {}
for k in RUNS:
    try:
        D[k] = load(k)
    except Exception as e:
        P('skip', k, e)
P('# energies 1e45 erg full sphere (ledger, change since t=0); T1 = gas T cell 1 [K]; rsurf = r(rho=1e-10) [Rsun]')
for k, (L, Tb, C) in D.items():
    tend = L['t'][-1] if L is not None else np.nan
    P('== %-18s t_end %.2f d' % (k, tend))
    for tt in [1, 2, 3, 5, 10, 20, 30]:
        if L is None or tt > L['t'][-1] + 0.02: continue
        i = np.argmin(np.abs(L['t'] - tt))
        j = np.argmin(np.abs(Tb['t'] - tt)) if len(Tb['t']) else None
        P('   t %5.2f dE_tot %+9.3f nontelesc %+9.3f boundary %+9.3f resets %+7.3f | T1 %.3e rsurf %6.2f'
          % (L['t'][i], L['dE'][i], L['nontel'][i], L['bnd'][i], L['reset'][i],
             Tb['T1'][j] if j is not None else np.nan, Tb['rsurf'][j] if j is not None else np.nan))
    if C is not None and len(C['t']):
        n = len(C['t'])
        P('   solves %d  mean iterations/solve %.1f  (max %d)  CPU in solves %.1f s over %.2f d -> %.1f s/day'
          % (n, C['niter'].mean(), C['niter'].max(), C['cpu'].sum(), C['t'][-1], C['cpu'].sum()/max(C['t'][-1], 1e-9)))

json.dump({k: {'ledger': {kk: v.tolist() for kk, v in L.items()} if L is not None else None,
               'tab': {kk: v.tolist() for kk, v in Tb.items()}} for k, (L, Tb, C) in D.items()},
          open('/raven/ptmp/jinma/athpp_fix_1016/fix.json', 'w'))

# ---------------- figure
COL = {'default': '#c0392b', 'converged ref': '#7f8c8d', 'fix': '#1f6fb4'}
def col(k):
    for c in COL:
        if k.endswith(c): return COL[c]
    return 'k'
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
for k, (L, Tb, C) in D.items():
    if L is None: continue
    ls = '-' if k.startswith('Bd1') else ('--' if k.startswith('tc4') else ':')
    lw = 2.2 if k.endswith('fix') else (4.5 if 'ref' in k else 1.4)
    a = ax[0, 0] if k.startswith('Bd1') else ax[0, 1]
    a.plot(L['t'], L['dE'], color=col(k), ls=ls, lw=lw, label=k)
    ax[1, 0].plot(L['t'], L['nontel'], color=col(k), ls=ls, lw=lw, label=k)
    if len(Tb['t']):
        ax[1, 1].semilogy(Tb['t'], Tb['T1'], color=col(k), ls=ls, lw=lw, label=k)
        ax[1, 2].plot(Tb['t'], Tb['rsurf'], color=col(k), ls=ls, lw=lw, label=k)
    if C is not None and len(C['t']):
        # iterations per solve (Jacobi sweeps for default; line-solve+sweep iterations for the fix), 50-solve running mean
        nw = 50
        if len(C['t']) > nw:
            rm = np.convolve(C['niter'], np.ones(nw)/nw, mode='valid'); tt = C['t'][nw-1:]
            cpd = C['cpu'].sum()/max(C['t'][-1], 1e-9)
            ax[0, 2].semilogy(tt, rm, color=col(k), ls=ls, lw=lw, label='%s (%.0f s CPU/day in solves)' % (k, cpd))
ax[0, 0].set_title('Closed-wall Bd1 (af1, taucell 4): dE_tot [1e45 erg]')
ax[0, 1].set_title("His BC (tc4: af0 / af1, taucell 4): dE_tot [1e45 erg]")
ax[0, 2].set_title('Iterations per solve (default: Jacobi sweeps; fix: line solve + 1 sweep,\n1 fix iteration ~ 6 sweeps of CPU), 50-solve mean', fontsize=9)
ax[1, 0].set_title('Non-telescoping rad. transport (energy destroyed) [1e45 erg]')
ax[1, 1].set_title('Gas T in cell 1 [K]')
ax[1, 2].set_title('r(rho = 1e-10 g/cc) [Rsun]')
for a in ax.flat:
    a.set_xlabel('t [d]'); a.grid(alpha=0.3); a.legend(fontsize=7)
for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
    a.set_xlim(0, 30)
ax[0, 0].set_ylim(-130, 40); ax[0, 1].set_ylim(-170, 40)
ax[0, 0].text(3.0, 25, 'fix and converged ref die at 2.70 / 2.68 d (base wall instability)', fontsize=8)
fig.suptitle('Athena++ implicit radiation: default Jacobi (nlimit 100, 1e-6) vs converged reference (1e-11) vs fix '
             '(sweep=2 x1-line solve + telescoping billing)', fontsize=11)
fig.tight_layout()
fig.savefig('/raven/ptmp/jinma/athpp_fix_1016/fix_user.png', dpi=110)
P('figure written')
