"""dhj round-1 progress page: reductions (cached per dump), figures, index.html.
Runs inside allocation 5274196 (CPU-only step) via refresh.sh. Reads run dirs read-only."""
import sys, os, glob, re, json, time, subprocess, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rochedhj_ana'); import dhjcs
MODE = sys.argv[1] if len(sys.argv) > 1 else 'full'
if MODE != 'html':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from multiprocessing import Pool

J = '/resnick/groups/carnegie_poc/jingze'
PAGE = '/resnick/home/jingze/ATHENAK/pages/dhj_r1'
IMG = PAGE + '/img'
CACHE = J + '/dhj_p1_page_cache'           # p1 production (11-01 ~15:30 PST)
OLDCACHE = J + '/dhj_r1b_page_cache'       # r1b/r1bo test runs, for the grey overlay
NOTE = 'production started 11-01 ~15:30 PST, setup final (see hydro_prep_1101)'
EOSF = J + '/w121r1_1101/anchor/eos_1x.txt'
NROT = 300
PROF = np.arange(-6.0, 2.01, 0.125)            # log10 p [bar]
LEV = {'1e-3': -3.0, '1e-1': -1.0}
NPROC = int(os.environ.get('NPROC', '14'))

R1 = J + '/r1planets_1101'
RUNS = [('w121h', 'WASP-121b (hydro)', J + '/w121r1_1101/runs/p1'),
        ('w76', 'WASP-76b', R1 + '/wasp76b/runs/p1'),
        ('toi1518', 'TOI-1518b', R1 + '/toi1518b/runs/p1'),
        ('w178', 'WASP-178b', R1 + '/wasp178b/runs/p1'),
        ('w189', 'WASP-189b', R1 + '/wasp189b/runs/p1')]
OLD = {'w121h': J + '/w121r1_1101/runs/r1b', 'w76': R1 + '/wasp76b/runs/r1b', 'toi1518': R1 + '/toi1518b/runs/r1b',
       'w178': R1 + '/wasp178b/runs/r1b', 'w189': R1 + '/wasp189b/runs/r1bo'}
OLDKEY = {'w189': 'w189o'}   # r1bo (outflow1 top, as p1) was cached under w189o
# WASP-121b MHD: any r1m* dir with dumps that is not stopped (or was relaunched after STOPPED.txt)
for d in sorted(glob.glob(J + '/w121r1_1101/runs/r1*m*')):
    logs = glob.glob(d + '/run.ia*.log'); st = d + '/STOPPED.txt'
    if not glob.glob(d + '/bin/*.bin') or not logs: continue
    if os.path.exists(st) and max(os.path.getmtime(l) for l in logs) <= os.path.getmtime(st) + 60: continue
    RUNS.append(('w121' + os.path.basename(d)[2:], 'WASP-121b (MHD, %s)' % os.path.basename(d), d))
BUDGET_KEYS = ['w121h', 'w76', 'toi1518', 'w178', 'w189']   # hydro p1 runs (conservation budgets section)
COL = {'w121h': '#1f4e79', 'w76': '#b5541c', 'toi1518': '#2e7d4f', 'w178': '#7b3f8c', 'w189': '#a3822a', 'w189o': '#5a6e7f'}

def deck(fn):
    d = {}
    for l in open(fn):
        m = re.match(r'^\s*(\w+)\s*=\s*([^#\s]+)', l)
        if m and m.group(1) not in d: d[m.group(1)] = m.group(2)
    return d

def hst(fn):
    lab = None
    for l in open(fn):
        if l.startswith('#  ['):
            lab = [s.split('=')[1] for s in l[1:].split()]; break
    rows = [l.split() for l in open(fn) if l.strip() and not l.startswith('#')]
    nc = len(lab) if lab else max(len(r) for r in rows)
    a = np.array([[float(x) for x in r] for r in rows if len(r) == nc], ndmin=2)   # 11-02: skip torn/partial lines
    t = a[:, 0]; _, idx = np.unique(t[::-1], return_index=True); idx = len(t) - 1 - idx   # last write wins
    a = a[np.sort(idx)]
    return {k: a[:, i] for i, k in enumerate(lab)}

# ------------------------------------------------------------------ per-dump reduction
_G = {}
def geom(raw):
    key = (raw['Nx1'] if 'Nx1' in raw else 0, len(raw['mb_geometry']), raw['nx2_mb'], raw['nx3_mb'], float(np.asarray(raw['mb_x1v'])[0, 0]))
    if key in _G: return _G[key]
    geo = dhjcs.cs_geometry(raw); g = np.asarray(raw['mb_geometry']); nmb = g.shape[0]; nj = raw['nx2_mb']; nk = raw['nx3_mb']
    A = np.zeros((nmb, nk, nj)); XYZ = np.zeros((nmb, nk, nj, 3))
    for m in range(nmb):
        xi = np.pi / 4 * (g[m, 2] + (g[m, 3] - g[m, 2]) * (np.arange(nj) + 0.5) / nj)
        eta = np.pi / 4 * (g[m, 4] + (g[m, 5] - g[m, 4]) * (np.arange(nk) + 0.5) / nk)
        XI, ETA = np.meshgrid(xi, eta); x = np.tan(XI); y = np.tan(ETA)
        A[m] = (1 + x * x) * (1 + y * y) / (1 + x * x + y * y) ** 1.5 * (np.pi / 4) ** 2 * (g[m, 3] - g[m, 2]) * (g[m, 5] - g[m, 4]) / (nj * nk)
        XYZ[m] = dhjcs.panel_map(m // (nmb // 6), XI, ETA)
    A /= A.sum()
    X, Y, Z = XYZ[..., 0], XYZ[..., 1], XYZ[..., 2]; rh = np.maximum(np.sqrt(X * X + Y * Y), 1e-12)
    G = dict(geo=geo, A=A, X=X, Y=Y, Z=Z, zx=-Y / rh, mx=-Z * X / rh,       # x-components of e_east, e_north
             lat=np.degrees(geo['lat']), lon=np.degrees(geo['lon']),
             term=np.abs(X) < 0.1, psi=np.degrees(np.arctan2(Z, -Y)))     # psi 0 = evening equator, +-180 morning
    _G[key] = G; return G

_EOS = None
def reduce_dump(args):
    path, P, want_fields = args
    global _EOS
    if _EOS is None: _EOS = dhjcs.EOS(EOSF)
    raw = dhjcs.bin_convert.read_binary(path); mb = raw['mb_data']; G = geom(raw)
    rho = np.asarray(mb['dens']); T, p = _EOS.invert(rho, np.asarray(mb['eint']))
    vr = np.asarray(mb['velx']); u, v = dhjcs.winds(raw, G['geo'])
    lp = np.log10(p) + 0.0
    A = G['A']; t = float(raw['time'])
    r = {'t': t, 'rot': t / P, 'nan': int((~np.isfinite(rho)).sum() + (~np.isfinite(T)).sum())}
    ul = dhjcs.level_interp(u, lp, PROF + 6)
    la, zm = dhjcs.zonal_mean(ul, G['geo']['lat'], dlat=6)
    eq = np.abs(la) < 15; ueq = np.nanmean(zm[:, eq], axis=1); ok = PROF <= 1.0
    k = int(np.nanargmax(np.where(ok & np.isfinite(ueq), ueq, -np.inf)))
    r['jet'] = float(ueq[k] / 100); r['jet_p'] = float(10 ** PROF[k])
    lv = [LEV[s] + 6 for s in LEV]
    Tl = dhjcs.level_interp(T, lp, lv); uL = dhjcs.level_interp(u, lp, lv); vL = dhjcs.level_interp(v, lp, lv)
    vrL = dhjcs.level_interp(vr, lp, lv)
    vx = vrL * G['X'] + uL * G['zx'] + vL * G['mx']                       # +x = away from the star = toward the transit observer
    def wm(f, msk):
        w = A * msk * np.isfinite(f); return float(np.nansum(np.where(np.isfinite(f), f, 0) * w) / max(w.sum(), 1e-300))
    day = G['X'] < 0; night = G['X'] > 0
    for i, s in enumerate(LEV):
        r['dT_' + s] = wm(Tl[i], day) - wm(Tl[i], night)
        r['Tday_' + s] = wm(Tl[i], day); r['Tnight_' + s] = wm(Tl[i], night)
        r['vterm_' + s] = wm(vx[i], G['term']) / 100
        r['vterm_ev_' + s] = wm(vx[i], G['term'] & (G['Y'] < 0)) / 100
        r['vterm_mo_' + s] = wm(vx[i], G['term'] & (G['Y'] > 0)) / 100
    if want_fields:
        F = {'la': la, 'zm': zm / 100, 'prof': PROF}
        for i, s in enumerate(LEV):
            lats, lons, m3 = dhjcs.latlon_bin(np.stack([Tl[i], uL[i] / 100, vL[i] / 100]), G['geo']['lat'], G['geo']['lon'], dlat=4, dlon=4)
            F['lats'], F['lons'] = lats, lons; F['T_' + s], F['u_' + s], F['v_' + s] = m3
            b = np.arange(-180, 181, 10); c = 0.5 * (b[1:] + b[:-1]); f = vx[i] / 100; tm = G['term'] & np.isfinite(f)
            F['psi'] = c
            F['ring_' + s] = np.array([np.nan if not (tm & (G['psi'] >= b[j]) & (G['psi'] < b[j + 1])).any() else
                                       np.sum((A * f)[tm & (G['psi'] >= b[j]) & (G['psi'] < b[j + 1])]) /
                                       np.sum(A[tm & (G['psi'] >= b[j]) & (G['psi'] < b[j + 1])]) for j in range(len(c))])
        r['_fields'] = F
    return os.path.basename(path), r

# ------------------------------------------------------------------ status from logs
def status(run, P):
    s = {'flags': []}
    logs = sorted(glob.glob(run + '/run.ia*.log'), key=os.path.getmtime)
    if not logs: s['flags'].append('no run log'); return s
    lg = logs[-1]; s['log'] = lg; s['age_min'] = (time.time() - os.path.getmtime(lg)) / 60
    tail = subprocess.run(['tail', '-n', '400', lg], capture_output=True, text=True).stdout.splitlines()
    rows = []
    for l in tail:
        m = re.match(r'elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)', l)
        if m: rows.append([float(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))])
    s['last_line'] = tail[-1].strip() if tail else ''
    if rows:
        a = np.array(rows); s['elapsed'] = a[-1, 0]; s['cycle'] = int(a[-1, 1]); s['time'] = a[-1, 2]; s['dt'] = a[-1, 3]
        n = min(len(a) - 1, 60)
        if n > 0:
            de = a[-1, 0] - a[-1 - n, 0]
            s['ms_cycle'] = 1e3 * de / (a[-1, 1] - a[-1 - n, 1])
            s['rot_per_h'] = (a[-1, 2] - a[-1 - n, 2]) / P / (de / 3600)
            s['eta_h'] = max(NROT * P - a[-1, 2], 0) / P / s['rot_per_h']
    pat = r'\bnan\b|\binf\b|abort|error|segmentation|killed|cancelled|terminat|non-converged'
    for f in [lg] + glob.glob(run + '/logs/*.out'):
        g = subprocess.run(['grep', '-ciE', pat, f], capture_output=True, text=True).stdout.strip()
        if g and g != '0': s['flags'].append('%s: %s line(s) match nan/abort/error' % (os.path.basename(f), g))
    if os.path.exists(run + '/STOPPED.txt'): s['flags'].append('STOPPED.txt present')
    if s['age_min'] > 15: s['flags'].append('log not written for %.0f min' % s['age_min'])
    # floor event counters (dhj.log, last row = last interval)
    try:
        hdr = None; last = None
        with open(run + '/dhj.log') as fh:
            for l in fh:
                if l.startswith('#  cycle'): hdr = l[1:].split()
                elif not l.startswith('#') and l.strip(): last = l
        s['floors'] = dict(zip(hdr, last.split()))
    except Exception:
        s['floors'] = {}
    return s

# ------------------------------------------------------------------ figures
if MODE != 'html': plt.rcParams.update({'font.size': 9.5, 'axes.titlesize': 10, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
                     'savefig.facecolor': 'white', 'axes.grid': True, 'grid.color': '#e3e3e3', 'grid.linewidth': 0.6,
                     'axes.spines.top': False, 'axes.spines.right': False})
C1, C2 = '#1f4e79', '#b5541c'

def fig_snapshot(key, name, F, r, fn):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.4))
    zm = F['zm']; lv = max(np.nanmax(np.abs(zm[F['prof'] <= 1.5])), 1.0); st = 10 ** np.floor(np.log10(lv / 2)) * (1 if lv / 10 ** np.floor(np.log10(lv / 2)) < 10 else 5); lv = np.ceil(lv / st) * st
    cf = ax[0, 0].contourf(F['la'], 10 ** F['prof'], zm, levels=np.linspace(-lv, lv, 21), cmap='RdBu_r', extend='both')
    ax[0, 0].contour(F['la'], 10 ** F['prof'], zm, levels=[0], colors='k', linewidths=0.6)
    ax[0, 0].set_yscale('log'); ax[0, 0].invert_yaxis(); ax[0, 0].set_xlabel('latitude [deg]'); ax[0, 0].set_ylabel('p [bar]')
    ax[0, 0].set_title('zonal-mean zonal wind u [m/s]'); fig.colorbar(cf, ax=ax[0, 0]); ax[0, 0].grid(False)
    a = ax[0, 1]
    a.plot(F['psi'], F['ring_1e-3'], c=C1, marker='o', ms=3, label='1e-3 bar')
    a.plot(F['psi'], F['ring_1e-1'], c=C2, marker='o', ms=3, ls='--', label='1e-1 bar')
    a.axhline(0, c='k', lw=0.6); a.set_xlim(-180, 180); a.set_xticks([-180, -90, 0, 90, 180])
    a.set_xticklabels(['morning\neq.', 'S pole', 'evening\neq.', 'N pole', 'morning\neq.'])
    a.set_ylabel('v toward observer [m/s]  (+ = blueshift)')
    a.set_title('terminator line-of-sight wind (|mu| < 0.1); 1e-3 bar mean %+.0f m/s' % r['vterm_1e-3']); a.legend(fontsize=8)
    for j, s in enumerate(['1e-3', '1e-1']):
        a = ax[1, j]; T = F['T_' + s]
        pm = a.pcolormesh(F['lons'], F['lats'], T, cmap='inferno', shading='nearest'); fig.colorbar(pm, ax=a, label='T [K]')
        sl = (slice(1, None, 3), slice(1, None, 3)); U = F['u_' + s][sl]; V = F['v_' + s][sl]
        sp = max(np.nanmax(np.hypot(U, V)), 1.0); LO, LA = np.meshgrid(F['lons'], F['lats'])
        q = a.quiver(LO[sl], LA[sl], U, V, color='#5ec8e6', scale=sp * 14, width=0.0032)
        a.quiverkey(q, 0.86, 1.03, round(sp, -2) if sp > 100 else sp, '%.0f m/s' % (round(sp, -2) if sp > 100 else sp), labelpos='W', fontproperties={'size': 8})
        for x in (-90, 90): a.axvline(x, c='w', lw=0.7, ls='--')
        a.set_xlim(-180, 180); a.set_ylim(-90, 90); a.set_xticks(range(-180, 181, 60)); a.set_yticks(range(-90, 91, 30))
        a.set_xlabel('longitude from substellar [deg] (+ east)'); a.set_ylabel('latitude [deg]'); a.grid(False)
        a.set_title('T and wind at %s bar' % s, loc='left')
    fig.suptitle('%s at %.1f rotations' % (name, r['rot']), fontsize=11.5, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(fn, dpi=90); plt.close(fig)

def runmean(t, y, w):
    if len(t) < 3: return y
    c = np.concatenate([[0], np.cumsum(y)]); out = np.empty_like(y)
    lo = np.searchsorted(t, t - w / 2); hi = np.searchsorted(t, t + w / 2, side='right')
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)

GREY = dict(c='#b9b9b9', lw=0.8, zorder=0)
def fig_series(key, name, H, scal, om, fn, Ho=None, so=None):
    fig, ax = plt.subplots(4, 2, figsize=(12, 12.4), sharex=True); ax = ax.ravel()
    rot = H['rot']; W = min(1.0, max(rot[-1] / 10, 0.02))   # running-mean window [rot]
    if Ho is not None:   # pre-fix run, same normalisation by its own t = 0 values
        ax[0].plot(Ho['rot'], (Ho['M'] - Ho['M'][0]) / Ho['M'][0], **GREY, label='earlier test run (r1b; W189: r1bo)')
        ax[1].plot(Ho['rot'], (Ho['E'] - Ho['E'][0]) / abs(Ho['E'][0]), **GREY)
        if Ho.get('Lrel') is not None: ax[2].plot(Ho['rot'], Ho['Lrel'] / (Ho['Iw'][0] * om), **GREY, label='earlier test run (r1b; W189: r1bo)')
        if Ho.get('Ttide') is not None: ax[3].plot(Ho['rot'], runmean(Ho['rot'], Ho['Ttide'], W) * 2 * np.pi / om / (Ho['Iw'][0] * om), **GREY)
        if Ho.get('Mfloor') is not None: ax[4].plot(Ho['rot'], runmean(Ho['rot'], Ho['Mfloor'], W), **GREY)
        if Ho.get('Mdot_top') is not None: ax[5].plot(Ho['rot'], runmean(Ho['rot'], Ho['Mdot_top'], W), **GREY)
        ax[6].plot(Ho['rot'], Ho['dt'], **GREY)
        if so: ax[7].plot([x['rot'] for x in so], [x['jet'] for x in so], **GREY, label='test-run jet')
    ax[0].plot(rot, (H['M'] - H['M'][0]) / H['M'][0], c=C1, lw=0.9); ax[0].set_ylabel('(M - M0) / M0'); ax[0].set_title('mass')
    ax[1].plot(rot, (H['E'] - H['E'][0]) / abs(H['E'][0]), c=C1, lw=0.9); ax[1].set_ylabel('(E - E0) / |E0|'); ax[1].set_title('total energy (hst tot-E)')
    if H.get('Lrel') is not None:
        ax[2].plot(rot, H['Lrel'] / (H['Iw'][0] * om), c=C1, lw=0.9)
    ax[2].set_ylabel('L_rel,z / (I_w0 Omega)'); ax[2].set_title('L_z of the winds (rotating frame)')
    if Ho is not None: ax[2].legend(fontsize=8)
    Pr = 2 * np.pi / om
    if H.get('Ttide') is not None:
        y = H['Ttide'] * Pr / (H['Iw'][0] * om)
        ax[3].plot(rot, y, c='#9db3c8', lw=0.5, label='per hst row'); ax[3].plot(rot, runmean(rot, y, W), c=C1, lw=1.2, label='%.2g-rot running mean' % W)
        ax[3].legend(fontsize=8)
    ax[3].axhline(0, c='k', lw=0.5)
    ax[3].set_ylabel('T_tide P / (I_w0 Omega)'); ax[3].set_title('tide torque about z (hst Ttide), per rotation')
    if H.get('Mfloor') is not None:
        ax[4].plot(rot, H['Mfloor'], c='#9db3c8', lw=0.5, label='per hst row')
        ax[4].plot(rot, runmean(rot, H['Mfloor'], W), c=C1, lw=1.2, label='%.2g-rot running mean' % W)
        ax[4].legend(fontsize=8)
    ax[4].set_ylabel('Mfloor [g/s]'); ax[4].set_title('mass added by floors (rate)')
    if H.get('Mdot_top') is not None:
        ax[5].plot(rot, H['Mdot_top'], c='#9db3c8', lw=0.5); ax[5].plot(rot, runmean(rot, H['Mdot_top'], W), c=C1, lw=1.2)
    ax[5].set_ylabel('Mdot_top [g/s]'); ax[5].set_title('mass flux through the top (hst Mdot_top)')
    ax[6].plot(rot, H['dt'], c=C1, lw=0.7); ax[6].set_ylabel('dt [s]'); ax[6].set_title('time step')
    if scal:
        sr = np.array([x['rot'] for x in scal])
        ax[7].plot(sr, [x['jet'] for x in scal], c=C1, marker='.', ms=3, lw=0.9, label='eq. jet max (|lat|<15, p<=10 bar)')
        ax[7].plot(sr, [x['vterm_1e-3'] for x in scal], c=C2, marker='.', ms=3, lw=0.9, label='terminator wind, 1e-3 bar')
        ax[7].legend(fontsize=8)
    ax[7].set_ylabel('m/s'); ax[7].set_title('equatorial jet and terminator wind (from dumps)')
    if Ho is not None and len(rot): 
        for a in ax: a.set_xlim(-0.02 * max(rot[-1], 1), max(rot[-1], 1) * 1.02)
    for a in ax[6:]: a.set_xlabel('rotations')
    fig.suptitle('%s: history' % name, fontsize=11.5, x=0.01, ha='left')
    fig.tight_layout(); fig.savefig(fn, dpi=90); plt.close(fig)

def fig_cross(rows, fn):
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
    for k, name, teq, sc in rows:
        last = sc[-10:]
        if not last: continue
        c = COL.get(k, '#444444'); R = last[-1]['rot']
        def pt(a, key, mk='o', lab=True, fill=True):
            y = np.array([x[key] for x in last]); m = y.mean()
            a.errorbar(teq, m, yerr=[[m - y.min()], [y.max() - m]], fmt=mk, c=c, mfc=c if fill else 'white', ms=6, capsize=2,
                       label=('%s (%.0f rot)' % (name, R)) if lab else None)
        pt(ax[0], 'jet'); pt(ax[1], 'dT_1e-3'); pt(ax[1], 'dT_1e-1', mk='s', lab=False, fill=False); pt(ax[2], 'vterm_1e-3')
    ax[0].set_title('equatorial jet max [m/s]'); ax[1].set_title('day - night T [K]: filled 1e-3 bar, open 1e-1 bar')
    ax[2].set_title('terminator wind at 1e-3 bar [m/s] (+ = toward observer)')
    for a in ax: a.set_xlabel('T_eq [K]')
    ax[2].axhline(0, c='k', lw=0.6)
    h, l = ax[0].get_legend_handles_labels(); fig.legend(h, l, loc='lower center', ncol=len(l), fontsize=8.5, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(fn, dpi=90); plt.close(fig)

def loadH(run, P):
    hf = [f for f in glob.glob(run + '/dhj.*.hst') if not f.endswith('user.hst')][0]
    h = hst(hf); u = hst(run + '/dhj.user.hst'); n = min(len(h['time']), len(u['time']))
    g = lambda k: u[k][:n] if k in u else None
    return {'rot': h['time'][:n] / P, 'M': h['mass'][:n], 'E': h['tot-E'][:n], 'dt': h['dt'][:n],
            'Lrel': g('Lrel'), 'Iw': g('Iw'), 'Mfloor': g('Mfloor'), 'Ttide': g('Ttide'), 'Mdot_top': g('Mdot_top')}

# ------------------------------------------------------------------ main
def main():
    t0 = time.time(); os.makedirs(IMG, exist_ok=True)
    info = {}; tasks = []
    for key, name, run in RUNS:
        d = deck(run + '/run.athinput'); om = float(d['omega']); P = 2 * np.pi / om
        cdir = '%s/%s' % (CACHE, key); os.makedirs(cdir, exist_ok=True)
        cf = cdir + '/scal.json'; sc = json.load(open(cf)) if os.path.exists(cf) else {}
        bins = sorted(glob.glob(run + '/bin/*.bin'))
        # skip a dump still being written (younger than 60 s)
        # pagefast-1102: only the newest few files can be younger than 60 s: stat those only (was ~4600 NFS stats, ~4 s)
        _pre = [b.rsplit('.', 2)[0] for b in bins]
        _tail = {i for i in range(len(bins)) if _pre[i:i + 5].count(_pre[i]) < 5}   # last 4 of each <prefix>.NNNNN.bin series
        bins = [b for i, b in enumerate(bins) if i not in _tail or time.time() - os.path.getmtime(b) > 60]
        latest = os.path.basename(bins[-1]) if bins else None
        fstate = cdir + '/fields_name.txt'
        have = open(fstate).read().strip() if os.path.exists(fstate) else ''
        for b in bins:
            bn = os.path.basename(b); wf = (bn == latest and have != latest)
            if bn not in sc or wf: tasks.append((key, (b, P, wf)))
        info[key] = dict(name=name, run=run, d=d, om=om, P=P, sc=sc, cf=cf, cdir=cdir, latest=latest)
    if MODE == 'html': tasks = []
    print('tasks: %d dumps to reduce' % len(tasks), flush=True)
    newsnap = set()
    if tasks:
        with Pool(min(NPROC, len(tasks))) as pool:
            for (key, _), (bn, r) in zip(tasks, pool.imap(reduce_dump, [t[1] for t in tasks], chunksize=1)):
                F = r.pop('_fields', None); info[key]['sc'][bn] = r
                if F is not None:
                    np.savez(info[key]['cdir'] + '/fields.npz', **F); open(info[key]['cdir'] + '/fields_name.txt', 'w').write(bn)
                    newsnap.add(key)
    if MODE != 'html':
        for key, I in info.items(): json.dump(I['sc'], open(I['cf'], 'w'))
    print('reductions done %.0f s' % (time.time() - t0), flush=True)

    cross = []
    for key, I in info.items():
        run = I['run']; d = I['d']; P = I['P']
        sc = [I['sc'][k] for k in sorted(I['sc'])]
        I['scl'] = sc
        H = loadH(run, P); I['H'] = H
        if MODE != 'html':
            Ho = so = None
            if key in OLD and os.path.exists(OLD[key] + '/dhj.user.hst'):
                try:
                    Ho = loadH(OLD[key], P); xm = max(H['rot'][-1], 1.0) * 1.02; m = Ho['rot'] <= xm
                    Ho = {k: (v[m] if v is not None else None) for k, v in Ho.items()}
                    of = '%s/%s/scal.json' % (OLDCACHE, OLDKEY.get(key, key))
                    if os.path.exists(of): so = [v for _, v in sorted(json.load(open(of)).items()) if v['rot'] <= xm]
                except Exception as e: print('old overlay failed', key, e)
            fig_series(key, I['name'], H, sc, I['om'], '%s/%s_series.png' % (IMG, key), Ho, so)
        snap = '%s/%s_snap.png' % (IMG, key)
        if MODE != 'html' and os.path.exists(I['cdir'] + '/fields.npz') and (key in newsnap or not os.path.exists(snap)):
            F = dict(np.load(I['cdir'] + '/fields.npz')); bn = open(I['cdir'] + '/fields_name.txt').read().strip()
            fig_snapshot(key, I['name'], F, I['sc'][bn], snap)
        I['st'] = status(run, P)
        cross.append((key, I['name'], float(d['Teq']), sc))
    if MODE != 'html': fig_cross(cross, IMG + '/cross_teq.png')
    # 11-02 USER: conservation budgets (mass, energy, L_z) of the hydro p1 runs -> budget.py
    if MODE != 'html':
        try:
            import budget
            bdir = CACHE + '/budget'; os.makedirs(bdir, exist_ok=True); Bs = []
            for key in BUDGET_KEYS:
                if key not in info: continue
                I = info[key]
                B = budget.analyse(key, I['name'], I['run'], I['om'], int(I['d']['nx1']), bdir)
                budget.fig_budget(I['name'], B, '%s/budget_%s.png' % (IMG, key)); Bs.append((key, I['name'], B))
            S = budget.summary(Bs)
            for key, name, B in Bs:   # implied top zonal wind if the residual L_z drift is the unbooked top flux
                r = S[key]; md = r['Mdot_top_mean']
                if 'L_res_spin' in r and md > 0:
                    r['Ttop_v'] = -r['L_res_spin'] * B['L']['scale'] / B['t'][-1] / (md * float(info[key]['d']['x1max'])) / 1e5
            budget.fig_summary(S, IMG + '/budget_summary.png')
            json.dump(S, open(bdir + '/summary.json', 'w'), indent=1)
            print('budgets done', flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); print('budgets FAILED', e, flush=True)
    build_html(info)
    tot = ' '.join('%s %.1f' % (k, I['st'].get('time', 0) / I['P']) for k, I in info.items())
    print('done %.0f s' % (time.time() - t0), flush=True)
    print('dhj_r1 READY %s rot (of %d)' % (tot, NROT), flush=True)

# ------------------------------------------------------------------ html
def rel(h):
    if h is None or not np.isfinite(h): return '-'
    if h < 1: return 'in ~%.0f min' % (60 * h)
    if h < 48: return 'in ~%.0f h' % h
    return 'in ~%.1f d' % (h / 24)

def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

CSS = """
:root{--bg:#fbfbf9;--fg:#1d2228;--muted:#5b636b;--line:#dfe1dc;--card:#ffffff;--accent:#1f4e79;--warn:#9a3412;--ok:#2e6b3f;
--mono:'IBM Plex Mono',ui-monospace,Menlo,Consolas,monospace;--sans:'Source Sans 3','Source Sans Pro',system-ui,-apple-system,'Segoe UI',sans-serif;color-scheme:light}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15181b;--fg:#e4e6e8;--muted:#9aa2aa;--line:#2c3237;--card:#1c2024;--accent:#8fb8de;--warn:#f0a070;--ok:#8fd19e;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#15181b;--fg:#e4e6e8;--muted:#9aa2aa;--line:#2c3237;--card:#1c2024;--accent:#8fb8de;--warn:#f0a070;--ok:#8fd19e;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font-family:var(--sans);font-size:16px;line-height:1.5}
.wrap{max-width:1100px;margin:0 auto;padding-block:24px 48px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
h1{font-size:1.6rem;margin:0 0 4px;font-weight:600;letter-spacing:-0.01em}
h2{font-size:1.2rem;margin:36px 0 8px;font-weight:600;border-bottom:1px solid var(--line);padding-bottom:4px}
h3{font-size:1.05rem;margin:20px 0 6px;font-weight:600}
p{margin:6px 0}
.sub{color:var(--muted);font-size:.95rem}
.tw{overflow-x:auto;-webkit-overflow-scrolling:touch;border:1px solid var(--line);border-radius:6px;background:var(--card)}
table{border-collapse:collapse;font-variant-numeric:tabular-nums;font-size:.9rem;width:100%}
th,td{padding:6px 10px;text-align:right;white-space:nowrap;border-bottom:1px solid var(--line)}
th{font-weight:600;color:var(--muted);font-size:.8rem;text-transform:none;background:transparent}
td:first-child,th:first-child{text-align:left}
tr:last-child td{border-bottom:none}
.ok{color:var(--ok)}.warn{color:var(--warn)}
figure{margin:12px 0 20px}
figure img{max-width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:4px;background:#fff}
figcaption{color:var(--muted);font-size:.88rem;margin-top:6px}
code,.mono{font-family:var(--mono);font-size:.82rem}
.log{font-family:var(--mono);font-size:.78rem;color:var(--muted);overflow-x:auto;white-space:pre;padding:6px 0}
ul.facts{margin:4px 0 8px;padding-left:20px}
ul.facts li{margin:2px 0}
nav a{color:var(--accent);margin-right:14px;text-decoration:none;white-space:nowrap}
"""

def build_html(info):
    els = [I['st'].get('elapsed', 0) for I in info.values()]
    upd = max(els) / 3600 if els else 0
    o = ['<title>Hot Jupiter Round 1</title>', '<meta name="claude:viewport-fit-sides" content="cover">',
         '<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>',
         '<link href="https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono&display=swap" rel="stylesheet">',
         '<style>' + CSS + '</style>', '<div class="wrap">',
         '<h1>Hot Jupiter Round 1</h1>',
         '<p class="sub">Deep-hot-Jupiter (dhj) AthenaK runs on Caltech: 1x solar metallicity, cubed sphere, host-star Roche tide, '
         '1 H200 per run, target %d rotations. <strong>Hydro production p1: %s.</strong> Time series show the stopped r1b test runs (WASP-189b: the outflow1-top test r1bo) as a faint grey overlay. Updated ~%.1f h into the runs (wall time of the longest-running log).</p>' % (NROT, NOTE, upd),
         '<nav>' + '<a href="#status">status</a><a href="#cross">cross-planet</a>' +
         ''.join('<a href="#%s">%s</a>' % (k, esc(I['name'])) for k, I in info.items()) + '<a href="#budget">conservation budgets</a><a href="#fehires">Fe I vs Seidel</a><a href="#obs">other observations</a></nav>',
         '<h2 id="status">Status</h2><div class="tw"><table><thead><tr><th>run</th><th>rotations</th><th>sim time</th><th>wall</th>'
         '<th>ms/cycle</th><th>rot/h</th><th>300 rot</th><th>dt now</th><th>floor events, last log interval<br>dens / eint / T</th><th>health</th></tr></thead><tbody>']
    for k, I in info.items():
        s = I['st']; rot = s.get('time', np.nan) / I['P']; fl = s.get('floors', {})
        hl = '<span class="ok">ok</span>' if not s['flags'] else '<span class="warn">' + esc('; '.join(s['flags'])) + '</span>'
        nn = sum(x.get('nan', 0) for x in I['scl'][-1:])
        if nn: hl += ' <span class="warn">NaN in last dump (%d)</span>' % nn
        o.append('<tr><td>%s</td><td>%.1f / %d</td><td>%.3g s = %.1f d</td><td>%.1f h</td><td>%.1f</td><td>%.2f</td><td>%s</td><td>%.2f s</td><td>%s / %s / %s</td><td>%s</td></tr>' % (
            esc(I['name']), rot, NROT, s.get('time', np.nan), s.get('time', np.nan) / 86400, s.get('elapsed', np.nan) / 3600,
            s.get('ms_cycle', np.nan), s.get('rot_per_h', np.nan), rel(s.get('eta_h')), s.get('dt', np.nan),
            fl.get('eos_dfloor', '-'), fl.get('eos_efloor', '-'), fl.get('eos_tfloor', '-'), hl))
    o.append('</tbody></table></div>')
    o.append('<p class="sub">ms/cycle, rot/h and the projection use the last ~6000 cycles of each log. Health flags: log lines matching '
             'nan/abort/error/terminate, STOPPED.txt, or no log write for 15 min. Floor events are the counters of the last row of '
             '<code>dhj.log</code> (cells hit per log interval).</p>')
    o.append('<h3>Last log line</h3>')
    for k, I in info.items():
        o.append('<div class="log">%-20s %s</div>' % (esc(I['name']), esc(I['st'].get('last_line', ''))))
    o.append('<h2 id="cross">Cross-planet</h2><figure><img src="img/cross_teq.png" alt="jet speed, day-night contrast and terminator wind versus Teq" loading="lazy">'
             '<figcaption>Mean over the last 10 dumps (5 rotations); bars span min to max. Runs are at different rotation counts '
             '(in the legend). Jet max: maximum of the zonal-mean u averaged over |lat| &lt; 15 deg, p &le; 10 bar. Day - night: '
             'area-weighted dayside minus nightside T on the isobar. Terminator wind: area-weighted line-of-sight velocity '
             '(+ = away from the star, toward a transit observer = blueshift) on the ring |cos(angle from substellar)| &lt; 0.1 at 1e-3 bar, '
             'unweighted by density or opacity. T_eq from each deck.</figcaption></figure>')
    o.append('<div class="tw"><table><thead><tr><th>run</th><th>T_eq</th><th>rot (dump)</th><th>jet max</th><th>at p</th>'
             '<th>dT d-n 1e-3</th><th>dT d-n 1e-1</th><th>term. wind 1e-3</th><th>morning / evening</th></tr></thead><tbody>')
    for k, I in info.items():
        if not I['scl']: continue
        r = I['scl'][-1]
        o.append('<tr><td>%s</td><td>%.0f K</td><td>%.1f</td><td>%.0f m/s</td><td>%.2g bar</td><td>%.0f K</td><td>%.0f K</td>'
                 '<td>%+.0f m/s</td><td>%+.0f / %+.0f m/s</td></tr>' % (
                     esc(I['name']), float(I['d']['Teq']), r['rot'], r['jet'], r['jet_p'], r['dT_1e-3'], r['dT_1e-1'],
                     r['vterm_1e-3'], r['vterm_mo_1e-3'], r['vterm_ev_1e-3']))
    o.append('</tbody></table></div>')
    for k, I in info.items():
        d = I['d']; H = I['H']; P = I['P']
        tbc = d.get('top_bc', 'hse (default)')
        wb = 'on' if d.get('roche_wb_tan', 'false').lower() == 'true' else 'off'
        met = re.sub(r'.*CK(\d+)X.*', r'\1x solar', d.get('ck_data_dir', '?'))
        mhd = 'MHD' if os.path.exists(I['run'] + '/dhj.mhd.hst') else 'hydro'
        sel = H['rot'] > H['rot'][-1] - 1
        mf = np.mean(H['Mfloor'][sel]) if H.get('Mfloor') is not None else np.nan
        o.append('<h2 id="%s">%s</h2>' % (k, esc(I['name'])))
        o.append('<ul class="facts">'
                 '<li>%s, %s; Roche tide %s (a = %.4g cm), well-balanced tangential tide (roche_wb_tan) %s; top boundary <code>%s</code>; top sponge %s; relax_drag %s; cs_seam_flux %s.</li>'
                 '<li>Rotation period %.4f d (omega = %s s<sup>-1</sup>, tidally locked); T_eq %.1f K; grid nx1 = %s, x1min = %s cm.</li>'
                 '<li>Sim time %.4g s = %.2f rotations (hst); dumps every 0.5 rotation, %d reduced.</li>'
                 '<li>Mfloor, mean over the last rotation: %.3g g/s (= %.2g of the mass per rotation). Mass change since t = 0: %+.2e; '
                 'energy change %+.2e of |E0|.</li></ul>' % (
                     mhd, met, 'on' if d.get('roche', 'false').lower() == 'true' else 'off', float(d.get('roche_a', 'nan')), wb, esc(tbc),
                     d.get('sponge_top', '-'), d.get('relax_drag', '-'), d.get('cs_seam_flux', '-'),
                     P / 86400, d['omega'], float(d['Teq']), d.get('nx1', '?'), d.get('x1min', '?'),
                     H['rot'][-1] * P, H['rot'][-1], len(I['scl']), mf, mf * P / H['M'][0],
                     (H['M'][-1] - H['M'][0]) / H['M'][0], (H['E'][-1] - H['E'][0]) / abs(H['E'][0])))
        if not os.path.exists('%s/%s_series.png' % (IMG, k)):
            o.append('<p class="sub">Figures for this run appear after the next in-allocation catch-up.</p>'); continue
        if os.path.exists('%s/%s_snap.png' % (IMG, k)):
            o.append('<figure><img src="img/%s_snap.png" alt="%s latest dump: zonal wind, terminator wind, T and wind maps" loading="lazy">'
                     '<figcaption>Latest dump. Top left: zonal-mean u(lat, p); black line u = 0. Top right: line-of-sight wind on the '
                     'terminator ring vs position along it. Bottom: T (colour) and horizontal wind (arrows, every 12 deg) on the isobar; '
                     'dashed lines are the terminators; substellar point at 0.</figcaption></figure>' % (k, esc(I['name'])))
        o.append('<figure><img src="img/%s_series.png" alt="%s history: mass, energy, L_z, Mfloor, dt, jet" loading="lazy">'
                 '<figcaption>From hst (every 0.01 rotation) except the last panel (dumps). Grey: the stopped r1b test run (WASP-189b: r1bo, outflow1 top), each normalised by its own t = 0 values. T_tide is the tide torque about z times one rotation period, in units of I_w0 Omega, so it compares directly with the L_z panel. Mfloor and Mdot_top are rates.</figcaption></figure>' % (k, esc(I['name'])))
    # 11-02 USER: conservation budgets (budget.py; figures and summary.json from the last full run)
    BSUM = CACHE + '/budget/summary.json'
    if os.path.exists(BSUM) and os.path.exists(IMG + '/budget_summary.png'):
        try:
            import budget
            o += budget.html(json.load(open(BSUM)), esc)
        except Exception as e:
            print('budget html failed', e)
    # fehook-1102 (USER): the fe_hires section is fe_section.py's block (written by ../fe_step.sh to logs/fe_section.html, images
    # img/fe_vs_rot.png + img/fehires_vs_seidel.png); until it exists, the previous static section (pre_fehook) is kept
    FES = PAGE + '/logs/fe_section.html'
    if os.path.exists(FES) and os.path.getsize(FES) > 0 and os.path.exists(IMG + '/fe_vs_rot.png'):
        o.append(open(FES).read().strip())
    else:
        FE = '/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires/fehires_vs_seidel'
        if os.path.exists(FE + '.png'):
            import shutil
            shutil.copyfile(FE + '.png', '%s/fehires_vs_seidel.png' % IMG)
            tk = open(FE + '.txt').read().strip() if os.path.exists(FE + '.txt') else ''
            o.append('<h2 id="fehires">Fe I terminator wind vs Seidel et al. 2026</h2>')
            o.append('<p><strong>%s</strong></p>' % esc(tk))
            o.append('<figure><img src="img/fehires_vs_seidel.png" alt="Fe I net blueshift vs T_eq: model vs Seidel 2026" loading="lazy">'
                     '<figcaption>Net Fe I blueshift (km/s, positive = blueshift) vs T_eq. Points with error bars: Seidel et al. 2026 (ESPRESSO/MAROON-X). '
                     'Model: viper\'s fe_hires pipeline (high-resolution Fe I transmission through the 3-D dumps, cross-correlation, Kp-marginalised offset), '
                     'filled = Fe removed above 1e-6 bar, open = no top cut; bars = 16-84 % range over Kp. Runs marked interim are not yet at 300 rotations.'
                     '</figcaption></figure>')
    # 11-02 USER: other observations (phase curves, emission spectra, time-resolved Fe I); built in obscheck_1102 (em3d GPU RT, tr_fe)
    OB = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/page'
    if os.path.exists(OB + '/obs_takeaways.json'):
        import shutil
        TKS = json.load(open(OB + '/obs_takeaways.json'))
        o.append('<h2 id="obs">Other observations: phase curves and time-resolved Fe</h2>')
        o.append('<p>%s</p>' % esc(TKS.get('_intro', '')))
        for fn in TKS.get('_order', []):
            if not os.path.exists(OB + '/' + fn): continue
            shutil.copyfile(OB + '/' + fn, IMG + '/' + fn)
            t = TKS.get(fn, {})
            o.append('<p><strong>%s</strong></p>' % esc(t.get('take', '')))
            o.append('<figure><img src="img/%s" alt="%s" loading="lazy"><figcaption>%s</figcaption></figure>' % (
                fn, esc(t.get('alt', fn)), esc(t.get('cap', ''))))
    o.append('<p class="sub">Generated by <code>pages/dhj_r1/refresh.sh</code> (tools/ana.py) inside allocation 5274196; dumps read '
             'read-only; EOS inversion with the 1x table (w121r1_1101/anchor/eos_1x.txt).</p></div>')
    tmp = PAGE + '/index.html.tmp%d' % os.getpid(); open(tmp, 'w').write('\n'.join(o) + '\n'); os.replace(tmp, PAGE + '/index.html')

if __name__ == '__main__':
    main()
