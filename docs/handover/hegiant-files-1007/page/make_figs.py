#!/usr/bin/env python3
"""Figures + numbers.json for the He giant 3-D results page (fig_rph, fig_hst, fig_profiles,
fig_slices, fig_shell). Reads the run directories in place (never writes there).

Shell means are solid-angle weighted (sin theta); R_ph = shell-mean tau_R = 2/3 from the top
(kappa_R from the run's Rosseland table at (<rho>, T_rad=(<E>/a)^1/4)), as check_link.py.

usage:
  python3 make_figs.py RUNDIR [RUNDIR ...] --out OUTDIR --table ROSS_TABLE --profile MESA_PROFILE
                       [--label NAME] [--vis ATHENAK/vis/python]
  - RUNDIR: one or more run directories IN TIME ORDER, each with bin/*.hydro_w.*.bin and
    bin/*.m1.*.bin dumps, hegiant.hydro.hst, hegiant.user.hst, optional run*.log and jobs.txt.
    Several dirs are JOINED: run i is used only for t < start time of run i+1 (its first hst
    row or first dump, whichever is earlier). This is how a run that changed the radial grid
    is followed (e.g. N445 to 10.42 d, then the remapped N897 continuation):
        make_figs.py .../n445_1g .../n897 --out page ...
    Each dump keeps its own radial grid; profile curves are drawn on their own r; history
    curves are concatenated; a dotted line marks each join.
  - all dumps of every run are used (shell means cached in OUTDIR/shellmeans_v2.npz, keyed by
    absolute path; delete it to recompute); the slice/shell figures use the latest dump.
  - Defaults may also come from the environment: HEGIANT_PAGE_OUT, HEGIANT_ROSS_TABLE,
    HEGIANT_MESA_PROFILE, ATHENAK_VIS (vis/python of an AthenaK checkout; default: the repo
    this file lives in, ../../../../vis/python).
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                       # noqa: E402
from scipy.interpolate import RectBivariateSpline     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('runs', nargs='+', help='run directories in time order')
ap.add_argument('--out', default=os.environ.get('HEGIANT_PAGE_OUT'))
ap.add_argument('--table', default=os.environ.get('HEGIANT_ROSS_TABLE',
                os.path.join(HERE, '..', 'rosseland_tops_hegiant_blend.txt')))
ap.add_argument('--profile', default=os.environ.get('HEGIANT_MESA_PROFILE',
                os.path.join(HERE, 'M2pt754_Porb100_profile29_cols.data')))
ap.add_argument('--label', default='He giant')
ap.add_argument('--vis', default=os.environ.get('ATHENAK_VIS',
                os.path.join(HERE, '..', '..', '..', '..', 'vis', 'python')))
args = ap.parse_args()
if not args.out:
    ap.error('--out (or HEGIANT_PAGE_OUT) is required')
os.makedirs(args.out, exist_ok=True)
sys.path.insert(0, os.path.abspath(args.vis))
import bin_convert as bc                              # noqa: E402

RS, AR, DAY, LSUN = 6.957e10, 7.5657332503e-15, 86400., 3.828e33
L_SURF = 1.62e4*LSUN
OUT, PROF, LAB = args.out, args.profile, args.label
RUNS = [os.path.abspath(d) for d in args.runs]
plt.rcParams.update({'font.size': 11, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
DPI = 120


def read_ross(fn, ld_ext=None):
    """he_star_m1 table format (copied from he_giant_1006/ic/make_ic_mlt_star_ge.py);
    ld_ext: HsReadOpacityTable's extension to low rho (log10 kappa continued in log rho with
    the first interval's slope clamped to [0, 1])"""
    with open(fn) as fh:
        for ln in fh:
            if ln.startswith('# ') and len(ln.split()) == 7 and ln.split()[1].isdigit():
                nT, nD, lT0, dlT, lD0, dlD = [float(v) for v in ln.split()[1:]]
    nT, nD = int(nT), int(nD)
    v = np.loadtxt(fn, comments='#').reshape(nT, nD)
    if ld_ext is not None and ld_ext < lD0 - 1e-9*dlD:
        nadd = int(np.ceil((lD0 - ld_ext)/dlD - 1e-9))
        sl = np.clip((v[:, 1] - v[:, 0])/dlD, 0.0, 1.0)
        ext = v[:, :1] - sl[:, None]*(nadd - np.arange(nadd))[None, :]*dlD
        v = np.concatenate([ext, v], axis=1)
        nD += nadd
        lD0 -= nadd*dlD
    return lT0 + dlT*np.arange(nT), lD0 + dlD*np.arange(nD), v


ty, tx, kv = read_ross(args.table, -21.0)
sk = RectBivariateSpline(ty, tx, kv, kx=1, ky=1)


def assemble(d, names):
    """full (nphi, ntheta, nr) arrays + cell centres/faces from the per-block data
    (one MeshBlock in x1, as in every He giant run)"""
    lg = np.array(d['mb_logical'])[:, :3]
    nx3, nx2, nx1 = d['mb_data'][names[0]][0].shape
    n2, n3 = lg[:, 1].max() + 1, lg[:, 2].max() + 1
    assert lg[:, 0].max() == 0, 'more than one MeshBlock in x1 is not supported'
    out = {k: np.empty((n3*nx3, n2*nx2, nx1)) for k in names}
    th, ph, thf, phf = np.empty(n2*nx2), np.empty(n3*nx3), np.empty(n2*nx2 + 1), np.empty(n3*nx3 + 1)
    for b in range(d['n_mbs']):
        j, k = lg[b, 1], lg[b, 2]
        for nm in names:
            out[nm][k*nx3:(k+1)*nx3, j*nx2:(j+1)*nx2, :] = d['mb_data'][nm][b]
        th[j*nx2:(j+1)*nx2] = d['mb_x2v'][b]
        thf[j*nx2:(j+1)*nx2 + 1] = d['mb_x2f'][b]
        ph[k*nx3:(k+1)*nx3] = d['mb_x3v'][b]
        phf[k*nx3:(k+1)*nx3 + 1] = d['mb_x3f'][b]
    return out, np.array(d['mb_x1v'][0]), np.array(d['mb_x1f'][0]), th, thf, ph, phf


def smean(f, w):
    return np.tensordot(w, f, axes=([0, 1], [0, 1]))


def rph_of(rho, E, r, rf):
    T = (np.maximum(E, 1e-30)/AR)**0.25
    k = 10**sk.ev(np.log10(T), np.log10(np.maximum(rho, 1e-21)))
    tau = np.cumsum((rho*k*np.diff(rf))[::-1])[::-1]
    i = np.where(tau >= 2.0/3.0)[0]
    return (r[i[-1]] if i.size else np.nan), T, tau


def dump_time(fn):
    """simulation time of a .bin dump from its header (no data read)"""
    with open(fn, 'rb') as f:
        for _ in range(8):
            ln = f.readline().decode('ascii', 'replace')
            m = re.search(r'time=\s*(\S+)', ln)
            if m:
                return float(m.group(1))
    return bc.read_binary(fn)['time']


def hst_load(d, name):
    fn = os.path.join(d, 'hegiant.%s.hst' % name)
    if not os.path.exists(fn):
        return None, None
    cols = None
    for ln in open(fn):
        if '[1]=' in ln:
            cols = {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', ln)}
            break
    return np.atleast_2d(np.loadtxt(fn)), cols


# ---------- segments: run i covers [t_start_i, t_start_{i+1})
seg = []
for d in RUNS:
    fh = sorted(glob.glob(os.path.join(d, 'bin', '*.hydro_w.*.bin')))
    fh = [f for f in fh if os.path.exists(f.replace('hydro_w', 'm1'))]
    td = [dump_time(f) for f in fh]
    hy, cy = hst_load(d, 'hydro')
    us, cu = hst_load(d, 'user')
    t0 = min([x for x in ([td[0]] if td else []) + ([hy[0, 0]] if hy is not None else [])]
             or [np.inf])
    seg.append(dict(dir=d, files=fh, tdump=td, hydro=hy, ch=cy, user=us, cu=cu, t0=t0))
for i, s in enumerate(seg):
    s['t1'] = seg[i + 1]['t0'] if i + 1 < len(seg) else np.inf
    # dump headers print time with 6 significant digits (rounding error <= half a unit in the
    # 6th digit): a dump of the earlier run counts as before the join if it is below t1 by
    # more than that; the later run owns a dump AT the join time if it wrote one
    tol = 0.5*10**(np.floor(np.log10(s['t1'])) - 5) if 0 < s['t1'] < np.inf else 0.0
    s['files'] = [f for f, t in zip(s['files'], s['tdump']) if t < s['t1'] - tol]
    print('segment %d: %s  t %.4f .. %s d, %d dumps' % (
        i, s['dir'], s['t0']/DAY, '%.4f' % (s['t1']/DAY) if np.isfinite(s['t1']) else 'end',
        len(s['files'])))
joins = [s['t0']/DAY for s in seg[1:]]
fh_all = [f for s in seg for f in s['files']]
assert fh_all, 'no dumps found'

# ---------- shell means of every dump (cached, per-dump radial grid)
cache = os.path.join(OUT, 'shellmeans_v2.npz')
C = dict(np.load(cache, allow_pickle=True)) if os.path.exists(cache) else {}
done = list(C.get('files', []))
recs = list(C.get('recs', []))
for fh in fh_all:
    if fh in done:
        continue
    a = bc.read_binary(fh)
    b = bc.read_binary(fh.replace('hydro_w', 'm1'))
    A, r, rf, th, thf, ph, phf = assemble(a, ['dens', 'velx'])
    B = assemble(b, ['m1_e'])[0]
    w = np.outer(np.diff(phf), np.diff(-np.cos(thf)))
    w /= w.sum()
    rho, vr, E = smean(A['dens'], w), smean(A['velx'], w), smean(B['m1_e'], w)
    vrms = np.sqrt(np.maximum(smean(A['velx']**2, w) - vr**2, 0))
    rph, T, tau = rph_of(rho, E, r, rf)
    j = np.where(rho >= 1e-12)[0]
    recs.append(dict(file=fh, t=a['time']/DAY, cycle=int(a['cycle']), nr=len(r), r=r,
                     rph=rph/RS, r12=(r[j[-1]]/RS if j.size else np.nan), rho=rho, vr=vr,
                     vrms=vrms, T=T))
    done.append(fh)
    del A, B
    np.savez(cache, files=np.array(done), recs=np.array(recs, dtype=object))
keep = set(fh_all)
recs = sorted([q for q in recs if q['file'] in keep], key=lambda q: q['t'])
t = np.array([q['t'] for q in recs])
rphs = np.array([q['rph'] for q in recs])
r12 = np.array([q['r12'] for q in recs])


def mark_joins(ax, label=True):
    for n, tj in enumerate(joins):
        ax.axvline(tj, color='k', ls=':', lw=1,
                   label=('run join (grid change)' if label and n == 0 else None))


# ---------- MESA profile
h = open(PROF).readlines()[5].split()
P = np.loadtxt(PROF, skiprows=6)
mR, mT, mD, mV = (P[:, h.index(c)] for c in ('logR', 'logT', 'logRho', 'velocity'))
mR = 10**mR

# ---------- fig 1: R_ph(t)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.axvspan(5, 10, color='0.88', label='scaffold ramp 5-10 d')
ax.plot(t, rphs, 'o-', ms=3, label=r'$R_{\rm ph}$ (shell-mean $\tau_R=2/3$)')
ax.plot(t, r12, 's-', ms=3, label=r'$r(\langle\rho\rangle=10^{-12})$')
tt = np.linspace(0, max(t.max(), 1), 50)
ax.plot(tt, 61.8 + 6.6e5*tt*DAY/RS, 'k--', label='MESA: 61.8 Rsun + 6.6 km/s t')
mark_joins(ax)
ax.set_xlabel('t [d]')
ax.set_ylabel(r'r [R$_\odot$]')
ax.set_title('%s: photosphere and outer density' % LAB)
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_rph.png'), dpi=DPI)
plt.close(fig)


# ---------- fig 2: history (concatenated over the segments)
def hst_cat(kind, name):
    tv, yv = [], []
    for s in seg:
        H, c = s[kind], s['cu' if kind == 'user' else 'ch']
        if H is None:
            continue
        m = H[:, 0] < s['t1']
        tv.append(H[m, 0])
        yv.append(sum(H[m, c[n]] for n in name.split('+')))
    return np.concatenate(tv)/DAY, np.concatenate(yv)


a0 = bc.read_binary(fh_all[0])
OMEGA = (np.cos(a0['x2min']) - np.cos(a0['x2max']))*(a0['x3max'] - a0['x3min'])
SC = 4*np.pi/OMEGA   # user.hst sums are over the wedge: scale to 4 pi
fig, axs = plt.subplots(3, 1, figsize=(7.5, 9), sharex=True)
ax = axs[0]
x, y = hst_cat('user', 'L_top')
ax.plot(x, y*SC/L_SURF, lw=0.8, label='L_top (r = 200 Rsun, lab) x 4pi/Omega')
x, y = hst_cat('user', 'L_bot')
ax.plot(x, y*SC/L_SURF, lw=0.8, label='L_bot (radiative, r_in) x 4pi/Omega')
ax.axhline(1, color='k', ls='--', lw=0.8)
ax.set_yscale('log')
ax.set_ylabel(r'L / L$_{\rm surf}$(MESA 1.62e4 L$_\odot$)')
ax.legend(fontsize=9)
ax = axs[1]
x, y = hst_cat('hydro', '1-KE+2-KE+3-KE')
ax.plot(x, y*SC, lw=0.8, label='total KE (hydro.hst 1-KE+2-KE+3-KE) x 4pi/Omega')
x, y = hst_cat('user', 'KEr_int')
ax.plot(x, y*SC, lw=0.8, label='radial KE (user.hst KEr_int) x 4pi/Omega')
ax.set_yscale('log')
ax.set_ylim(1e40, None)
ax.set_ylabel('KE [erg]')
ax.legend(fontsize=9)
ax = axs[2]
x, y = hst_cat('hydro', 'mass')
ax.plot(x, (y - y[0])/y[0], lw=0.8, label='(M - M0)/M0 (hydro.hst mass)')
ax.set_ylabel(r'$\Delta M / M_0$')
ax.legend(fontsize=9, title='user.hst Mdot_top = 0 at every row (top wall closed)',
          title_fontsize=8)
ax.set_xlabel('t [d]')
for ax in axs:
    ax.axvspan(5, 10, color='0.9')
    mark_joins(ax, label=False)
    ax.grid(alpha=0.3)
axs[0].set_title('%s: history (grey = scaffold ramp, dotted = run join)' % LAB)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_hst.png'), dpi=DPI)
plt.close(fig)

# ---------- fig 3: profiles (each on its own radial grid)
sel = sorted(set(int(np.argmin(abs(t - x))) for x in [0, 2.5, 5] + joins + [t.max()]))
fig, axs = plt.subplots(2, 2, figsize=(11, 8))
for n, i in enumerate(sel):
    q, c = recs[i], 'C%d' % n
    rR = q['r']/RS
    lab = 't = %.2f d (N%d)' % (q['t'], q['nr'])
    axs[0, 0].plot(rR, np.log10(np.maximum(q['rho'], 1e-30)), c, label=lab)
    axs[0, 1].plot(rR, np.log10(q['T']), c, label=lab)
    axs[1, 0].plot(rR, q['vr']/1e5, c, label=lab)
    axs[1, 1].plot(rR, q['vrms']/1e5, c, label=lab)
axs[0, 0].plot(mR, mD, 'k--', lw=1, label='MESA')
axs[0, 1].plot(mR, mT, 'k--', lw=1, label='MESA')
axs[1, 0].plot(mR, mV/1e5, 'k--', lw=1, label='MESA v')
axs[0, 0].set_ylabel(r'log $\langle\rho\rangle$ [g cm$^{-3}$]')
axs[0, 1].set_ylabel(r'log $T_{\rm rad}$ [K]')
axs[1, 0].set_ylabel(r'$\langle v_r\rangle$ [km/s]')
axs[1, 1].set_ylabel(r'rms $(v_r-\langle v_r\rangle)$ [km/s]')
axs[0, 0].set_ylim(-17, -2)
axs[0, 1].set_ylim(3.4, 6.7)
for ax in axs.flat:
    ax.set_xlim(3, 100)
    ax.set_xlabel(r'r [R$_\odot$]')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig.suptitle('%s: shell means (solid-angle weighted)' % LAB)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_profiles.png'), dpi=DPI)
plt.close(fig)

# ---------- latest dump, full arrays
fl = recs[-1]['file']
a, b = bc.read_binary(fl), bc.read_binary(fl.replace('hydro_w', 'm1'))
A, r, rf, th, thf, ph, phf = assemble(a, ['dens', 'velx'])
B = assemble(b, ['m1_e'])[0]
tl = a['time']/DAY
k0 = len(ph)//2
rho_m = 0.5*(A['dens'][k0 - 1] + A['dens'][k0])
vr_m = 0.5*(A['velx'][k0 - 1] + A['velx'][k0])
T_m = (np.maximum(0.5*(B['m1_e'][k0 - 1] + B['m1_e'][k0]), 1e-30)/AR)**0.25
RF, TF = np.meshgrid(rf/RS, thf)
X, Z = RF*np.sin(TF), RF*np.cos(TF)
panels = [(np.log10(rho_m), dict(cmap='viridis', vmin=-14, vmax=-3), r'log $\rho$'),
          (vr_m/1e5, dict(cmap='RdBu_r', vmin=-30, vmax=30), r'$v_r$ [km/s] ($\pm$30 clip)'),
          (np.log10(T_m), dict(cmap='inferno', vmin=3.6, vmax=6.6), r'log $T_{\rm rad}$')]
rmax = 80
fig, axs = plt.subplots(2, 3, figsize=(13, 9), gridspec_kw=dict(height_ratios=[1.6, 1]))
for n, (f, kw, lab) in enumerate(panels):
    for row, (r0, r1) in enumerate([(3, rmax), (50, 75)]):
        ax = axs[row, n]
        ir = (rf[:-1]/RS < r1)
        nr = ir.sum()
        im = ax.pcolormesh(X[:, :nr + 1], Z[:, :nr + 1], f[:, :nr], shading='flat', **kw)
        ax.set_aspect('equal')
        if row == 1:
            ax.set_xlim(r0*np.sin(np.pi/3) - 1, r1)
            ax.set_ylim(-r1*0.5, r1*0.5)
        ax.set_xlabel(r'x [R$_\odot$]')
        ax.set_ylabel(r'z [R$_\odot$]')
        fig.colorbar(im, ax=ax, shrink=0.8, label=lab)
    axs[0, n].set_title(lab.split(' (')[0] + (', r < %d Rsun' % rmax))
    axs[1, n].set_title('zoom 50-75 Rsun')
fig.suptitle(r'%s: meridional slice at $\phi$ = %.0f deg, t = %.2f d' %
             (LAB, np.degrees(phf[k0]), tl))
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_slices.png'), dpi=DPI)
plt.close(fig)

# ---------- fig 5: theta-phi shells
w = np.outer(np.diff(phf), np.diff(-np.cos(thf)))
w /= w.sum()
rho_s = smean(A['dens'], w)
E_s = smean(B['m1_e'], w)
rph_l = rph_of(rho_s, E_s, r, rf)[0]/RS
shells = [20.0, rph_l - 2.0]
PH, THm = np.meshgrid(np.degrees(phf), 90 - np.degrees(thf), indexing='ij')
fig, axs = plt.subplots(2, 2, figsize=(11, 9))
for row, rs in enumerate(shells):
    i = int(np.argmin(abs(r/RS - rs)))
    v = A['velx'][:, :, i]/1e5
    vm = np.percentile(abs(v), 99)
    im = axs[row, 0].pcolormesh(PH, THm, v, cmap='RdBu_r', vmin=-vm, vmax=vm)
    fig.colorbar(im, ax=axs[row, 0], label=r'$v_r$ [km/s]')
    dr = A['dens'][:, :, i]/np.sum(w*A['dens'][:, :, i])
    im = axs[row, 1].pcolormesh(PH, THm, dr, cmap='PuOr_r', vmin=0.5, vmax=1.5)
    fig.colorbar(im, ax=axs[row, 1], label=r'$\rho/\langle\rho\rangle$')
    for ax in axs[row]:
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_ylabel('latitude [deg]')
        ax.set_title('r = %.1f Rsun' % (r[i]/RS))
fig.suptitle('%s: shells at t = %.2f d (R_ph = %.1f Rsun)' % (LAB, tl, rph_l))
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig_shell.png'), dpi=DPI)
plt.close(fig)


# ---------- numbers.json
def job_status(d):
    fn = os.path.join(d, 'jobs.txt')
    if not os.path.exists(fn):
        return {}
    ids = open(fn).read().split()
    try:
        q = subprocess.run(['squeue', '-h', '-j', ','.join(ids), '-o', '%i %T %M'],
                           capture_output=True, text=True).stdout.split('\n')
    except OSError:
        return {i: 'squeue unavailable' for i in ids}
    st = {z.split()[0]: ' '.join(z.split()[1:]) for z in q if z.strip()}
    return {i: st.get(i, 'not in queue (finished)') for i in ids}


links = []
for s in seg:
    for lg in sorted(glob.glob(os.path.join(s['dir'], 'run*.log')), key=os.path.getmtime):
        txt = open(lg, errors='replace').read()
        c = re.findall(r'^elapsed=(\S+) cycle=(\d+) time=(\S+) dt=(\S+)', txt, re.M)
        e = np.array([[float(x) for x in z] for z in c]) if c else np.zeros((0, 4))
        n = len(e)//10
        links.append(dict(log=os.path.relpath(lg, os.path.dirname(s['dir'])),
                          FATAL=txt.count('FATAL'),
                          NONCONV=len(re.findall('NON-CONVERGED after', txt)),
                          cycle_last=int(e[-1, 1]) if len(e) else None,
                          t_last_d=float(e[-1, 2]/DAY) if len(e) else None,
                          s_per_cycle=(float((e[-1, 0] - e[n, 0])/max(e[-1, 1] - e[n, 1], 1))
                                       if len(e) > 2 else None),
                          finished=('Terminating on' in txt or 'wall time limit' in txt.lower()
                                    or 'time limit reached' in txt.lower())))
num = dict(runs=RUNS, joins_t_d=[round(x, 4) for x in joins],
           segments=[dict(dir=s['dir'], t0_d=round(s['t0']/DAY, 4), ndumps=len(s['files']))
                     for s in seg],
           latest_dump=fl, latest_dump_t_d=round(tl, 4), latest_dump_cycle=int(a['cycle']),
           latest_dump_nr=len(r),
           latest_log_t_d=links[-1]['t_last_d'] if links else None,
           latest_log_cycle=links[-1]['cycle_last'] if links else None,
           Rph_table=[dict(t_d=round(q['t'], 3), nr=int(q['nr']), Rph_Rsun=round(float(q['rph']), 3),
                           **{'r_rho1e-12_Rsun': round(float(q['r12']), 3)}) for q in recs],
           max_abs_mean_vr_kms=[dict(t_d=round(q['t'], 3),
                                     vr=round(float(abs(q['vr']).max()/1e5), 3),
                                     at_Rsun=round(float(q['r'][np.argmax(abs(q['vr']))]/RS), 2))
                                for q in recs],
           link_status={k: v for s in seg for k, v in job_status(s['dir']).items()},
           links=links, NONCONV_total=int(sum(x['NONCONV'] for x in links)),
           note='shell means solid-angle weighted; R_ph from shell-mean tau_R=2/3 (T_rad, blend table)')
json.dump(num, open(os.path.join(OUT, 'numbers.json'), 'w'), indent=1)
print('latest t %.3f d (N%d), R_ph %.3f Rsun' % (tl, len(r), rph_l))
