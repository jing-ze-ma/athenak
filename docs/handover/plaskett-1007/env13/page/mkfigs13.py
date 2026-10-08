#!/usr/bin/env python3
"""Plaskett accretor page figures, env13 version (adapted from viper's mkfigs_plaskett.py,
origin/accretor-1006:docs/handover/plaskett-1007/page/mkfigs_plaskett.py; the original is kept
next to this file as mkfigs_plaskett_orig.py).

Changes against the original:
  * measuring radius: MR/JR/Menv/Jenv refer to problem/r_meas when that key is set (env13:
    9.615675, above the bulged photosphere), else to R_acc.  The radius actually used is
    taken from the pgen's "measuring face ... at r = X" line in run.log when present.
    j is normalised by j_K(r_meas); the ballistic j and the star-rotation j are evaluated at
    r_meas too.  For r_meas = R_acc (q12) this reproduces the old numbers exactly.
  * maps: r_meas drawn as a dashed blue circle when it differs from R_acc.
  * envelope edge / breathing (new): from every bin dump, the first radius above 8.5 Rsun
    where the theta-mean density drops below edge_frac * env_rho_ph, per phi; compared with
    the first dump and with the Roche equipotential through (R_acc, phi 90) (= photosphere).
    Cached per arm in edge_<arm>.npz.  fig_edge_cmp.png, numbers edge_*.
  * no-data / early-run handling: an arm without a usable history writes
    numbers_<arm>.json with no_data = true and is left out of the figures.
  * --extra-logs: slurm logs outside the run dir (FATAL count, min dt).
  * maps / edge are skipped when their output is newer than every bin dump (--force redoes).
  * default --vis is ./vis (bin_convert.py from accretor-1007: the rt-integration copy in
    /u/jma20/ATHENAK/athenak/vis/python cannot read the version-1.2 bin files).

usage: mkfigs13.py --run DIR [DIR ...] [--names ...] --out DIR [--vis DIR] [--cmp-only]
"""
import argparse
import glob
import json
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

G, MSUN, RSUN, KB, MH = 6.674e-8, 1.989e33, 6.957e10, 1.381e-16, 1.6726e-24
GU = G*MSUN/(RSUN*1e10)          # G Msun in (km/s)^2 Rsun (as the pgen)
COLS = ('#2a78d6', '#eb6834', '#2f9e44', '#8e44ad', '#888888')
HERE = os.path.dirname(os.path.abspath(__file__))
EDGE_RMIN = 8.5                  # edge search starts here (below every photosphere)


# ---------------------------------------------------------------- input parameters
def read_input_text(run, inp):
    if inp:
        return open(inp).read(), inp
    rst = sorted(glob.glob(os.path.join(run, 'rst', '*.rst'))
                 + glob.glob(os.path.join(run, '*.rst')))
    for fn in rst:
        with open(fn, 'rb') as f:
            head = f.read(200000).decode('latin-1')
        k = head.find('<par_end>')
        if k > 0:
            return head[:k], fn
    ath = sorted(glob.glob(os.path.join(run, '*.athinput')))
    if ath:
        return open(ath[0]).read(), ath[0]
    return None, None


def parse_input(txt):
    par, blk = {}, None
    for line in txt.splitlines():
        line = line.split('#')[0].strip()
        m = re.match(r'<(\S+)>', line)
        if m:
            blk = m.group(1)
            continue
        if '=' in line and blk:
            k, v = line.split('=', 1)
            par[(blk, k.strip())] = v.strip()
    return par


def rmeas_from_log(run):
    fn = os.path.join(run, 'run.log')
    if os.path.isfile(fn):
        with open(fn, errors='replace') as f:
            for line in f:
                m = re.search(r'measuring face \(MR.*?\) at r = ([0-9.eE+-]+)', line)
                if m:
                    return float(m.group(1))
    return None


class Binary:
    """Roche geometry + the ballistic L1 stream, same launch as the pgen
    (from L1 toward the accretor at the donor sound speed)."""

    def __init__(self, par, rmeas_log=None):
        def g(k, d):
            return float(par.get(('problem', k), d))
        self.ma, self.md = g('m_acc', 6.24), g('m_don', 1.69)
        self.a = g('a_sep', 30.3)
        self.racc = g('r_acc', 4.06)
        self.rmeas = rmeas_log if rmeas_log else g('r_meas', self.racc)
        self.tdon, self.mudon = g('t_don', 6250.0), g('mu_don', 1.27)
        self.rin = float(par[('mesh', 'x1min')])
        self.rout = float(par[('mesh', 'x1max')])
        self.nx1 = int(par[('mesh', 'nx1')])
        self.stream = par.get(('problem', 'stream'), 'true').lower() == 'true'
        self.spin = g('spin', 1.0)
        self.rho_ph = g('env_rho_ph', np.nan)
        self.top_mode = par.get(('problem', 'env_top_mode'), 'default (sphere)')
        self.stream_width = par.get(('problem', 'stream_width'), 'default (c_s/Omega)')
        self.om = np.sqrt(GU*(self.ma + self.md)/self.a**3)     # code 1/time
        self.porb = 2*np.pi/self.om
        self.aom = self.a*self.om                               # km/s
        self.jk = np.sqrt(GU*self.ma*self.rmeas)                # j_Kep(r_meas)
        self.jk_acc = np.sqrt(GU*self.ma*self.racc)             # j_Kep(R_acc)
        self.jcor = self.spin*self.om*self.rmeas**2             # rigid rotation at r_meas
        mu = self.md/(self.ma + self.md)
        self.mu = mu
        xa, xd = -mu, 1 - mu
        eps = np.sqrt(KB*self.tdon/(self.mudon*MH))/1e5/self.aom

        def gx(x):
            return (-(1 - mu)*(x - xa)/abs(x - xa)**3 - mu*(x - xd)/abs(x - xd)**3 + x)
        xl1 = brentq(gx, xa + 1e-3, xd - 1e-3)

        def rhs(t, s):
            x, y, vx, vy = s
            r1, r2 = np.hypot(x - xa, y), np.hypot(x - xd, y)
            return [vx, vy, -(1 - mu)*(x - xa)/r1**3 - mu*(x - xd)/r2**3 + x + 2*vy,
                    -(1 - mu)*y/r1**3 - mu*y/r2**3 + y - 2*vx]

        sol = solve_ivp(rhs, [0, 3], [xl1 - 1e-4, 0, -eps, 0], rtol=1e-11, atol=1e-13,
                        max_step=1e-3)
        x, y, vx, vy = sol.y
        self.bx, self.by = (x - xa)*self.a, y*self.a       # accretor at 0, donor at +x
        rb = np.hypot(self.bx, self.by)
        # inertial specific AM at the first inward crossing of r_meas (code units)
        jd = (x - xa)*vy - y*vx + ((x - xa)**2 + y**2)
        i = np.argmax(rb < self.rmeas)
        self.jbal = np.nan
        if rb[i] < self.rmeas and i > 0:
            f = (rb[i - 1] - self.rmeas)/(rb[i - 1] - rb[i])
            self.jbal = (jd[i - 1] + f*(jd[i] - jd[i - 1]))*self.a*self.aom
            ph = np.degrees(np.arctan2(self.by[i - 1], self.bx[i - 1]))
            print(f'  ballistic at r_meas {self.rmeas:.4f}: phi {ph:.3f} deg, '
                  f'j/j_K {self.jbal/self.jk:.4f}')

    def roche(self, r, phi):
        """Roche potential (synchronous frame, equator), accretor at 0, donor at +x."""
        x, y = r*np.cos(phi), r*np.sin(phi)
        xcm = self.a*self.md/(self.ma + self.md)
        return (-GU*self.ma/r - GU*self.md/np.hypot(x - self.a, y)
                - 0.5*self.om**2*((x - xcm)**2 + y**2))

    def photosphere(self, phis):
        """Equipotential through (R_acc, phi 90): the env13 photosphere (spin 1 only)."""
        p0 = self.roche(self.racc, np.pi/2)
        return np.array([brentq(lambda r: self.roche(r, p) - p0, 0.8*self.racc,
                                1.3*self.racc) for p in phis])


# ---------------------------------------------------------------- history
def read_hst(fn):
    lab = None
    with open(fn) as f:
        for line in f:
            if line.startswith('#') and '[1]=' in line:
                lab = re.findall(r'\[\d+\]=(\S+)', line)
    if lab is None:
        return None
    try:
        h = np.loadtxt(fn, ndmin=2)
    except ValueError:            # a half-written last line while the run is going
        rows = []
        with open(fn) as f:
            for line in f:
                if not line.startswith('#'):
                    v = line.split()
                    if len(v) == len(lab):
                        try:
                            rows.append([float(x) for x in v])
                        except ValueError:
                            pass
        h = np.array(rows, ndmin=2)
    if h.size == 0 or h.shape[0] < 2:
        return None
    d = {k: h[:, i] for i, k in enumerate(lab)}
    # restarts append and restart the cumulative integrals at zero: keep the newest
    # segment where segments overlap and add the previous value at each join
    t = d['time']
    brk = [0] + [i for i in range(1, len(t)) if t[i] <= t[i - 1]] + [len(t)]
    if len(brk) > 2:
        keep, segs = [], [(brk[k], brk[k + 1]) for k in range(len(brk) - 1)]
        for k, (s0, s1) in enumerate(segs):
            tn = t[segs[k + 1][0]] if k + 1 < len(segs) else np.inf
            keep.append(np.arange(s0, s1)[t[s0:s1] < tn])
        cum = [k for k in lab if not k.startswith('d') and k[0] in 'MJ'
               and k not in ('Mdom', 'Jdom', 'Menv', 'Jenv')]
        for c in cum:
            v = d[c].copy()
            for k in range(1, len(segs)):
                s0 = segs[k][0]
                p0, p1 = segs[k - 1]     # the whole previous segment (already offset)
                off = np.interp(t[s0], t[p0:p1], v[p0:p1])
                v[segs[k][0]:segs[k][1]] += off
            d[c] = v
        idx = np.concatenate(keep)
        d = {k: v[idx] for k, v in d.items()}
    return d


def fatal_and_dt(run, extra=()):
    nf, dts = 0, []
    logs = [f for f in glob.glob(os.path.join(run, '*'))
            if re.search(r'(\.log|\.out|\.err)$', f) and os.path.isfile(f)]
    logs += [f for f in extra if os.path.isfile(f)]
    for fn in logs:
        with open(fn, errors='replace') as f:
            for line in f:
                if re.search(r'FATAL ERROR|### FATAL', line):
                    nf += 1
                m = re.search(r'\bdt=([0-9.eE+-]+)', line)
                if m:
                    try:
                        dts.append(float(m.group(1)))
                    except ValueError:
                        pass
    return nf, (min(dts) if dts else np.nan), [os.path.basename(x) for x in logs]


# ---------------------------------------------------------------- maps
def load(fn, vis):
    if vis not in sys.path:
        sys.path.insert(0, vis)
    import bin_convert as bc
    fd = bc.read_binary(fn)
    nb = fd['n_mbs']
    if fd['Nx1'] != fd['nx1_out_mb']:
        sys.exit('maps assume one MeshBlock in x1 (meshblock/nx1 = mesh/nx1)')
    if 'mb_x1f' in fd:
        rf, pfs = fd['mb_x1f'][0], [fd['mb_x3f'][b] for b in range(nb)]
    else:   # older bin_convert: faces from the centres
        r = fd['mb_x1v'][0]
        rf = np.concatenate([[r[0] - 0.5*(r[1] - r[0])], 0.5*(r[1:] + r[:-1]),
                             [r[-1] + 0.5*(r[-1] - r[-2])]])
        pfs = []
        for b in range(nb):
            p = fd['mb_x3v'][b]
            dp = p[1] - p[0]
            pfs.append(np.concatenate([p - 0.5*dp, [p[-1] + 0.5*dp]]))
    order = np.argsort([p[0] for p in pfs])
    pf = np.concatenate([pfs[b][:-1] for b in order] + [pfs[order[-1]][-1:]])
    d = np.concatenate([fd['mb_data']['dens'][b] for b in order], axis=0).mean(axis=1)
    return fd['time'], rf, pf, d


def circle(ax, r, **kw):
    th = np.linspace(0, 2*np.pi, 400)
    ax.plot(r*np.cos(th), r*np.sin(th), **kw)


def facemap(ax, dump, B, vis, rmax, vmin, vmax, title=''):
    t, rf, pf, d = dump
    R, P = np.meshgrid(rf, pf)
    m = ax.pcolormesh(R*np.cos(P), R*np.sin(P), np.log10(np.maximum(d, 1e-30)),
                      cmap='magma', vmin=vmin, vmax=vmax, shading='flat', rasterized=True)
    th = np.linspace(0, 2*np.pi, 400)
    rin = rf[0]
    ax.fill(rin*np.cos(th), rin*np.sin(th), color='#cfd3d8', zorder=3)
    ax.plot(rin*np.cos(th), rin*np.sin(th), color='black', lw=1.0, zorder=4)
    ax.text(0, 0, 'inner boundary\n(not simulated)', ha='center', va='center',
            fontsize=8, color='0.25', zorder=5, clip_on=True)
    circle(ax, B.racc, color='#2a78d6', lw=1.0, zorder=4)
    if abs(B.rmeas - B.racc) > 1e-6:
        circle(ax, B.rmeas, color='#7fb3ff', lw=0.9, ls='--', zorder=4)
    rb = np.hypot(B.bx, B.by)
    inside = (rb <= rf[-1]) & (rb >= rin)
    # first passage only: stop where the orbit leaves the grid inward
    ilast = np.argmax(rb < rin) if np.any(rb < rin) else len(rb)
    inside[ilast:] = False
    ax.plot(np.where(inside, B.bx, np.nan), np.where(inside, B.by, np.nan),
            color='white', lw=0.9, ls='--', zorder=4)
    if rmax[0] == 'full':
        ax.set_xlim(-rf[-1], rf[-1])
        ax.set_ylim(-rf[-1], rf[-1])
    else:
        ax.set_xlim(*rmax[0])
        ax.set_ylim(*rmax[1])
    ax.set_aspect('equal')
    ax.set_title(f'{title}t = {t/B.porb:.3f} orbit', fontsize=10)
    ax.set_xlabel('x [R$_\\odot$] (donor to +x)')
    return m


def bin_time(fn):
    with open(fn, 'rb') as f:
        m = re.search(rb'time=([0-9.eE+-]+)', f.read(400))
    return float(m.group(1)) if m else np.nan


def bin_list(run):
    return sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))


def _src_sig(run):
    return dict(run=os.path.realpath(run),
                bins=[[os.path.basename(f), os.path.getmtime(f)] for f in bin_list(run)])


def up_to_date(out_fn, run, force):
    """True if out_fn exists, was made from this run's current bin dumps (sidecar
    <out_fn>.src.json) and is newer than this script."""
    side = out_fn + '.src.json'
    if force or not os.path.isfile(out_fn) or not os.path.isfile(side):
        return False
    try:
        if json.load(open(side)) != _src_sig(run):
            return False
    except Exception:
        return False
    return os.path.getmtime(out_fn) > os.path.getmtime(os.path.abspath(__file__))


def mark_made(out_fn, run):
    with open(out_fn + '.src.json', 'w') as f:
        json.dump(_src_sig(run), f)


def fig_maps(arm, run, B, vis, out, vmin, vmax, force=False):
    fns = bin_list(run)
    fn_out = os.path.join(out, f'fig_{arm}_maps.png')
    if not fns:
        print(f'  {arm}: no bin dumps, maps skipped')
        return None
    if up_to_date(fn_out, run, force):
        print(f'  {arm}: maps up to date, kept {fn_out}')
        return fn_out
    # drop a dump that nearly repeats the next one (the final dump at tlim)
    tb = np.array([bin_time(f) for f in fns])
    if len(fns) > 2:
        gap = np.median(np.diff(tb))
        fns = [f for k, f in enumerate(fns)
               if k == len(fns) - 1 or tb[k + 1] - tb[k] > 0.25*gap]
    pick = np.unique(np.linspace(0, len(fns) - 1, 4).round().astype(int))
    if len(fns) > 4 and pick[0] == 0:
        pick = np.unique(np.linspace(1, len(fns) - 1, 4).round().astype(int))
    dumps = [load(fns[i], vis) for i in pick]
    n = len(dumps)
    # zoom box: the accretor surface on the stream side (impact region of the
    # ballistic orbit), +-0.45 R_acc around the point where it crosses R_acc
    rb = np.hypot(B.bx, B.by)
    i = np.argmax(rb < B.racc) if np.any(rb < B.racc) else np.argmin(rb)
    cx, cy, w = B.bx[i], B.by[i], 0.45*B.racc
    zoom = ((cx - w, cx + w), (cy - w, cy + w))
    fig, axs = plt.subplots(2, n, figsize=(3.9*n + 1.2, 7.6), constrained_layout=True,
                            squeeze=False)
    for k, dmp in enumerate(dumps):
        m = facemap(axs[0, k], dmp, B, vis, ('full',), vmin, vmax, f'{arm}\n')
        axs[0, k].plot([zoom[0][0], zoom[0][1], zoom[0][1], zoom[0][0], zoom[0][0]],
                       [zoom[1][0], zoom[1][0], zoom[1][1], zoom[1][1], zoom[1][0]],
                       color='cyan', lw=0.8, zorder=6)
        facemap(axs[1, k], dmp, B, vis, zoom, vmin, vmax, 'zoom, ')
    for r in range(2):
        axs[r, 0].set_ylabel('y [R$_\\odot$]')
    fig.colorbar(m, ax=axs, shrink=0.8,
                 label='log$_{10}$ $\\rho$ / $\\rho_{stream}$ (θ-mean)')
    fig.savefig(fn_out, dpi=100)
    plt.close(fig)
    mark_made(fn_out, run)
    print('  wrote', fn_out)
    return fn_out


# ---------------------------------------------------------------- envelope edge
def edge_radius(rf, d, thr, rmin=EDGE_RMIN):
    """Per phi row: first radius above rmin where d drops below thr (log-interpolated
    between cell centres); nan if it never does."""
    rc = 0.5*(rf[1:] + rf[:-1])
    i0 = np.searchsorted(rc, rmin)
    ld, lt = np.log(np.maximum(d[:, i0:], 1e-300)), np.log(thr)
    below = ld < lt
    has = below.any(axis=1) & ~below[:, 0]
    k = np.argmax(below, axis=1)
    out = np.full(d.shape[0], np.nan)
    j = np.where(has)[0]
    kk = k[j] + i0
    l0, l1 = ld[j, k[j] - 1], ld[j, k[j]]
    f = (l0 - lt)/(l0 - l1)
    out[j] = rc[kk - 1] + f*(rc[kk] - rc[kk - 1])
    return out


def edge_data(arm, run, B, vis, out, frac, force=False):
    fns = bin_list(run)
    cache = os.path.join(out, f'edge_{arm}.npz')
    if not fns or not np.isfinite(B.rho_ph):
        return None
    if up_to_date(cache, run, force):
        z = np.load(cache)
        if float(z['frac']) == frac and len(z['t']) == len(fns):
            print(f'  {arm}: edge cache up to date')
            return dict(z)
    ts, edges, pc = [], [], None
    thr = frac*B.rho_ph
    for fn in fns:
        t, rf, pf, d = load(fn, vis)
        pc = 0.5*(pf[1:] + pf[:-1])
        ts.append(t)
        edges.append(edge_radius(rf, d, thr))
    z = dict(t=np.array(ts), phi=pc, edge=np.array(edges), frac=frac, thr=thr,
             rph=B.photosphere(pc[::16]), phi_rph=pc[::16])
    np.savez(cache, **z)
    mark_made(cache, run)
    print('  wrote', cache)
    return z


STREAM_SECTOR = (-35.0, 50.0)    # deg: stream arrival + impact zone, left out of the edge metrics


def edge_metrics(z, B):
    """Breathing numbers from the edge data, away from the stream sector (STREAM_SECTOR):
    change against the first dump (median and max over phi; time series and its max), and
    the offset above the photosphere equipotential at phi 90 and 180 (initial, final)."""
    ph = np.degrees(np.angle(np.exp(1j*z['phi'])))
    keep = ~((ph > STREAM_SECTOR[0]) & (ph < STREAM_SECTOR[1]))
    e = z['edge'][:, keep]
    de = np.abs(e - e[0][None, :])
    med_t = np.nanmedian(de, axis=1)
    max_t = np.nanmax(de, axis=1)
    rph = np.interp(z['phi'], z['phi_rph'], z['rph'], period=2*np.pi)
    off = z['edge'] - rph[None, :]

    def at(i, p):
        s = np.abs(np.angle(np.exp(1j*(z['phi'] - np.radians(p))))) < np.radians(5)
        v = off[i, s]
        return float(np.nanmedian(v)) if np.isfinite(v).any() else None
    return dict(edge_dr_median_max=float(np.nanmax(med_t)),
                edge_dr_max_max=float(np.nanmax(max_t)),
                edge_dr_median_final=float(med_t[-1]),
                edge_minus_rph_phi90_initial=at(0, 90), edge_minus_rph_phi90_final=at(-1, 90),
                edge_minus_rph_phi180_initial=at(0, 180),
                edge_minus_rph_phi180_final=at(-1, 180),
                edge_minus_rph_phi0_initial=at(0, 0),
                edge_frac_of_rho_ph=float(z['frac'])), med_t, max_t


def fig_edge(arms, out):
    have = [(a, z, B) for a, z, B in arms if z is not None]
    if not have:
        print('  no edge data, fig_edge_cmp skipped')
        return None
    n = len(have)
    fig = plt.figure(figsize=(4.3*n + 1, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(2, n)
    cmap = plt.get_cmap('viridis')
    for k, (arm, z, B) in enumerate(have):
        ax = fig.add_subplot(gs[0, k])
        nt = len(z['t'])
        ph = np.degrees(np.angle(np.exp(1j*z['phi'])))
        o = np.argsort(ph)
        for i in range(nt):
            ax.plot(ph[o], z['edge'][i][o], lw=0.8, color=cmap(i/max(nt - 1, 1)),
                    label='first dump' if i == 0 else ('last dump' if i == nt - 1 else None))
        pr = np.degrees(np.angle(np.exp(1j*z['phi_rph'])))
        o = np.argsort(pr)
        ax.plot(pr[o], z['rph'][o], color='k', ls='--', lw=1.2,
                label='photosphere (Roche eqp. through R$_{acc}$, φ 90°)')
        ax.axhline(B.racc, color='#2a78d6', lw=0.8)
        if abs(B.rmeas - B.racc) > 1e-6:
            ax.axhline(B.rmeas, color='#7fb3ff', lw=0.9, ls='--', label='r$_{meas}$')
        ax.axvspan(*STREAM_SECTOR, color='0.85', zorder=0, lw=0)
        ax.set_xlim(-180, 180)
        ax.set_xticks(range(-180, 181, 90))
        top = np.nanpercentile(z['edge'], 99.5)
        ax.set_ylim(min(B.racc, np.nanmin(z['rph'])) - 0.05, max(top, B.rmeas) + 0.1)
        ax.set_xlabel('φ [deg] (donor at 0)')
        ax.set_title(f'{arm}: edge ρ = {z["frac"]:g} ρ$_{{ph}}$, t = 0 … '
                     f'{z["t"][-1]/B.porb:.3f} orb', fontsize=10)
        if k == 0:
            ax.set_ylabel('envelope edge radius [R$_\\odot$]')
        ax.legend(frameon=False, fontsize=8, loc='upper center')
    ax = fig.add_subplot(gs[1, :])
    for k, (arm, z, B) in enumerate(have):
        _, med_t, max_t = edge_metrics(z, B)
        c = COLS[k % len(COLS)]
        ax.plot(z['t']/B.porb, med_t, color=c, lw=2, marker='o', ms=3,
                label=f'{arm}: median over φ')
        ax.plot(z['t']/B.porb, max_t, color=c, lw=1, ls='--', label=f'{arm}: max over φ')
    ax.set_yscale('symlog', linthresh=1e-3)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel('|r$_{edge}$(t) − r$_{edge}$(first dump)| [R$_\\odot$]')
    ax.set_title(f'envelope breathing: edge displacement against the first dump, '
                 f'outside the stream sector φ {STREAM_SECTOR[0]:.0f}…{STREAM_SECTOR[1]:.0f}°')
    ax.legend(frameon=False, fontsize=9, ncol=2)
    fn = os.path.join(out, 'fig_edge_cmp.png')
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    print('  wrote', fn)
    return fn


# ---------------------------------------------------------------- budget
def mlab(B):
    return 'r$_{meas}$' if abs(B.rmeas - B.racc) > 1e-6 else 'R$_{acc}$'


def jref_lines(ax, B, col='0.4', tag='', x0=0.02):
    for v, lab, ls in ((B.jbal, 'ballistic', ':'), (B.jcor, 'star rotation', '-.')):
        if np.isfinite(v):
            ax.axhline(v/B.jk, color=col, ls=ls, lw=1.1)
            ax.text(x0, v/B.jk, f' {tag}{lab} {v/B.jk:.3f}', va='bottom', fontsize=8,
                    color=col, transform=ax.get_yaxis_transform())


def jmean_last(t, M, J, tw):
    sel = t >= t[-1] - tw
    if sel.sum() < 2:
        return np.nan
    dm = M[-1] - M[sel][0]
    return (J[-1] - J[sel][0])/dm if abs(dm) > 0 else np.nan


def fig_budget(arm, h, B, out):
    t = h['time']/B.porb
    ml = mlab(B)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    ax = axs[0]
    ax.plot(t, h['Min'], color=COLS[0], lw=2, label='stream inflow (window)')
    ax.plot(t, h['MR'], color=COLS[1], lw=2, label=f'net mass in through {ml} = {B.rmeas:.3f}')
    ax.plot(t, h['Mout'], color=COLS[2], lw=1.5, label='outflow through r$_{out}$')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel('cumulative mass [ρ$_{stream}$ R$_\\odot^3$]')
    ax.set_title(arm)
    ax.legend(frameon=False, fontsize=9)
    ax = axs[1]
    for M, J, c, lab in ((h['Min'], h['Jin'], COLS[0], 'brought in by the stream'),
                         (h['MR'], h['JR'], COLS[1], f'gas through {ml}')):
        ok = np.abs(M) > 1e-2*max(np.abs(M).max(), 1e-30)
        if ok.any() and np.abs(M).max() > 0:
            ax.plot(t[ok], J[ok]/M[ok]/B.jk, color=c, lw=2, label=lab + ' (cumulative)')
    jref_lines(ax, B)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel(f'inertial specific AM / j$_{{Kep}}$({ml})')
    ax.legend(frameon=False, fontsize=9)
    fn = os.path.join(out, f'fig_{arm}_budget.png')
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('  wrote', fn)


def fig_cmp(arms, out):
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.4), constrained_layout=True)
    for k, (arm, h, B) in enumerate(arms):
        t = h['time']/B.porb
        c = COLS[k % len(COLS)]
        ml = mlab(B)
        axs[0].plot(t, h['MR'], color=c, lw=2, label=f'{arm}: through {ml} ({B.rmeas:.2f})')
        axs[0].plot(t, h['Min'], color=c, lw=1, ls='--', label=f'{arm}: stream inflow')
        if h['Min'][-1] > 0:
            axs[1].plot(t, h['MR']/h['Min'][-1], color=c, lw=2, label=arm)
        ok = np.abs(h['MR']) > 1e-2*max(np.abs(h['MR']).max(), 1e-30)
        if ok.any() and np.abs(h['MR']).max() > 0:
            axs[2].plot(t[ok], h['JR'][ok]/h['MR'][ok]/B.jk, color=c, lw=2,
                        label=f'{arm} (/ j$_K$({B.rmeas:.2f}))')
    seen = set()
    for k, (arm, h, B) in enumerate(arms):
        key = round(B.rmeas, 4)
        if key in seen:
            continue
        seen.add(key)
        jref_lines(axs[2], B, col=COLS[k % len(COLS)] if len(arms) > 1 else '0.4',
                   tag=f'r {B.rmeas:.2f}: ', x0=0.02 + 0.5*(len(seen) - 1))
    lo, hi = axs[2].get_ylim()
    if lo < -1 or hi > 2:      # zero crossings of the cumulative MR give spikes
        axs[2].set_ylim(max(lo, -1), min(hi, 2))
        axs[2].text(0.98, 0.02, 'clipped to [-1, 2]; spikes = MR crossing 0', ha='right',
                    fontsize=8, color='0.4', transform=axs[2].transAxes)
    axs[2].axhline(1, color='0.75', lw=0.8)
    axs[2].text(0.02, 1, ' Keplerian', va='bottom', fontsize=9, color='0.5',
                transform=axs[2].get_yaxis_transform())
    for ax in axs:
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set_xlabel('time [orbits]')
        ax.legend(frameon=False, fontsize=8)
    axs[0].set_ylabel('cumulative mass [ρ$_{stream}$ R$_\\odot^3$]')
    axs[1].set_ylabel('net mass through measuring face / final stream inflow')
    axs[2].set_ylabel('j of gas through the measuring face (cumulative) / j$_{Kep}$')
    axs[0].set_title('accreted mass')
    axs[1].set_title('accreted mass, per unit stream mass')
    axs[2].set_title('specific AM through the measuring face')
    fn = os.path.join(out, 'fig_plaskett_cmp.png')
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('  wrote', fn)


def fig_envonly(runs, out, names):
    """Envelope-only tests (no stream): net mass through the measuring face."""
    fig, ax = plt.subplots(1, 1, figsize=(7.5, 4.2), constrained_layout=True)
    n = 0
    for k, (run, nm) in enumerate(zip(runs, names)):
        hf = sorted(glob.glob(os.path.join(run, '*.user.hst')))
        txt, _ = read_input_text(run, None)
        if not hf or txt is None:
            continue
        h = read_hst(hf[0])
        if h is None:
            continue
        par = parse_input(txt)
        porb = 2*np.pi/np.sqrt(GU*(float(par[('problem', 'm_acc')])
                                   + float(par[('problem', 'm_don')]))
                               / float(par[('problem', 'a_sep')])**3)
        rho_ph = float(par.get(('problem', 'env_rho_ph'), 1.0))
        ax.plot(h['time']/porb, h['MR']/rho_ph, color=COLS[k % len(COLS)], lw=2,
                label=f'{nm} (r$_{{meas}}$ {rmeas_from_log(run) or float(par[("problem", "r_acc")]):.3f})')
        n += 1
    if not n:
        plt.close(fig)
        return None
    ax.set_yscale('symlog', linthresh=1e-6)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel('net mass in through the measuring face / ρ$_{ph}$ [R$_\\odot^3$]')
    ax.set_title('envelope only (no stream): breathing through the measuring face')
    ax.legend(frameon=False, fontsize=9)
    fn = os.path.join(out, 'fig_envonly.png')
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('  wrote', fn)
    return fn


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', nargs='+', required=True, help='arm run directories')
    ap.add_argument('--names', nargs='*', help='arm names (default: dir basenames)')
    ap.add_argument('--vis', default=os.path.join(HERE, 'vis'),
                    help='dir with a bin_convert.py that reads version-1.2 bins')
    ap.add_argument('--out', required=True)
    ap.add_argument('--input', help='input file (default: from the run directory)')
    ap.add_argument('--extra-logs', nargs='*', default=[],
                    help='logs outside the run dirs (e.g. slurm out), scanned for every arm '
                         'whose name appears in the file name')
    ap.add_argument('--vmin', type=float, default=-7.0)
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--edge-frac', type=float, default=1e-2,
                    help='edge threshold in units of env_rho_ph')
    ap.add_argument('--no-maps', action='store_true')
    ap.add_argument('--no-edge', action='store_true')
    ap.add_argument('--force', action='store_true', help='redo maps/edge even if up to date')
    ap.add_argument('--envonly', nargs='*', default=[],
                    help='envelope-only test run dirs for fig_envonly.png (name=dir)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    names = a.names or [os.path.basename(os.path.normpath(r)) for r in a.run]
    arms, edges = [], []
    for arm, run in zip(names, a.run):
        txt, src = read_input_text(run, a.input)
        hf = sorted(glob.glob(os.path.join(run, '*.user.hst')))
        h = read_hst(hf[0]) if hf else None
        fn = os.path.join(a.out, f'numbers_{arm}.json')
        if txt is None or h is None:
            why = 'no input/restart' if txt is None else 'no usable user history yet'
            print(f'{arm}: NO DATA ({why}) in {run}')
            with open(fn, 'w') as f:
                json.dump(dict(arm=arm, no_data=True, reason=why, run=run), f, indent=1)
            continue
        print(f'{arm}: input from {src}')
        B = Binary(parse_input(txt), rmeas_from_log(run))
        print(f'  M_a {B.ma} M_d {B.md} a {B.a} R_acc {B.racc:.4f} r_meas {B.rmeas:.4f} '
              f'P_orb {B.porb:.6f} code, j_K(r_meas) {B.jk:.1f}, stream {B.stream}')
        maps = None
        if not a.no_maps:
            maps = fig_maps(arm, run, B, a.vis, a.out, a.vmin, a.vmax, a.force)
        z = None if a.no_edge else edge_data(arm, run, B, a.vis, a.out, a.edge_frac, a.force)
        edges.append((arm, z, B))
        fig_budget(arm, h, B, a.out)
        extra = [x for x in a.extra_logs if arm in os.path.basename(x)]
        nf, dtl, logs = fatal_and_dt(run, extra)
        hd = read_hst(glob.glob(os.path.join(run, '*.hydro.hst'))[0]) \
            if glob.glob(os.path.join(run, '*.hydro.hst')) else None
        # the last history row is the shortened step onto tlim: leave it out
        hdt = h['dt'][:-1] if len(h['dt']) > 2 else h['dt']
        hddt = (hd['dt'][:-1] if len(hd['dt']) > 2 else hd['dt']) if hd else None
        dtmin = np.nanmin([dtl, hdt.min(), hddt.min() if hd else np.nan])
        jl = jmean_last(h['time'], h['MR'], h['JR'], 0.1*B.porb)
        num = dict(arm=arm, run=run, no_data=False, maps_png=bool(maps),
                   t_final_orbits=float(h['time'][-1]/B.porb),
                   t_final_code=float(h['time'][-1]), P_orb_code=B.porb,
                   R_acc=B.racc, r_meas=B.rmeas, nx1=B.nx1, spin=B.spin,
                   env_rho_ph=B.rho_ph, env_top_mode=B.top_mode,
                   stream_width=B.stream_width, stream=B.stream,
                   M_stream_in_cum=float(h['Min'][-1]),
                   M_through_Racc_cum=float(h['MR'][-1]),
                   M_through_rmeas_min=float(h['MR'].min()),
                   M_through_rmeas_max=float(h['MR'].max()),
                   accreted_fraction=float(h['MR'][-1]/h['Min'][-1])
                   if h['Min'][-1] != 0 else None,
                   M_out_cum=float(h['Mout'][-1]),
                   M_env_initial=float(h['Menv'][0]), M_env_final=float(h['Menv'][-1]),
                   M_dom_initial=float(h['Mdom'][0]), M_dom_final=float(h['Mdom'][-1]),
                   j_through_Racc_cum_over_jK=float(h['JR'][-1]/h['MR'][-1]/B.jk)
                   if h['MR'][-1] != 0 else None,
                   j_through_Racc_last0p1orb_over_jK=float(jl/B.jk)
                   if np.isfinite(jl) else None,
                   j_through_cum_abs=float(h['JR'][-1]/h['MR'][-1])
                   if h['MR'][-1] != 0 else None,
                   j_stream_over_jK=float(h['Jin'][-1]/h['Min'][-1]/B.jk)
                   if h['Min'][-1] != 0 else None,
                   j_ballistic_over_jK=float(B.jbal/B.jk) if np.isfinite(B.jbal) else None,
                   j_Kep_Racc_code=B.jk, j_Kep_Racc_true_code=B.jk_acc,
                   j_star_rotation_Racc_over_jK=float(B.jcor/B.jk),
                   fatal_count=nf, logs_scanned=logs, dt_min_code=float(dtmin),
                   mass_unit='rho_stream Rsun^3', j_unit='Rsun km/s (inertial)',
                   note='"Racc" keys refer to the measuring face r_meas (= R_acc when '
                        'problem/r_meas is unset); j_K is j_Kep(r_meas)')
        if z is not None:
            num.update(edge_metrics(z, B)[0])
            num['edge_n_dumps'] = int(len(z['t']))
        with open(fn, 'w') as f:
            json.dump(num, f, indent=1)
        print('  wrote', fn)
        arms.append((arm, h, B))
    if arms:
        fig_cmp(arms, a.out)
    if not a.no_edge:
        fig_edge(edges, a.out)
    if a.envonly:
        nm = [s.split('=', 1)[0] for s in a.envonly]
        dr = [s.split('=', 1)[1] for s in a.envonly]
        fig_envonly(dr, a.out, nm)


if __name__ == '__main__':
    main()
