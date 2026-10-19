#!/usr/bin/env python3
"""Plaskett accretor page figures (ry_per_accretor, problem/inner = envelope).

Per arm (one --run directory each):
  fig_<arm>_maps.png    face-on theta-mean density, 4 times across the run, full domain
                        (top) and a zoom on the accretor surface (bottom); ballistic
                        stream (white dashed, clipped to the grid), inner boundary (grey
                        disk), R_acc (blue circle)
  fig_<arm>_budget.png  cumulative stream inflow / mass through R_acc / outflow, and the
                        specific AM of the stream and of the gas through R_acc
  numbers_<arm>.json    final t, cumulative inflow, mass through R_acc, mean j over the
                        last 0.1 orbit, FATAL count, min dt
and, for all arms together, fig_plaskett_cmp.png.

Binary parameters (masses, separation, R_acc, donor T, mu) are read from the input
embedded in the first restart file (includes command-line overrides), else from the
*.athinput in the run directory, else from --input.  History columns are found by their
labels in <run>/ryper.user.hst (RyPerHistEnv order: Mdom Jdom Menv Jenv Min Mout Mwal
Jwal Jstr Jout MR JR Jin, then the rates).  Units: Rsun, km/s, time 6.957e5 s.
Deps: numpy, scipy, matplotlib, and <athenak>/vis/python/bin_convert.py (--vis).

usage: mkfigs_plaskett.py --run DIR [DIR ...] --vis <athenak>/vis/python --out DIR
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
COLS = ('#2a78d6', '#eb6834', '#2f9e44', '#8e44ad')


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
    sys.exit(f'no restart or *.athinput in {run}; pass --input')


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


class Binary:
    """Roche geometry + the ballistic L1 stream, same launch as the pgen
    (from L1 toward the accretor at the donor sound speed)."""

    def __init__(self, par):
        def g(k, d):
            return float(par.get(('problem', k), d))
        self.ma, self.md = g('m_acc', 6.24), g('m_don', 1.69)
        self.a = g('a_sep', 30.3)
        self.racc = g('r_acc', 4.06)
        self.tdon, self.mudon = g('t_don', 6250.0), g('mu_don', 1.27)
        self.rin = float(par[('mesh', 'x1min')])
        self.rout = float(par[('mesh', 'x1max')])
        self.stream = par.get(('problem', 'stream'), 'true').lower() == 'true'
        self.spin = g('spin', 1.0)
        self.om = np.sqrt(GU*(self.ma + self.md)/self.a**3)     # code 1/time
        self.porb = 2*np.pi/self.om
        self.aom = self.a*self.om                               # km/s
        self.jk = np.sqrt(GU*self.ma*self.racc)                 # j_Kep(R_acc)
        self.jcor = self.spin*self.om*self.racc**2              # star's rotation at R_acc
        mu = self.md/(self.ma + self.md)
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
        # inertial specific AM at the first inward crossing of R_acc (code units)
        jd = (x - xa)*vy - y*vx + ((x - xa)**2 + y**2)
        i = np.argmax(rb < self.racc)
        self.jbal = np.nan
        if rb[i] < self.racc and i > 0:
            f = (rb[i - 1] - self.racc)/(rb[i - 1] - rb[i])
            self.jbal = (jd[i - 1] + f*(jd[i] - jd[i - 1]))*self.a*self.aom
            ph = np.degrees(np.arctan2(self.by[i - 1], self.bx[i - 1]))
            print(f'  ballistic at R_acc: phi {ph:.3f} deg, '
                  f'j/j_K {self.jbal/self.jk:.4f}')


# ---------------------------------------------------------------- history
def read_hst(fn):
    lab = None
    with open(fn) as f:
        for line in f:
            if line.startswith('#') and '[1]=' in line:
                lab = re.findall(r'\[\d+\]=(\S+)', line)
    h = np.loadtxt(fn, ndmin=2)
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


def fatal_and_dt(run):
    nf, dts = 0, []
    logs = [f for f in glob.glob(os.path.join(run, '*'))
            if re.search(r'(\.log|\.out|\.err)$', f) and os.path.isfile(f)]
    for fn in logs:
        with open(fn, errors='replace') as f:
            for line in f:
                if 'FATAL' in line:
                    nf += 1
                m = re.search(r'\bdt=([0-9.eE+-]+)', line)
                if m:
                    dts.append(float(m.group(1)))
    return nf, (min(dts) if dts else np.nan), [os.path.basename(x) for x in logs]


# ---------------------------------------------------------------- maps
def load(fn, vis):
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
    ax.plot(B.racc*np.cos(th), B.racc*np.sin(th), color='#2a78d6', lw=1.0, zorder=4)
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


def fig_maps(arm, run, B, vis, out, vmin, vmax):
    fns = sorted(glob.glob(os.path.join(run, 'bin', '*.hydro_w.*.bin')))
    if not fns:
        print(f'  {arm}: no bin dumps, maps skipped')
        return
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
    # ballistic orbit), +-0.35 R_acc around the point where it crosses R_acc
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
    fn = os.path.join(out, f'fig_{arm}_maps.png')
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    print('  wrote', fn)


# ---------------------------------------------------------------- budget
def jref_lines(ax, B):
    for v, lab, ls in ((B.jbal, 'ballistic', ':'), (B.jcor, 'star rotation', '-.')):
        if np.isfinite(v):
            ax.axhline(v/B.jk, color='0.4', ls=ls, lw=1.1)
            ax.text(0.02, v/B.jk, f' {lab} {v/B.jk:.3f}', va='bottom', fontsize=9,
                    color='0.3', transform=ax.get_yaxis_transform())


def jmean_last(t, M, J, tw):
    sel = t >= t[-1] - tw
    if sel.sum() < 2:
        return np.nan
    dm = M[-1] - M[sel][0]
    return (J[-1] - J[sel][0])/dm if abs(dm) > 0 else np.nan


def fig_budget(arm, h, B, out):
    t = h['time']/B.porb
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    ax = axs[0]
    ax.plot(t, h['Min'], color=COLS[0], lw=2, label='stream inflow (window)')
    ax.plot(t, h['MR'], color=COLS[1], lw=2, label='net mass in through R$_{acc}$')
    ax.plot(t, h['Mout'], color=COLS[2], lw=1.5, label='outflow through r$_{out}$')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel('cumulative mass [ρ$_{stream}$ R$_\\odot^3$]')
    ax.set_title(arm)
    ax.legend(frameon=False, fontsize=9)
    ax = axs[1]
    for M, J, c, lab in ((h['Min'], h['Jin'], COLS[0], 'brought in by the stream'),
                         (h['MR'], h['JR'], COLS[1], 'gas through R$_{acc}$')):
        ok = np.abs(M) > 1e-2*max(np.abs(M).max(), 1e-30)
        if ok.any():
            ax.plot(t[ok], J[ok]/M[ok]/B.jk, color=c, lw=2, label=lab + ' (cumulative)')
    jref_lines(ax, B)
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('time [orbits]')
    ax.set_ylabel('inertial specific AM / j$_{Kep}$(R$_{acc}$)')
    ax.legend(frameon=False, fontsize=9)
    fn = os.path.join(out, f'fig_{arm}_budget.png')
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('  wrote', fn)


def fig_cmp(arms, out):
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for k, (arm, h, B) in enumerate(arms):
        t = h['time']/B.porb
        c = COLS[k % len(COLS)]
        axs[0].plot(t, h['MR'], color=c, lw=2, label=f'{arm}: through R$_{{acc}}$')
        axs[0].plot(t, h['Min'], color=c, lw=1, ls='--', label=f'{arm}: stream inflow')
        ok = np.abs(h['MR']) > 1e-2*max(np.abs(h['MR']).max(), 1e-30)
        if ok.any():
            axs[1].plot(t[ok], h['JR'][ok]/h['MR'][ok]/B.jk, color=c, lw=2, label=arm)
    B = arms[0][2]
    jref_lines(axs[1], B)
    axs[1].axhline(1, color='0.75', lw=0.8)
    axs[1].text(0.02, 1, ' Keplerian', va='bottom', fontsize=9, color='0.5',
                transform=axs[1].get_yaxis_transform())
    for ax in axs:
        ax.axhline(0, color='0.6', lw=0.8)
        ax.set_xlabel('time [orbits]')
        ax.legend(frameon=False, fontsize=9)
    axs[0].set_ylabel('cumulative mass [ρ$_{stream}$ R$_\\odot^3$]')
    axs[1].set_ylabel('j of gas through R$_{acc}$ (cumulative) / j$_{Kep}$(R$_{acc}$)')
    axs[0].set_title('Plaskett: accreted mass')
    axs[1].set_title('Plaskett: specific AM through R$_{acc}$')
    fn = os.path.join(out, 'fig_plaskett_cmp.png')
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print('  wrote', fn)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', nargs='+', required=True, help='arm run directories')
    ap.add_argument('--names', nargs='*', help='arm names (default: dir basenames)')
    ap.add_argument('--vis', required=True, help='<athenak>/vis/python')
    ap.add_argument('--out', required=True)
    ap.add_argument('--input', help='input file (default: from the run directory)')
    ap.add_argument('--vmin', type=float, default=-7.0)
    ap.add_argument('--vmax', type=float, default=3.0)
    ap.add_argument('--no-maps', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    names = a.names or [os.path.basename(os.path.normpath(r)) for r in a.run]
    arms = []
    for arm, run in zip(names, a.run):
        txt, src = read_input_text(run, a.input)
        print(f'{arm}: input from {src}')
        B = Binary(parse_input(txt))
        print(f'  M_a {B.ma} M_d {B.md} a {B.a} R_acc {B.racc:.4f} P_orb {B.porb:.6f} '
              f'code, j_K {B.jk:.1f}, stream {B.stream}')
        hf = sorted(glob.glob(os.path.join(run, '*.user.hst')))
        if not hf:
            sys.exit(f'no *.user.hst in {run}')
        h = read_hst(hf[0])
        if not a.no_maps:
            fig_maps(arm, run, B, a.vis, a.out, a.vmin, a.vmax)
        fig_budget(arm, h, B, a.out)
        nf, dtl, logs = fatal_and_dt(run)
        hd = read_hst(glob.glob(os.path.join(run, '*.hydro.hst'))[0]) \
            if glob.glob(os.path.join(run, '*.hydro.hst')) else {}
        dtmin = np.nanmin([dtl, h['dt'].min(), hd['dt'].min() if hd else np.nan])
        jl = jmean_last(h['time'], h['MR'], h['JR'], 0.1*B.porb)
        num = dict(arm=arm, t_final_orbits=float(h['time'][-1]/B.porb),
                   t_final_code=float(h['time'][-1]), P_orb_code=B.porb,
                   M_stream_in_cum=float(h['Min'][-1]),
                   M_through_Racc_cum=float(h['MR'][-1]),
                   M_out_cum=float(h['Mout'][-1]),
                   j_through_Racc_cum_over_jK=float(h['JR'][-1]/h['MR'][-1]/B.jk)
                   if h['MR'][-1] != 0 else None,
                   j_through_Racc_last0p1orb_over_jK=float(jl/B.jk)
                   if np.isfinite(jl) else None,
                   j_stream_over_jK=float(h['Jin'][-1]/h['Min'][-1]/B.jk)
                   if h['Min'][-1] != 0 else None,
                   j_ballistic_over_jK=float(B.jbal/B.jk), j_Kep_Racc_code=B.jk,
                   j_star_rotation_Racc_over_jK=float(B.jcor/B.jk),
                   fatal_count=nf, logs_scanned=logs, dt_min_code=float(dtmin),
                   mass_unit='rho_stream Rsun^3', j_unit='Rsun km/s (inertial)')
        fn = os.path.join(a.out, f'numbers_{arm}.json')
        with open(fn, 'w') as f:
            json.dump(num, f, indent=1)
        print('  wrote', fn)
        arms.append((arm, h, B))
    fig_cmp(arms, a.out)


if __name__ == '__main__':
    main()
