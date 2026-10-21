#!/usr/bin/env python3
"""rsgE page (DeltaAI): movies from the reduced dumps on NATIVE cubed-sphere cells (READ ONLY on the run; reads only
OUT/data/red_*.npz with ver >= 2).  usage: mk_movies.py <red_dir> <out_dir>
Movies (OUT/publish/movies/), one frame per dump, colour limits pooled over ALL dumps (frames/<movie>/limits.json; if a new
dump moves a limit by more than LIM_TOL of its range every frame of that movie is re-rendered, else only new frames):
  rsgE_eq.mp4     equatorial plane, whole star |x|,|y| < 1.3 R_ph: log rho, v_r, T_gas            (2400 px wide)
  rsgE_inner.mp4  equatorial plane, inner zoom r < 0.25 R_ph: dlog rho, v_r, dT_gas/T_gas          (2400 px wide)
  rsgE_phot.mp4   photosphere band 0.9-1.2 R_ph, eq ring unwrapped (phi, r): log rho, v_r, dT_gas/T_gas  (1920 x 1080)
  rsgE_map.mp4    Mollweide, native cells, shells 0.9 and 0.98 R_ph: v_r - <v_r>, dT_gas/T_gas     (1920 px wide)
H.264 yuv420p, CRF (default 18), FPS (6), +faststart; if a file exceeds MAXMB (14) it is re-encoded with CRF + 3, + 6, ...
Also OUT/img/fig_eq_latest.jpg (latest full-star frame) and OUT/img/movies.json (limits, times) for build_page.py."""
import os
import sys
import glob
import json
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import csplot as cs  # noqa: E402

RED, OUT = sys.argv[1], sys.argv[2]
RPH0 = float(os.environ.get('RPH0', '5.80019e13'))
RMAX = float(os.environ.get('MOVIE_RMAX', '1.3'))
ZR = float(os.environ.get('ZOOM_R', '0.25'))
PB = (0.9, 1.2)
FPS = os.environ.get('FPS', '6')
CRF = int(os.environ.get('CRF', '18'))
MAXMB = float(os.environ.get('MAXMB', '14'))
TOL = float(os.environ.get('LIM_TOL', '0.15'))
VRPOOL = float(os.environ.get('VRPOOL', '1.0'))
DAY = 86400.
FF = os.environ.get('FFMPEG') or shutil.which('ffmpeg')
FD = os.path.join(os.path.abspath(OUT), 'frames')
MV = os.path.join(OUT, 'publish', 'movies')
MAPSH = ('090', '098')
for d in (MV, os.path.join(OUT, 'img')):
    os.makedirs(d, exist_ok=True)
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#0b0d10', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
                     'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'font.size': 13})

reds = []
for f in sorted(glob.glob(os.path.join(RED, 'red_*.npz'))):
    try:
        z = np.load(f)
        if 'ver' in z.files and int(z['ver']) >= 2:
            reds.append(f)
    except Exception:
        pass
if not reds:
    print('mk_movies: no reduced dumps with ver >= 2'); sys.exit(0)
idx = [os.path.basename(f)[4:9] for f in reds]


def load(f):
    d = np.load(f)
    rf = np.asarray(d['rf'], float) / RPH0
    ns = int(d['ns'])
    lr = np.asarray(d['eqa_lrho'], float); Tg = np.asarray(d['eqa_Tg'], float)
    D = dict(t=float(d['t']) / DAY, rf=rf, ns=ns, phi=np.asarray(d['eqa_phi'], float), lr=lr, Tg=Tg,
             dlr=lr - np.log10(np.asarray(d['rho'][:ns], float))[None, :], vr=np.asarray(d['eqa_vr'], float),
             dT=100 * (Tg / np.asarray(d['Tgas'][:ns], float)[None, :] - 1), fl=np.asarray(d['eqa_Tfl']))
    D['rc'] = 0.5 * (rf[1:ns + 1] + rf[:ns])
    for t in MAPSH:
        v = np.asarray(d['sh_vr_' + t], float)
        D['mvr_' + t] = v - cs.shell_mean(v)
        D['mdT_' + t] = 100 * np.asarray(d['sh_dT_' + t], float)
        D['mfl_' + t] = np.asarray(d['sh_Tfl_' + t])
    return D


def sub(a, m=None):
    a = np.asarray(a)[m] if m is not None else np.asarray(a)
    a = a.ravel()
    return a[np.isfinite(a)][::5]


# ---------------- per-movie pooled limits ----------------
POOLS = {'eq': {'lr': [], 'vr': [], 'Tg': []}, 'inner': {'dlr': [], 'vr': [], 'dT': []},
         'phot': {'lr': [], 'vr': [], 'dT': []}, 'map': {k + '_' + t: [] for k in ('mvr', 'mdT') for t in MAPSH}}
for f in reds:
    D = load(f)
    rc = D['rc']; nofl = ~D['fl']
    full = rc < RMAX; inn = rc < ZR; ph = (rc > PB[0]) & (rc < PB[1])
    POOLS['eq']['lr'].append(sub(D['lr'][:, full]))
    POOLS['eq']['vr'].append(sub(D['vr'][:, rc < VRPOOL]))
    POOLS['eq']['Tg'].append(sub(D['Tg'][:, full], nofl[:, full]))
    for k in ('dlr', 'vr', 'dT'):
        POOLS['inner'][k].append(sub(D[k][:, inn]))
    for k in ('lr', 'vr', 'dT'):
        POOLS['phot'][k].append(sub(D[k][:, ph], nofl[:, ph] if k == 'dT' else None))
    for t in MAPSH:
        POOLS['map']['mvr_' + t].append(sub(D['mvr_' + t]))
        POOLS['map']['mdT_' + t].append(sub(D['mdT_' + t], ~D['mfl_' + t]))


def P(a, lo, hi):
    a = np.concatenate(a)
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))]


def S(a, floor):
    return max(floor, float(np.percentile(np.abs(np.concatenate(a)), 99.5)))


TLO = 2000 * 1.003
LIM = {'eq': dict(lr=P(POOLS['eq']['lr'], 0.5, 99.5), vr=S(POOLS['eq']['vr'], 1e-3),
                  Tg=[TLO, max(TLO * 1.01, P(POOLS['eq']['Tg'], 0, 99.9)[1])]),
       'inner': dict(dlr=S(POOLS['inner']['dlr'], 1e-6), vr=S(POOLS['inner']['vr'], 1e-3), dT=S(POOLS['inner']['dT'], 1e-5)),
       'phot': dict(lr=P(POOLS['phot']['lr'], 0.5, 99.5), vr=S(POOLS['phot']['vr'], 1e-3),
                    dT=S(POOLS['phot']['dT'], 1e-4)),
       'map': {k: S(v, 1e-5) for k, v in POOLS['map'].items()}}
del POOLS


def moved(new, old):
    for k, v in new.items():
        o = old.get(k)
        if o is None:
            return True
        if isinstance(v, list):
            span = abs(v[1] - v[0]) or 1
            if abs(v[0] - o[0]) > TOL * span or abs(v[1] - o[1]) > TOL * span:
                return True
        elif abs(v - o) > TOL * abs(v):
            return True
    return False


REDO = {}
for name, lim in LIM.items():
    d = os.path.join(FD, name); os.makedirs(d, exist_ok=True)
    lf = os.path.join(d, 'limits.json')
    old = json.load(open(lf)) if os.path.exists(lf) else None
    REDO[name] = old is None or moved(lim, old)
    if REDO[name]:
        for g in glob.glob(os.path.join(d, '*.png')):
            os.remove(g)
        json.dump(lim, open(lf, 'w'), indent=1)
    else:
        LIM[name] = old
print('mk_movies limits', json.dumps(LIM), 'redo', REDO, flush=True)


def cbar(fig, m, a, lab, **kw):
    cb = fig.colorbar(m, ax=a, label=lab, extend='both', **kw)
    (cb.ax.xaxis if kw.get('orientation') == 'horizontal' else cb.ax.yaxis).set_major_locator(MaxNLocator(6))
    return cb


def frame_disc(D, fn, which):
    L = LIM[which]
    fig, ax = plt.subplots(1, 3, figsize=(24, 8.6), constrained_layout=True)
    if which == 'eq':
        i0, i1 = cs.rrange(D['rf'], 0, RMAX * 1.42); ext = RMAX
        pans = [(D['lr'], 'viridis', cs.norm_of('lin', L['lr']), 'log ρ [g cm⁻³]', None),
                (D['vr'], 'RdBu_r', cs.norm_of('sym', L['vr']), 'v_r [km/s]', None),
                (D['Tg'], cs.cmap_floor('inferno'), cs.norm_of('log', L['Tg']), 'T_gas [K] (grey: 2000 K floor)', D['fl'])]
        ttl = 'equatorial plane z = 0, whole star'
    else:
        i0, i1 = cs.rrange(D['rf'], 0, ZR * 1.42); ext = ZR
        pans = [(D['dlr'], 'RdBu_r', cs.norm_of('sym', L['dlr']), 'log ρ − ⟨log ρ⟩_shell', None),
                (D['vr'], 'RdBu_r', cs.norm_of('sym', L['vr']), 'v_r [km/s]', None),
                (D['dT'], 'PuOr_r', cs.norm_of('sym', L['dT']), 'δT_gas / ⟨T_gas⟩_shell [%]', None)]
        ttl = 'equatorial plane z = 0, inner zoom r < %.2f R_ph (dashed: MESA CZ base)' % ZR
    i1 = min(i1, D['ns'])
    for a, (q, cm, nm, lab, fl) in zip(ax, pans):
        m = cs.ring_xy(a, D['phi'], D['rf'], cs.with_floor(q, fl, nm), i0, i1, cm, nm)
        cbar(fig, m, a, lab, shrink=0.92)
        a.add_patch(plt.Circle((0, 0), 1.0, fill=False, ec='w', lw=0.6, ls=':'))
        if which == 'inner':
            a.add_patch(plt.Circle((0, 0), 0.0494, fill=False, ec='w', lw=0.6, ls='--'))
        a.set_xlim(-ext, ext); a.set_ylim(-ext, ext); a.set_aspect('equal')
        a.set_xlabel('x / R_ph,IC'); a.set_ylabel('y / R_ph,IC')
    fig.suptitle('RSG E (DeltaAI) — %s — native cells — t = %6.2f d' % (ttl, D['t']), fontsize=16)
    fig.savefig(fn, dpi=100)
    plt.close(fig)


def frame_phot(D, fn):
    L = LIM['phot']
    i0, i1 = cs.rrange(D['rf'], *PB); i1 = min(i1, D['ns'])
    fig, ax = plt.subplots(3, 1, figsize=(19.2, 10.8), constrained_layout=True, sharex=True)
    pans = [(D['lr'], 'viridis', cs.norm_of('lin', L['lr']), 'log ρ [g cm⁻³]', None),
            (D['vr'], 'RdBu_r', cs.norm_of('sym', L['vr']), 'v_r [km/s]', None),
            (D['dT'], cs.cmap_floor('PuOr_r'), cs.norm_of('sym', L['dT']), 'δT_gas / ⟨T_gas⟩ [%] (grey: floor)', D['fl'])]
    for a, (q, cm, nm, lab, fl) in zip(ax, pans):
        m = cs.ring_unwrap(a, D['phi'], D['rf'], cs.with_floor(q, fl, nm), i0, i1, cm, nm)
        cbar(fig, m, a, lab, pad=0.01)
        a.set_ylabel('r / R_ph,IC'); a.set_ylim(PB)
        for e in (-45, 45, 135):
            a.axvline(e, color='w', lw=0.4, ls=':')
    ax[-1].set_xlim(-135, 225); ax[-1].set_xlabel('φ [deg] in the equatorial plane (dotted: panel edges)')
    fig.suptitle('RSG E (DeltaAI) — photosphere band %.1f–%.1f R_ph,IC unwrapped, native cells (%d radial) — t = %6.2f d'
                 % (PB[0], PB[1], i1 - i0, D['t']), fontsize=15)
    fig.savefig(fn, dpi=100)
    plt.close(fig)


def frame_map(D, fn):
    L = LIM['map']
    fig, ax = plt.subplots(2, 2, figsize=(19.2, 10.4), constrained_layout=True)
    for r, t in enumerate(MAPSH):
        rr = int(t) / 100
        nv = cs.norm_of('sym', L['mvr_' + t]); nt = cs.norm_of('sym', L['mdT_' + t])
        m = cs.moll(ax[r, 0], D['mvr_' + t], 'RdBu_r', nv, W=900, cache=FD)
        cbar(fig, m, ax[r, 0], 'v_r − ⟨v_r⟩ at %.2f R_ph,IC [km/s]' % rr, orientation='horizontal', shrink=0.6, pad=0.02)
        m = cs.moll(ax[r, 1], D['mdT_' + t], cs.cmap_floor('PuOr_r'), nt, W=900, cache=FD, floor=D['mfl_' + t])
        cbar(fig, m, ax[r, 1], 'δT_gas / ⟨T_gas⟩ at %.2f R_ph,IC [%%]' % rr, orientation='horizontal', shrink=0.6, pad=0.02)
    fig.suptitle('RSG E (DeltaAI) — Mollweide, native 0.70° cells (graticule 45°: cube edges at lon ±45°, ±135°) — t = %6.2f d'
                 % D['t'], fontsize=15)
    fig.savefig(fn, dpi=100)
    plt.close(fig)


RENDER = {'eq': lambda D, fn: frame_disc(D, fn, 'eq'), 'inner': lambda D, fn: frame_disc(D, fn, 'inner'),
          'phot': frame_phot, 'map': frame_map}
nnew = {k: 0 for k in RENDER}
for f, i in zip(reds, idx):
    todo = [k for k in RENDER if not os.path.exists(os.path.join(FD, k, i + '.png'))]
    if not todo:
        continue
    D = load(f)
    for k in todo:
        RENDER[k](D, os.path.join(FD, k, i + '.png'))
        nnew[k] += 1
print('mk_movies: new frames', nnew, 'of', len(reds), flush=True)

from PIL import Image  # noqa: E402
Image.open(os.path.join(FD, 'eq', idx[-1] + '.png')).convert('RGB').save(os.path.join(OUT, 'img', 'fig_eq_latest.jpg'), quality=90)
D0, D1 = load(reds[0]), load(reds[-1])
json.dump(dict(t_first=D0['t'], t_last=D1['t'], n=len(reds), lim=LIM, zr=ZR, pb=PB, rmax=RMAX, mapsh=MAPSH),
          open(os.path.join(OUT, 'img', 'movies.json'), 'w'))


def encode(name, fname):
    out = os.path.join(MV, fname)
    stamp = os.path.join(FD, name + '.stamp')
    key = '%d %s %s %d' % (len(reds), idx[-1], FPS, CRF)
    if nnew[name] == 0 and os.path.exists(out) and os.path.exists(stamp) and open(stamp).read() == key:
        return
    lst = os.path.join(FD, name + '.txt')
    with open(lst, 'w') as fl:
        for i in idx:
            fl.write("file '%s'\nduration %.6f\n" % (os.path.join(FD, name, i + '.png'), 1.0 / float(FPS)))
        fl.write("file '%s'\n" % os.path.join(FD, name, idx[-1] + '.png'))
    crf = CRF
    while True:
        cmd = [FF, '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', lst,
               '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2,fps=%s' % FPS, '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
               '-crf', str(crf), '-preset', 'slow', '-movflags', '+faststart', out + '.tmp.mp4']
        subprocess.run(cmd, check=True)
        mb = os.path.getsize(out + '.tmp.mp4') / 1e6
        if mb <= MAXMB or crf >= 40:
            break
        crf += 3
    os.replace(out + '.tmp.mp4', out)
    open(stamp, 'w').write(key)
    print('wrote', out, '%.2f MB crf %d' % (mb, crf), flush=True)


for name, fname in (('eq', 'rsgE_eq.mp4'), ('inner', 'rsgE_inner.mp4'), ('phot', 'rsgE_phot.mp4'), ('map', 'rsgE_map.mp4')):
    encode(name, fname)
old = os.path.join(MV, 'rsgE_map090.mp4')        # superseded by rsgE_map.mp4
if os.path.exists(old):
    os.remove(old)
