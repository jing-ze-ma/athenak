#!/usr/bin/env python3
"""rsgE page (DeltaAI): movies from the reduced dumps on NATIVE cubed-sphere cells (READ ONLY on the run; reads only
OUT/data/red_*.npz with ver >= 2).  usage: mk_movies.py <red_dir> <out_dir>
Movies (OUT/publish/movies/), one frame per dump, colour limits pooled over ALL dumps (frames/<movie>/limits.json; if a new
dump moves a limit by more than LIM_TOL of its range every frame of that movie is re-rendered, else only new frames):
  rsgE_eq.mp4     equatorial plane, whole star |x|,|y| < 1.3 R_ph: log rho, v_r, T_gas            (2400 px wide)
  rsgE_inner.mp4  equatorial plane, inner zoom r < 0.25 R_ph: dlog rho, v_r, dT_gas/T_gas          (2400 px wide)
  rsgE_phot.mp4   photosphere band 0.9-1.2 R_ph, eq ring unwrapped (phi, r): log rho, v_r, dT_gas/T_gas  (1920 x 1080)
  rsgE_map.mp4    Mollweide, native cells, shells 0.9 and 0.96 R_ph: v_r - <v_r>, dT_gas/T_gas     (1920 px wide)
H.264 yuv420p, CRF (default 18), FPS (6), +faststart; if a file exceeds MAXMB (14) it is re-encoded with CRF + 3, + 6, ...
Also OUT/img/fig_eq_latest.jpg (latest full-star frame) and OUT/img/movies.json (limits, times) for build_page.py."""
import os
import sys
import glob
import json
import shutil
import subprocess
import math
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import csplot as cs  # noqa: E402

RED, OUT = sys.argv[1], sys.argv[2]
RPH0 = float(os.environ.get('RPH0', '5.80019e13'))
RMAX = float(os.environ.get('MOVIE_RMAX', '2.5'))   # 10-27 USER: CSM out to 2.5 R_ph (sponge 2.72; slices to the edge, VER 6)
ZR = float(os.environ.get('ZOOM_R', '0.25'))
PB = (0.9, 1.2)
FPS = os.environ.get('FPS', '6')
CRF = int(os.environ.get('CRF', '18'))
MAXMB = float(os.environ.get('MAXMB', '15'))   # artifact per-file cap 16 MB
DPI = int(os.environ.get('MOVIE_DPI', '160'))   # 10-23 USER: frames were 100 dpi (coarser than native near the surface)
SHARD = os.environ.get('SHARD')                    # 10-23 'k/n': render every n-th frame only, no encode (parallel renders)
TOL = float(os.environ.get('LIM_TOL', '0.15'))
VRPOOL = float(os.environ.get('VRPOOL', '1.0'))
DAY = 86400.
FF = os.environ.get('FFMPEG') or shutil.which('ffmpeg')
FD = os.path.join(os.path.abspath(OUT), 'frames')
TGREY = os.environ.get('T_FLOOR_GREY', '1') == '1'   # greytest-1102: 0 = continuous inferno down to the floor (no grey mask) in the T_gas panel
MV = os.path.join(OUT, 'publish', 'movies')
MAPSH = ('090', '096')   # USER 10-22: 0.98 sat on the photospheric T front (ring artefact)
for d in (MV, os.path.join(OUT, 'img')):
    os.makedirs(d, exist_ok=True)
plt.rcParams.update({'figure.facecolor': '#0b0d10', 'axes.facecolor': '#0b0d10', 'savefig.facecolor': '#0b0d10',
                     'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
                     'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'font.size': 13})

# pagefast-1102: the ver >= 2 check is cached per file (name, size, mtime_ns) in frames/reds_ver.json (was: every npz opened by
# every pass/shard: ~3900 zip opens each on NFS)
_VC = os.path.join(FD, 'reds_ver.json')
try:
    _vc = json.load(open(_VC))
except Exception:
    _vc = {}
reds = []; _vn = {}; _chg = False
for f in sorted(glob.glob(os.path.join(RED, 'red_*.npz'))):
    try:
        st_ = os.stat(f); k_ = os.path.basename(f); sig_ = [st_.st_size, st_.st_mtime_ns]
        if k_ in _vc and _vc[k_][:2] == sig_:
            ok_ = _vc[k_][2]
        else:
            z = np.load(f)
            ok_ = bool('ver' in z.files and int(z['ver']) >= 2); _chg = True
        _vn[k_] = sig_ + [ok_]
        if ok_:
            reds.append(f)
    except Exception:
        pass
if _chg or len(_vn) != len(_vc):
    try:
        os.makedirs(FD, exist_ok=True)
        json.dump(_vn, open(_VC + '.%d.tmp' % os.getpid(), 'w')); os.replace(_VC + '.%d.tmp' % os.getpid(), _VC)
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
def _pool_one(f):
    """pagefast-1102: the per-dump pool samples (same code as the old serial loop), run in worker processes"""
    out = {'eq': {}, 'inner': {}, 'phot': {}, 'map': {}}
    D = load(f)
    rc = D['rc']; nofl = ~D['fl']
    full = rc < RMAX; inn = rc < ZR; ph = (rc > PB[0]) & (rc < PB[1])
    out['eq']['lr'] = sub(D['lr'][:, full])
    out['eq']['vr'] = sub(D['vr'][:, rc < VRPOOL])
    out['eq']['Tg'] = sub(D['Tg'][:, full], nofl[:, full])
    for k in ('dlr', 'vr', 'dT'):
        out['inner'][k] = sub(D[k][:, inn])
    for k in ('lr', 'vr', 'dT'):
        out['phot'][k] = sub(D[k][:, ph], nofl[:, ph] if k == 'dT' else None)
    for t in MAPSH:
        out['map']['mvr_' + t] = sub(D['mvr_' + t])
        out['map']['mdT_' + t] = sub(D['mdT_' + t], ~D['mfl_' + t])
    return out


POOLS = {'eq': {'lr': [], 'vr': [], 'Tg': []}, 'inner': {'dlr': [], 'vr': [], 'dT': []},
         'phot': {'lr': [], 'vr': [], 'dT': []}, 'map': {k + '_' + t: [] for k in ('mvr', 'mdT') for t in MAPSH}}
# pagefast-1102: LIMITS_FROM=1 (shard + final passes of catchup_job.sh, after a LIMITS_ONLY pass of the same data) reads
# frames/<movie>/limits.json + frames/redo.json instead of re-pooling every dump (was 1 + NSH + 1 full passes over ~3900 npz)
_FROMFILE = bool(os.environ.get('LIMITS_FROM')) and all(os.path.exists(os.path.join(FD, n_, 'limits.json')) for n_ in POOLS) \
    and os.path.exists(os.path.join(FD, 'redo.json'))
if not _FROMFILE:
    import multiprocessing as _mp
    _np_ = max(1, min(int(os.environ.get('POOL_PROCS', '16')), len(reds)))
    with _mp.get_context('fork').Pool(_np_) as _pl:
        for o_ in _pl.imap(_pool_one, reds, chunksize=4):
            for m_, dd_ in o_.items():
                for k_, v_ in dd_.items():
                    POOLS[m_][k_].append(v_)



def P(a, lo, hi):
    a = np.concatenate(a)
    return [float(np.percentile(a, lo)), float(np.percentile(a, hi))]


def S(a, floor):
    return max(floor, float(np.percentile(np.abs(np.concatenate(a)), 99.5)))


if not _FROMFILE:
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
        elif abs(v) > (1 + TOL) * abs(o) or abs(v) < 0.5 * abs(o):   # 10-23: scalar amplitudes redo only when exceeded (or halved)
            return True
    return False


REDO = {}
if _FROMFILE:
    LIM = {n_: json.load(open(os.path.join(FD, n_, 'limits.json'))) for n_ in ('eq', 'inner', 'phot', 'map')}
    REDO = json.load(open(os.path.join(FD, 'redo.json')))
for name, lim in ({} if _FROMFILE else LIM).items():
    d = os.path.join(FD, name); os.makedirs(d, exist_ok=True)
    lf = os.path.join(d, 'limits.json')
    old = json.load(open(lf)) if os.path.exists(lf) else None
    REDO[name] = old is None or moved(lim, old)
    if REDO[name]:
        lim = {k: (v * 1.25 if not isinstance(v, list) else v) for k, v in lim.items()}   # 10-23: headroom -> rare full re-renders
        LIM[name] = lim
        for g in glob.glob(os.path.join(d, '*.png')):
            os.remove(g)
        json.dump(lim, open(lf, 'w'), indent=1)
    else:
        LIM[name] = old
print('mk_movies limits', json.dumps(LIM), 'redo', REDO, flush=True)
if os.environ.get('LIMITS_ONLY'):
    # pagefast-1102: record the redo flags (the final pass re-encodes a movie in full after a limits redo: with shards its own
    # new-frame count is 0, so the old append test missed it) and the missing frames as one task list for balanced shards
    json.dump(REDO, open(os.path.join(FD, 'redo.json.tmp'), 'w')); os.replace(os.path.join(FD, 'redo.json.tmp'), os.path.join(FD, 'redo.json'))
    _have = {n_: set(os.listdir(os.path.join(FD, n_))) for n_ in LIM}
    _tasks = [[n__, i_, n_] for n__, i_ in enumerate(idx) for n_ in ('eq', 'inner', 'phot', 'map') if i_ + '.png' not in _have[n_]]
    json.dump(_tasks, open(os.path.join(FD, 'tasks.json.tmp'), 'w')); os.replace(os.path.join(FD, 'tasks.json.tmp'), os.path.join(FD, 'tasks.json'))
    print('mk_movies: %d frames to render' % len(_tasks), flush=True)
    sys.exit(0)


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
                (D['Tg'], cs.cmap_floor('inferno') if TGREY else plt.get_cmap('inferno'), cs.norm_of('log', L['Tg']), 'T_gas [K] (grey: 2000 K floor)' if TGREY else 'T_gas [K]', D['fl'] if TGREY else None)]
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
        if which == 'inner':
            a.add_patch(plt.Circle((0, 0), 0.0494, fill=False, ec='w', lw=0.6, ls='--'))
        a.set_xlim(-ext, ext); a.set_ylim(-ext, ext); a.set_aspect('equal')
        a.set_xlabel('x / R_ph,IC'); a.set_ylabel('y / R_ph,IC')
    fig.suptitle('RSG E (Caltech) — %s — native cells — t = %6.2f d' % (ttl, D['t']), fontsize=16)
    fig.savefig(fn, dpi=DPI)
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
    fig.suptitle('RSG E (Caltech) — photosphere band %.1f–%.1f R_ph,IC unwrapped, native cells (%d radial) — t = %6.2f d'
                 % (PB[0], PB[1], i1 - i0, D['t']), fontsize=15)
    fig.savefig(fn, dpi=DPI)
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
    fig.suptitle('RSG E (Caltech) — Mollweide, native 0.70° cells (graticule 45°: cube edges at lon ±45°, ±135°) — t = %6.2f d'
                 % D['t'], fontsize=15)
    fig.savefig(fn, dpi=DPI)
    plt.close(fig)


RENDER = {'eq': lambda D, fn: frame_disc(D, fn, 'eq'), 'inner': lambda D, fn: frame_disc(D, fn, 'inner'),
          'phot': frame_phot, 'map': frame_map}
nnew = {k: 0 for k in RENDER}
_TJ = os.path.join(FD, 'tasks.json')
if _FROMFILE and os.path.exists(_TJ):
    _k, _n = (int(x) for x in SHARD.split('/')) if SHARD else (0, 1)
    _byd = {}
    for n__, i_, k_ in json.load(open(_TJ))[_k::_n]:
        _byd.setdefault((n__, i_), []).append(k_)
    for (n__, i_), ks_ in sorted(_byd.items()):
        f = os.path.join(RED, 'red_%s.npz' % i_)
        todo = [k for k in ks_ if not os.path.exists(os.path.join(FD, k, i_ + '.png'))]
        if not todo or not os.path.exists(f):
            continue
        D = load(f)
        for k in todo:
            RENDER[k](D, os.path.join(FD, k, i_ + '.png'))
            nnew[k] += 1
else:
    for n_, (f, i) in enumerate(zip(reds, idx)):
        if SHARD and n_ % int(SHARD.split('/')[1]) != int(SHARD.split('/')[0]):
            continue
        todo = [k for k in RENDER if not os.path.exists(os.path.join(FD, k, i + '.png'))]
        if not todo:
            continue
        D = load(f)
        for k in todo:
            RENDER[k](D, os.path.join(FD, k, i + '.png'))
            nnew[k] += 1

print('mk_movies: new frames', nnew, 'of', len(reds), flush=True)
if SHARD:
    sys.exit(0)

from PIL import Image  # noqa: E402
Image.open(os.path.join(FD, 'eq', idx[-1] + '.png')).convert('RGB').save(os.path.join(OUT, 'img', 'fig_eq_latest.jpg'), quality=90)
D0, D1 = load(reds[0]), load(reds[-1])
json.dump(dict(t_first=D0['t'], t_last=D1['t'], n=len(reds), lim=LIM, zr=ZR, pb=PB, rmax=RMAX, mapsh=MAPSH),
          open(os.path.join(OUT, 'img', 'movies.json'), 'w'))
# 10-24 USER "make the 2d videos as fast as the 3d vis video": unless FPS is set, play the 2-D movies at the 3-D movie's pace
# (simulated days per second of video, from vis3d/<run>/out/rsgE_3d.json), clamped to 6..60 fps
J3D = os.environ.get('MATCH3D', '/resnick/groups/carnegie_poc/jingze/vis3d/rsgE_glob/out/rsgE_3d.json')
if 'FPS' not in os.environ and os.path.exists(J3D):
    try:
        i3 = json.load(open(J3D))
        pace = (i3['t1_d'] - i3['t0_d']) / (i3['nframes'] / float(i3['fps']))      # days per video second
        FPS = '%.2f' % max(6.0, min(60.0, len(reds) / ((float(D1['t']) - float(D0['t'])) / pace)))   # reduction t is in days
    except Exception as e:
        print('MATCH3D skipped:', e, flush=True)
print('mk_movies fps', FPS, flush=True)


PREFIX = 'rsgE'
def _enc(lst_frames, fps, crf, dst, threads='8'):
    """movfix-1102 (USER): numbered-symlink image sequence + -framerate (as cs3d.py): every frame exactly once.  The old concat
    list (1/fps durations on PNG inputs, rounded to the image demuxer's 1/25 s) dropped ~58 % of the dumps and doubled others.
    Same size filter, codec, pix_fmt, crf, preset and faststart as before."""
    tdir = dst + '.in'
    shutil.rmtree(tdir, ignore_errors=True); os.makedirs(tdir)
    for k, f in enumerate(lst_frames):
        os.symlink(f, os.path.join(tdir, '%06d.png' % k))
    cmd = [FF, '-y', '-loglevel', 'error', '-framerate', str(fps), '-i', os.path.join(tdir, '%06d.png'),
           '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2', '-c:v', 'libx264', '-threads', threads, '-pix_fmt', 'yuv420p',
           '-crf', str(crf), '-preset', 'slow', '-movflags', '+faststart', dst]
    try:
        subprocess.run(cmd, check=True)
    finally:
        shutil.rmtree(tdir, ignore_errors=True)


PARTMB = float(os.environ.get('PARTMB', '14'))      # movfix-1102 (USER): no movie file above ~15 MB -> parts of <= 14 MB
SEGN = int(os.environ.get('SEGN', '200'))           # frames per encoded segment (segments are joined into parts with -c copy)
CHANGED = []


def _join(SD, segs, dst):
    cl = dst + '.join.txt'
    with open(cl, 'w') as fl:
        for g in segs:
            fl.write("file '%s'\n" % os.path.join(SD, g['f']))
    subprocess.run([FF, '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', cl, '-c', 'copy',
                    '-movflags', '+faststart', dst + '.tmp.mp4'], check=True)
    os.replace(dst + '.tmp.mp4', dst); os.remove(cl)


def _tday(i):
    try:
        return float(np.load(os.path.join(RED, 'red_%s.npz' % i))['t']) / DAY
    except Exception:
        return float('nan')


def encode(name, fname):
    """movfix-1102: the movie is a list of PARTS (<prefix>_<name>_pNN.mp4, oldest first), each <= PARTMB, each the -c copy join of
    segments of <= SEGN frames (FD/<name>_segs/*.mp4).  New dumps: new segment(s) go into the newest part while it stays
    <= PARTMB, else a new part starts; closed parts are never touched again.  Everything is re-encoded (segments in
    parallel) only when the colour limits were redone, DPI/CRF changed, the 3-D-matched fps drifted > 10 %, or the frame list
    changed.  crf = the movie's current crf (frames/<name>.crf), never raised any more (the parts carry the size)."""
    pre = fname[:-4]
    SD = os.path.join(FD, name + '_segs'); os.makedirs(SD, exist_ok=True)
    pj = os.path.join(SD, 'parts.json')
    st = json.load(open(pj)) if os.path.exists(pj) else None
    frames = [os.path.join(FD, name, i + '.png') for i in idx]
    cf = os.path.join(FD, name + '.crf')
    crf = max(CRF, int(open(cf).read())) if os.path.exists(cf) else CRF
    full = (st is None or st.get('dpi') != DPI or st.get('crf0') != CRF or st.get('crf') != crf
            or abs(float(st['fps']) / float(FPS) - 1) > 0.10 or REDO.get(name)
            or any(not os.path.exists(os.path.join(SD, g['f'])) for p in st['parts'] for g in p['segs'])
            or any(not os.path.exists(os.path.join(MV, p['file'])) for p in st['parts']))
    if not full:
        n0 = st['n']
        if idx[:n0] != st['idx'][:n0] or nnew[name] > len(idx) - n0:
            full = True
        elif len(idx) == n0:
            return
    tag = '%d' % int(time.time())
    if full:
        fps = FPS; n0 = 0
        for p in (st or {}).get('parts', []):
            for g in p['segs']:
                try:
                    os.remove(os.path.join(SD, g['f']))
                except OSError:
                    pass
            try:
                os.remove(os.path.join(MV, p['file']))
            except OSError:
                pass
        st = dict(fps=fps, crf=crf, crf0=CRF, dpi=DPI, n=0, idx=[], parts=[])
    fps = st['fps']
    cuts = list(range(n0, len(idx), SEGN)) + [len(idx)]
    new = [dict(f='%s_%03d_%s.mp4' % (name, k, tag), a=cuts[k], b=cuts[k + 1]) for k in range(len(cuts) - 1)]
    with ThreadPoolExecutor(int(os.environ.get('SEG_PAR', '4'))) as ex2:
        list(ex2.map(lambda g: _enc(frames[g['a']:g['b']], fps, crf, os.path.join(SD, g['f'])), new))
    touched = set()
    for g in new:
        mb = os.path.getsize(os.path.join(SD, g['f'])) / 1e6
        seg = dict(f=g['f'], n=g['b'] - g['a'], mb=mb, i0=idx[g['a']], i1=idx[g['b'] - 1])
        if st['parts'] and st['parts'][-1]['mb'] + mb <= PARTMB:
            p = st['parts'][-1]
        else:
            p = dict(file='%s_p%02d.mp4' % (pre, len(st['parts']) + 1), segs=[], mb=0.0)
            st['parts'].append(p)
        p['segs'].append(seg); p['mb'] += mb; touched.add(p['file'])
    for p in st['parts']:
        if p['file'] in touched:
            _join(SD, p['segs'], os.path.join(MV, p['file']))
            p['mb'] = os.path.getsize(os.path.join(MV, p['file'])) / 1e6
            p['n'] = sum(g['n'] for g in p['segs']); p['i0'] = p['segs'][0]['i0']; p['i1'] = p['segs'][-1]['i1']
            p['t0'], p['t1'] = _tday(p['i0']), _tday(p['i1'])
            CHANGED.append(p['file'])
    st['n'] = len(idx); st['idx'] = list(idx)
    json.dump(st, open(pj + '.tmp', 'w')); os.replace(pj + '.tmp', pj)
    for leg in (os.path.join(MV, fname), os.path.join(SD, 'segs.json')):   # the single-file movie and its old state
        if os.path.exists(leg):
            os.remove(leg)
    print('%s %s: %d new frames in %d segment(s) -> parts %s (%s MB), crf %d, fps %s' % (
        'rebuilt' if full else 'appended', pre, len(idx) - n0, len(new), ','.join(sorted(touched)),
        ', '.join('%.1f' % p['mb'] for p in st['parts']), crf, fps), flush=True)


import time  # noqa: E402
MOVIES = (('eq', PREFIX + '_eq.mp4'), ('inner', PREFIX + '_inner.mp4'), ('phot', PREFIX + '_phot.mp4'), ('map', PREFIX + '_map.mp4'))
with ThreadPoolExecutor(4) as ex:
    list(ex.map(lambda a: encode(*a), MOVIES))
# movfix-1102: the parts of every movie for build_page.py (newest last) + the files changed by this run (READY list)
MJp = json.load(open(os.path.join(OUT, 'img', 'movies.json')))
MJp['parts'] = {}
for name, fname in MOVIES:
    s_ = json.load(open(os.path.join(FD, name + '_segs', 'parts.json')))
    MJp['parts'][fname] = [dict(file=p['file'], n=p['n'], i0=p['i0'], i1=p['i1'], t0=p['t0'], t1=p['t1'], mb=p['mb']) for p in s_['parts']]
json.dump(MJp, open(os.path.join(OUT, 'img', 'movies.json'), 'w'))
with open(os.path.join(FD, 'changed_2d.txt'), 'w') as fl:
    fl.write(''.join('movies/%s\n' % f for f in sorted(set(CHANGED))))
