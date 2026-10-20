"""ACC3D page bundle 10-20 (portable): shared settings + readers. READ ONLY on the run.
Everything comes from the environment (defaults = Raven, run12). upd_page.sh exports them.
  RUN     run parent dir holding static/ and stream/ (each: ryper.{user,hydro,trc}.hst, run.log, bin/ryper.{hydro_w,m1}.NNNNN.bin)
  OUT     work dir of the page (data caches, img/, movie frames, publish/)
  VISPY   optional dir with an AthenaK vis/python bin_convert.py to use instead of the bundled binread.py
  EOSTBL  the run's general-EOS table dump (<hydro>/eos_table_dump; md5 EOSMD5)
  OPAC    Rosseland opacity table of the deck (for the optical depth of the accreted gas)
  SITE, LAYOUT, BIN_LABEL, DECK_LABEL, RUN_LABEL   text on the page / in the figures
Grid: TRUE stretched faces from the format-1.2 per-MeshBlock arrays mb_x1f/mb_x2f/mb_x3f (never linspace)."""
import os
import re
import sys
import glob
import json
import hashlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E = os.environ.get
RUN = E('RUN', '/raven/ptmp/jinma/acc3d_1017/prod/run12')
OUT = E('OUT', '/raven/ptmp/jinma/acc3d_1017/page12')
EOSTBL = E('EOSTBL', '/raven/ptmp/jinma/acc3d_1017/mk/eosdump/eos_x0.70_y0.28_logd-18.txt')
EOSMD5 = E('EOSMD5', 'e548d9c2777b8b63f979c024ab2dd564')
OPAC = E('OPAC', os.path.join(HERE, '..', 'acc3d-prod-1018', 'rosseland_mesakap_gs98_x0.70_z0.02.txt'))
SITE = E('SITE', 'Raven')
LAYOUT = E('LAYOUT', '2 nodes x 4 A100, 8 MeshBlocks 420x58x64')
BIN_LABEL = E('BIN_LABEL', 'prodmerge-1019b @ 13206541 (A100 nofma md5 29caa1ab)')
DECK_LABEL = E('DECK_LABEL', 'acc3d_prod_v11amwtr')
RUN_LABEL = E('RUN_LABEL', 'fresh tracer run (run12)')
FORCE = E('FORCE', '') == '1'
NP = int(E('NP', '4'))
TS = float(E('TS', '0.05'))              # static phase end = stream switch-on (code time)
PORB = 0.473169                          # code time (3.81 d)
DAYS = 8.0521                            # days per code time
MDOT = float(E('STREAM_MDOT', '254.843'))   # nozzle rate, code mass per code time (1e-4 Msun/yr)
DL1, RIN, RACC, RMEAS = 16.2266, 6.3, 9.0013, 9.572888
OM = 13.28                               # Omega_orb, code
ARAD, TU = 3.0981247e-09, 120.272355     # a_rad (code), K per code temperature unit
RHO_U = 1.7205983e25 / 6.957e10 ** 3     # g/cc
E_U = 1.0e10                             # (km/s)^2 -> erg/g
DIRS = [d for d in (os.path.join(RUN, 'static'), os.path.join(RUN, 'stream')) if os.path.isdir(d)]

if E('VISPY'):
    sys.path.insert(0, E('VISPY'))
    from bin_convert import read_binary  # noqa: E402,F401
else:
    sys.path.insert(0, HERE)
    from binread import read_binary  # noqa: E402,F401
sys.path.insert(0, HERE)

DARK = {'figure.facecolor': '#0b0d10', 'axes.facecolor': '#12161b', 'savefig.facecolor': '#0b0d10',
        'text.color': '#e8ecf1', 'axes.labelcolor': '#e8ecf1', 'xtick.color': '#a3adb8',
        'ytick.color': '#a3adb8', 'axes.edgecolor': '#3a424c', 'axes.titlesize': 11,
        'font.size': 10, 'legend.fontsize': 8, 'legend.facecolor': '#12161b',
        'legend.edgecolor': '#3a424c', 'grid.color': '#262c33'}
C1, C2, C3, C4, C5 = '#ffb347', '#4cc3ff', '#7dff8a', '#c77dff', '#ff5c7a'


def mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update(DARK)
    return plt


def md5(fn, n=None):
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()[:n] if n else h.hexdigest()


# ---------------------------------------------------------------- history files
def readhst(f):
    lab = None
    for ln in open(f):
        if ln.startswith('#  ['):
            lab = [x.split('=')[1] for x in ln[1:].split()]
            break
    d = np.atleast_2d(np.loadtxt(f, comments='#'))
    return {n: d[:, i] for i, n in enumerate(lab)}


def hstcat(name):
    """static + stream hst chained in time. A restart appends in place and repeats rows from the rst time on:
    the LAST occurrence of each time wins (rows of a later segment replace earlier rows with t >= its first t)."""
    rows, lab = [], None
    for d in DIRS:
        f = f'{d}/ryper.{name}.hst'
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            continue
        c = readhst(f)
        lab = lab or list(c)
        rows.append(np.stack([c[k] for k in lab], axis=1))
    if not rows:
        return None
    A = np.concatenate(rows)
    keep = np.ones(len(A), bool)
    tmin = np.inf
    for i in range(len(A) - 1, -1, -1):        # walk backwards: drop rows not earlier than a later row
        if A[i, 0] >= tmin:
            keep[i] = False
        else:
            tmin = A[i, 0]
    A = A[keep]
    return {k: A[:, i] for i, k in enumerate(lab)}


def has_trc():
    return any(os.path.exists(f'{d}/ryper.trc.hst') for d in DIRS)


# Cumulative registers: pre-acctrc-1019c binaries restarted them from 0 in every link / phase (chain them);
# 1019c and later (13206541) carry them in the rst, so they are continuous. HST_CHAIN=auto: chain only without trc.hst.
CHAIN = {'auto': not has_trc(), '1': True, '0': False}[E('HST_CHAIN', 'auto')]


def cum(c):
    out = np.array(c, float).copy()
    if not CHAIN:
        return out
    off = 0.0
    for i in range(1, len(c)):
        if abs(c[i]) < 0.5 * abs(c[i - 1]) and abs(c[i - 1] - c[i]) > 1e-6 * (abs(c[i - 1]) + 1e-30):
            off += c[i - 1]
        out[i] = c[i] + off
    return out


def rate(t, c, w=5):
    cc = cum(c)
    if len(t) <= w:
        return t[:0], cc[:0]
    tm = 0.5 * (t[w:] + t[:-w])
    return tm, (cc[w:] - cc[:-w]) / np.maximum(t[w:] - t[:-w], 1e-30)


# ---------------------------------------------------------------- dumps
_TC = None


def bin_time(f):
    """time stamp from the text pre-header (cached by path, size, mtime)."""
    global _TC
    cf = os.path.join(OUT, 'data', 'bintimes.json')
    if _TC is None:
        try:
            _TC = json.load(open(cf))
        except Exception:
            _TC = {}
    st = os.stat(f)
    key = f'{f}|{st.st_size}|{int(st.st_mtime)}'
    if key not in _TC:
        with open(f, 'rb') as fp:
            head = fp.read(600).decode('latin-1')
        _TC[key] = float(re.search(r'time=([0-9.eE+-]+)', head).group(1))
        os.makedirs(os.path.dirname(cf), exist_ok=True)
        try:
            json.dump(_TC, open(cf + '.tmp', 'w')); os.replace(cf + '.tmp', cf)
        except Exception:
            pass
    return _TC[key]


def dumps(kind='hydro_w', settle=120):
    """time-ordered [(t, path, phase)] over static + stream; files younger than `settle` s are skipped (being written).
    Same time in both phases: the stream copy wins."""
    out = {}
    for d in DIRS:
        ph = os.path.basename(d)
        for f in sorted(glob.glob(f'{d}/bin/ryper.{kind}.*.bin')):
            if time_now() - os.path.getmtime(f) < settle:
                continue
            out[round(bin_time(f), 6)] = (bin_time(f), f, ph)
    return [out[k] for k in sorted(out)]


def time_now():
    import time
    return time.time()


def m1_for(th, m1list=None):
    """m1 dump at the same time as a hydro dump (|dt| < 1e-4), or None."""
    m1list = dumps('m1') if m1list is None else m1list
    for t, f, _ in m1list:
        if abs(t - th) < 1e-4:
            return f
    return None


def dnum(f):
    return int(f.split('.')[-2])


def load(f, vars=None):
    """one .bin dump on the true grid: dict with x1f x2f x3f (faces), x1v x2v x3v (centres from the file), Time,
    and every variable as (nphi, ntheta, nr) float32 arrays. Uniform-level grids only (as all production runs)."""
    d = read_binary(f)
    if 'mb_x1f' not in d:
        raise RuntimeError(f'{f}: bin format 1.1 has no true face arrays; the stretched r grid cannot be rebuilt')
    lv = d['mb_logical']
    if len(set(lv[:, 3])) != 1:
        raise RuntimeError('mesh refinement in the dump: not supported')
    b1, b2, b3 = d['nx1_out_mb'], d['nx2_out_mb'], d['nx3_out_mb']
    n1, n2, n3 = (lv[:, 0].max() + 1) * b1, (lv[:, 1].max() + 1) * b2, (lv[:, 2].max() + 1) * b3
    x = {}
    for dd, (n, b, col) in enumerate(((n1, b1, 0), (n2, b2, 1), (n3, b3, 2))):
        xf = np.full(n + 1, np.nan); xv = np.full(n, np.nan)
        for m in range(len(lv)):
            o = lv[m, col] * b
            xf[o:o + b + 1] = d[f'mb_x{dd + 1}f'][m]; xv[o:o + b] = d[f'mb_x{dd + 1}v'][m]
        x[f'x{dd + 1}f'] = xf; x[f'x{dd + 1}v'] = xv
    out = dict(x, Time=float(d['time']), cycle=int(d['cycle']), var_names=d['var_names'])
    for v in (vars or d['var_names']):
        if v not in d['mb_data']:
            continue
        a = np.empty((n3, n2, n1), np.float32)
        for m in range(len(lv)):
            i, j, k = lv[m, 0] * b1, lv[m, 1] * b2, lv[m, 2] * b3
            a[k:k + b3, j:j + b2, i:i + b1] = d['mb_data'][v][m]
        out[v] = a
    dr = np.diff(out['x1f'])
    if n1 > 8 and np.ptp(dr) < 1e-6 * dr.mean():
        print(f'WARNING {f}: radial faces are uniform -- expected the stretched production grid', flush=True)
    return out


def dvol(h):
    rf, tf, pf = h['x1f'], h['x2f'], h['x3f']
    return ((rf[1:] ** 3 - rf[:-1] ** 3) / 3)[None, None, :] * (np.cos(tf[:-1]) - np.cos(tf[1:]))[None, :, None] \
        * np.diff(pf)[:, None, None]


_TG = None


def tgas():
    global _TG
    if _TG is None:
        from tgas import TGas
        _TG = TGas()
    return _TG


def thermo(h, m1f=None):
    """log10 T_gas [K] and log10 P_gas [code] from the EOS table (float64), and P_rad/P_gas with P_rad = E/3
    from the m1 dump (the m1 dump carries E and F only, no VET tensor), when m1f is given."""
    TG = tgas()
    rho = h['dens'].astype(np.float64)
    lT = TG.logT(rho, h['eint'].astype(np.float64))
    lPg = TG.logPgas(rho, lT)
    o = dict(lT=lT, lPg=lPg)
    if m1f:
        mm = load(m1f, ['m1_e'])
        assert abs(mm['Time'] - h['Time']) < 1e-4, (m1f, mm['Time'], h['Time'])
        o['m1_t'] = mm['Time']
        o['lPr'] = np.log10(np.maximum(mm['m1_e'].astype(np.float64), 1e-300) / 3.0)
        o['lbeta'] = o['lPr'] - lPg
        o['lTrad'] = 0.25 * np.log10(np.maximum(mm['m1_e'].astype(np.float64), 1e-300) / ARAD * TU ** 4)
    return o


def stream_dirs_label():
    return ', '.join(os.path.relpath(d, RUN) for d in DIRS)


def nrm(key, lo, hi):
    """pcolormesh colour-scale kwargs; P_rad/P_gas is diverging about 1 (log 0) with unequal ends."""
    if key == 'lbeta':
        from matplotlib.colors import TwoSlopeNorm
        return dict(norm=TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi))
    return dict(vmin=lo, vmax=hi)
