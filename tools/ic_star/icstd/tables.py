"""Opacity (ext2) and general-EOS tables by ONE recipe for any composition, cached by a
composition hash (tools/ic_star).

OPACITY = the AG Car E1 'ext2' recipe (Viper lbv_1008/agcar/tables_ext/scripts/build_ext.py
fergR, 10-09), generalised only in the composition inputs (= the bsg_newsetup_1015
build_ext_bsg.py generalisation):
  * HIGH T: LANL TOPS/ATOMIC Rosseland + Planck for the mixture (number fractions; GS98 metals
    scaled to Z via data/stellar_opac/planck_tools/mixspec.py, or a custom mixture string),
    T 0.0005-1 keV (48 T), 'hi' file rho 1e-14..1 (71 points), optional 'lo' file rho
    1e-21..1e-14 (36 points; without it the extension starts at log rho -14 in every TOPS row,
    the bsg_newsetup 'ext2nolo' case);
  * LOW T: Ferguson et al. 2005 GS98 Rosseland (*.tron) and Planck (*.tpon), linear in X between
    the bracketing X tables and linear in log Z between the bracketing Z tables (a weight of
    exactly 0 skips that axis: X 0.36 Z 0.02 = the AG Car (0.35, 0.5) X-interpolation, X 0.7
    Z 0.008 = the BSG log-Z interpolation, bit for bit); inside the Ferguson box the
    convert_ferguson.py bilinear (log T, log R) interpolation;
  * BLEND linear in log T of log kappa over log T 4.0-4.2;
  * LOW-RHO EXTENSION below each source's valid floor: Rosseland = Saha e-scattering +
    absorption (rho/rho_e)^s, Planck = kP_e (rho/rho_e)^s, s = local slope over 0.6 dex clamped
    [0, 1], boxcar-smoothed over +-0.1 dex in log T (saha_es.py);
  * grid log T 2.6..8.0 (0.025), log rho -21..0 (0.05); printed '%.5f'.
  NOT included (and listed so): the AG Car E1 table's 'bitwise patch' (log T >= 4.2 and log rho
  >= -14 copied from the older plain-TOPS table) and the 'Orion' low-rho Planck stitch
  (planck_ext2orion, an X 0.36 Z 0.02 product).

EOS = the code's general-EOS table dump (<hydro>/eos_table_dump) from a 1-cell shock_tube deck
with the composition keys (eos_xh, eos_yhe, eos_a_metal) and the table range keys; the
dump of the production binary is what the IC builder interpolates (bicubic spline).
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
PT = os.path.join(REPO, 'data', 'stellar_opac', 'planck_tools')
if PT not in sys.path:
    sys.path.insert(0, PT)
RECIPE_VERSION = 'ext2-v1 (build_ext.py fergR 10-09, generalised 10-15)'

KEV_K = 1.1604518e7
LT0, DLT, NT = 2.6, 0.025, 217
LD0, DLD, ND = -21.0, 0.05, 421
LB0, LB1 = 4.0, 4.2
SLOPE_DEX = 0.6
lT = np.round(LT0 + DLT*np.arange(NT), 6)
lD = np.round(LD0 + DLD*np.arange(ND), 6)

TOPS_HI_ARGS = ('0.0005', '1', '1e-14', '1', '71')
TOPS_LO_ARGS = ('0.0005', '1', '1e-21', '1e-14', '36')


def md5(fn):
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def md5_data(fn):
    """md5 of the non-comment lines (numbers only: header text excluded)"""
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for ln in f:
            if not ln.startswith(b'#'):
                h.update(ln)
    return h.hexdigest()


# --------------------------------------------------------------------------- composition
def mixture_string(comp):
    """TOPS mixture (number fractions) for the composition block"""
    if comp.get('tops_mixture'):
        return comp['tops_mixture'].strip()
    if comp.get('mixture', 'GS98').upper() != 'GS98':
        raise ValueError('mixture must be GS98 or give tops_mixture explicitly')
    from mixspec import mixture
    s, _ = mixture(comp['X'], comp['Z'])
    return s


def comp_hash(comp, mix):
    key = dict(X=round(comp['X'], 10), Z=round(comp['Z'], 10), mix=mix,
               ferg=comp.get('ferguson_xz'), recipe=RECIPE_VERSION)
    if comp.get('saha_mixture'):
        key['saha'] = comp['saha_mixture'].strip()
    return hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]


def eos_hash(comp, eosr):
    key = dict(X=comp['X'], Y=comp['Y'], a_metal=comp.get('a_metal', 16.0), r=eosr)
    return hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------- TOPS raw
def tops_raw(cache, h, mix, spec_src, log):
    """hi / lo raw TOPS files in the cache (copied from spec paths, or fetched)."""
    d = os.path.join(cache, 'raw', 'tops')
    os.makedirs(d, exist_ok=True)
    out = {}
    for part, args in (('hi', TOPS_HI_ARGS), ('lo', TOPS_LO_ARGS)):
        fn = os.path.join(d, '%s_%s.dat' % (h, part))
        src = (spec_src or {}).get(part)
        if not os.path.exists(fn) and src:
            first = open(src).readline()
            if mix.split()[:4] != first.split('mixture:')[1].split()[:4]:
                log('  WARNING: %s mixture header differs from the spec mixture' % src)
            shutil.copyfile(src, fn)
            log('  TOPS %s: copied %s -> %s' % (part, src, fn))
        if not os.path.exists(fn):
            log('  TOPS %s: fetching from LANL (fetch_tops.py %s)' % (part, ' '.join(args)))
            env = dict(os.environ)
            env.setdefault('SSL_CERT_FILE', '/etc/ssl/ca-bundle.pem')
            try:
                subprocess.run([sys.executable, os.path.join(PT, 'fetch_tops.py'), mix] +
                               list(args) + [fn + '.tmp'], check=True, env=env, timeout=1800,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                os.rename(fn + '.tmp', fn)
            except Exception as e:                       # noqa: BLE001
                log('  TOPS %s: FETCH FAILED (%s)' % (part, str(e)[:200]))
        out[part] = fn if os.path.exists(fn) else None
    if out['hi'] is None:
        raise RuntimeError('no TOPS hi file for this mixture (LANL down?): give '
                           'opacity.tops_hi in the spec or retry later')
    return out


# --------------------------------------------------------------------------- ext2 build
class Ext2:
    def __init__(self, ferg_dir, tops, mixfile, X, Z, log=print):
        from .saha_es import Saha
        self.saha = Saha(mixfile)
        self.ferg_dir = ferg_dir
        self.tops = tops
        self.X, self.Z = X, Z
        self.log = log

    def kes_row(self, t, lds):
        lds = np.asarray(lds, float)
        return np.log10(self.saha.kes(t, 10**lds))

    @staticmethod
    def slope(ld_e, xs, ys):
        y0 = np.interp(ld_e, xs, ys)
        y1 = np.interp(ld_e + SLOPE_DEX, xs, ys)
        return (y1 - y0)/SLOPE_DEX

    @staticmethod
    def smooth(lts, s, half=0.1):
        s = np.asarray(s, float)
        out = s.copy()
        for i, t in enumerate(lts):
            m = (abs(lts - t) <= half + 1e-9) & np.isfinite(s)
            if m.any() and np.isfinite(s[i]):
                out[i] = s[m].mean()
        return out

    def extend_row(self, t, xs, kR, kP, ld_e, sPo=None, sRo=None):
        R = np.interp(lD, xs, kR)
        P = np.interp(lD, xs, kP)
        lo = lD < ld_e - 1e-9
        kRe = np.interp(ld_e, xs, kR)
        kPe = np.interp(ld_e, xs, kP)
        sP = float(np.clip(self.slope(ld_e, xs, kP), 0.0, 1.0)) if sPo is None else sPo
        P[lo] = kPe + sP*(lD[lo] - ld_e)
        xe = np.array([ld_e, ld_e + SLOPE_DEX])
        ke = self.kes_row(t, np.concatenate([xe, lD[lo]]))
        kRxe = np.array([kRe, np.interp(ld_e + SLOPE_DEX, xs, kR)])
        abs_x = 10**kRxe - 10**ke[:2]
        if abs_x[0] > 0.0:
            sR = (np.log10(abs_x[1]) - np.log10(abs_x[0]))/SLOPE_DEX if abs_x[1] > 0 else 1.0
            sR = float(np.clip(sR, 0.0, 1.0)) if sRo is None else sRo
            R[lo] = np.log10(10**ke[2:] + abs_x[0]*10**(sR*(lD[lo] - ld_e)))
        else:
            sR = np.nan
            R[lo] = kRe + ke[2:] - ke[0]
        return R, P, sP, sR

    def tops_field(self):
        hi = np.loadtxt(self.tops['hi'])
        if self.tops.get('lo'):
            lo = np.loadtxt(self.tops['lo'])
            d = np.vstack([lo, hi[hi[:, 1] > 1.0001e-14]])
        else:
            d = hi
        T = np.unique(d[:, 0])
        Rh = np.unique(d[:, 1])
        kR = np.full((len(T), len(Rh)), np.nan)
        kP = np.full((len(T), len(Rh)), np.nan)
        ti = {v: i for i, v in enumerate(T)}
        ri = {v: i for i, v in enumerate(Rh)}
        for t, r, a, b in d:
            kR[ti[t], ri[r]] = a
            kP[ti[t], ri[r]] = b
        assert np.isfinite(kR).all()
        lts, lds = np.log10(T*KEV_K), np.log10(Rh)
        j14 = int(np.argmin(abs(lds + 14.0)))
        rowsR, rowsP, info, rows = [], [], [], []
        for i in range(len(T)):
            j = 0
            while j + 1 < len(lds) and kR[i, j+1] == kR[i, 0] and kP[i, j+1] == kP[i, 0]:
                j += 1
            jf = j
            je = min(jf, j14)
            _, _, sP, sR = self.extend_row(T[i]*KEV_K, lds[je:], np.log10(kR[i, je:]),
                                           np.log10(kP[i, je:]), lds[je])
            if jf > j14:
                sP = float(np.clip(self.slope(lds[jf], lds[jf:], np.log10(kP[i, jf:])), 0, 1))
            rows.append((jf, je, sP, sR))
        sPs = self.smooth(lts, [r[2] for r in rows])
        sRs = self.smooth(lts, [r[3] for r in rows])
        for i in range(len(T)):
            jf, je, sP0, sR0 = rows[i]
            R, P, sP, sR = self.extend_row(T[i]*KEV_K, lds[je:], np.log10(kR[i, je:]),
                                           np.log10(kP[i, je:]), lds[je], sPs[i], sRs[i])
            rowsR.append(R)
            rowsP.append(P)
            info.append((lts[i], lds[jf], sR, sP))
        return lts, np.array(rowsR), np.array(rowsP), info

    # ---- Ferguson
    @staticmethod
    def _tag_val(tag, isx):
        if isx and tag == '10':
            return 1.0
        return float('0.' + tag) if tag != '0' else 0.0

    def _ferg_nodes(self, kind):
        d = os.path.join(self.ferg_dir, 'ross' if kind == 'R' else 'planck')
        pre = 'g98.' if kind == 'R' else 'g98.pl.'
        suf = '.tron' if kind == 'R' else '.tpon'
        nodes = {}
        for f in os.listdir(d):
            if f.startswith(pre) and f.endswith(suf):
                xs, zs = f[len(pre):-len(suf)].split('.')
                nodes[(self._tag_val(xs, True), self._tag_val(zs, False))] = os.path.join(d, f)
        return nodes

    def _bracket(self, vals, v, logz=False):
        vals = sorted(set(vals))
        if any(abs(v - x) <= 1e-12 for x in vals):
            x = min(vals, key=lambda x: abs(x - v))
            return x, x, 0.0
        lo = max([x for x in vals if x < v])
        hi = min([x for x in vals if x > v])
        if logz and lo > 0.0:
            w = (np.log10(v) - np.log10(lo))/(np.log10(hi) - np.log10(lo))
        else:
            w = (v - lo)/(hi - lo)
        return lo, hi, w

    def ferg(self, kind):
        from compare import read_ferg
        nodes = self._ferg_nodes(kind)
        nR, nP = self._ferg_nodes('R'), self._ferg_nodes('P')
        common = set(nR) & set(nP)
        x0, x1, wx = self._bracket([k[0] for k in common], self.X)
        zc = [k[1] for k in common if k[0] in (x0, x1)]
        z0, z1, wz = self._bracket(zc, self.Z, logz=True)

        def rd(x, z):
            return read_ferg(nodes[(x, z)])
        t, r, F00 = rd(x0, z0)
        if wx == 0.0 and wz == 0.0:
            F = F00
        elif wz == 0.0:
            t2, r2, F10 = rd(x1, z0)
            assert np.allclose(t, t2) and np.allclose(r, r2)
            F = (1 - wx)*F00 + wx*F10
        elif wx == 0.0:
            t2, r2, F01 = rd(x0, z1)
            assert np.allclose(t, t2) and np.allclose(r, r2)
            F = (1 - wz)*F00 + wz*F01
        else:
            F10, F01, F11 = rd(x1, z0)[2], rd(x0, z1)[2], rd(x1, z1)[2]
            F = (1 - wz)*((1 - wx)*F00 + wx*F10) + wz*((1 - wx)*F01 + wx*F11)
        note = ('Ferguson+2005 GS98 %s: X %g..%g (w %.4f linear), Z %g..%g (w %.4f linear in '
                'log Z)' % ('Rosseland' if kind == 'R' else 'Planck', x0, x1, wx, z0, z1, wz))
        return t, r, F, note

    def lowT_field(self):
        tA, rA, KA, noteR = self.ferg('R')
        tF, rF, KF, noteP = self.ferg('P')

        def rows(ts, rs, K, which):
            sl = []
            for i, t in enumerate(ts):
                xs = rs + 3*t - 18.0
                _, _, sP, sR = self.extend_row(10**t, xs, K[i], K[i], xs[0])
                sl.append(sR if which == 'R' else sP)
            sl = self.smooth(ts, sl)
            out = []
            for i, t in enumerate(ts):
                xs = rs + 3*t - 18.0
                if which == 'R':
                    R, _, _, _ = self.extend_row(10**t, xs, K[i], K[i], xs[0], None, sl[i])
                    out.append(R)
                else:
                    _, P, _, _ = self.extend_row(10**t, xs, K[i], K[i], xs[0], sl[i], None)
                    out.append(P)
            return np.array(out)
        return (tA, rows(tA, rA, KA, 'R'), tF, rows(tF, rF, KF, 'P'), noteR, noteP,
                (rA[0], rA[-1], rF[0], rF[-1], tA[0], tA[-1], tF[0], tF[-1]),
                (tA, rA, KA), (tF, rF, KF))

    @staticmethod
    def tinterp(ts, F, lt):
        out = np.empty((len(lt), F.shape[1]))
        for j in range(F.shape[1]):
            out[:, j] = np.interp(lt, ts, F[:, j])
        return out

    def build(self, outdir, tag, hdr):
        from convert_tops import bilin
        lts, TR, TP, info = self.tops_field()
        KR_t, KP_t = self.tinterp(lts, TR, lT), self.tinterp(lts, TP, lT)
        tA, AR, tF, FP, noteR, noteP, rng, fR, fP = self.lowT_field()
        KR_l, KP_l = self.tinterp(tA, AR, lT), self.tinterp(tF, FP, lT)
        TT, DD = np.meshgrid(lT, lD, indexing='ij')
        RRg = DD - 3*TT + 18.0
        for K_, (ts, rs, F) in ((KR_l, fR), (KP_l, fP)):
            m = (TT >= ts[0]) & (TT <= ts[-1]) & (RRg >= rs[0]) & (RRg <= rs[-1])
            K_[m] = bilin(ts, rs, F, TT, RRg)[m]
        w = np.clip((lT - LB0)/(LB1 - LB0), 0.0, 1.0)[:, None]
        KR = (1 - w)*KR_l + w*KR_t
        KP = (1 - w)*KP_l + w*KP_t
        sR = np.array([['%.5f\n' % v for v in row] for row in KR], dtype=object)
        sP = np.array([['%.5f\n' % v for v in row] for row in KP], dtype=object)
        os.makedirs(outdir, exist_ok=True)
        fr = os.path.join(outdir, 'rosseland_ext2_%s.txt' % tag)
        fp = os.path.join(outdir, 'planck_ext2_%s.txt' % tag)
        common = hdr + [
            'TOPS hi %s md5 %s; lo %s' % (self.tops['hi'], md5(self.tops['hi']),
                                          ('%s md5 %s' % (self.tops['lo'], md5(self.tops['lo'])))
                                          if self.tops.get('lo') else
                                          'NONE (extension from log rho -14 in every TOPS row)'),
            'BLEND linear in log T of log10 kappa over log T %.2f..%.2f (Ferguson below, TOPS '
            'above); LOW-RHO EXTENSION below each source floor: Rosseland = Saha e-scattering + '
            'absorption x (rho/rho_e)^s; Planck = kP_e (rho/rho_e)^s; s = local log slope over '
            '0.6 dex clamped [0,1], boxcar-averaged +-0.1 dex in log T' % (LB0, LB1)]
        write(fr, sR, ['AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on '
                       '(log10 T[K], log10 rho[g/cm3]); INCLUDES electron scattering'] +
              common + ['LOW T: ' + noteR, 'Read with he_opac_table (he_star_m1).'])
        write(fp, sP, ['AthenaK stellar Planck opacity table, log10 kappa_P [cm^2/g] on '
                       '(log10 T[K], log10 rho[g/cm3]); ABSORPTION ONLY (no scattering)'] +
              common + ['LOW T: ' + noteP, 'Read with he_planck_table (he_star_m1).'])
        np.savez(os.path.join(outdir, 'build_fields.npz'), lT=lT, lD=lD, KR=KR, KP=KP,
                 tops_lt=np.array([a for a, b, c, d in info]),
                 tops_floor=np.array([b for a, b, c, d in info]),
                 ferg_lt=np.array([rng[4], rng[5]]), ferg_lr=np.array([rng[0], rng[1]]))
        return fr, fp


def write(fn, S, comments):
    with open(fn, 'w') as f:
        for c in comments:
            f.write('# %s\n' % c)
        f.write('# grid: nT nD lTmin dlT lDmin dlD\n')
        f.write('# %d %d %g %g %g %g\n' % (NT, ND, LT0, DLT, LD0, DLD))
        f.write('# then nT*nD rows, T slowest: log10 kappa [cm^2/g]\n')
        for i in range(NT):
            f.write(''.join(S[i]))


def opacity_tables(spec, cache, log=print):
    """-> dict(ross, planck, hash, mixture, notes); builds into cache/tables/<hash>/ once."""
    comp, op = spec['composition'], spec.get('opacity', {})
    if op.get('rosseland'):                 # explicit table (e.g. a regression)
        return dict(ross=op['rosseland'], planck=op.get('planck'), hash=None,
                    mixture=None, built=False, fields=None,
                    note='explicit tables from the spec (no build)')
    mix = mixture_string(comp)
    h = comp_hash(comp, mix)
    tag = 'x%g_z%g_%s' % (comp['X'], comp['Z'], h)
    outdir = os.path.join(cache, 'tables', 'opac_' + h)
    fr = os.path.join(outdir, 'rosseland_ext2_%s.txt' % tag)
    fp = os.path.join(outdir, 'planck_ext2_%s.txt' % tag)
    if os.path.exists(fr) and os.path.exists(fp):
        log('  opacity: cached %s' % outdir)
    else:
        os.makedirs(outdir, exist_ok=True)
        mixfile = os.path.join(outdir, 'mixture.txt')
        open(mixfile, 'w').write(mix + '\n')
        if comp.get('saha_mixture'):     # the star's own mixture for the Saha-es extension
            mixfile = os.path.join(outdir, 'mixture_saha.txt')   # (when TOPS is interim)
            open(mixfile, 'w').write(comp['saha_mixture'].strip() + '\n')
        tops = tops_raw(cache, h, mix, dict(hi=op.get('tops_hi'), lo=op.get('tops_lo')), log)
        ferg = op.get('ferguson_dir', os.path.join(cache, 'raw', 'ferguson05'))
        log('  opacity: building ext2 for X %g Z %g (hash %s)' % (comp['X'], comp['Z'], h))
        b = Ext2(ferg, tops, mixfile, comp['X'], comp['Z'], log)
        hdr = ['tools/ic_star ext2 table, recipe %s; composition X %g Y %g Z %g, mixture %s '
               '(TOPS number fractions in mixture.txt)' % (RECIPE_VERSION, comp['X'],
                                                           comp['Y'], comp['Z'],
                                                           comp.get('mixture', 'GS98'))]
        fr, fp = b.build(outdir, tag, hdr)
        json.dump(dict(X=comp['X'], Y=comp['Y'], Z=comp['Z'], mixture=mix, tops=tops,
                       ferguson=ferg, recipe=RECIPE_VERSION), open(os.path.join(
                           outdir, 'provenance.json'), 'w'), indent=1)
    return dict(ross=fr, planck=fp, hash=h, mixture=mix, built=True,
                fields=os.path.join(outdir, 'build_fields.npz'),
                note='ext2 recipe, cache %s' % outdir)


# --------------------------------------------------------------------------- EOS dump
EOS_DECK = """<comment>
problem = tools/ic_star general-EOS table dump (X {X} Y {Y} a_metal {am})
<job>
basename = dump
<mesh>
nghost = 2
nx1 = 16
x1min = -0.5
x1max = 0.5
ix1_bc = outflow
ox1_bc = outflow
nx2 = 1
x2min = -0.5
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic
nx3 = 1
x3min = -0.5
x3max = 0.5
ix3_bc = periodic
ox3_bc = periodic
<meshblock>
nx1 = 16
nx2 = 1
nx3 = 1
<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.3
nlim = 0
tlim = 1.0
<units>
length_cgs = 1.0
mass_cgs = 1.0
time_cgs = 1.0
mu = 1.0
<hydro>
eos = general
general_eos = table
gamma = 1.6666666667
eos_h2 = false
eos_ionization = true
eos_metal_ionization = false
eos_radiation = false
eos_xh = {X}
eos_yhe = {Y}
eos_a_metal = {am}
eos_logd_min = {ldmin}
eos_logd_max = {ldmax}
eos_logt_min = {ltmin}
eos_logt_max = {ltmax}
eos_dlog = {dlog}
eos_table_dump = {out}
reconstruct = plm
rsolver = hllc
<problem>
pgen_name = shock_tube
shock_dir = 1
xshock = 0.0
dl = 1.0e-6
pl = 1.0e6
ul = 0.0
vl = 0.0
wl = 0.0
dr = 1.0e-6
pr = 1.0e6
ur = 0.0
vr = 0.0
wr = 0.0
"""


def eos_dump(spec, cache, log=print):
    comp, e = spec['composition'], spec.get('eos', {})
    if e.get('dump'):
        return dict(dump=e['dump'], hash=None, note='explicit dump from the spec')
    rng = dict(ldmin=e.get('logd_min', -21.0), ldmax=e.get('logd_max', -1.0),
               ltmin=e.get('logt_min', 3.4), ltmax=e.get('logt_max', 6.8),
               dlog=e.get('dlog', 0.02))
    h = eos_hash(comp, rng)
    d = os.path.join(cache, 'tables', 'eos_' + h)
    out = os.path.join(d, 'eos_dump.txt')
    if os.path.exists(out):
        log('  EOS: cached %s' % out)
        return dict(dump=out, hash=h, note='cache', deck=os.path.join(d, 'dump.athinput'))
    os.makedirs(d, exist_ok=True)
    deck = os.path.join(d, 'dump.athinput')
    am = comp.get('a_metal', 16.0)
    open(deck, 'w').write(EOS_DECK.format(X=repr(comp['X']), Y=repr(comp['Y']), am=repr(am),
                                          out=out, **{k: repr(v) for k, v in rng.items()}))
    binary = e.get('binary', spec['_cfg']['eos_binary'])
    log('  EOS: dumping with %s (md5 %s)' % (binary, md5(binary)))
    env = dict(os.environ, OMP_NUM_THREADS='1')
    r = subprocess.run([binary, '-i', deck], cwd=d, env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, timeout=1800)
    open(os.path.join(d, 'dump.log'), 'wb').write(r.stdout)
    if not os.path.exists(out):
        raise RuntimeError('EOS dump failed, see %s/dump.log' % d)
    json.dump(dict(binary=binary, binary_md5=md5(binary), comp=comp, range=rng),
              open(os.path.join(d, 'provenance.json'), 'w'), indent=1, default=str)
    return dict(dump=out, hash=h, note='built', deck=deck)
