"""MESA-consistent low-T Rosseland opacity (tools/ic_star, 10-18; spec opacity.mesa_kap).

Reproduces MESA's own radiative opacity for a star built with MESA r15140 kap defaults
(kap/defaults/kap.defaults, md5 e0f6eacc: kap_file_prefix 'gs98', kap_lowT_prefix 'lowT_fa05_gs98',
kap_blend_logT_lower/upper_bdy 3.80/3.88, linear interpolation in X and Z) from the kap_data
MESA ships (kap/kap_data.tar.xz, git-LFS sha256 95ca3578...), and splices it into the recipe
table below a join in log T:
  log kappa = MESA(T, rho)                         for log T <= join_lo,
            = linear in log T to the recipe table  for join_lo < log T < join_hi,
            = the recipe table (TOPS)              above join_hi;
outside MESA's log R = log rho - 3 log T + 18 range ([-8, 1]) the recipe table is kept, blended
over 0.5 dex of log R at the low-R edge.  Type2 (C/O-enhanced) tables are NOT used: for the He
giant X 9.5e-4 (Type2 full off at X >= 1e-3) and dZ = Z - Zbase ~0 make MESA's frac_Type2 ~0;
specify only for stars where that holds (the build logs the Type2 check inputs).
Within MESA's own track this reproduces the MESA profile's opacity column to 4-10 %
(RESULTS sect. 11)."""
import hashlib
import json
import os
import tarfile

import numpy as np
from scipy.interpolate import RectBivariateSpline

XS = [0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98]
VERSION = 'mesakap-v1 (r15140 defaults: gs98 + lowT_fa05_gs98, blend 3.80-3.88 cos, lin X Z)'


def md5(fn):
    h = hashlib.md5()
    with open(fn, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def kap_dir(tarball, cache):
    """extract MESA's kap_data once, keyed by the tarball md5"""
    h = md5(tarball)
    d = os.path.join(cache, 'raw', 'mesa_kap', 'kap_data_' + h[:12])
    if not os.path.isdir(os.path.join(d, 'kap_data')):
        os.makedirs(d, exist_ok=True)
        with tarfile.open(tarball) as t:
            t.extractall(d)
    return os.path.join(d, 'kap_data'), h


def _read(fn):
    L = open(fn).read().splitlines()
    h = L[2].split()
    lR = np.array([float(v) for v in L[5].split()])
    rows = [ln.split() for ln in L[7:] if ln.strip()]
    lT = np.array([float(r[0]) for r in rows])
    v = np.array([[float(x) for x in r[1:]] for r in rows])
    return float(h[2]), float(h[3]), lT, lR, v


class _Tab:
    def __init__(self, fn):
        self.X, self.Z, self.lT, self.lR, v = _read(fn)
        self.sp = RectBivariateSpline(self.lT, self.lR, v, kx=3, ky=3)

    def __call__(self, lT, lR):
        return self.sp.ev(np.clip(lT, self.lT[0], self.lT[-1]),
                          np.clip(lR, self.lR[0], self.lR[-1]))


def _fname(kd, prefix, z, x):
    xs = '0.0' if x == 0.0 else '%g' % x
    for zz in ('%g' % z if z else '0.0', '%.4f' % z, '%.3f' % z, '%.2f' % z):
        f = os.path.join(kd, '%s_z%s_x%s.data' % (prefix, zz, xs))
        if os.path.exists(f):
            return f
    raise FileNotFoundError('MESA kap table %s z %g x %g' % (prefix, z, x))


def _zs(kd, prefix):
    out = set()
    for f in os.listdir(kd):
        if f.startswith(prefix + '_z') and f.endswith('_x0.0.data'):
            out.add(float(f[len(prefix) + 2:].split('_x')[0]))
    return sorted(out)


class MesaKap:
    def __init__(self, kd, X, Z, hi='gs98', lo='lowT_fa05_gs98', blo=3.80, bhi=3.88):
        self.blo, self.bhi, self.hi, self.lo = blo, bhi, hi, lo
        self.parts, self.files = {}, []
        for pre in (hi, lo):
            ix = max(0, min(int(np.searchsorted(XS, X)) - 1, len(XS) - 2))
            zs = _zs(kd, pre)
            iz = max(0, min(int(np.searchsorted(zs, Z)) - 1, len(zs) - 2))
            x0, x1, z0, z1 = XS[ix], XS[ix+1], zs[iz], zs[iz+1]
            ax, az = (X - x0)/(x1 - x0), (Z - z0)/(z1 - z0)
            ws = []
            for z, wz in ((z0, 1 - az), (z1, az)):
                for x, wx in ((x0, 1 - ax), (x1, ax)):
                    fn = _fname(kd, pre, z, x)
                    self.files.append(os.path.basename(fn))
                    ws.append((wx*wz, _Tab(fn)))
            self.parts[pre] = ws
        t = self.parts[lo][0][1]
        self.lr_lo = max(self.parts[hi][0][1].lR[0], t.lR[0])
        self.lr_hi = min(self.parts[hi][0][1].lR[-1], t.lR[-1])

    def logkap(self, lT, lD):
        lR = lD - 3*lT + 18
        kh = sum(w*t(lT, lR) for w, t in self.parts[self.hi])
        kl = sum(w*t(lT, lR) for w, t in self.parts[self.lo])
        a = np.clip((lT - self.blo)/(self.bhi - self.blo), 0, 1)
        a = 0.5*(1 - np.cos(np.pi*a))
        return a*kh + (1 - a)*kl, lR


def splice(base_fn, out_fn, kd, mk, log, hdr):
    """the recipe table base_fn with MESA's kappa spliced in below mk['join']"""
    L = open(base_fn).read().splitlines()
    com = [ln for ln in L if ln.startswith('#')]
    g = [ln for ln in com if len(ln.split()) == 7 and ln.split()[1].isdigit()][0].split()[1:]
    nT, nD = int(g[0]), int(g[1])
    lT0, dlT, lD0, dlD = [float(x) for x in g[2:]]
    v = np.array([float(x) for x in L if x.strip() and not x.startswith('#')]).reshape(nT, nD)
    T = lT0 + dlT*np.arange(nT)
    D = lD0 + dlD*np.arange(nD)
    M = MesaKap(kd, mk['X'], mk['Z'], mk.get('hi', 'gs98'), mk.get('lo', 'lowT_fa05_gs98'),
                mk.get('blend_lo', 3.80), mk.get('blend_hi', 3.88))
    j0, j1 = mk.get('join', [4.0, 4.2])
    LT, LD = np.meshgrid(T, D, indexing='ij')
    km, lR = M.logkap(LT, LD)
    wT = np.clip((j1 - LT)/(j1 - j0), 0, 1)              # 1 below join_lo, 0 above join_hi
    wR = np.clip((lR - M.lr_lo)/0.5, 0, 1)*(lR <= M.lr_hi)
    w = wT*wR
    out = w*km + (1 - w)*v
    with open(out_fn, 'w') as f:
        for c in hdr:
            f.write('# %s\n' % c)
        for c in com:
            if 'grid' in c or (len(c.split()) == 7 and c.split()[1].isdigit()) or 'then nT' in c:
                f.write(c + '\n')
        for x in out.ravel():
            f.write('%.5f\n' % x)
    d = out - v
    log('  mesa_kap: spliced MESA kappa (%s) below log T %.3g (blend to %.3g), log R >= %.2f '
        '(0.5 dex edge); max |dlog kappa| vs the recipe %.3f, %d of %d nodes changed'
        % (', '.join(M.files), j0, j1, M.lr_lo, np.abs(d).max(), int((np.abs(d) > 1e-6).sum()),
           d.size))
    return out_fn


def derive(base, spec, cache, log):
    """opacity.mesa_kap -> the spliced Rosseland table (cache keyed by base md5, kap_data md5,
    X, Z, join, blend, prefixes)"""
    mk = dict(spec['opacity']['mesa_kap'])
    comp = spec['composition']
    mk.setdefault('X', comp['X'])
    mk.setdefault('Z', comp['Z'])
    kd, kh = kap_dir(mk['kap_data'], cache)
    key = dict(base=md5(base['ross']), kap_data=kh, X=mk['X'], Z=mk['Z'],
               join=mk.get('join', [4.0, 4.2]), blend=[mk.get('blend_lo', 3.80),
                                                      mk.get('blend_hi', 3.88)],
               hi=mk.get('hi', 'gs98'), lo=mk.get('lo', 'lowT_fa05_gs98'), version=VERSION)
    h = hashlib.sha1(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]
    od = os.path.join(cache, 'tables', 'opac_mk_' + h)
    fn = os.path.join(od, 'rosseland_mesakap_%s.txt' % h)
    if os.path.exists(fn):
        log('  mesa_kap: cached %s' % fn)
    else:
        os.makedirs(od, exist_ok=True)
        splice(base['ross'], fn, kd, mk, log,
               ['AthenaK stellar Rosseland opacity table, log10 kappa_R [cm^2/g] on (log10 T[K], '
                'log10 rho[g/cm3]); INCLUDES electron scattering',
                'tools/ic_star %s; MESA kap_data %s (md5 %s); X %g Z %g; join log T %s; base '
                '%s (md5 %s)' % (VERSION, mk['kap_data'], kh, mk['X'], mk['Z'],
                                 key['join'], base['ross'], key['base']),
                'Read with he_opac_table (he_star_m1).'])
        json.dump(key, open(os.path.join(od, 'provenance.json'), 'w'), indent=1)
    return dict(base, ross=fn, mesa_kap=dict(key, hash=h),
                note=base.get('note', '') + '; mesa_kap splice %s' % od)
