"""INTERIM APPROXIMATION (not the paper's observable): Fe I opacity-weighted terminator line-of-sight
velocity from the dhj p1 hydro dumps, for comparison with Seidel et al. 2026 (Fig. 1 net Fe I blueshift).
No line RT, no CCF, no transit-phase geometry. Reads dumps read-only. + = toward the transit observer (blueshift).
Weight per cell (terminator annulus |cos(angle from substellar)| < 0.1, cell volume dV):
  w = n_FeI * exp(-E_eff/kT) / U_FeI(T) * dV,  n_FeI = A_Fe n_H * f_gas(T,p) / (1 + S_Fe(T)/n_e),
  n_e from a Saha charge balance (Na K Ca Al Mg Fe Si H, thermal only), Fe gas cut below the Fe condensation curve
  (Visscher+2010: 1e4/T = 5.44 - 0.48 log p[bar]). E_eff = 1.5 eV. Velocity = wind (rotating frame, incl. v_r)
  + synchronous solid-body rotation Omega r cos(lat)."""
import sys, os, glob, re, json, numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/rochedhj_ana'); import dhjcs
from multiprocessing import Pool
J = '/resnick/groups/carnegie_poc/jingze'; R1 = J + '/r1planets_1101'
EOSF = J + '/w121r1_1101/anchor/eos_1x.txt'
RUNS = [('WASP-76b', R1 + '/wasp76b/runs/p1', 2229, 6.95, 0.36),
        ('WASP-121b', J + '/w121r1_1101/runs/p1', 2358, 4.75, 0.30),
        ('WASP-178b', R1 + '/wasp178b/runs/p1', 2398, 4.79, 1.19),
        ('TOI-1518b', R1 + '/toi1518b/runs/p1', 2491, 3.82, 0.50),
        ('WASP-189b', R1 + '/wasp189b/runs/p1', 2638, 1.05, 1.59)]
NROT = 5.0; XH = 0.7381; MH = 1.6726e-24; KEV = 8.617333e-5
ELEM = [  # A (by number rel. H, Asplund+09), chi [eV], g = 2 U+/U0
    (1.74e-6, 5.139, 1.0), (1.07e-7, 4.341, 1.0), (2.19e-6, 6.113, 4.0), (2.82e-6, 5.986, 0.33),
    (3.98e-5, 7.646, 4.0), (3.16e-5, 7.902, 3.2), (3.24e-5, 8.152, 1.2), (1.0, 13.598, 1.0)]
AFE = 3.16e-5; EEFF = 1.5

def saha_S(T, chi, g): return 2.415e15 * T ** 1.5 * g * np.exp(-chi / (KEV * T))

def n_e(nH, T):
    S = [(a, saha_S(T, c, g)) for a, c, g in ELEM]
    lo = np.full(nH.shape, -20.0); hi = np.log10(nH * 1.2 + 1.0)
    for _ in range(60):
        mid = 0.5 * (lo + hi); ne = 10 ** mid
        f = sum(a * nH * s / (s + ne) for a, s in S) - ne
        lo = np.where(f > 0, mid, lo); hi = np.where(f > 0, hi, mid)
    return 10 ** (0.5 * (lo + hi))

def deck(fn):
    d = {}
    for l in open(fn):
        m = re.match(r'^\s*(\w+)\s*=\s*([^#\s]+)', l)
        if m and m.group(1) not in d: d[m.group(1)] = m.group(2)
    return d

_G = {}
def geom(raw):
    key = (len(raw['mb_geometry']), raw['nx2_mb'], raw['nx3_mb'])
    if key in _G: return _G[key]
    geo = dhjcs.cs_geometry(raw); g = np.asarray(raw['mb_geometry']); nmb = g.shape[0]; nj = raw['nx2_mb']; nk = raw['nx3_mb']
    A = np.zeros((nmb, nk, nj)); XYZ = np.zeros((nmb, nk, nj, 3))
    for m in range(nmb):
        xi = np.pi / 4 * (g[m, 2] + (g[m, 3] - g[m, 2]) * (np.arange(nj) + 0.5) / nj)
        eta = np.pi / 4 * (g[m, 4] + (g[m, 5] - g[m, 4]) * (np.arange(nk) + 0.5) / nk)
        XI, ETA = np.meshgrid(xi, eta); x = np.tan(XI); y = np.tan(ETA)
        A[m] = (1 + x * x) * (1 + y * y) / (1 + x * x + y * y) ** 1.5 * (np.pi / 4) ** 2 * (g[m, 3] - g[m, 2]) * (g[m, 5] - g[m, 4]) / (nj * nk)
        XYZ[m] = dhjcs.panel_map(m // (nmb // 6), XI, ETA)
    A /= A.sum(); X, Y, Z = XYZ[..., 0], XYZ[..., 1], XYZ[..., 2]; rh = np.maximum(np.sqrt(X * X + Y * Y), 1e-12)
    G = dict(geo=geo, A=A, X=X, Y=Y, zx=-Y / rh, mx=-Z * X / rh, coslat=np.cos(geo['lat']))
    _G[key] = G; return G

_EOS = None
def reduce_dump(args):
    path, om = args
    global _EOS
    if _EOS is None: _EOS = dhjcs.EOS(EOSF)
    raw = dhjcs.bin_convert.read_binary(path); mb = raw['mb_data']; G = geom(raw)
    rho = np.asarray(mb['dens']); T, p = _EOS.invert(rho, np.asarray(mb['eint'])); lp = np.log10(p) - 6.0   # bar
    r = np.asarray(raw['mb_x1v'])[:, None, None, :]; dr = np.gradient(np.asarray(raw['mb_x1v']), axis=-1)[:, None, None, :]
    u, v = dhjcs.winds(raw, G['geo']); vr = np.asarray(mb['velx'])
    urot = om * r * G['coslat'][..., None]
    X = G['X'][..., None]; zx = G['zx'][..., None]; mx = G['mx'][..., None]
    vw = vr * X + u * zx + v * mx; vrot = urot * zx; vx = (vw + vrot) / 1e5      # km/s, + = toward observer
    term = (np.abs(G['X']) < 0.1)[..., None] & np.ones_like(rho, bool)
    mo = term & (G['Y'] > 0)[..., None]; ev = term & (G['Y'] < 0)[..., None]
    nH = XH * rho / MH; ok = np.isfinite(T) & np.isfinite(p) & (T > 300)
    Tc = np.where(ok, T, 1000.0); ne = n_e(np.where(ok, nH, 1.0), Tc)
    fI = 1.0 / (1.0 + saha_S(Tc, 7.902, 3.2) / ne)
    Tcond = 1e4 / (5.44 - 0.48 * np.where(ok, lp, 0.0)); fgas = (Tc > Tcond).astype(float)
    U = 25.0 * (Tc / 2000.0) ** 0.45
    kap = AFE * nH * fgas * fI * np.exp(-EEFF / (KEV * Tc)) / U
    dV = G['A'][..., None] * r * r * dr
    out = {'t': float(raw['time'])}
    def wmean(f, w, m):
        w = np.where(m & ok & np.isfinite(f), w, 0.0); s = w.sum(); return float((f * w).sum() / s) if s > 0 else np.nan
    win = {'w43': (lp >= -4) & (lp <= -3), 'w52': (lp >= -5) & (lp <= -2), 'w31': (lp >= -3) & (lp <= -1)}
    W = kap * dV
    for k, msk in win.items():
        out['v_' + k] = wmean(vx, W, term & msk); out['vmo_' + k] = wmean(vx, W, mo & msk); out['vev_' + k] = wmean(vx, W, ev & msk)
        out['vwind_' + k] = wmean(vw / 1e5, W, term & msk)
        out['T_' + k] = wmean(T, W, term & msk)
        out['frot_mo_' + k] = wmean(vrot / 1e5, W, mo & msk)
    # density-weighted, no Fe factor, 1e-4..1e-3
    out['vrho_w43'] = wmean(vx, rho * dV, term & win['w43'])
    # Fe-weighted log p median over the whole terminator column set (formation-pressure proxy, no window)
    w = np.where(term & ok, W, 0.0).ravel(); o = np.argsort(lp.ravel()); c = np.cumsum(w[o]); c /= c[-1]
    out['lp16'], out['lp50'], out['lp84'] = [float(lp.ravel()[o][np.searchsorted(c, q)]) for q in (0.16, 0.5, 0.84)]
    out['v_all'] = wmean(vx, W, term)
    # crude page number for reference: unweighted area mean at 1e-3 bar, wind only
    vl = dhjcs.level_interp(vw / 1e5, lp + 6.0, [3.0])[0]; m = np.abs(G['X']) < 0.1
    out['v_crude'] = float(np.nansum(np.where(m & np.isfinite(vl), vl * G['A'], 0)) / np.sum(np.where(m & np.isfinite(vl), G['A'], 0)))
    return path, out

def main():
    tasks = []; meta = {}
    for name, run, teq, vs, es in RUNS:
        d = deck(run + '/run.athinput'); om = float(d['omega']); P = 2 * np.pi / om
        bins = sorted(glob.glob(run + '/bin/*.bin'))
        import time; bins = [b for b in bins if time.time() - os.path.getmtime(b) > 120]
        nb = int(round(NROT / (float(55077.0) / 55077.0 * 0.5))) + 1   # dumps every 0.5 rot -> 11
        sel = bins[-nb:]; meta[name] = dict(P=P, teq=teq, vs=vs, es=es, files=sel, teq_deck=float(d['Teq']))
        tasks += [(b, om) for b in sel]
    with Pool(int(os.environ.get('NPROC', '12'))) as pool: res = dict(pool.map(reduce_dump, tasks, chunksize=1))
    tab = {}
    for name, M in meta.items():
        rs = [res[f] for f in M['files']]
        keys = [k for k in rs[0] if k != 't']
        A = {k: np.array([x[k] for x in rs]) for k in keys}
        rots = np.array([x['t'] for x in rs]) / M['P']
        tab[name] = dict(M, files=[os.path.basename(f) for f in M['files']], rot0=float(rots.min()), rot1=float(rots.max()),
                         mean={k: float(np.nanmean(a)) for k, a in A.items()}, sd={k: float(np.nanstd(a)) for k, a in A.items()},
                         per=[dict(x, rot=x['t'] / M['P']) for x in rs])
    json.dump(tab, open(os.path.dirname(os.path.abspath(__file__)) + '/fewind.json', 'w'), indent=1)
    for name, M in tab.items():
        m, s = M['mean'], M['sd']
        print('%-10s rot %.1f-%.1f  v43 %+.2f+-%.2f  mo %+.2f ev %+.2f  wind-only %+.2f  v52 %+.2f v31 %+.2f vrho %+.2f vall %+.2f crude %+.2f  T43 %.0f  lp %.2f/%.2f/%.2f  Seidel %.2f+-%.2f'
              % (name, M['rot0'], M['rot1'], m['v_w43'], s['v_w43'], m['vmo_w43'], m['vev_w43'], m['vwind_w43'], m['v_w52'], m['v_w31'],
                 m['vrho_w43'], m['v_all'], m['v_crude'], m['T_w43'], m['lp16'], m['lp50'], m['lp84'], M['vs'], M['es']))

if __name__ == '__main__': main()
