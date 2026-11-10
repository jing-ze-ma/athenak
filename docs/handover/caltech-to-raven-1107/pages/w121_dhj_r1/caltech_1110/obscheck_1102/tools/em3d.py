"""em3d.py: 3-D thermal-emission ray tracer (GPU, CuPy RawKernel; CPU fallback = numba, slow) for
AthenaK cubed-sphere dhj dumps -> planet flux vs orbital phase and wavelength.

Geometry (fe3d.py convention, substellar point = -x, rotation about +z): the unit vector toward
the observer at orbital phase phi (0 = mid-transit, 0.5 = mid-eclipse) is
    n_obs = (cos 2 pi phi, -sin 2 pi phi, 0)   [planet frame, tidally locked]
At phi < 0.5 the observer sees the hemisphere east (downstream of the rotation) of the
substellar point, so an eastward hotspot gives a peak before eclipse.

Rays: parallel to n_obs on an (impact parameter b, position angle theta) grid over the projected
disc out to the top of the grid; each ray is cut at every radial face crossing (both sides of
the tangent point) and every segment split into NSUB sub-segments; T, p, rho are interpolated
linearly in radius between cell centres of the NEAREST column (as fe3d).  Rays below the inner
boundary end on it with I += B(T_bottom) e^-tau (it is at ~400 bar: never visible).
Formal solution per (ray, wavelength bin, g-point): I = sum B(T_s) e^-tau_s (1 - e^-dtau_s),
no scattering (thermal emission only; Rayleigh/e- scattering left out, no clouds, no reflected
light).  Opacity: premixed correlated-k table ln kappa[T, p, w, g] (cm^2/g, from opac.py).

  python em3d.py --planet P.json --eos EOS --table tables/opac_R200.npz --out OUT dump.bin [key=val]
keys: nph (40), nb (96), nth (64), nsub (3), tag
Output OUT/em_<tag>.npz: phases, lam_um, Fp[phase, w] (erg/s/cm/sr-integrated = sum I dA, i.e.
planet flux x d^2) and Iray maps at phi 0 and 0.5.
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, '/resnick/groups/carnegie_poc/jingze/fewind_1102/fehires_fast/code')
import csio  # noqa: E402

H, C, KB = 6.62607015e-27, 2.99792458e10, 1.380649e-16

KSRC = r'''
extern "C" __global__ void trace(const float* lk, const int* iT, const float* fT, const int* iP,
    const float* fP, const float* rds, const float* Ts, const float* Tbot, const double* area,
    const double* lam, const double* wg, int nray, int ns, int nP, int NW, int NG,
    double* F, double* Iw, int storeI)
{
  long t = (long)blockIdx.x*blockDim.x + threadIdx.x;
  int NWG = NW*NG;
  if (t >= (long)nray*NWG) return;
  int ray = t / NWG; int q = t % NWG; int w = q / NG; int g = q % NG;
  double L = lam[w];
  double c1 = 2.0*6.62607015e-27*2.99792458e10*2.99792458e10/(L*L*L*L*L);
  double c2 = 6.62607015e-27*2.99792458e10/(L*1.380649e-16);
  double tau = 0.0, I = 0.0;
  long o0 = (long)ray*ns;
  for (int s = 0; s < ns; s++) {
    long o = o0 + s;
    float rd = rds[o];
    if (rd <= 0.f) continue;
    int a = iT[o], b = iP[o];
    float ft = fT[o], fp = fP[o];
    long base = ((long)a*nP + b)*NWG + q;
    float l = (1.f-ft)*((1.f-fp)*lk[base] + fp*lk[base+NWG])
            + ft*((1.f-fp)*lk[base+(long)nP*NWG] + fp*lk[base+(long)(nP+1)*NWG]);
    double dt = exp((double)l)*rd;
    double B = c1/expm1(c2/Ts[o]);
    I += B*exp(-tau)*(-expm1(-dt));
    tau += dt;
    if (tau > 50.0) break;
  }
  if (tau <= 50.0 && Tbot[ray] > 0.f) I += c1/expm1(c2/Tbot[ray])*exp(-tau);
  atomicAdd(&F[w], wg[g]*area[ray]*I);
  if (storeI) atomicAdd(&Iw[(long)ray*NW + w], wg[g]*I);
}
'''


def _mk_cpu():
    from numba import njit, prange
    import math

    @njit(parallel=True, cache=True)
    def tr(lk, iT, fT, iP, fP, rds, Ts, Tbot, area, lam, wg, nP, NW, NG, F, Iw, storeI):
        nray, ns = rds.shape
        NWG = NW*NG
        Fr = np.zeros((nray, NW))
        for ray in prange(nray):
            for q in range(NWG):
                w = q//NG; g = q % NG
                L = lam[w]
                c1 = 2.0*H*C*C/L**5
                c2 = H*C/(L*KB)
                tau = 0.0; I = 0.0
                for s in range(ns):
                    rd = rds[ray, s]
                    if rd <= 0: continue
                    a = iT[ray, s]; b = iP[ray, s]; ft = fT[ray, s]; fp = fP[ray, s]
                    l = (1-ft)*((1-fp)*lk[a, b, q] + fp*lk[a, b+1, q]) + ft*((1-fp)*lk[a+1, b, q] + fp*lk[a+1, b+1, q])
                    dt = math.exp(l)*rd
                    I += c1/math.expm1(c2/Ts[ray, s])*math.exp(-tau)*(-math.expm1(-dt))
                    tau += dt
                    if tau > 50.0: break
                if tau <= 50.0 and Tbot[ray] > 0:
                    I += c1/math.expm1(c2/Tbot[ray])*math.exp(-tau)
                Fr[ray, w] += wg[g]*I
        for ray in range(nray):
            for w in range(NW):
                F[w] += area[ray]*Fr[ray, w]
                if storeI: Iw[ray, w] += Fr[ray, w]
    return tr


_TR = []


def _trace_cpu(*a):
    if not _TR:
        _TR.append(_mk_cpu())
    _TR[0](*a)


def load_model(dump, eosfile):
    raw = csio.read_dump(dump)
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], float)
    T, p = csio.EOS(eosfile).invert(rho, np.asarray(mb['eint'], float))
    rf, rc = csio.radii(raw)
    X, _, _, _ = csio.geometry(raw)
    nr = rho.shape[-1]
    ncol = X[..., 0].size
    return dict(rho=rho.reshape(ncol, nr), T=T.reshape(ncol, nr), lp=np.log10(p.reshape(ncol, nr)/1e6),
                rf=rf, rc=rc, X=X.reshape(ncol, 3), time=raw['time'])


def colmap(X, nlat=720, nlon=1440):
    from scipy.spatial import cKDTree
    lat = (np.arange(nlat) + 0.5)/nlat*np.pi - np.pi/2
    lon = (np.arange(nlon) + 0.5)/nlon*2*np.pi - np.pi
    LA, LO = np.meshgrid(lat, lon, indexing='ij')
    d = np.stack([np.cos(LA)*np.cos(LO), np.cos(LA)*np.sin(LO), np.sin(LA)], -1).reshape(-1, 3)
    _, i = cKDTree(X).query(d)
    return i.reshape(nlat, nlon).astype(np.int32)


def ray_grid(rf, nb, nth):
    """impact parameters: 40 % of nodes uniform in [0, r_in], 45 % in [r_in, r_in + 0.25 (r_top-r_in)]
    (the photosphere/limb), 15 % above; midpoint rule, weights b db dtheta."""
    r0, r1 = rf[0], rf[-1]
    n1, n2 = int(0.4*nb), int(0.45*nb)
    n3 = nb - n1 - n2
    e = np.concatenate([np.linspace(0, r0, n1 + 1)[:-1], np.linspace(r0, r0 + 0.25*(r1 - r0), n2 + 1)[:-1],
                        np.linspace(r0 + 0.25*(r1 - r0), r1, n3 + 1)])
    b = 0.5*(e[1:] + e[:-1])
    db = np.diff(e)
    th = (np.arange(nth) + 0.5)/nth*2*np.pi
    B, TH = np.meshgrid(b, th, indexing='ij')
    area = (B*np.repeat(db[:, None], nth, 1)*(2*np.pi/nth)).ravel()
    return B.ravel(), TH.ravel(), area


class Tracer:
    def __init__(self, model, table, nb=96, nth=64, nsub=3, xp=None):
        if os.environ.get('EM3D_CPU'):
            cp = np
        else:
            import cupy as cp
        self.cp = cp
        self.m = model
        tb = np.load(table)
        self.lam = tb['lam_um']*1e-4
        self.wg = tb['wg'].astype(np.float64)
        self.Tg = tb['T'].astype(np.float64)
        self.lpg = np.log10(tb['p_bar']).astype(np.float64)
        lk = tb['lnk'].astype(np.float32)          # (nT, nP, NW, NG)
        self.nT, self.nP, self.NW, self.NG = lk.shape
        self.lk = cp.asarray(lk.reshape(self.nT, self.nP, -1))
        self.b, self.th, self.area = ray_grid(model['rf'], nb, nth)
        self.nsub = nsub
        self.cmap = cp.asarray(colmap(model['X']))
        self.gT = cp.asarray(model['T'], dtype=cp.float64)
        self.glp = cp.asarray(model['lp'], dtype=cp.float64)
        self.grho = cp.asarray(model['rho'], dtype=cp.float64)
        self.kern = None if cp is np else cp.RawKernel(KSRC, 'trace')

    def samples(self, phi):
        """per (ray, sample): table indices/weights, rho ds, T; ordered from the observer inward."""
        cp = self.cp
        m = self.m
        rf = cp.asarray(m['rf'])
        rc = cp.asarray(m['rc'])
        nr = len(m['rc'])
        a = 2*np.pi*phi
        nob = np.array([np.cos(a), -np.sin(a), 0.0])
        e1 = np.array([np.sin(a), np.cos(a), 0.0])
        e2 = np.array([0.0, 0.0, 1.0])
        b = cp.asarray(self.b)[:, None]
        th = cp.asarray(self.th)[:, None]
        # face crossings s_k = sqrt(rf_k^2 - b^2) (0 where rf_k < b)
        sk = cp.sqrt(cp.maximum(rf[None, :]**2 - b**2, 0.0))            # (nray, nr+1)
        # near side segments from the top: k = nr-1 .. 0, s from s_{k+1} down to s_k
        sub = (cp.arange(self.nsub) + 0.5)/self.nsub
        hi = sk[:, 1:][:, ::-1]
        lo = sk[:, :-1][:, ::-1]                                         # (nray, nr)
        s_near = hi[:, :, None] - (hi - lo)[:, :, None]*sub[None, None, :]
        ds_near = cp.repeat(((hi - lo)/self.nsub)[:, :, None], self.nsub, 2)
        hit = (b[:, 0] < m['rf'][0])
        # far side: k = 0 .. nr-1, s from -s_k to -s_{k+1}; absent for rays hitting the inner boundary
        lo2 = sk[:, :-1]
        hi2 = sk[:, 1:]
        s_far = -(lo2[:, :, None] + (hi2 - lo2)[:, :, None]*sub[None, None, :])
        ds_far = cp.repeat(((hi2 - lo2)/self.nsub)[:, :, None], self.nsub, 2)
        ds_far = ds_far*(~hit)[:, None, None]
        s = cp.concatenate([s_near.reshape(len(self.b), -1), s_far.reshape(len(self.b), -1)], 1)
        ds = cp.concatenate([ds_near.reshape(len(self.b), -1), ds_far.reshape(len(self.b), -1)], 1)
        P = (b*cp.cos(th))[..., None]*cp.asarray(e1) + (b*cp.sin(th))[..., None]*cp.asarray(e2) \
            + s[..., None]*cp.asarray(nob)
        r = cp.sqrt((P**2).sum(-1))
        lat = cp.arcsin(cp.clip(P[..., 2]/cp.maximum(r, 1.0), -1, 1))
        lon = cp.arctan2(P[..., 1], P[..., 0])
        nlat, nlon = self.cmap.shape
        il = cp.clip(((lat + np.pi/2)/np.pi*nlat).astype(cp.int32), 0, nlat - 1)
        io = cp.clip(((lon + np.pi)/(2*np.pi)*nlon).astype(cp.int32), 0, nlon - 1)
        col = self.cmap[il, io]
        fi = cp.interp(r.ravel(), rc, cp.arange(nr, dtype=cp.float64)).reshape(r.shape)
        i0 = cp.clip(cp.floor(fi).astype(cp.int32), 0, nr - 2)
        f = fi - i0
        def ip(A):
            return A[col, i0]*(1 - f) + A[col, i0 + 1]*f
        T = ip(self.gT)
        lp = ip(self.glp)
        rho = cp.exp(ip(cp.log(self.grho)))
        Tg = cp.asarray(self.Tg)
        jT = cp.clip(cp.searchsorted(Tg, T) - 1, 0, self.nT - 2)
        fT = cp.clip((T - Tg[jT])/(Tg[jT + 1] - Tg[jT]), 0, 1)
        lpg = self.lpg
        x = (lp - lpg[0])/(lpg[1] - lpg[0])
        jP = cp.clip(cp.floor(x).astype(cp.int32), 0, self.nP - 2)
        fP = cp.clip(x - jP, 0, 1)
        rds = rho*ds
        self._colir = (col, cp.where(f > 0.5, i0 + 1, i0).astype(cp.int32))
        Tbot = cp.where(hit, self.gT[col[:, self.nsub*nr - 1], 0], 0.0)   # last near-side sample
        return (jT.astype(cp.int32), fT.astype(cp.float32), jP.astype(cp.int32), fP.astype(cp.float32),
                rds.astype(cp.float32), T.astype(cp.float32), Tbot.astype(cp.float32))

    def flux(self, phi, storeI=False):
        cp = self.cp
        jT, fT, jP, fP, rds, T, Tbot = self.samples(phi)
        nray, ns = rds.shape
        F = cp.zeros(self.NW, cp.float64)
        Iw = cp.zeros((nray, self.NW) if storeI else (1, 1), cp.float64)
        n = nray*self.NW*self.NG
        bs = 256
        if cp is np:
            _trace_cpu(self.lk, jT, fT, jP, fP, rds, T, Tbot, self.area, self.lam, self.wg, self.nP, self.NW,
                       self.NG, F, Iw, storeI)
            return F, (Iw if storeI else None)
        self.kern(((n + bs - 1)//bs,), (bs,),
                  (self.lk, jT, fT, jP, fP, rds, T, Tbot, cp.asarray(self.area), cp.asarray(self.lam),
                   cp.asarray(self.wg), np.int32(nray), np.int32(ns), np.int32(self.nP),
                   np.int32(self.NW), np.int32(self.NG), F, Iw, np.int32(1 if storeI else 0)))
        return cp.asnumpy(F), (cp.asnumpy(Iw) if storeI else None)


# ======================================================================================= clouds
# Post-processed condensate clouds (no radiative feedback), the Parmentier+2016 / Roman & Rauscher
# recipe: a species condenses in every cell whose T is below its condensation curve T_c(p)
# (Visscher+2010 / Wakeford+2017 fits, [M/H] = 0); in each column the cloud base is the deepest
# condensing cell and the cloud is kept from the base up to EXT decades of pressure above it
# (EXT = inf: every condensing cell).  Condensate mass = the full condensable mass of the limiting
# element (1x solar, Asplund 2009), uniform in the cloud.  Opacity per gram of condensate from the
# pRT 3.4 Mie tables (crystalline, monodisperse radius R): 'abs' = absorption only (lower bound:
# the thermal-emission tracer has no scattering), 'ext' = absorption + (1 - g) scattering
# treated as absorption (delta-scaled extinction; upper bound for blocking).
CLOUD_DIR = '/resnick/groups/carnegie_poc/jingze/obscheck_1102/prt_data/input_data/opacities/continuum/clouds/'
COND = {   # 1e4/T_c = a + b log10 p[bar] + c (log10 p)^2   ([M/H] = 0)
    'Mg2SiO4': (5.89, -0.37, 0.0, 2.05e-3, 'Mg2SiO4(s)_crystalline_000'),     # Visscher+10; X: Mg/2 limited
    'MgSiO3': (6.26, -0.35, 0.0, 9.2e-4, 'MgSiO3(s)_crystalline_000'),        # Visscher+10; X: Si left after Mg2SiO4
    'Fe': (5.44, -0.48, 0.0, 1.29e-3, 'Fe(s)_crystalline_000'),               # Visscher+10
    'Al2O3': (5.0139, -0.21794, 2.2636e-3, 1.05e-4, 'Al2O3(s)_crystalline_000'),  # Wakeford+17
    'CaTiO3': (5.125, -0.277, 0.0, 8.9e-6, 'CaTiO3(s)_crystalline_000'),      # Wakeford+17; X: Ti limited
}
SETS = {'A': ['Al2O3', 'CaTiO3'],                                  # high-T condensates only
        'C': ['Al2O3', 'CaTiO3', 'Mg2SiO4', 'MgSiO3'],              # no Fe (nucleation-limited)
        'B': ['Al2O3', 'CaTiO3', 'Mg2SiO4', 'MgSiO3', 'Fe']}        # all


def tcond(sp, lp):
    a, b, c = COND[sp][:3]
    return 1e4/(a + b*lp + c*lp*lp)


def cloud_kappa(sp, R_cm, mode, lam_um):
    import glob, h5py
    fn = glob.glob(CLOUD_DIR + COND[sp][4] + '/*/*Mie.R39*.h5')[0]
    f = h5py.File(fn, 'r')
    r = f['particles_radii'][:]
    i = int(np.argmin(abs(np.log(r) - np.log(R_cm))))
    lw = np.log(1e4/f['wavenumbers'][:])
    k = f['absorption_opacities'][i].astype(float)
    if mode == 'ext':
        k = k + (1 - f['asymmetry_parameters'][i])*f['scattering_opacities'][i]
    o = np.argsort(lw)
    return np.exp(np.interp(np.log(lam_um), lw[o], np.log(np.maximum(k[o], 1e-30))))


def cloud_mass(m, species, ext):
    """(ncol, nr, ns) condensate mass fraction per cell for one case."""
    T, lp = m['T'], m['lp']
    out = np.zeros(T.shape + (len(species),), np.float32)
    for j, sp in enumerate(species):
        c = T < tcond(sp, lp)
        if np.isfinite(ext):
            lpb = np.where(c, lp, -np.inf).max(1)                       # base: deepest condensing cell
            c &= lp >= (lpb - ext)[:, None]
        out[..., j] = COND[sp][3]*c
    return out


KSRC_CL = KSRC.replace('void trace(', 'void trace_cl(').replace(
    'double* F, double* Iw, int storeI)', 'double* F, double* Iw, int storeI, const float* cw, const float* kc, int ncs)').replace(
    '    double dt = exp((double)l)*rd;',
    '    double kcl = 0.0;\n    for (int c = 0; c < ncs; c++) kcl += (double)cw[o*ncs + c]*kc[c*NW + w];\n'
    '    double dt = (exp((double)l) + kcl)*rd;')


def main_clouds(a, kv):
    """clouds=CASES.json: list of {name, set, ext (decades or 'inf'), R_um, mode}; one dump, all cases."""
    import cupy as cp
    cases = json.load(open(kv['clouds']))
    t0 = time.time()
    m = load_model(a.dump, a.eos)
    tr = Tracer(m, a.table, int(kv.get('nb', 96)), int(kv.get('nth', 64)), int(kv.get('nsub', 3)))
    kern = cp.RawKernel(KSRC_CL, 'trace_cl')
    nph = int(kv.get('nph', 20))
    phases = np.arange(nph)/nph
    lam = tr.lam*1e4
    CW, KC, frac = [], [], []
    for c in cases:
        sp = SETS[c['set']]
        ext = float(c['ext'])
        M = cloud_mass(m, sp, ext)*float(c.get('f', 1.0))
        CW.append(cp.asarray(M))
        KC.append(cp.asarray(np.array([cloud_kappa(s_, c['R_um']*1e-4, c['mode'], lam) for s_ in sp], np.float32)))
        frac.append([float((M[..., j] > 0).mean()) for j in range(len(sp))])
    cp.cuda.Device().synchronize()
    t1 = time.time()
    Fp = np.zeros((len(cases), nph, tr.NW))
    area, lamc, wg = cp.asarray(tr.area), cp.asarray(tr.lam), cp.asarray(tr.wg)
    for i, ph in enumerate(phases):
        jT, fT, jP, fP, rds, T, Tbot = tr.samples(ph)
        col, ir = tr._colir
        nray, ns = rds.shape
        n = nray*tr.NW*tr.NG
        for k, c in enumerate(cases):
            cw = cp.ascontiguousarray(CW[k][col, ir]).astype(cp.float32)          # (nray, ns, ncs)
            F = cp.zeros(tr.NW, cp.float64)
            Iw = cp.zeros((1, 1), cp.float64)
            kern(((n + 255)//256,), (256,),
                 (tr.lk, jT, fT, jP, fP, rds, T, Tbot, area, lamc, wg, np.int32(nray), np.int32(ns), np.int32(tr.nP),
                  np.int32(tr.NW), np.int32(tr.NG), F, Iw, np.int32(0), cw, KC[k], np.int32(KC[k].shape[0])))
            Fp[k, i] = cp.asnumpy(F)
    cp.cuda.Device().synchronize()
    t2 = time.time()
    tag = kv.get('tag', os.path.basename(a.dump).rsplit('.', 1)[0])
    np.savez(os.path.join(a.out, 'emcl_%s.npz' % tag), phases=phases, lam_um=lam, Fp=Fp, cases=json.dumps(cases),
             cell_frac=json.dumps(frac), rot=m['time'], table=a.table, dump=a.dump, timing=[t1 - t0, t2 - t1])
    print('%s clouds: %d cases x %d phases: setup %.1f s, GPU %.1f s (%.1f s/case)' % (
        tag, len(cases), nph, t1 - t0, t2 - t1, (t2 - t1)/len(cases)), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eos', required=True)
    ap.add_argument('--table', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('dump')
    ap.add_argument('keys', nargs='*')
    a = ap.parse_args()
    kv = dict(k.split('=', 1) for k in a.keys)
    if 'clouds' in kv:
        return main_clouds(a, kv)
    t0 = time.time()
    m = load_model(a.dump, a.eos)
    if 'iso' in kv:                     # validation: isothermal T on the model's rho, p
        m['T'][:] = float(kv['iso'])
    t1 = time.time()
    tr = Tracer(m, a.table, int(kv.get('nb', 96)), int(kv.get('nth', 64)), int(kv.get('nsub', 3)))
    nph = int(kv.get('nph', 40))
    phases = np.arange(nph)/nph
    if os.environ.get('EM3D_CPU'):
        cp = None
    else:
        import cupy as cp
        cp.cuda.Device().synchronize()
    t2 = time.time()
    Fp = np.zeros((nph, tr.NW))
    Imaps = {}
    for i, ph in enumerate(phases):
        st = abs(ph - 0.5) < 1e-9 or ph == 0.0
        Fp[i], I = tr.flux(ph, storeI=st)
        if st:
            Imaps['I_%03d' % int(round(ph*1000))] = I.astype(np.float32)
    if cp is not None:
        cp.cuda.Device().synchronize()
    t3 = time.time()
    tag = kv.get('tag', os.path.basename(a.dump).rsplit('.', 1)[0])
    os.makedirs(a.out, exist_ok=True)
    np.savez(os.path.join(a.out, 'em_%s.npz' % tag), phases=phases, lam_um=tr.lam*1e4, Fp=Fp, b=tr.b, th=tr.th,
             area=tr.area, rot=m['time'], table=a.table, dump=a.dump, kv=str(kv),
             timing=[t1 - t0, t2 - t1, t3 - t2], **Imaps)
    print('%s: load %.1f s, setup %.1f s, %d phases x %d rays x %d bins x %d g on GPU %.1f s (%.2f s/phase)' % (
        tag, t1 - t0, t2 - t1, nph, len(tr.b), tr.NW, tr.NG, t3 - t2, (t3 - t2)/nph), flush=True)


if __name__ == '__main__':
    main()
