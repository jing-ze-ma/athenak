"""fe3d.py: 3-D Fe I high-resolution transmission through an AthenaK WASP-121b dump.

Per exposure (orbital phase) the planet is rotated (tidally locked, prograde; observer geometry
as kpvsys/kp.py = obs.winds_one rotated with phase), rays are cast on a (position angle theta,
impact parameter b) annulus grid, each ray is sampled uniformly in path length s, and
tau(lambda) = sum_lines tau_l + tau_cont is integrated with the local LOS velocity Doppler-shifting
every line (planet frame: winds + rotation; + = redshift).  Line opacity: every Kurucz gfall Fe I
line in 380-790 nm (vacuum) with exact per-sample Gaussian (thermal) cores and Lorentz wings
(radiative + vdW damping); per-line E_low dependence through an E_low grid (0.05 eV, linear
interpolation).  Continuum: H- bf/ff + H2/He/H Rayleigh at 9 nodes, interpolated in ln lambda.
Excess absorption per ray = exp(-tau_c) (1 - exp(-tau_line)), summed over rays in front of the
stellar disc (uniform disc, partial overlap per ray), normalised by pi R_s^2.

  python fe3d.py --planet planets/wasp121b.json --eos EOS_TABLE --out OUTDIR \
         [--deck run.athinput] DUMP.bin [key=val ...]          -> OUTDIR/spec_<tag>.npz
  (DUMP: AthenaK cubed-sphere .bin with dens, velx-z, eint; hydro_w or mhd_w_bcc)
  env NPROC (default 16): worker processes, one exposure each.
fehires_fast (11-02): pcuts=1e-6,0 traces once for several cuts -> spec_<tag>_pc6.npz,
  spec_<tag>_pc0.npz (or tags=a,b); gpu=1 (or env FEHIRES_GPU=1) runs the rays on a GPU (fe3d_gpu.py,
  CuPy); default = the CPU Pool path.
keys: pcut (bar; default planet pcut_bar = 1e-6), dth (deg, 2), bsub (b points per radial
      cell, 2), ds (cm, 5e7), nph (planet nph = 31), afe (planet a_fe = 3.2e-5),
      tcond (0 = no condensation), static (0/1: isothermal static test), rot (0/1 solid body only
      for the static test), tag, phases (comma list overrides nph), Tiso (2500)
"""
import argparse
import os
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fe_core as F  # noqa: E402
import csio          # noqa: E402
import planet        # noqa: E402
NPROC = int(os.environ.get('NPROC', '16'))
P = {}                                     # planet.load() dict, set by configure()


def configure(cfg):
    """system geometry (a/R*, b, R*, Rp, Omega, P, limb darkening) from planet.load()."""
    global AOR, BIMP, KRP, PH14, RS
    P.clear()
    P.update(cfg)
    F.configure(cfg)
    AOR, BIMP, KRP, PH14, RS = cfg['a_over_Rs'], cfg['b'], cfg['KRP'], cfg['PH14'], cfg['RS']


configure(planet.load())
NV = 150                                   # line velocity grid +-NV pixels (+-45 km/s)
VG = F.DV*np.arange(-NV, NV + 1)
EG = np.arange(0.0, 6.2001, 0.05)          # E_low grid [eV]
LNODE = np.exp(np.linspace(np.log(F.LAM0), np.log(F.LAM1), 9))
NPIX = int(np.log(F.LAM1/F.LAM0)/F.DLN) + 1
TAUL_MIN = 1e-4                            # per-ray line pruning (peak tau); key taumin
G = {}
from numba import njit  # noqa: E402
import math  # noqa: E402


@njit(cache=True)
def _deposit(keep, lc, gr, gw, ie, fe, iq, fq, Gm, Am, Bm, tau, nv):
    for k in keep:
        e = ie[k]
        f = fe[k]
        c = lc[k]
        a = gr[k]
        b = gw[k]
        base = iq[k] - nv
        q = fq[k]
        for j in range(2*nv + 1):
            t = c*((Gm[e, j]*(1 - f) + Gm[e + 1, j]*f) + a*(Am[e, j]*(1 - f) + Am[e + 1, j]*f)
                   + b*(Bm[e, j]*(1 - f) + Bm[e + 1, j]*f))
            if t > 1e-14:
                tau[base + j] += t*(1 - q)
                tau[base + j + 1] += t*q


@njit(cache=True)
def _absorb(tau, lo, hi, nodei, nodew, ltauc, tmpl, row, area):
    W = 0.0
    for p in range(lo, hi):
        t = tau[p]
        if t != 0.0:
            if t > 1e-10:
                i = nodei[p]
                lt = ltauc[i]*(1 - nodew[p]) + ltauc[i + 1]*nodew[p]
                ex = math.exp(-math.exp(lt))*(-math.expm1(-t))*area
                row[p] += ex
                W += ex*tmpl[p]
            tau[p] = 0.0
    return W


@njit(cache=True)
def _absorbN(taus, lo, hi, nodei, nodew, ltauc, tmpl, rows, area, W):
    """_absorb for several cuts at once (taus, rows: (ncut, .)); the continuum factor is shared;
    per cut the arithmetic and the summation order are those of _absorb (bitwise identical)."""
    nc = taus.shape[0]
    for p in range(lo, hi):
        have = False
        cf = 0.0
        for c in range(nc):
            t = taus[c, p]
            if t != 0.0:
                if t > 1e-10:
                    if not have:
                        i = nodei[p]
                        lt = ltauc[i]*(1 - nodew[p]) + ltauc[i + 1]*nodew[p]
                        cf = math.exp(-math.exp(lt))
                        have = True
                    ex = cf*(-math.expm1(-t))*area
                    rows[c, p] += ex
                    W[c] += ex*tmpl[p]
                taus[c, p] = 0.0


def cut_list(kv):
    """pressure cuts [bar]: key pcuts=a,b,.. (one trace for all), else pcut (one cut)."""
    if 'pcuts' in kv:
        return [float(x) for x in kv['pcuts'].split(',')]
    return [float(kv.get('pcut', P['pcut_bar']))]


def cut_tags(kv, base):
    """output tags: tags=a,b,.. explicit; one cut -> tag (default dump basename, as before);
    several cuts -> <tag>_pc0 for 0, <tag>_pc<n> for 10^-n bar (pc6 = 1e-6), else _pc<value>."""
    cuts = cut_list(kv)
    tag = kv.get('tag', base.rsplit('.', 1)[0])
    if 'tags' in kv:
        t = kv['tags'].split(',')
        assert len(t) == len(cuts), 'tags and pcuts differ in length'
        return t
    if len(cuts) == 1:
        return [tag]
    out = []
    for c in cuts:
        n = -np.log10(c) if c > 0 else 0.0
        out.append('%s_pc%s' % (tag, '0' if c == 0 else ('%d' % round(n) if abs(n - round(n)) < 1e-9
                                                            else ('%g' % c))))
    return out


def setup_lines():
    lam, gf, El = F.lines()
    gr, gw = F.line_damping()
    m = (lam > F.LAM0*(1 + 50e5/F.CL)) & (lam < F.LAM1*(1 - 50e5/F.CL)) & (El < EG[-1])
    lam, gf, El, gr, gw = lam[m], gf[m], El[m], gr[m], gw[m]
    q = (np.log(lam) - np.log(F.LAM0))/F.DLN
    iq = np.floor(q).astype(np.int64)
    fe = El/0.05
    ie = np.minimum(fe.astype(int), len(EG) - 2)
    G.update(lc=F.SIGF*gf*lam, iq=iq, fq=q - iq, ie=ie, fe=fe - ie, lam=lam,
             gr=gr*lam/(4*np.pi**2), gw=gw*lam/(4*np.pi**2), nl=len(lam))
    G['tmpl'] = np.maximum(F.template_hp(P['LSF']), 0)


def load(dump, eosfile, kv, deck=None):
    comp = F.Comp()
    afe = float(kv.get('afe', F.A_FE_SUN))
    tcond = float(kv.get('tcond', 0.0))
    pcuts = [c*1e6 for c in cut_list(kv)]             # Fe removed at p < pcut [bar]
    raw = csio.read_dump(dump)
    mb = raw['mb_data']
    eos = csio.EOS(eosfile)
    rho = np.asarray(mb['dens'], float)
    T, p = eos.invert(rho, np.asarray(mb['eint'], float))
    rf, rc = csio.radii(raw, deck)
    X, EXI, EET, _ = csio.geometry(raw)
    ni = rho.shape[-1]
    if int(kv.get('static', 0)):
        # 1-D isothermal hydrostatic atmosphere on the same grid: T = Tiso, p(r) from the
        # column-median p at the first cell above 0.1 bar, g = G M / r^2, v = 0 (or solid body)
        Tiso = float(kv.get('Tiso', 2500.0))
        mu = 2.33
        gm = 6.674e-8*P['MP']
        pm = np.median(p.reshape(-1, ni), axis=0)
        j = int(np.argmax(pm < 1e5))
        lp_ = np.log(pm[j]) - gm*mu*F.MU_/(F.KB*Tiso)*(1/rc[j] - 1/rc)
        pp = np.exp(lp_)
        T = np.broadcast_to(np.full(ni, Tiso), rho.shape).copy()
        p = np.broadcast_to(pp, rho.shape).copy()
        rho = np.empty_like(p)
        # rho from the run composition at (p, T): iterate mu
        r_ = p.reshape(-1, ni)[0]*mu*F.MU_/(F.KB*Tiso)
        for _ in range(6):
            c = comp(r_, np.full(ni, Tiso))
            r_ = p.reshape(-1, ni)[0]*c['mu']*F.MU_/(F.KB*Tiso)
        rho[...] = r_
        vr = v2 = v3 = np.zeros_like(rho)
        omega = P['OMEGA'] if int(kv.get('rot', 0)) else 0.0
    else:
        vr, v2, v3 = [np.asarray(mb[k], float) for k in ('velx', 'vely', 'velz')]
        omega = P['OMEGA']
    ncol = X[..., 0].size
    pf = p.reshape(ncol, ni)
    i0 = max(int(np.min(np.argmax(pf < 1e6, axis=1))) - 1, 0)
    sl = np.s_[..., i0:]
    Vx = (vr*X[..., 0:1] + v2*EXI[..., 0:1] + v3*EET[..., 0:1])[sl]
    Vy = (vr*X[..., 1:2] + v2*EXI[..., 1:2] + v3*EET[..., 1:2])[sl]
    rcl = rc[i0:]
    Vx = Vx - omega*rcl*X[..., 1:2]
    Vy = Vy + omega*rcl*X[..., 0:1]
    rho_, T_ = rho[sl].reshape(ncol, -1), T[sl].reshape(ncol, -1)
    nr = rho_.shape[1]
    lw = np.empty((len(pcuts), ncol, nr), np.float32)
    nP = np.empty((ncol, nr), np.float32)
    lac = np.empty((ncol, nr, len(LNODE)), np.float32)
    ch = max(1, 200000//nr)
    for s in range(0, ncol, ch):
        sp = F.species(rho_[s:s+ch], T_[s:s+ch], comp, afe=afe, tcond=tcond)
        for ic, pcut in enumerate(pcuts):          # one Fe weight per cut, same values as before
            lw[ic, s:s+ch] = np.log(np.maximum(sp['fe1_over_U'], 1e-300))
            lw[ic, s:s+ch][pf[s:s+ch, i0:] < pcut] = -700.0
        nP[s:s+ch] = sp['nH'] + sp['nH2'] + sp['nHe']
        lac[s:s+ch] = np.log(np.maximum(F.alpha_cont(LNODE, T_[s:s+ch], sp), 1e-300))
    from scipy.spatial import cKDTree
    G.update(rc=rcl, rf=rf[i0:], r1=rf[-1], tree=cKDTree(X.reshape(ncol, 3)), lw=lw, nP=nP,
             lac=lac, T=T_.astype(np.float32), lp=np.log10(pf[:, i0:]/1e6).astype(np.float32),
             Vx=Vx.reshape(ncol, -1).astype(np.float32), Vy=Vy.reshape(ncol, -1).astype(np.float32),
             rot=raw['time']/P['PROT'], i0=i0)
    # impact parameters: bsub points per radial cell (cell-centred sub-intervals)
    bsub = int(kv.get('bsub', 2))
    rfl = G['rf']
    bb = np.concatenate([rfl[k] + (rfl[k+1] - rfl[k])*(np.arange(bsub) + 0.5)/bsub
                         for k in range(len(rfl) - 1)])
    db = np.concatenate([np.full(bsub, (rfl[k+1] - rfl[k])/bsub) for k in range(len(rfl) - 1)])
    G.update(bb=bb, db=db)


def ray_one(r, col, ds, ob, area, rows, taumin):
    """one ray, all cuts: deposit line tau, add excess absorption * area into rows[cut]; returns
    (template-weighted absorption W[cut], log p_tan).  Per cut the arithmetic is the one-cut code's
    (the cut only changes the Fe weight lw); the cut-independent work is done once."""
    rc = G['rc']
    fi = np.interp(r, rc, np.arange(len(rc)))
    i0 = np.clip(fi.astype(int), 0, len(rc) - 2)
    f = fi - i0
    c = col

    def ip(A):
        return A[c, i0]*(1 - f) + A[c, i0 + 1]*f
    lp = ip(G['lp'])
    lac = G['lac'][c, i0]*(1 - f)[:, None] + G['lac'][c, i0 + 1]*f[:, None]
    tauc = (np.exp(lac)*ds).sum(0)
    lpt = lp.max()
    LW = G['lw']
    nc = LW.shape[0]
    Wc = np.zeros(nc)
    if tauc.min() > 15:
        return Wc, lpt
    lws = [ip(LW[k]) for k in range(nc)]
    act = [lw.max() >= -600 for lw in lws]
    if not any(act):
        return Wc, lpt
    T = ip(G['T'])
    v = -(ip(G['Vx'])*ob[0] + ip(G['Vy'])*ob[1])
    nP = ip(G['nP'])*(T/1e4)**0.3
    bth = np.sqrt(2*F.KB*T/F.M_FE)
    dv = VG[None, :] - v[:, None]
    x = dv/bth[:, None]
    Phi = np.exp(-x*x)/(np.sqrt(np.pi)*bth[:, None])
    Kw = (1 - np.exp(-(x/2.5)**4))/np.maximum(dv*dv, 1.0)
    Ez = np.exp(-np.outer(F.EV/(F.KB*T), EG))                   # (S, nE), w applied per cut
    ie, fe = G['ie'], G['fe']
    taus = G['tau']
    lo, hi = NPIX, 0
    for k in range(nc):
        if not act[k]:
            continue
        lw = lws[k]
        m = lw > lw.max() - 25
        w = np.exp(lw[m])*ds
        Bz = Ez[m]*w[:, None]
        Phm, Kwm = Phi[m], Kw[m]
        Gm = Bz.T @ Phm
        Am = Bz.T @ Kwm
        Bm = (Bz*nP[m][:, None]).T @ Kwm
        gmax = Gm.max(1)
        tpk = G['lc']*(gmax[ie]*(1 - fe) + gmax[ie + 1]*fe)
        keep = np.flatnonzero(tpk > taumin)
        if keep.size == 0:
            continue
        _deposit(keep, G['lc'], G['gr'], G['gw'], G['ie'], G['fe'], G['iq'], G['fq'], Gm, Am, Bm,
                 taus[k], NV)
        lo = min(lo, max(int(G['iq'][keep].min()) - NV - 1, 0))
        hi = max(hi, min(int(G['iq'][keep].max()) + NV + 3, NPIX))
    if hi > lo:
        _absorbN(taus, lo, hi, G['nodei'], G['nodew'], np.log(np.maximum(tauc, 1e-300)),
                 G['tmpl'], rows, area, Wc)
    return Wc, lpt


def one_phase(args):
    ph, kv = args
    t0 = time.time()
    a = 2*np.pi*ph
    ob = np.array([np.cos(a), -np.sin(a), 0.0])
    e1 = np.array([np.sin(a), np.cos(a), 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    dth = float(kv.get('dth', 2.0))
    taumin = float(kv.get('taumin', TAUL_MIN))
    th = np.deg2rad(np.arange(0, 360, dth) + dth/2)
    ds = float(kv.get('ds', 5e7))
    bb, db = G['bb'], G['db']
    r1 = G['r1']
    xs, ys = AOR*np.sin(a), BIMP*np.cos(a)
    u1, u2 = P['limb_darkening']
    ldn = 1 - u1/3 - u2/6                  # disc-integrated intensity / (pi R*^2 I_centre)
    nc = G['lw'].shape[0]
    spec = np.zeros((nc, 2, NPIX))
    G['tau'] = np.zeros((nc, NPIX + 2))
    lgl = np.log(F.LAM0) + F.DLN*np.arange(NPIX)
    lgn = np.log(LNODE)
    G['nodei'] = np.clip(np.searchsorted(lgn, lgl) - 1, 0, len(lgn) - 2)
    G['nodew'] = (lgl - lgn[G['nodei']])/(lgn[G['nodei'] + 1] - lgn[G['nodei']])
    Wm = np.zeros((nc, len(th), len(bb)))
    Pt = np.zeros((len(th), len(bb)))
    for it, t_ in enumerate(th):
        limb = 0 if np.cos(t_) > 0 else 1
        sx = xs + bb*np.cos(t_)/RS
        sy = ys + bb*np.sin(t_)/RS
        onstar = (sx**2 + sy**2) < 1.0
        if u1 != 0.0 or u2 != 0.0:         # quadratic limb darkening, I/I_centre per ray
            mu1 = 1 - np.sqrt(np.clip(1 - sx**2 - sy**2, 0, 1))
            ldw = (1 - u1*mu1 - u2*mu1**2)/ldn
        else:                              # uniform disc (the 10-07 baseline)
            ldw = None
        for ib, b_ in enumerate(bb):
            if not onstar[ib]:
                continue
            smax = np.sqrt(max(r1**2 - b_**2, 0.0))
            ns = max(int(2*smax/ds), 4)
            s = (np.arange(ns) + 0.5)/ns*2*smax - smax
            P0 = b_*(np.cos(t_)*e1 + np.sin(t_)*e2)
            pts = P0[None, :] + s[:, None]*ob[None, :]
            r = np.linalg.norm(pts, axis=1)
            ok = (r >= G['rc'][0]) & (r <= r1)
            if ok.sum() < 2:
                continue
            pts, r = pts[ok], r[ok]
            _, col = G['tree'].query(pts/r[:, None])
            area = b_*db[ib]*np.deg2rad(dth)
            if ldw is not None:
                area = area*ldw[ib]
            Wm[:, it, ib], Pt[it, ib] = ray_one(r, col, 2*smax/ns, ob, area, spec[:, limb],
                                                taumin)
    norm = np.pi*RS**2
    return dict(spec=(spec/norm).astype(np.float32), W=Wm/norm, ptan=Pt, th=th,
                sec=time.time() - t0)


def run(dump, eosfile, outdir, kv, deck=None):
    from multiprocessing import Pool
    setup_lines()
    t0 = time.time()
    load(dump, eosfile, kv, deck)
    if 'phases' in kv:
        phases = np.array([float(x) for x in kv['phases'].split(',')])
    elif P['phase_range']:
        phases = np.linspace(P['phase_range'][0], P['phase_range'][1],
                             int(kv.get('nph', P['nph'])))
    else:
        phases = np.linspace(-PH14, PH14, int(kv.get('nph', P['nph'])))
    base = os.path.basename(dump)
    print('%s %s: %s rot %.2f loaded in %.0f s: %d b, i0 %d, lines %d' % (
        P['name'], base, dump, G['rot'], time.time() - t0, len(G['bb']), G['i0'], G['nl']),
        flush=True)
    cuts, tags = cut_list(kv), cut_tags(kv, base)
    pcs = kv['pcuts'].split(',') if 'pcuts' in kv else [kv.get('pcut', repr(cuts[0]))]
    if int(kv.get('gpu', os.environ.get('FEHIRES_GPU', '0'))):
        sys.modules['fe3d'] = sys.modules[__name__]       # fe3d_gpu shares this module's G, P
        import fe3d_gpu
        if kv.get('prep', os.environ.get('FEHIRES_PREP', 'cpu')) == 'gpu':   # 11-03: prep on the GPU
            out = fe3d_gpu.run_phases_gpuprep(phases, kv)
        else:
            out = fe3d_gpu.run_phases(phases, kv)
    else:
        with Pool(min(NPROC, len(phases))) as pool:
            out = pool.map(one_phase, [(ph, kv) for ph in phases], chunksize=1)
    os.makedirs(outdir, exist_ok=True)
    for ic, (pc, tag) in enumerate(zip(cuts, tags)):
        kvc = {k: v for k, v in kv.items() if k not in ('pcuts', 'tags', 'pcut', 'tag')}
        kvc.update(pcut=pcs[ic], tag=tag)                 # as the one-cut runs record it
        np.savez(os.path.join(outdir, 'spec_%s.npz' % tag), phases=phases,
                 spec=np.array([o['spec'][ic] for o in out]),
                 W=np.array([o['W'][ic] for o in out]), ptan=np.array([o['ptan'] for o in out]),
                 th=out[0]['th'], bb=G['bb'], db=G['db'], rot=G['rot'], sec=[o['sec'] for o in out],
                 kv=str(kvc), planet=str({k: P[k] for k in planet.DEFAULTS}), dump=dump)
    print('%s done in %.0f s (per phase %.0f s)' % (','.join(tags), time.time() - t0,
                                                    np.mean([o['sec'] for o in out])), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--planet', help='planet JSON (default: WASP-121b values)')
    ap.add_argument('--eos', required=True, help="the run's eos_table*.txt (general EOS dump)")
    ap.add_argument('--out', required=True, help='output directory')
    ap.add_argument('--deck', help='run .athinput (only for version-1.1 dumps)')
    ap.add_argument('dump')
    ap.add_argument('keys', nargs='*', help='key=val overrides (see module doc)')
    a = ap.parse_args()
    configure(planet.load(a.planet))
    run(a.dump, a.eos, a.out, dict(k.split('=', 1) for k in a.keys), a.deck)
