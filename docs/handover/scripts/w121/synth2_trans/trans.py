"""trans.py: WASP-121b GCM transmission spectra of Fe I and Na D through the transit.

Chords through the 3-D rot-300 state (dhjcs panel frames, geometry of synth.py 'winds'),
per-cell FastChem equilibrium (chem.py tables), Kurucz gfall Fe I / Na I lines with
Voigt profiles (thermal + radiative + van der Waals, first-order-in-a Voigt), continuum
H- bf/ff (John 1988) + Rayleigh (H2, H, He).  Each cell's opacity is Doppler-shifted by its
inertial LOS velocity (wind + solid-body rotation).  Output per epoch and limb: excess
absorption (transit-depth units) in velocity windows around the Fe mask lines and Na D.

  python trans.py run <run> <tag> [key=val ...]   -> out/<run>_<tag>.npz
  python trans.py rhotest <run>                    -> out/<run>_rhotest.txt  (validation 7b)

keys: wind=1 rot=1 ext=0 (isothermal hydrostatic extension above x1max up to REXT)
      pext=1 (extension pressure factor) fe=1 (Fe abundance factor) chem=cond|gas
      iso=0 (T [K] of a static isothermal test atmosphere; 0 = the GCM state)
      nep=15 (epochs from 1st to 4th contact) ep=all|mid  npa=180 nproc=15 ld=1
"""
import os
import sys
import time
import numpy as np
import scipy.sparse as sps
from scipy.special import dawsn, wofz
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts')
import dhjcs   # noqa: E402

W = '/viper/ptmp2/jinma/w121prod_0929/synth2_trans/'
D = W + 'data/'
h, c, kB, amu = 6.62607015e-27, 2.99792458e10, 1.380649e-16, 1.66053907e-24
PIE2MC = 0.0265400      # pi e^2 / (m_e c) [cm^2 s^-1]
RJ, RSUN = 7.1492e9, 6.957e10
RP, RS = 1.742*RJ, 1.461*RSUN
AOR, BIMP = 3.7844, 0.10          # a/R*, impact parameter (Sing+2024 / Bourrier+2020)
U1, U2 = 0.35, 0.28              # quadratic limb darkening, approximate V band 6600 K
OMEGA = 5.704026e-05
GM = 1167.8749*1.126575e10**2    # grav * ap^2 (1x input; point mass)
RUNS = {
    '1x': dict(bin='/viper/ptmp2/jinma/w121prod_0929/w1x/bin/dhj.hydro_w.%05d.bin',
               eos='/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt',
               R0=1.126575e10, R1=1.631037e10, nr=76,
               CSTR=[-0.379629, -1.797216, 7.873416, -79.711522, 347.885804, -753.943230,
                     784.221097, -314.370119]),
    '10x': dict(bin='/viper/ptmp2/jinma/w121prod_0929/w10x/bin/dhj.hydro_w.%05d.bin',
                eos='/viper/ptmp2/jinma/w121prod_0929/synth_rot300_10x/diag/eos_table_10x.txt',
                R0=1.116018e10, R1=1.530921e10, nr=74,
                CSTR=[-0.090789, -4.624280, 43.684697, -282.659406, 932.446804, -1642.475346,
                      1459.195482, -514.767866])}
DUMP = 150
REXT = 2.2e10          # outer radius of the (optional) isothermal extension
DU = 0.5e5             # velocity pixel [cm/s] (R = 600 000)
WFE, VP = 40e5, 50e5   # Fe window half width, profile half extent
WNA = 100e5            # Na window half width
UW = np.arange(-WFE, WFE + 1, DU)
UP = np.arange(-VP, VP + 1, DU)
UNA = np.arange(-WNA, WNA + 1, DU)
DE = 0.05              # eV, excitation-energy group spacing
EV = 1.602176634e-12
LC = np.linspace(0.37, 0.80, 12)   # um, continuum reference wavelengths
MFE, MNA = 55.845*amu, 22.98977*amu


# ------------------------------------------------------------------ data
def gfall(fn, l0=380.0, l1=790.0):
    """Kurucz gfall (air wavelengths nm > 200): lam[nm], loggf, Elow[eV], logGr, logGw"""
    out = []
    for ln in open(fn):
        try:
            lam = float(ln[0:11])
            if lam < l0 or lam > l1:
                continue
            e1, e2 = abs(float(ln[24:36])), abs(float(ln[52:64]))
            out.append((lam, float(ln[11:18]), min(e1, e2)*1.239841984e-4,
                        float(ln[80:86]), float(ln[92:98])))
        except ValueError:
            continue
    return np.array(out)


def partfn(sym):
    """Barklem & Collet 2016 atomic partition function -> callable Q(T)"""
    Tg = None
    for ln in open(D + 'table8.dat'):
        if ln.startswith('#') and 'T [K]' in ln:
            Tg = np.array(ln.split()[3:], float)
        elif ln.split() and ln.split()[0] == sym:
            q = np.array(ln.split()[1:], float)
            break
    return lambda T: np.interp(np.log(T), np.log(Tg), q)


QFE, QNA = partfn('Fe_I'), partfn('Na_I')


def select_lines(nmask=300, sep=20e5, Tref=2500.0):
    L = gfall(D + 'gf2600.all')
    S = 10**L[:, 1]*L[:, 0]*np.exp(-L[:, 2]*EV/(kB*Tref))
    o = np.argsort(-S)
    mask = []
    for i in o:
        if all(abs(L[i, 0] - L[j, 0])/L[j, 0]*c > sep for j in mask):
            mask.append(i)
        if len(mask) == nmask:
            break
    mask = np.array(sorted(mask, key=lambda i: L[i, 0]))
    blend = S >= 0.01*S[mask].min()
    return L, mask, blend


def build_M(L, mask, blend):
    """sparse map: group profiles tau_g(UP) (ng*nup) -> window tau (nmask*nuw)"""
    bl = np.where(blend)[0]
    ng = int(np.ceil(L[bl, 2].max()/DE)) + 2
    rows, cols, vals = [], [], []
    nuw, nup = len(UW), len(UP)
    for w, iw in enumerate(mask):
        lw = L[iw, 0]
        du = (L[bl, 0] - lw)/lw*c
        sel = bl[np.abs(du) < WFE + VP]
        du = (L[sel, 0] - lw)/lw*c
        for l, d in zip(sel, du):
            up = UW - d
            ok = np.abs(up) < VP - DU
            fk = (up[ok] + VP)/DU
            k0 = np.floor(fk).astype(int)
            fk -= k0
            fg = L[l, 2]/DE
            g0 = int(fg)
            fg -= g0
            coef = 10**L[l, 1]*L[l, 0]*1e-7
            ri = w*nuw + np.where(ok)[0]
            for dg, wg in ((0, 1 - fg), (1, fg)):
                for dk, wk in ((0, 1 - fk), (1, fk)):
                    rows.append(ri)
                    cols.append((g0 + dg)*nup + k0 + dk)
                    vals.append(coef*wg*wk)
    M = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(len(mask)*nuw, ng*nup))
    return M, ng


def nad():
    L = gfall(D + 'gf1100.all', 588.0, 590.0)
    L = L[np.argsort(-L[:, 1])][:2]
    return L[np.argsort(L[:, 0])]     # D2 588.995, D1 589.592


# ------------------------------------------------------------------ physics
def voigt_v(u, b, a):
    """normalised Voigt profile in velocity (s/cm): H(a,x)/(sqrt(pi) b), first order in a"""
    x = u/b
    H = np.exp(-x*x) + a*(2/np.sqrt(np.pi))*(2*x*dawsn(x) - 1)
    big = (a > 0.02).ravel()
    if big.any():               # exact Faddeeva where the damping is not small
        xb = np.broadcast_to(x, H.shape)[big]
        ab = np.broadcast_to(a, H.shape)[big]
        H[big] = wofz(xb + 1j*ab).real
    return H/(np.sqrt(np.pi)*b)


def cont_alpha(lam_um, T, nHm, nH, ne, nH2, nHe):
    """continuum absorption coefficient [cm^-1], shape (nlam, nsamp)"""
    lam = lam_um[:, None]
    # H- bound-free, John (1988)
    l0 = 1.6419
    x = np.clip(1/lam - 1/l0, 0, None)
    C = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]
    f = sum(Cn*x**(n/2) for n, Cn in enumerate(C))
    sbf = 1e-18*lam**3*x**1.5*f
    stim = 1 - np.exp(-h*c/(lam*1e-4*kB*T[None, :]))
    a = nHm[None, :]*sbf*stim
    # H- free-free, John (1988), lambda > 0.3645 um  [cm^4 dyn^-1] * n_H * P_e
    FF = np.array([[0, 2483.346, -3449.889, 2200.040, -696.271, 88.283],
                   [0, 285.827, -1158.382, 2427.719, -1841.400, 444.517],
                   [0, -2054.291, 8746.523, -13651.105, 8624.970, -1863.864],
                   [0, 2827.776, -11485.632, 16755.524, -10051.530, 2095.288],
                   [0, -1341.537, 5303.609, -7510.494, 4400.067, -901.788],
                   [0, 208.952, -812.939, 1132.738, -655.020, 128.608]])
    th = 5040.0/T[None, :]
    kff = 0
    for n in range(6):
        A_, B_, C_, D_, E_, F_ = FF[n]
        kff = kff + th**((n + 2)/2)*(lam**2*A_ + B_ + C_/lam + D_/lam**2 + E_/lam**3
                                      + F_/lam**4)
    a = a + 1e-29*kff*nH[None, :]*ne[None, :]*kB*T[None, :]
    # Rayleigh (lambda in Angstrom)
    lA = lam*1e4
    a = a + nH2[None, :]*(8.14e-13/lA**4 + 1.28e-6/lA**6 + 1.61/lA**8)
    a = a + nH[None, :]*(5.799e-13/lA**4 + 1.422e-6/lA**6 + 2.784/lA**8)
    a = a + nHe[None, :]*5.484e-14/lA**4*(1 + 2.44e5/lA**2)**2
    return a


class Chem:
    def __init__(self, run, kind):
        z = np.load(W + 'chem_%s.npz' % run)
        self.LT, self.LP = z['LT'], z['LP']
        self.t = {s: z['%s_%s' % (kind, s)] for s in ('Fe', 'Na', 'e-', 'H1-', 'H', 'H2', 'He')}
        self.mu = z['gas_mu']

    def _ip(self, A, lT, lP):
        x = np.clip((lT - self.LT[0])/(self.LT[1] - self.LT[0]), 0, len(self.LT) - 1.001)
        y = np.clip((lP - self.LP[0])/(self.LP[1] - self.LP[0]), 0, len(self.LP) - 1.001)
        i, j = x.astype(int), y.astype(int)
        fx, fy = x - i, y - j
        return (A[i, j]*(1 - fx)*(1 - fy) + A[i + 1, j]*fx*(1 - fy) + A[i, j + 1]*(1 - fx)*fy
                + A[i + 1, j + 1]*fx*fy)

    def __call__(self, T, lpbar):
        n = 10**lpbar*1e6/(kB*T)
        lT = np.log10(T)
        return {s: n*10**self._ip(A, lT, lpbar) for s, A in self.t.items()}

    def mmw(self, T, lpbar):
        return self._ip(self.mu, np.log10(T), lpbar)


# ------------------------------------------------------------------ state
def geometry(raw):
    g = np.asarray(raw['mb_geometry'])
    nmb, nj, nk = g.shape[0], raw['nx2_mb'], raw['nx3_mb']
    per = nmb//6
    X = np.zeros((nmb, nk, nj, 3))
    EXI = np.zeros_like(X)
    EET = np.zeros_like(X)
    for m in range(nmb):
        p = m//per
        xi = np.pi/4*(g[m, 2] + (g[m, 3] - g[m, 2])*(np.arange(nj) + 0.5)/nj)
        eta = np.pi/4*(g[m, 4] + (g[m, 5] - g[m, 4])*(np.arange(nk) + 0.5)/nk)
        XI, ETA = np.meshgrid(xi, eta)
        X[m] = dhjcs.panel_map(p, XI, ETA)
        d = 1e-6
        e = dhjcs.panel_map(p, XI + d, ETA) - dhjcs.panel_map(p, XI - d, ETA)
        EXI[m] = e/np.linalg.norm(e, axis=-1)[..., None]
        e = dhjcs.panel_map(p, XI, ETA + d) - dhjcs.panel_map(p, XI, ETA - d)
        EET[m] = e/np.linalg.norm(e, axis=-1)[..., None]
    return X, EXI, EET


def rgrid(R):
    xi = (np.arange(R['nr']) + 0.5)/R['nr']
    u = xi.copy()
    xk = xi.copy()
    for ck in R['CSTR']:
        u += ck*xk*(1 - xi)
        xk *= xi
    return R['R0'] + (R['R1'] - R['R0'])*u


def load_state(run, opt, chem):
    R = RUNS[run]
    raw = dhjcs.bin_convert.read_binary(R['bin'] % DUMP)
    X, EXI, EET = geometry(raw)
    mb = raw['mb_data']
    rho = np.asarray(mb['dens'], float)
    eos = dhjcs.EOS(R['eos'])
    eint = np.asarray(mb['eint'], float)
    T, p = np.empty_like(rho), np.empty_like(rho)
    for m in range(rho.shape[0]):         # per MeshBlock: keeps the inversion RSS small
        T[m], p[m] = eos.invert(rho[m], eint[m])
    vr, v2, v3 = (np.asarray(mb[k], float) for k in ('velx', 'vely', 'velz'))
    rc = rgrid(R)
    V = (vr[..., None]*X[:, :, :, None, :] + v2[..., None]*EXI[:, :, :, None, :]
         + v3[..., None]*EET[:, :, :, None, :])
    ncol = X[..., 0].size
    Xf = X.reshape(ncol, 3)
    S = dict(X=Xf, rc=rc, T=T.reshape(ncol, -1), lp=np.log10(p.reshape(ncol, -1)/1e6),
             lrho=np.log(rho.reshape(ncol, -1)),
             vx=V[..., 0].reshape(ncol, -1), vy=V[..., 1].reshape(ncol, -1),
             time=raw['time'], nr=len(rc))
    if opt['iso'] > 0:          # static isothermal test atmosphere, 100 bar at rc[0]
        Ti = float(opt['iso'])
        mu = chem.mmw(np.array([Ti]), np.array([-4.0]))[0]
        lp = 2.0 - GM*mu*amu/(kB*Ti)*(1/rc[0] - 1/rc)/np.log(10)
        S['T'][:] = Ti
        S['lp'][:] = lp[None, :]
        S['vx'][:] = 0
        S['vy'][:] = 0
    if opt['ext']:
        re = np.arange(rc[-1] + 0.5*(rc[-1] - rc[-2]), REXT + 1, 4e7)
        Tt, lpt = S['T'][:, -1], S['lp'][:, -1]
        mu = chem.mmw(Tt, lpt)
        lpe = (lpt[:, None] + np.log10(opt['pext'])
               - (GM*mu*amu/(kB*Tt))[:, None]*(1/rc[-1] - 1/re[None, :])/np.log(10))
        S['rc'] = np.concatenate([rc, re])
        S['lp'] = np.concatenate([S['lp'], lpe], 1)
        for k in ('T', 'vx', 'vy', 'lrho'):
            S[k] = np.concatenate([S[k], np.repeat(S[k][:, -1:], len(re), 1)], 1)
        S['lrho'][:, len(rc):] += np.log(10)*(lpe - lpt[:, None])
    S['vx'] *= opt['wind']
    S['vy'] *= opt['wind']
    return S


# ------------------------------------------------------------------ chords
def chord_samples(S, tree, ob, e1, e2, t_, bb, s):
    """samples of all chords at position angle t_: returns chord index, and per sample the
    column, radial index/fraction (linear in r between cell centres, as synth.py)"""
    rc = S['rc']
    nr = len(rc)
    P0 = bb[:, None]*(np.cos(t_)*e1 + np.sin(t_)*e2)[None, :]
    pts = P0[:, None, :] + s[None, :, None]*ob[None, None, :]
    r = np.linalg.norm(pts, axis=-1)
    ok = (r >= rc[0]) & (r <= rc[-1])
    ci = np.nonzero(ok)[0]
    pts, r = pts[ok], r[ok]
    _, col = tree.query(pts/r[:, None])
    fi = np.interp(r, rc, np.arange(nr))
    i0 = np.clip(fi.astype(int), 0, nr - 2)
    f = fi - i0
    return ci, pts, col, i0, f


def ip(A, col, i0, f):
    return A[col, i0]*(1 - f) + A[col, i0 + 1]*f


def epochs(nep):
    a4 = np.degrees(np.arcsin(np.sqrt((1 + RP/RS*1.0)**2 - BIMP**2)/AOR))
    return np.linspace(-a4, a4, nep)


def stellar_weight(alpha, bb, t_, ld=1):
    """I(mu)/F_* per unit sky area [R*^-2 -> cm^-2 via RS]: chords at (b, t_) at epoch alpha"""
    a = np.deg2rad(alpha)
    Xs = AOR*np.sin(a) + bb*np.cos(t_)/RS
    Ys = BIMP + bb*np.sin(t_)/RS
    d2 = Xs**2 + Ys**2
    mu = np.sqrt(np.clip(1 - d2, 0, 1))
    u1, u2 = (U1, U2) if ld else (0.0, 0.0)
    I = (1 - u1*(1 - mu) - u2*(1 - mu)**2)*(d2 < 1)
    return I/(np.pi*RS**2*(1 - u1/3 - u2/6))


def run_epoch(args):
    run, opt, alpha = args
    chem = Chem(run, opt['chem'])
    S = load_state(run, opt, chem)
    from scipy.spatial import cKDTree
    tree = cKDTree(S['X'])
    L, mask, blend = select_lines()
    M, ng = build_M(L, mask, blend)
    MT = M.T.tocsr()
    NA = nad()
    Lm = L[mask]
    lgr = np.median(np.where(Lm[:, 3] != 0, Lm[:, 3], 8.0))
    lgw = np.median(np.where(Lm[:, 4] != 0, Lm[:, 4], -7.5))
    lamw_um = Lm[:, 0]*1e-3
    nw, nuw, nup, nuna = len(mask), len(UW), len(UP), len(UNA)
    Eg = np.arange(ng)*DE*EV
    rc = S['rc']
    th = np.deg2rad(np.arange(0, 360, 360.0/opt['npa']) + 180.0/opt['npa'])
    dth = 2*np.pi/opt['npa']
    bb = np.linspace(RUNS[run]['R0'], RUNS[run]['R1'], 120)
    db = bb[1] - bb[0]
    if opt['ext']:
        bb = np.concatenate([bb, np.arange(bb[-1] + db, rc[-1], db)])
    Rmax = rc[-1]
    s = np.arange(-Rmax, Rmax + 1, 2*RUNS[run]['R1']/800)
    ds = s[1] - s[0]
    a = np.deg2rad(alpha)
    ob = np.array([np.cos(a), -np.sin(a), 0.0])
    e1 = np.array([np.sin(a), np.cos(a), 0.0])
    e2 = np.array([0.0, 0.0, 1.0])
    lpb = np.arange(-14, 3.01, 0.1)
    out = dict(Afe=np.zeros((2, nw, nuw)), Afe_c=np.zeros((2, nw)),
               Ana=np.zeros((2, 2, nuna)), Ana_c=np.zeros((2, 2)),
               CFfe=np.zeros((2, len(lpb))), CFna=np.zeros((2, len(lpb))),
               CFfe_top=np.zeros(2), CFna_top=np.zeros(2),
               area=np.zeros(2), top_lp=[], EWb=np.zeros((2, len(bb))), CEb=np.zeros((2, len(bb))), bb=bb)
    fe_scale = opt['fe']
    mref = np.ones(nw) if 'mw' not in opt else opt['mw']
    for t_ in th:
        limb = 0 if np.cos(t_) > 0 else 1        # 0 morning (leading, +e1), 1 evening
        ci, pts, col, i0, f = chord_samples(S, tree, ob, e1, e2, t_, bb, s)
        if len(ci) == 0:
            continue
        T = ip(S['T'], col, i0, f)
        lp = ip(S['lp'], col, i0, f)
        vob = ((ip(S['vx'], col, i0, f) + opt['rot']*OMEGA*(-pts[:, 1]))*ob[0]
               + (ip(S['vy'], col, i0, f) + opt['rot']*OMEGA*pts[:, 0])*ob[1])
        n = chem(T, lp)
        wgt = stellar_weight(alpha, bb, t_, opt['ld'])*bb*db*dth
        out['area'][limb] += wgt.sum()
        # continuum per chord at LC, then interpolated to the windows
        ac = cont_alpha(LC, T, n['H1-'], n['H'], n['e-'], n['H2'], n['He'])*ds
        tc = np.zeros((len(bb), len(LC)))
        for k in range(len(LC)):
            tc[:, k] = np.bincount(ci, ac[k], minlength=len(bb))
        # Fe
        b_fe = np.sqrt(2*kB*T/MFE)
        npert = n['H'] + n['H2'] + 0.4*n['He']
        Gam = 10**lgr + 10**lgw*npert*(T/1e4)**0.3
        afe = Gam*550e-7/(4*np.pi*b_fe)
        Afe = PIE2MC*ds*fe_scale*n['Fe']/QFE(T)
        b_na = np.sqrt(2*kB*T/MNA)
        Gna = 10**NA[0, 3] + 10**NA[0, 4]*npert*(T/1e4)**0.3
        ana = Gna*589.3e-7/(4*np.pi*b_na)
        Ana = PIE2MC*ds*n['Na']/QNA(T)          # gf (= g_low f) and lambda per line below
        gfl_na = 10**NA[:, 1]*NA[:, 0]*1e-7
        rs = np.linalg.norm(pts, axis=1)
        cuts = np.searchsorted(ci, np.arange(len(bb) + 1))
        live = np.nonzero((wgt > 0) & (np.diff(cuts) >= 3))[0]
        if len(live) == 0:
            continue
        live_all = live
        for live in np.array_split(live_all, int(np.ceil(len(live_all)/24))):
            tg = np.zeros((len(live), ng*nup))
            tna = np.zeros((len(live), nuna))
            keep = []
            for ii, kc in enumerate(live):
                sl = slice(cuts[kc], cuts[kc + 1])
                Wg = Afe[sl][None, :]*np.exp(-Eg[:, None]/(kB*T[sl][None, :]))
                G = voigt_v(UP[None, :] + vob[sl][:, None], b_fe[sl][:, None], afe[sl][:, None])
                tg[ii] = (Wg @ G).ravel()
                Gn = voigt_v(UNA[None, :] + vob[sl][:, None], b_na[sl][:, None], ana[sl][:, None])
                tna[ii] = Ana[sl] @ Gn
                keep.append((sl, Wg, G, Gn))
            # Fe windows
            tw = (M @ tg.T).T.reshape(len(live), nw, nuw)
            tcw = np.exp(np.array([np.interp(lamw_um, LC, np.log(np.maximum(tc[kc], 1e-300)))
                                   for kc in live]))
            trans = np.exp(-(tw + tcw[:, :, None]))
            wl = wgt[live]
            out['Afe'][limb] += np.einsum('c,cwu->wu', wl, 1 - trans)
            out['Afe_c'][limb] += np.einsum('c,cw->w', wl, 1 - np.exp(-tcw))
            out['EWb'][0, live] += wl*np.einsum('w,cwu->c', mref, np.exp(-tcw)[:, :, None] - trans)
            out['CEb'][0, live] += wl*np.einsum('w,cw->c', mref, np.exp(-tcw) - trans[:, :, nuw//2])
            # Na D (gf lambda of each line; profile shared)
            tcn = np.exp(np.array([np.interp(0.5893, LC, np.log(np.maximum(tc[kc], 1e-300)))
                                   for kc in live]))
            tq = tna[:, None, :]*gfl_na[None, :, None] + tcn[:, None, None]
            out['Ana'][limb] += np.einsum('c,cqu->qu', wl, 1 - np.exp(-tq))
            out['Ana_c'][limb] += (wl @ (1 - np.exp(-tcn)))
            out['EWb'][1, live] += wl*(np.exp(-tcn)[:, None, None] - np.exp(-tq)).sum((1, 2))
            out['CEb'][1, live] += wl*(np.exp(-tcn) - np.exp(-tq[:, 0, nuna//2]))
            # contribution functions: marginal line absorption d(1-e^-tau)/d(sample) =
            # dtau_line(sample) e^-tau, summed over the window pixels (Fe: mask weights)
            yv = (mref[None, :, None]*trans).reshape(len(live), -1)*wl[:, None]
            Y = (MT @ yv.T).T.reshape(len(live), ng, nup)
            yn = np.einsum('q,cqu->cu', gfl_na, np.exp(-tq))*wl[:, None]
            for ii, (sl, Wg, G, Gn) in enumerate(keep):
                cf = (Wg.T*(G @ Y[ii].T)).sum(1)
                ib = np.clip(((lp[sl] - lpb[0])/0.1).astype(int), 0, len(lpb) - 1)
                top = rs[sl] > RUNS[run]['R1']
                out['CFfe'][limb] += np.bincount(ib, cf, len(lpb))
                out['CFfe_top'][limb] += cf[top].sum()
                cfn = Ana[sl]*(Gn @ yn[ii])
                out['CFna'][limb] += np.bincount(ib, cfn, len(lpb))
                out['CFna_top'][limb] += cfn[top].sum()
    out['lpb'] = lpb
    return alpha, out


def parse(argv):
    opt = dict(wind=1.0, rot=1.0, ext=0, pext=1.0, fe=1.0, chem='cond', iso=0.0, nep=15,
               ep='all', npa=180, nproc=15, ld=1)
    for a in argv:
        k, v = a.split('=')
        opt[k] = v if k in ('chem', 'ep') else float(v)
    for k in ('ext', 'nep', 'npa', 'nproc', 'ld'):
        opt[k] = int(opt[k])
    return opt


def main_run(run, tag, argv):
    from multiprocessing import Pool
    opt = parse(argv)
    if 'mw' in opt:
        del opt['mw']
    tfile = W + 'out/%s_template.npz' % run
    if os.path.exists(tfile) and tag != 'template':
        opt['mw'] = np.load(tfile)['mw']
    al = epochs(opt['nep']) if opt['ep'] == 'all' else np.array([0.0])
    t0 = time.time()
    with Pool(min(opt['nproc'], len(al))) as P:
        res = P.map(run_epoch, [(run, opt, x) for x in al])
    res.sort(key=lambda z: z[0])
    L, mask, _ = select_lines()
    NA = nad()
    o = {k: np.array([r[1][k] for r in res]) for k in res[0][1] if k not in ('lpb', 'top_lp')}
    o['lpb'] = res[0][1]['lpb']
    o['alpha'] = al
    o['lam_mask'] = L[mask, 0]
    o['lam_na'] = NA[:, 0]
    o['UW'], o['UNA'] = UW, UNA
    o['opt'] = repr({k: v for k, v in opt.items() if k != 'mw'})
    if tag == 'template':
        # mask weights: windless, rotationless line depth at the line centre, all epochs
        dep = (o['Afe'][:, :, :, len(UW)//2] - o['Afe_c'])[:, :, :].sum((0, 1))
        o['mw'] = dep/dep.max()
    os.makedirs(W + 'out', exist_ok=True)
    np.savez(W + 'out/%s_%s.npz' % (run, tag), **o)
    print('%s %s done in %.0f s, opt %s' % (run, tag, time.time() - t0, o['opt']))


def rhotest(run):
    """validation 7b: the old synth.py rho ds window average through chord_samples"""
    from scipy.spatial import cKDTree
    opt = parse([])
    S = load_state(run, opt, None)
    tree = cKDTree(S['X'])
    R = RUNS[run]
    th = np.deg2rad(np.arange(0, 360, 2.0) + 1.0)
    bb = np.linspace(S['rc'][0], S['rc'][-1], 120)
    s = np.linspace(-R['R1'], R['R1'], 801)
    lines = []
    for alpha in (0.0, -8.5, 8.5):
        a = np.deg2rad(alpha)
        ob = np.array([np.cos(a), -np.sin(a), 0.0])
        e1 = np.array([np.sin(a), np.cos(a), 0.0])
        e2 = np.array([0.0, 0.0, 1.0])
        acc = {}
        for t_ in th:
            limb = 'morning' if np.cos(t_) > 0 else 'evening'
            ci, pts, col, i0, f = chord_samples(S, tree, ob, e1, e2, t_, bb, s)
            rr = np.exp(ip(S['lrho'], col, i0, f))
            lpp = ip(S['lp'], col, i0, f)
            vobs = ((ip(S['vx'], col, i0, f) + OMEGA*(-pts[:, 1]))*ob[0]
                    + (ip(S['vy'], col, i0, f) + OMEGA*pts[:, 0])*ob[1])
            cnt = np.bincount(ci, minlength=len(bb))
            lpt = np.full(len(bb), -99.0)
            np.maximum.at(lpt, ci, lpp)
            m = (lpp >= -5) & (lpp <= -3)
            good = cnt[ci] >= 3
            ww = rr*(s[1] - s[0])*bb[ci]
            rsel = ((lpt >= -5) & (lpt <= -3))[ci] & good
            for wt, sel in (('ray', rsel), ('win', m & good)):
                q = acc.setdefault((limb, wt), np.zeros(2))
                q += [ww[sel].sum(), (ww[sel]*vobs[sel]).sum()]
        for k in sorted(acc):
            lines.append('alpha %+5.1f  Fe(1e-5..1e-3) %-8s w=%s  RV %+7.2f' % (
                alpha, k[0], k[1], -acc[k][1]/acc[k][0]/1e5))
        for wt in ('ray', 'win'):
            Wt = sum(acc[k][0] for k in acc if k[1] == wt)
            lines.append('alpha %+5.1f  Fe(1e-5..1e-3) BOTH     w=%s  RV %+7.2f' % (
                alpha, wt, -sum(acc[k][1] for k in acc if k[1] == wt)/Wt/1e5))
    os.makedirs(W + 'out', exist_ok=True)
    open(W + 'out/%s_rhotest.txt' % run, 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def top(run):
    """pressure and T of the top cell (centre, just below x1max) at the terminators"""
    S = load_state(run, parse([]), None)
    X = S['X']
    lat = np.degrees(np.arcsin(X[:, 2]))
    lon = np.degrees(np.arctan2(-X[:, 1], -X[:, 0]))
    for nm, l0 in (('morning (lon -90)', -90), ('evening (lon +90)', 90)):
        for lb in (30, 90):
            m = (abs(lon - l0) < 10) & (abs(lat) < lb)
            q = np.percentile(S['lp'][m, -1], [0, 16, 50, 84, 100])
            print('%s %s |lat|<%d: top-cell log10 p [bar] min/16/50/84/max %s ; T median %.0f'
                  ' K, fraction of T=200 K floor cells %.2f' % (
                      run, nm, lb, ' '.join('%.2f' % x for x in q),
                      np.median(S['T'][m, -1]), (S['T'][m, -1] < 201).mean()))


if __name__ == '__main__':
    if sys.argv[1] == 'top':
        top(sys.argv[2])
    if sys.argv[1] == 'run':
        main_run(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif sys.argv[1] == 'rhotest':
        rhotest(sys.argv[2])
