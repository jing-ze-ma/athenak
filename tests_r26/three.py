"""tests_r26: 3-D shell diagnostics for the r22 settled wedges.  READ-ONLY on data.

One dump at a time; only 1-D profiles are kept.
"""
import os, sys
import numpy as np
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc  # noqa: E402

RS = 2.3717e11
TURN = 4705.0
LW = 4.0776e37
LSTAR = 2.3066e38
OMFRAC = 0.17678
ARAD = 7.5657332503e-15
CLIGHT = 2.99792458e10
SIGSB = 5.670374419e-5
KB = 1.380649e-16
MH = 1.6605390e-24
GM = None
CPOLY = [0.992525, -0.404821, -0.372255, -1.499759]
RAD_LO, RAD_HI = -11.0, -10.0
VCEIL = 1.0e8
EOSTBL = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
          'eos_table_box_w5.txt')
OPAC = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt'
BASE = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r22'


# ---------------------------------------------------------------- EOS table (e and p)
def _load_tbl(fn):
    with open(fn) as fh:
        ln = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in ln[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    return dict(nx=nx, ny=ny, xmin=xmin, dx=dx, ymin=ymin, dy=dy,
                E=d[:, 0].reshape(ny, nx), P=d[:, 1].reshape(ny, nx))


TB = _load_tbl(EOSTBL)


def _bil(F, ix, iy, u, v):
    return ((1-u)*(1-v)*F[iy, ix] + u*(1-v)*F[iy, ix+1]
            + (1-u)*v*F[iy+1, ix] + u*v*F[iy+1, ix+1])


def _idx(x, y):
    gx = (x - TB['xmin'])/TB['dx']
    gy = (y - TB['ymin'])/TB['dy']
    ix = np.clip(np.floor(gx).astype(np.int32), 0, TB['nx']-2)
    iy = np.clip(np.floor(gy).astype(np.int32), 0, TB['ny']-2)
    return ix, iy, np.clip(gx-ix, -3, 4), np.clip(gy-iy, -3, 4)


def weight(lrho):
    s = (lrho - RAD_LO)/(RAD_HI - RAD_LO)
    return np.where(s <= 0, 0.0, np.where(s >= 1, 1.0, s*s*(3-2*s)))


def temp_of(rho, eint):
    x = np.log10(np.maximum(rho, 1e-300))
    w = weight(x)
    lo = np.full(rho.shape, 2.5)
    hi = np.full(rho.shape, 7.5)
    for _ in range(50):
        mid = 0.5*(lo+hi)
        ix, iy, u, v = _idx(x, mid)
        f = rho*10.0**_bil(TB['E'], ix, iy, u, v) + w*ARAD*(10.0**mid)**4 - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    return 10.0**(0.5*(lo+hi))


def pgas_of(rho, T):
    x = np.log10(np.maximum(rho, 1e-300))
    ix, iy, u, v = _idx(x, np.log10(np.maximum(T, 1e-30)))
    return rho*10.0**_bil(TB['P'], ix, iy, u, v)


# ---------------------------------------------------------------- opacity
def _load_opac(fn):
    hdr, vals = [], []
    for ln in open(fn):
        (hdr if ln.startswith('#') else vals).append(ln)
    g = [h for h in hdr if h.strip().startswith('# 2')]
    nT, nD, lTmin, dlT, lDmin, dlD = [float(v) for v in g[-1].strip('# \n').split()]
    k = np.array([float(v) for v in vals]).reshape(int(nT), int(nD))
    return k, lTmin, dlT, int(nT), lDmin, dlD, int(nD)


OP = _load_opac(OPAC)


def kappa(T, rho):
    k, lTmin, dlT, nT, lDmin, dlD, nD = OP
    x = np.clip((np.log10(T) - lTmin)/dlT, 0, nT - 1.000001)
    y = np.clip((np.log10(np.maximum(rho, 1e-300)) - lDmin)/dlD, 0, nD - 1.000001)
    i0, j0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - i0, y - j0
    return 10.0**((1-fx)*(1-fy)*k[i0, j0] + fx*(1-fy)*k[i0+1, j0]
                  + (1-fx)*fy*k[i0, j0+1] + fx*fy*k[i0+1, j0+1])


# ---------------------------------------------------------------- grid / loader
def stretch(xi):
    u = xi.copy()
    xik = xi.copy()
    for kk in range(1, 5):
        u = u + CPOLY[kk-1]*xik*(1.0-xi)
        xik = xik*xi
    return u


def load(run, idx, varlist=('dens', 'velx', 'vely', 'velz', 'eint')):
    d = bc.read_binary(os.path.join(BASE, run, 'bin', 'he4.hydro_w.%05d.bin' % idx))
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    nj, nk = d['nx2_mb'], d['nx3_mb']
    out = {}
    for nm in varlist:
        g = np.empty((n3, n2, n1), dtype=np.float64)
        src = d['mb_data'][nm]
        for m in range(d['n_mbs']):
            _, lx2, lx3, _ = d['mb_logical'][m]
            g[lx3*nk:(lx3+1)*nk, lx2*nj:(lx2+1)*nj, :] = src[m]
        out[nm] = g
    rf = d['x1min'] + (d['x1max']-d['x1min'])*stretch(np.linspace(0.0, 1.0, n1+1))
    rl, rr = rf[:-1], rf[1:]
    q = rl/rr
    rc = 0.25*(q*q+1.0)/((1./3.)*(q*q+q+1.))*(rr+rl)
    dr = rr-rl
    thf = np.linspace(d['x2min'], d['x2max'], n2+1)
    wj = np.cos(thf[:-1]) - np.cos(thf[1:])
    dphi = (d['x3max']-d['x3min'])/n3
    dV = ((rr**3-rl**3)/3.0)[None, None, :]*wj[None, :, None]*dphi
    meta = dict(t=d['time'], rc=rc, dr=dr, rl=rl, rr=rr, wj=wj, dphi=dphi,
                dV=dV, n1=n1, n2=n2, n3=n3)
    del d
    return meta, out


def shmean(F, wj):
    """solid-angle weighted horizontal mean -> (n1,)"""
    return np.tensordot(wj, F.sum(axis=0), axes=([0], [0]))/(wj.sum()*F.shape[0])


def profile(run, idx, gmass=None):
    """all shell diagnostics for one dump."""
    M, g = load(run, idx)
    rc, dr, wj, dV = M['rc'], M['dr'], M['wj'], M['dV']
    rho, vr, vt, vp, ei = (g['dens'], g['velx'], g['vely'], g['velz'], g['eint'])
    del g
    n3 = rho.shape[0]
    T = np.empty_like(rho)
    for a in range(0, rho.shape[2], 16):
        b = min(a+16, rho.shape[2])
        T[:, :, a:b] = temp_of(rho[:, :, a:b], ei[:, :, a:b])
    lr = np.log10(np.maximum(rho, 1e-300))
    w = weight(lr)
    pg = pgas_of(rho, T)
    pr = w*ARAD*T**4/3.0
    p = pg + pr
    mu = KB*T/(MH*np.maximum(pg/rho, 1e-30))
    beta = pg/np.maximum(p, 1e-300)
    gam1 = beta*(5./3.) + (1-beta)*(4./3.)
    cs = np.sqrt(gam1*p/np.maximum(rho, 1e-300))
    kap = kappa(T, rho)
    kr = kap*rho
    # specific entropy proxy (gas ideal with table mu + radiation)
    sg = (KB/(mu*MH))*np.log(np.maximum(T, 1e-30)**1.5/np.maximum(rho, 1e-300))
    sr = 4.0*w*ARAD*T**3/(3.0*np.maximum(rho, 1e-300))
    s = sg + sr
    # ---- shell means
    dm = rho*dV
    W = dV.sum(axis=(0, 1))
    Msh = dm.sum(axis=(0, 1))
    rm = dm.sum(axis=(0, 1))/W                     # = volume-weighted <rho>
    def vw(F):   # volume weighted
        return (F*dV).sum(axis=(0, 1))/W
    def mw(F):   # mass weighted
        return (F*dm).sum(axis=(0, 1))/Msh
    Tm, pgm, prm, pm_ = vw(T), vw(pg), vw(pr), vw(p)
    sm = mw(s)
    v1m = mw(vr)
    vr_rms = np.sqrt(mw(vr**2))
    vr_rms_p = np.sqrt(mw((vr - v1m[None, None, :])**2))
    vh_rms = np.sqrt(mw(vt**2 + vp**2))
    csm = mw(cs)
    mach = np.sqrt(mw((vr**2 + vt**2 + vp**2)/cs**2))
    # fluxes
    dTdr = np.gradient(T, rc, axis=2)
    Fcol = -(4.0*ARAD*CLIGHT/3.0)*T**3/kr*dTdr
    Fdiff = shmean(Fcol, wj)
    Tmm = shmean(T, wj)
    rmm = shmean(rho, wj)
    Fmean = -(4.0*ARAD*CLIGHT/3.0)*Tmm**3/(kappa(Tmm, rmm)*rmm)*np.gradient(Tmm, rc)
    # enthalpy + kinetic flux (mass-flux-weighted fluctuations)
    h = (ei + p)/rho
    hbar = mw(h)
    Fenth = mw(rho*vr*(h - hbar[None, None, :]))*0.0 + \
        (rho*vr*(h - hbar[None, None, :])*dV).sum(axis=(0, 1))/W
    Fkin = (0.5*rho*vr*(vr**2+vt**2+vp**2)*dV).sum(axis=(0, 1))/W
    Fadv = (rho*vr*h*dV).sum(axis=(0, 1))/W
    # porosity
    ikrm = shmean(1.0/kr, wj)
    krm = shmean(kr, wj)
    P1 = ikrm*krm
    rrms = np.sqrt(shmean((rho/rmm - 1.0)**2, wj))
    trms = np.sqrt(shmean((T/Tmm - 1.0)**2, wj))
    # tau from the top, with real kappa (shell mean of kappa*rho)
    tau = np.cumsum((krm*dr)[::-1])[::-1]
    # ceiling census
    vmag = np.sqrt(vr**2+vt**2+vp**2)
    nceil = (vmag > 0.99*VCEIL).sum(axis=(0, 1))
    nceil5 = (vmag > 0.5*VCEIL).sum(axis=(0, 1))
    ncell = rho.shape[0]*rho.shape[1]
    # dt census
    sig = 0.3*dr[None, None, :]/(np.abs(vr) + cs)
    sigt = 0.3*(rc[None, None, :]*(np.pi/2/rho.shape[1]))/(np.abs(vt) + cs)
    sigp = 0.3*(rc[None, None, :]*np.sin(np.pi/2)*(np.pi/2/n3))/(np.abs(vp) + cs)
    out = dict(t=M['t'], r=rc, dr=dr, rho=rm, rho_h=rmm, T=Tm, T_h=Tmm, pg=pgm, pr=prm,
               p=pm_, beta=pgm/np.maximum(pm_, 1e-300), s=sm, v1=v1m,
               vr_rms=vr_rms, vr_rms_p=vr_rms_p, vh_rms=vh_rms, cs=csm, mach=mach,
               menc=np.cumsum(Msh), Msh=Msh, Mtot=Msh.sum(),
               Fdiff=Fdiff, Fmean=Fmean, Fenth=Fenth, Fkin=Fkin, Fadv=Fadv,
               P1=P1, rrms=rrms, trms=trms, tau=tau, kr=krm, kap=shmean(kap, wj),
               nceil=nceil, nceil5=nceil5, ncell=ncell,
               dtmin=np.array([sig.min(), sigt.min(), sigp.min()]),
               dtargmin=np.array([np.unravel_index(sig.argmin(), sig.shape)[2],
                                  np.unravel_index(sigt.argmin(), sigt.shape)[2],
                                  np.unravel_index(sigp.argmin(), sigp.shape)[2]]),
               Tmax=T.max(axis=(0, 1)), rhomin=rho.min(axis=(0, 1)),
               rhomax=rho.max(axis=(0, 1)),
               Lsh=4*np.pi*rc**2*OMFRAC)
    return out


def avg(profs):
    """time average of a list of profile dicts (element-wise on arrays)."""
    o = {}
    for k, v in profs[0].items():
        if isinstance(v, np.ndarray):
            o[k] = np.mean([p[k] for p in profs], axis=0)
        else:
            o[k] = np.mean([p[k] for p in profs])
    return o
