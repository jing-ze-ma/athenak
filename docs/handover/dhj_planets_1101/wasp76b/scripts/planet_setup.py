#!/usr/bin/env python3
"""planet_setup.py -- planet, orbit and radial-grid keys for the deep hot Jupiter pgen
(src/pgen/deep_hot_jupiter_rt.cpp) from literature parameters.

What the pgen needs and how it defines it
  <problem>/ap    radius where p = 250 bar (get_init_eos_arr: p0 = 250 bar at z = 0,
                  z = r - ap); we put the inner wall there: <mesh>/x1min = ap.
  <problem>/grav  g at ap; with grav_point_mass = true, g(r) = grav (ap/r)^2, i.e.
                  GM_p = grav ap^2.
  <problem>/Teq   T_eff,* sqrt(R_*/(2a)) (zero albedo, full redistribution; the pgen
                  uses T_irr = sqrt2 Teq, the substellar flux, and T_int(Teq)).
  <problem>/omega 2 pi / P_rot, P_rot = P_orb (tidally locked).
The observed radius R_p is the area-equivalent radius of the terminator on the p_ref
isobar (p_ref = the transit radius of the adopted band, transit.py in the setup
directory).  Integrating hydrostatic equilibrium with point-mass gravity through the IC
profile T(p) and the general-EOS table, from p_ref inward to 250 bar, gives ap and grav:
    d(1/r)/d ln p = (p/rho)/(G M_p)          (exact for g = G M_p / r^2).
With problem/rot_potential the isobars follow Phi = G M (1/ap - 1/r) - Omega^2 (r
sin theta)^2 / 2 (the pgen's TotPotAt): the IC column is the POLAR one and the equator
bulges; the polar radius of p_ref is chosen so that <r(theta)^2> over the terminator
(a great circle through both poles for a synchronous rotator) equals R_p^2.
<mesh>/x1max (default rule "rce"): 3 % of (r(p_top) - ap) inside r(p_top) of the IC
column, so no cell of the fresh start lies above the profile's top (nothing starts on a
floor); hotter / colder bounding columns (T x f above 0.01 bar) are reported.  Rule
"scaled" instead scales the reference planet's MEASURED 3-D extent in potential depth,
    G M (1/ap - 1/x1max) = s_top * G M_ref (1/ap_ref - 1/x1max_ref),
    s_top = Phi(p_top) / Phi_ref(p_top),  Phi(p) = int_p^{250 bar} dp'/rho.
Radial grids (default basis "ic", on the IC column, H = p/(rho g)):
  sparc : >= 3 cells per H where p > 1e-6 bar, free above; smallest even nx1.
  prod  : plan A of nx1time_0924/design256.py (1e-9:2.5, 1e-8:3, 1e-7:4, 1e-6..300 bar:5
          cells per H), 8 coefficients at nx1 = 256.
  dt ~ CFL 0.3 min(dr/(c_s + 1 km/s)) (radial; c_s = sqrt(Gamma_1 p/rho) of the column).
Basis "template" (the regression path) maps the reference planet's 3-D columns
(radres_0924/meas.npz) at equal pressure and applies sparc_0925/design_c8.py and
design256.py unchanged; with the reference planet it reproduces both grids.

Usage
  planet_setup.py --regress                      reference planet round trip (gates)
  planet_setup.py --regress --grids --grid_basis template     (also the grids)
  planet_setup.py --mp 1.170 --rp 1.7420 --pref 1.8e-4 --teff 6628 --rstar 1.461 \
      --a_au 0.02571 --porb 1.27492504 --profile ic.txt --eos eos_table.txt \
      --grids --save setup.npz --envfile grid.env      (WASP-121b, Sing+2024)
  eos_table.txt / eos_grid.txt: the general-EOS dumps (deepconv_0925 dump pgen).
Units: M_p [M_J], R_p [R_J] (equatorial Jupiter 7.1492e9 cm), p [bar], a [AU],
P [d], R_* [R_sun].  Needs numpy, scipy.
"""
import argparse
import os
import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.optimize import minimize

G_CGS = 6.67430e-8
MJ = 1.26686534e23/G_CGS       # IAU 2015 nominal GM_J / G
RJ = 7.1492e9                  # equatorial Jupiter radius (IAU 2015 nominal)
RSUN = 6.957e10
AU = 1.495978707e13
DAY = 86400.0
BAR = 1.0e6
P0 = 250.0*BAR                 # the pgen's anchor pressure at r = ap (hard-coded there)

# the reference (current) dhj planet: sparc_0925/sparc.athinput and its IC
REF = dict(grav=942.0, ap=9.44e9, x1max=2.0556e10, teq=2500.0, omega=2.06e-5,
           profile='/viper/ptmp2/jinma/ckrce_0925/ckrce_nq2_B.txt',
           npz='/viper/ptmp2/jinma/ckrce_0925/ckrce_nq2_B.npz',
           eos='/viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt',
           sparc_nx1=66,
           sparc_c=[-0.248682, -1.290656, -1.933528, 17.707345, -66.701191, 89.639593,
                    -30.359120, -11.805760],
           prod_c=[0.418017, 4.818585, -84.625150, 397.754052, -955.667782, 1270.732651,
                   -887.750373, 254.765233])
TEMPLATE = '/viper/ptmp2/jinma/radres_0924/meas.npz'


# ------------------------------------------------------------------------------ EOS
class EOSTable:
    """rho(p, T) from the code's general-EOS table dump (log rho x log T grid, column p),
    by bisection in log rho -- as ckrce_0925/rce.py."""

    def __init__(self, fn):
        with open(fn) as fh:
            fh.readline()
            fh.readline()
            nx, ny, x0, dx, y0, dy = [float(v) for v in fh.readline()[1:].split()]
        nx, ny = int(nx), int(ny)
        tab = np.loadtxt(fn, comments='#')
        self.lr = x0 + dx*np.arange(nx)
        lt = y0 + dy*np.arange(ny)
        self.lpr = RGI((self.lr, lt), tab[:, 1].reshape(ny, nx).T, bounds_error=False,
                       fill_value=None)
        # Gamma_1 = chi_rho + p chi_T^2/(rho T c_v) from the EOSTable::Eval dump next to
        # it (eos_grid.txt: lrho lT e p chi_rho chi_T c_v mu; as dclib.EOS)
        g = np.loadtxt(os.path.join(os.path.dirname(fn), 'eos_grid.txt'), comments='#',
                       usecols=(0, 1, 3, 4, 5, 6))
        glr, glt = np.unique(g[:, 0]), np.unique(g[:, 1])
        rho, T = 10**g[:, 0], 10**g[:, 1]
        g1 = g[:, 3] + g[:, 2]*g[:, 4]**2/(rho*T*g[:, 5])
        self.g1 = RGI((glr, glt), g1.reshape(len(glt), len(glr)).T, bounds_error=False,
                      fill_value=None)

    def gamma1(self, rho, T):
        return self.g1(np.stack([np.log10(rho), np.log10(T)], -1))

    def rho(self, p, T):
        lp = np.log10(np.atleast_1d(p)).ravel()
        lt = np.log10(np.atleast_1d(T)).ravel()*np.ones_like(lp)
        lo = np.full(lp.shape, self.lr[0] - 1)
        hi = np.full(lp.shape, self.lr[-1])
        for _ in range(45):
            mid = 0.5*(lo + hi)
            up = mid + self.lpr(np.stack([mid, lt], -1)) > lp
            hi = np.where(up, mid, hi)
            lo = np.where(up, lo, mid)
        return 10**(0.5*(lo + hi))


def read_profile(fn):
    a = np.loadtxt(fn, comments='#')
    o = np.argsort(a[:, 0])
    return a[o, 0], a[o, 1]              # log10 p [barye] ascending, T [K]


class Column:
    """Phi(p) = int_p^{p0} dp/rho = G M (1/ap - 1/r(p)) on a fine ln p grid."""

    def __init__(self, profile, eos, ptop=1e-9, n=20001, tfac=1.0):
        lpt, Tt = read_profile(profile)
        if tfac != 1.0:
            # a hotter / colder column: T x tfac above 0.01 bar, tapering (log p) to the
            # unchanged profile at 0.1 bar -- used only to bound x1max
            w = np.clip((np.log10(0.1*BAR) - lpt)/1.0, 0.0, 1.0)
            Tt = Tt*(1.0 + (tfac - 1.0)*w)
        if lpt[-1] < np.log10(P0):
            sys.exit('profile %s must reach 250 bar (as problem/ic_profile)' % profile)
        # above the table's top row T is held (the pgen's lookup clamps the same way)
        self.lp = np.linspace(min(np.log10(ptop*BAR), lpt[0]), lpt[-1], n)
        self.T = np.interp(self.lp, lpt, Tt)
        p = 10**self.lp
        self.rho = eos.rho(p, self.T)
        self.pr = p/self.rho                                  # p/rho [cm^2/s^2]
        self.cs = np.sqrt(eos.gamma1(self.rho, self.T)*self.pr)
        f = self.pr*np.log(10.0)                              # dPhi/dlog10 p
        cum = np.concatenate([[0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(self.lp))])
        self.phi = np.interp(np.log10(P0), self.lp, cum) - cum   # 0 at 250 bar

    def phi_at(self, pbar):
        return np.interp(np.log10(np.asarray(pbar)*BAR), self.lp, self.phi)

    def pr_at(self, pbar):
        return np.interp(np.log10(np.asarray(pbar)*BAR), self.lp, self.pr)

    def T_at(self, pbar):
        return np.interp(np.log10(np.asarray(pbar)*BAR), self.lp, self.T)


def anchor(col, gm, rp, pref):
    """ap, grav from R_p at p_ref."""
    inv_ap = 1.0/rp + col.phi_at(pref)/gm
    ap = 1.0/inv_ap
    return ap, gm/ap**2


def limb_radii(gm, ap, r_pole, om2, n=181):
    """Radius of the equipotential through r_pole (pole) at colatitude theta, with the
    pgen's rot_potential: Phi = GM(1/ap - 1/r) - om2 (r sin theta)^2 / 2 (TotPotAt)."""
    th = np.linspace(0.0, 0.5*np.pi, n)
    phi0 = gm*(1.0/ap - 1.0/r_pole)
    r = np.full(n, r_pole)
    for _ in range(50):
        f = gm*(1.0/ap - 1.0/r) - 0.5*om2*(r*np.sin(th))**2 - phi0
        df = gm/r**2 - om2*r*np.sin(th)**2
        r = r - f/df
    return th, r


def anchor_limb(col, gm, rp, pref, om2):
    """ap, grav and the polar / equatorial radius of the p_ref isobar such that the
    area-equivalent radius of the terminator (a great circle through both poles for a
    synchronous rotator seen at transit; sin^2 theta uniform along it),
    R_eff^2 = <r(theta)^2>, equals the observed R_p.  om2 = 0 gives r_pole = R_p."""
    r_pole = rp
    for _ in range(30):
        ap, grav = anchor(col, gm, r_pole, pref)
        th, r = limb_radii(gm, ap, r_pole, om2)
        reff = np.sqrt(np.trapezoid(r**2, th)/th[-1])
        r_pole *= rp/reff
    ap, grav = anchor(col, gm, r_pole, pref)
    th, r = limb_radii(gm, ap, r_pole, om2)
    return ap, grav, r_pole, r[-1], np.sqrt(np.trapezoid(r**2, th)/th[-1])


def r_of_p(col, gm, ap, pbar):
    return 1.0/(1.0/ap - col.phi_at(pbar)/gm)


def code_ic_fill(col, eos, gm, ap, x1min, x1max, pref, N=10000):
    """Emulate get_init_eos_arr (forward Euler in ln p on N uniform z steps from ap,
    zmax = 1.1 (x1max - x1min), point-mass g at z(n)) and return r(p_ref)."""
    dz = 1.1*(x1max - x1min)/N
    z = dz*np.arange(N)
    lnp = np.empty(N)
    lnp[0] = np.log(P0)
    lpc = col.lp*np.log(10.0)
    for n in range(N - 1):
        p = np.exp(lnp[n])
        T = np.interp(lnp[n], lpc, col.T)
        rho = eos.rho(p, T)[0]
        lnp[n + 1] = lnp[n] - gm/(ap + z[n])**2*dz*rho/p
        if lnp[n + 1] < np.log(pref*BAR) - 1.0:
            break
    m = n + 2
    return ap + np.interp(-np.log(pref*BAR), -lnp[:m], z[:m])


def p_of_r(col, gm, ap, r):
    """pressure [bar] at radius r on a column anchored at ap (inverse of r_of_p)."""
    phi = gm*(1.0/ap - 1.0/np.asarray(r))
    lp = np.interp(-phi, -col.phi, col.lp)
    # beyond the column top: isothermal continuation, ln p falls by dPhi/(p/rho)_top
    lp = np.where(phi > col.phi[0], col.lp[0] - (phi - col.phi[0])/col.pr[0]/np.log(10.0),
                  lp)
    return 10**lp/BAR


def x1max_rce(col, gm, ap, ptop, margin, cols_bound=()):
    """x1max = ap + (1 - margin)(r(ptop) - ap) on the IC column itself, so that every
    cell of the fresh start lies below the profile's p_top (nothing starts on a floor).
    Returns x1max and p(x1max) on each bounding column (hot/cold variants)."""
    rt = r_of_p(col, gm, ap, ptop)
    x1max = ap + (1.0 - margin)*(rt - ap)
    return x1max, [p_of_r(c, gm, ap, x1max) for c in cols_bound]


def x1max_scaled(col, colref, gm, ap, gmref, apref, x1maxref, ptop):
    s_top = col.phi_at(ptop)/colref.phi_at(ptop)
    dphi = s_top*gmref*(1.0/apref - 1.0/x1maxref)
    return 1.0/(1.0/ap - dphi/gm), s_top


# --------------------------------------------------------------- grid design (verbatim
# recipes of design_lit.py / design_c8.py / design256.py on a mapped template)
TODAY = np.array([-0.068392, -2.191487, 2.464818, -1.366698])
A320 = np.array([0.792150, -4.844414, 4.828881, -0.707535])
A288 = np.array(REF['prod_c'])
C0 = np.array([-0.254526, -0.100163, -15.776475, 86.273189, -231.670141, 283.094593,
               -131.562381, 4.460242])
XS = np.linspace(0, 1, 801)
REG = [('top 1e-9..1e-6', 1e-9, 1e-6), ('upper 1e-6..1e-3', 1e-6, 1e-3),
       ('phot/jet-base 1e-3..1e-1', 1e-3, 1e-1), ('jet 1e-1..10', 1e-1, 10),
       ('deep 10..300', 10, 300)]
PLAN_A = [(1e-9, 2.5), (1e-8, 3), (1e-7, 4), (1e-6, 5), (300, 5)]


def poly_u(c, xi):
    u = xi.copy()
    for k, ck in enumerate(c, start=1):
        u = u + ck*xi**k*(1 - xi)
    return u


class Design:
    def __init__(self, rc, R0, R1):
        self.rc, self.R0, self.R1 = rc, R0, R1

    def need_of(self, c, dens):
        u = poly_u(c, XS)
        du = np.gradient(u, XS)
        return np.max(np.interp(self.R0 + (self.R1 - self.R0)*u, self.rc, dens)*du
                      * (self.R1 - self.R0))

    def fit_poly(self, dens, ncoef, starts):
        R0, R1, rc = self.R0, self.R1, self.rc

        def worst(c):
            u = poly_u(c, XS)
            du = np.gradient(u, XS)
            if du.min() <= 0.02:
                return 1e6*(1 + 0.02 - du.min())
            q = np.interp(R0 + (R1 - R0)*u, rc, dens)*du*(R1 - R0)
            return np.log(np.sum(q**16)/len(q))/16
        best = None
        for x0 in starts:
            x0 = np.array(list(x0[:ncoef]) + [0]*(ncoef - len(x0)))
            for _ in range(3):
                r = minimize(worst, x0, method='Nelder-Mead',
                             options={'maxiter': 40000, 'maxfev': 40000, 'xatol': 1e-6,
                                      'fatol': 1e-9})
                x0 = r.x
            if best is None or r.fun < best.fun:
                best = r
        return best.x

    def grid(self, c, n):
        xi = np.linspace(0, 1, n + 1)
        re = self.R0 + (self.R1 - self.R0)*poly_u(c, xi)
        dr = np.diff(re)
        q = re[:-1]/re[1:]
        rcn = 0.25*(q*q + 1)/((q*q + q + 1)/3)*(re[1:] + re[:-1])
        return re, rcn, dr

    def code_faces(self, c, n):
        """x1f exactly as the Mesh builds it (LeftEdgeX + ApplyRStretch)."""
        R0, R1 = self.R0, self.R1
        f = np.empty(n + 1)
        for i in range(n + 1):
            x = float(i)/float(n)
            r = (x*R1 - x*R0) - (0.5*R1 - 0.5*R0) + (0.5*R0 + 0.5*R1)
            xi = (r - R0)/(R1 - R0)
            u, xik = xi, xi
            for k in range(len(c)):
                u += c[k]*xik*(1.0 - xi)
                xik *= xi
            f[i] = R0 + (R1 - R0)*u
        return f


def cph_table(rc_t, rc, dr, P, H, CZ, sub=slice(None)):
    nH = np.exp(np.array([np.interp(rc, rc_t, x) for x in np.log(H[sub])]))
    pb = np.exp(np.array([np.interp(rc, rc_t, x) for x in np.log(P[sub])]))/BAR
    x = nH/dr[None, :]
    out = {}
    for nm, lo, hi in REG:
        for side, m in (('day', CZ[sub] > 0.5), ('night', CZ[sub] < -0.5)):
            sel = (pb >= lo) & (pb < hi) & m[:, None]
            out[(nm, side)] = ((np.percentile(x[sel], 10), np.median(x[sel])) if sel.any()
                               else (np.nan, np.nan))
    return out


def fmt_table(t):
    return ' '.join('%-24s' % ('%.1f/%.1f | %.1f/%.1f'
                               % (t[(r[0], 'day')] + t[(r[0], 'night')])) for r in REG)


class Template:
    """radres_0924 columns mapped to (gm, ap, x1max, col) at equal pressure."""

    def __init__(self, fn, colref, col, gmref, apref, gm, ap, x1max, tags):
        Z = np.load(fn, allow_pickle=True)
        rc0, re0 = Z['rc'], Z['re']
        self.R0, self.R1 = ap, x1max
        self.rc = ap + (x1max - ap)*(rc0 - re0[0])/(re0[-1] - re0[0])
        d = {k: [] for k in ('P', 'Hp', 'Hhse', 'CZ', 'S', 'F')}
        dphi_old = gmref*(1.0/apref - 1.0/rc0)
        lp_hi = np.log10(P0/BAR)
        for t in tags:
            P = Z[t + '_p'].astype(float)
            pb = P/BAR
            pbc = np.clip(pb, 10**col.lp[0]/BAR, 10**lp_hi*0.999)
            # near the anchor Phi -> 0 on both columns: use the local ratio q there
            s = np.where(pb >= 150.0, col.pr_at(pbc)/colref.pr_at(pbc),
                         col.phi_at(pbc)/np.maximum(colref.phi_at(pbc), 1e-300))
            q = col.pr_at(pbc)/colref.pr_at(pbc)
            rnew = 1.0/(1.0/ap - s*dphi_old[None, :]/gm)
            jac = np.gradient(rnew, axis=1)/np.gradient(rc0)[None, :]
            gold = gmref/rc0[None, :]**2
            gnew = gm/rnew**2
            Hp = Z[t + '_Hp'].astype(float)*jac
            Hh = Z[t + '_Hhse'].astype(float)*q*gold/gnew
            S = (np.abs(Z[t + '_w'].astype(float))
                 + Z[t + '_fast'].astype(float))*np.sqrt(q)
            F = Z[t + '_fast'].astype(float)*np.sqrt(q)
            # resample every column on the common grid rc
            for nm, v in (('P', P), ('Hp', Hp), ('Hhse', Hh), ('S', S), ('F', F)):
                d[nm].append(np.array([np.exp(np.interp(self.rc, rnew[i], np.log(v[i])))
                                       for i in range(v.shape[0])]))
            d['CZ'].append(Z[t + '_cosz'].astype(float))
        self.d = {k: np.concatenate(v) for k, v in d.items()}
        self.ntag = len(tags)


def design_sparc(tp, verbose=True):
    """sparc_0925/design_lit.py (4c, Hmed/3) -> design_c8.py (8c, day p10 H/2.8)."""
    P, H, CZ, F = tp.d['P'], tp.d['Hp'], tp.d['CZ'], tp.d['F']
    ds = Design(tp.rc, tp.R0, tp.R1)
    fmax = F.max(0)
    low = np.median(P, 0)/BAR > 1e-6
    cS = ds.fit_poly(np.where(low, 3.0/np.median(H, 0), 0.0), 4,
                     [TODAY, A320, np.zeros(4)])
    day = CZ > 0.5
    dens = np.where(np.median(P[day], 0)/BAR > 1e-6,
                    2.8/np.percentile(H[day], 10, axis=0), 0.0)
    c8 = ds.fit_poly(dens, 8, [C0, A288, np.zeros(8), cS])
    rows, best = [], None
    for n in range(40, 121, 2):
        re, rc, dr = ds.grid(c8, n)
        t = cph_table(tp.rc, rc, dr, P, H, CZ)
        ok = min(t[(REG[i][0], 'day')][0] for i in (1, 3, 4)) >= 2.8
        dt = 0.184*np.min(dr/np.interp(rc, tp.rc, fmax))
        rows.append((n, ok, dr.min(), dr.max(), dt, t))
        if ok and best is None:
            best = n
        if best is not None and n >= best + 4:
            break
    n = best + (best % 2)
    if verbose:
        print('## SPARC-like grid (design_c8 recipe): need nx1 = %.1f -> nx1 = %d'
              % (ds.need_of(c8, dens), n))
        print('   c = ' + ', '.join('%.6f' % v for v in c8))
        print('%-6s %4s %9s %9s %7s ' % ('', 'nx1', 'dr_min', 'dr_max', 'dt_est')
              + ' '.join('%-24s' % r[0] for r in REG))
        for (m, ok, a, b, dt, t) in rows:
            print('%-6s %4d %9.3g %9.3g %7.2f %s' % ('ok' if ok else '', m, a, b, dt,
                                                     fmt_table(t)))
    f = ds.code_faces(c8, n)
    u = poly_u(c8, np.linspace(0, 1, 20001))
    mono = np.gradient(u, np.linspace(0, 1, 20001)).min()
    return dict(nx1=n, c=c8, faces=f, dudxi_min=mono,
                dt=[r[4] for r in rows if r[0] == n][0])


def design_prod(tp, n=256, verbose=True, dt_h_ref=10.45, cfl=0.3):
    """nx1time_0924/design256.py plan A, 8 coefficients at nx1 = n."""
    P, H, CZ, S = tp.d['P'], tp.d['Hhse'], tp.d['CZ'], tp.d['S']
    ds = Design(tp.rc, tp.R0, tp.R1)
    lp = np.log10([a for a, b in PLAN_A])
    nv = np.array([b for a, b in PLAN_A], float)
    dens = np.percentile(np.interp(np.log10(P/BAR), lp, nv)/H, 95, axis=0)
    c8 = ds.fit_poly(dens, 8, [TODAY, A320])
    re, rc, dr = ds.grid(c8, n)
    smax = S.max(0)
    d1 = cfl*np.min(dr/np.interp(rc, tp.rc, smax))
    f = n/ds.need_of(c8, dens)
    t = cph_table(tp.rc, rc, dr, P, H, CZ, sub=slice(0, None, 7))
    if verbose:
        print('## production grid (design256 plan A, 8c) at nx1 = %d: f = %.2f of the '
              'plan' % (n, f))
        print('   c = ' + ', '.join('%.6f' % v for v in c8))
        print('   dr_min %.3g dr_max %.3g  dt_radial %.2f s  (1024 horizontal limit, '
              'scaled: %.2f s)' % (dr.min(), dr.max(), d1, dt_h_ref))
        print('   cells/H p10/median day | night: ' + fmt_table(t))
    return dict(nx1=n, c=c8, f=f, dt=min(d1, dt_h_ref), dt_radial=d1,
                faces=ds.code_faces(c8, n))


# ---------------------------------------------------------------------------- main
def design_ic(col, gm, ap, x1max, cols_bound=(), nprod=256, vwind=1.0e5, cfl=0.3,
              verbose=True):
    """Grids designed on the IC column itself (the state every column starts from):
    sparc: >= 3 cells per H for p > 1e-6 bar, free above (the user's rule); smallest even
    nx1.  prod: plan A (1e-9:2.5, 1e-8:3, 1e-7:4, 1e-6..300 bar:5 cells per H) at nprod.
    H = p/(rho g).  dt = cfl min(dr/(c_s + v_wind)) with c_s = sqrt(Gamma_1 p/rho) of the
    column and v_wind an assumed peak RADIAL speed (the horizontal limit, ~100 s at C32,
    is not included)."""
    rc = np.linspace(ap, x1max, 1201)
    ds = Design(rc, ap, x1max)

    def colq(c):
        pb = p_of_r(c, gm, ap, rc)
        lp = np.log10(np.clip(pb*BAR, 10**c.lp[0], 10**c.lp[-1]))
        H = np.interp(lp, c.lp, c.pr)/(gm/rc**2)
        return pb, H, np.interp(lp, c.lp, c.cs)
    pb, H, cs = colq(col)
    dens = np.where(pb > 1e-6, 3.0/H, 0.0)
    c4 = ds.fit_poly(dens, 4, [TODAY, A320, np.zeros(4)])
    c8 = ds.fit_poly(dens, 8, [C0, A288, np.zeros(8), c4])
    nneed = ds.need_of(c8, dens)
    ns = int(np.ceil(nneed))
    ns += ns % 2
    lpA = np.log10([x for x, y in PLAN_A])
    nvA = np.array([y for x, y in PLAN_A], float)
    densA = np.interp(np.log10(pb), lpA, nvA)/H
    cA = ds.fit_poly(densA, 8, [TODAY, A320, A288])
    fA = nprod/ds.need_of(cA, densA)
    res = {}
    for nm, c, n in (('sparc', c8, ns), ('prod', cA, nprod)):
        re, rcn, dr = ds.grid(c, n)
        dt = cfl*np.min(dr/(np.interp(rcn, rc, cs) + vwind))
        rows = []
        for cc, lab in [(col, 'IC')] + [(b, 'bound %d' % i)
                                        for i, b in enumerate(cols_bound)]:
            pbb, Hb, _ = colq(cc)
            x = np.interp(rcn, rc, Hb)/dr
            pc = np.interp(rcn, rc, pbb)
            rows.append((lab, [(np.min(x[(pc >= lo) & (pc < hi)])
                                if np.any((pc >= lo) & (pc < hi)) else np.nan)
                               for _, lo, hi in REG]))
        res[nm] = dict(nx1=n, c=c, dt=dt, dr=dr, faces=ds.code_faces(c, n), rows=rows)
        if verbose:
            why = (' (need %.1f for 3/H at p > 1e-6 bar)' % nneed if nm == 'sparc'
                   else ' (f = %.2f of plan A)' % fA)
            print('## %s grid on the IC column: nx1 = %d%s  dr %.3g..%.3g cm  '
                  'dt ~ %.2f s (c_s + %.0f km/s, CFL %.1f)'
                  % (nm, n, why, dr.min(), dr.max(), dt, vwind/1e5, cfl))
            print('   c = ' + ', '.join('%.6f' % v for v in c))
            print('   min cells/H per band  '
                  + ' '.join('%-10s' % r[0][:10] for r in REG))
            for lab, v in rows:
                print('   %-20s ' % lab + ' '.join('%-10.2f' % x for x in v))
    return res['sparc'], res['prod']


def derive(a, col, colref, eos, verbose=True):
    gm = a.mp*MJ*G_CGS
    rp = a.rp*RJ
    om2 = (2*np.pi/(a.porb*DAY))**2 if (a.porb and a.rot_potential) else 0.0
    ap, grav, r_pole, r_eq, r_eff = anchor_limb(col, gm, rp, a.pref, om2)
    gmref = REF['grav']*REF['ap']**2
    x1max_s, s_top = x1max_scaled(col, colref, gm, ap, gmref, REF['ap'], REF['x1max'],
                                  a.ptop)
    bounds = []
    if a.x1max_rule == 'rce':
        cb = [Column(a.profile, eos, a.ptop, tfac=f) for f in a.tbound] \
            if getattr(a, 'profile', None) and a.tbound else []
        x1max, pb = x1max_rce(col, gm, ap, a.ptop, a.x1max_margin, cb)
        bounds = list(zip(a.tbound, pb, [r_of_p(c, gm, ap, a.ptop) for c in cb]))
    else:
        x1max = x1max_s
    teq = a.teff*np.sqrt(a.rstar*RSUN/(2*a.a_au*AU)) if a.teff else None
    omega = 2*np.pi/(a.porb*DAY) if a.porb else None
    r_ic = code_ic_fill(col, eos, gm, ap, ap, x1max, a.pref)     # the polar column
    out = dict(gm=gm, rp=rp, ap=ap, grav=grav, x1min=ap, x1max=x1max, s_top=s_top,
               r_pole=r_pole, r_eq=r_eq, r_eff=r_eff,
               x1max_scaled=x1max_s, p_x1max=p_of_r(col, gm, ap, x1max),
               teq=teq, omega=omega, r_ic=r_ic, r_top1d=r_of_p(col, gm, ap, a.ptop))
    if verbose:
        print('# M_p = %.4f M_J (GM = %.6e), R_p = %.4f R_J = %.6e cm at p_ref = %.3e bar'
              % (a.mp, gm, a.rp, rp, a.pref))
        if teq:
            print('# Teq = Teff sqrt(R*/2a) = %.1f K   omega = 2pi/P = %.6e s^-1 '
                  '(P = %.1f s)' % (teq, omega, a.porb*DAY))
        print('p_ref isobar (rot_potential, omega^2 = %.3e): pole %.6e, equator %.6e, '
              'terminator area-equivalent %.6e cm (= R_p)' % (om2, r_pole, r_eq, r_eff))
        print('ap = x1min = %.6e cm   grav = %.3f cm/s^2   g(R_p) = %.3f'
              % (ap, grav, gm/rp**2))
        print('x1max = %.6e cm (rule %s): p(x1max) on the IC column = %.3e bar; '
              '1-D r(%.0e bar) = %.6e; scaled-extent rule would give %.6e (s_top %.4f)'
              % (x1max, a.x1max_rule, out['p_x1max'], a.ptop, out['r_top1d'], x1max_s,
                 s_top))
        for f, pb, rt in bounds:
            print('   bound: T x %.2f above 0.01 bar -> r(p_top) = %.5e, '
                  'p(x1max) = %.3e bar' % (f, rt, pb))
        for pb in (100, 10, 1, 0.1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9):
            print('   r(%7.0e bar) = %.5e cm  T = %6.0f K  H = %.3e cm' %
                  (pb, r_of_p(col, gm, ap, pb), col.T_at(pb),
                   col.pr_at(pb)/(gm/r_of_p(col, gm, ap, pb)**2)))
        print('code IC fill (get_init_eos_arr emulation, polar column): r(p_ref) = %.6e '
              'cm, vs r_pole %+.3e (%+.4f %%)'
              % (r_ic, r_ic - r_pole, 100*(r_ic/r_pole - 1)))
    return out


def regress(a):
    """Round trip on the reference planet: R_p at 1 mbar from the RCE's own radii."""
    eos = EOSTable(REF['eos'])
    col = Column(REF['profile'], eos, a.ptop)
    Z = np.load(REF['npz'])
    rf = Z['rf']
    rcc = 0.5*(rf[1:] + rf[:-1])
    ok = True
    for pref in (1e-3, 1e-2, 1.0):
        rp = np.interp(np.log10(pref*BAR), Z['lpc'][::-1], rcc[::-1])
        gm = REF['grav']*REF['ap']**2
        a.mp, a.rp, a.pref = gm/G_CGS/MJ, rp/RJ, pref
        a.x1max_rule = 'scaled'
        a.porb = None
        o = derive(a, col, col, eos, verbose=False)
        e1 = o['x1min']/REF['ap'] - 1
        e2 = o['grav']/REF['grav'] - 1
        e3 = o['x1max']/REF['x1max'] - 1
        print('regress p_ref %.0e bar: R_p %.6e -> x1min %.6e (%+.2e) grav %.3f (%+.2e) '
              'x1max %.6e (%+.2e) IC-fill r(p_ref) %+.2e'
              % (pref, rp, o['x1min'], e1, o['grav'], e2, o['x1max'], e3,
                 o['r_ic']/rp - 1))
        ok &= abs(e1) < 1e-3 and abs(e2) < 2e-3 and abs(e3) < 1e-3
    print('Teq(old) check: omega %.4e -> P %.4e s' % (REF['omega'], 2*np.pi/REF['omega']))
    if a.grids:
        tp = Template(a.template, col, col, gm, REF['ap'], gm, REF['ap'], REF['x1max'],
                      ['s4'])
        sp = design_sparc(tp)
        print('regress sparc: nx1 %d (ref %d), max |c - c_ref| = %.2e'
              % (sp['nx1'], REF['sparc_nx1'], np.max(np.abs(sp['c'] - REF['sparc_c']))))
        tp = Template(a.template, col, col, gm, REF['ap'], gm, REF['ap'], REF['x1max'],
                      ['s0', 's1', 's2', 's3', 's4'])
        pr = design_prod(tp)
        print('regress prod: max |c - c_ref| = %.2e'
              % np.max(np.abs(pr['c'] - REF['prod_c'])))
    print('REGRESSION', 'PASS' if ok else 'FAIL')


def main():
    ap_ = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap_.add_argument('--regress', action='store_true')
    ap_.add_argument('--mp', type=float, help='M_p [M_J]')
    ap_.add_argument('--rp', type=float, help='R_p [R_J] at p_ref')
    ap_.add_argument('--pref', type=float, default=1e-3, help='p_ref [bar]')
    ap_.add_argument('--ptop', type=float, default=1e-9, help='p_top [bar]')
    ap_.add_argument('--teff', type=float, help='host T_eff [K]')
    ap_.add_argument('--rstar', type=float, help='R_* [R_sun]')
    ap_.add_argument('--a_au', type=float, help='a [AU]')
    ap_.add_argument('--porb', type=float, help='P_orb = P_rot [d]')
    ap_.add_argument('--no_rot_potential', dest='rot_potential', action='store_false',
                     help='anchor R_p on a spherical isobar '
                          '(problem/rot_potential false)')
    ap_.add_argument('--profile', help='IC profile (log10 p [barye], T [K])')
    ap_.add_argument('--eos', default=REF['eos'], help='general-EOS table dump')
    ap_.add_argument('--x1max_rule', choices=['rce', 'scaled'], default='rce',
                     help='rce: just inside r(p_top) of the IC column (default); scaled: '
                          'the reference planet\'s measured 3-D extent, scaled')
    ap_.add_argument('--x1max_margin', type=float, default=0.03,
                     help='rce rule: fraction of (r(p_top) - ap) kept inside')
    ap_.add_argument('--tbound', type=float, nargs='*', default=[0.6, 1.3],
                     help='T factors (above 0.01 bar) of the bounding cold/hot columns')
    ap_.add_argument('--grids', action='store_true', help='also design the radial grids')
    ap_.add_argument('--grid_basis', choices=['ic', 'template'], default='ic',
                     help='ic: design on the IC column (default); template: the '
                          'reference planet\'s 3-D columns mapped at equal pressure '
                          '(regression path)')
    ap_.add_argument('--template', default=TEMPLATE)
    ap_.add_argument('--nprod', type=int, default=256)
    ap_.add_argument('--save', help='write the derived numbers to this .npz')
    ap_.add_argument('--envfile',
                     help='write the grid overrides (G_SPARC, G_PROD, G) here')
    a = ap_.parse_args()
    if a.regress:
        return regress(a)
    eos = EOSTable(a.eos)
    col = Column(a.profile, eos, a.ptop)
    colref = Column(REF['profile'], EOSTable(REF['eos']), a.ptop)
    o = derive(a, col, colref, eos)
    if a.grids and a.grid_basis == 'ic':
        cb = [Column(a.profile, eos, a.ptop, tfac=f) for f in a.tbound]
        sp, pr = design_ic(col, o['gm'], o['ap'], o['x1max'], cb, a.nprod)
    elif a.grids:
        gmref = REF['grav']*REF['ap']**2
        tp = Template(a.template, colref, col, gmref, REF['ap'], o['gm'], o['ap'],
                      o['x1max'], ['s4'])
        sp = design_sparc(tp)
        tp = Template(a.template, colref, col, gmref, REF['ap'], o['gm'], o['ap'],
                      o['x1max'], ['s0', 's1', 's2', 's3', 's4'])
        # horizontal limit of the 1024 grid: scale 10.45 s by r and 1/sound speed
        pl = np.logspace(-6, 2, 50)
        qm = np.median(col.pr_at(pl)/colref.pr_at(pl))
        pr = design_prod(tp, a.nprod, dt_h_ref=10.45*o['ap']/REF['ap']/np.sqrt(qm))
    if a.grids:
        o.update(sparc_nx1=sp['nx1'], sparc_c=sp['c'], sparc_dt=sp['dt'],
                 sparc_faces=sp['faces'], prod_nx1=pr['nx1'], prod_c=pr['c'],
                 prod_dt=pr['dt'], prod_faces=pr['faces'])
        env = []
        for nm, g in (('sparc', sp), ('prod', pr)):
            keys = ('mesh/nx1=%d meshblock/nx1=%d mesh/x1min=%.6e mesh/x1max=%.6e '
                    'mesh/use_grid_stretch_r_poly=true '
                    % (g['nx1'], g['nx1'], o['x1min'], o['x1max'])
                    + ' '.join('mesh/f_stretch_r_c%d=%.6f' % (k + 1, v)
                               for k, v in enumerate(g['c'])))
            print('%s keys: %s' % (nm, keys))
            env.append('G_%s="%s"' % (nm.upper(), keys))
        if a.envfile:
            with open(a.envfile, 'w') as fh:
                fh.write('# planet_setup.py grid overrides, basis %s (sparc: >= 3 '
                         'cells/H at p > 1e-6 bar; prod: plan A, 8c)\n' % a.grid_basis
                         + '\n'.join(env)
                         + '\nG="$G_SPARC"\n')
    teq = '%.1f' % o['teq'] if o['teq'] else '-'
    om = '%.6e' % o['omega'] if o['omega'] else '-'
    print('<problem> keys: ap = %.6e  grav = %.4f  Teq = %s  omega = %s   <mesh> x1min = '
          '%.6e  x1max = %.6e' % (o['ap'], o['grav'], teq, om, o['x1min'], o['x1max']))
    if a.save:
        np.savez(a.save, **{k: np.asarray(v) for k, v in o.items() if v is not None})


if __name__ == '__main__':
    main()
