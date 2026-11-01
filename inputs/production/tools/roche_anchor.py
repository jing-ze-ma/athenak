#!/usr/bin/env python3
"""roche_anchor.py -- ap / x1min / grav / roche_a for deep_hot_jupiter_rt with the host-star
Roche tide (problem/roche) from literature parameters, for a planet that has never been run.

Sibling of planet_setup.py (same anchoring idea, extended to the full Roche potential).

What the pgen does (src/pgen/deep_hot_jupiter_rt.cpp, roche-dhj-1031 + w121hroche-1101)
  * one 1-D column: p(z) integrated in POINT-MASS gravity g = G M_p/(ap + z)^2 through the
    IC table T(p) and the EOS, p = 250 bar at z = 0 (r = ap); below the table bottom T(p)
    follows the adiabat of the bottom row (problem/ic_adiabat_ext, default on), above the
    table top T is held.  In potential form Phi_col(p) = int_p^{250 bar} dp/rho.
  * every cell samples that column at z_eff(Phi_tot): GM (1/ap - 1/(ap + z_eff)) = Phi_tot,
        Phi_tot = GM (1/ap - 1/r) - Omega^2 R_perp^2/2 + TidePot(r, mu),
        TidePot = -G M_*/d + G M_*/a + G M_* r mu/a^2,  d = |x - a e_s|,
    with G M_* = Omega^2 a^3 - G M_p (Kepler, as the pgen derives it from omega and
    roche_a), mu = cosine of the angle to the substellar point.
  So the p_ref isobar is the equipotential Phi_tot = Phi_col(p_ref).
Anchor: with G M_p fixed (M_p or g(R_p) from the literature), find ap such that the
area-equivalent radius of the TERMINATOR ring (the great circle mu = 0 through both poles,
as seen at mid transit), sqrt(<r(psi)^2>) uniform in the ring angle psi, equals R_p.
Then x1min = ap, grav = G M_p/ap^2 (the pgen's g at ap), roche_a = a.
Also reported: the p_ref surface in the four directions, L1, the inner-wall pressure and
temperature range (the wall sphere r = x1min cuts equipotentials; with the tide the
substellar/antistellar wall is much deeper than the pole), and whether it leaves the
IC table (then it is on the adiabat extension), and the deck keys.

EOS: a dhjcs-format general-EOS dump (log10 rho x log10 T grid; rows log10(e/rho)
log10(p/rho) mu), e.g. w121prod_0930/ana_3x_rot254/eos_3x.txt, made for any composition by
eos_dump_dhj (ana_3x_rot254/eos_dump_dhj.cpp: eos_dump_dhj out.txt X Y [M/H]).
The adiabat below the table uses de = (p/rho^2) d rho on that dump.

Usage (WASP-121b, the w121hroche_1101 deck):
  roche_anchor.py --mp 1.170 --rp 1.742 --pref 1.8e-4 --a_au 0.02571 --porb 1.27492504 \\
      --teff 6628 --rstar 1.461 --mstar 1.358 --profile ic_w121_3x_f5.txt \\
      --eos eos_3x.txt --x1max 2.1e10
  --no_roche: the rot_potential-only anchor (planet_setup.py's anchor_limb, for comparison).
Units: M_p [M_J], R_p [R_J] (7.1492e9 cm), g [cm s^-2], p [bar], a [AU], P [d],
R_* [R_sun], M_* [M_sun].  Needs numpy only.
"""
import argparse
import sys
import numpy as np

G_CGS = 6.67430e-8
GMJ = 1.26686534e23            # IAU 2015 nominal G M_J
RJ = 7.1492e9
RSUN = 6.957e10
GMSUN = 1.3271244e26           # IAU 2015 nominal G M_sun
AU = 1.495978707e13
DAY = 86400.0
BAR = 1.0e6
P0 = 250.0*BAR                 # the pgen's anchor pressure at r = ap (hard-coded there)


# ------------------------------------------------------------------------------ EOS
class EOSDump:
    """general-EOS dump in the dhjcs format: rho(p, T) and the adiabat dT/dp."""

    def __init__(self, fn):
        with open(fn) as fh:
            fh.readline()
            fh.readline()
            nx, ny, x0, dx, y0, dy = [float(v) for v in fh.readline()[1:].split()]
        self.nx, self.ny = int(nx), int(ny)
        self.x0, self.dx, self.y0, self.dy = x0, dx, y0, dy
        d = np.loadtxt(fn, comments='#', usecols=(0, 1))
        self.le = d[:, 0].reshape(self.ny, self.nx)      # log10 e/rho  [T row, rho col]
        self.lpr = d[:, 1].reshape(self.ny, self.nx)     # log10 p/rho
        self.ld = x0 + dx*np.arange(self.nx)

    def _bil(self, f, lr, lt):
        x = (lr - self.x0)/self.dx
        y = (lt - self.y0)/self.dy
        ix = np.clip(np.floor(x).astype(int), 0, self.nx - 2)
        iy = np.clip(np.floor(y).astype(int), 0, self.ny - 2)
        fx, fy = x - ix, y - iy
        return ((f[iy, ix]*(1 - fx) + f[iy, ix + 1]*fx)*(1 - fy)
                + (f[iy + 1, ix]*(1 - fx) + f[iy + 1, ix + 1]*fx)*fy)

    def rho(self, p, T):
        """rho(p, T) by bisection in log rho (log p = log rho + log p/rho is monotone)."""
        lp = np.log10(np.atleast_1d(np.asarray(p, float)))
        lt = np.log10(np.atleast_1d(np.asarray(T, float)))*np.ones_like(lp)
        lo = np.full(lp.shape, self.ld[0])
        hi = np.full(lp.shape, self.ld[-1])
        for _ in range(55):
            mid = 0.5*(lo + hi)
            up = mid + self._bil(self.lpr, mid, lt) > lp
            hi = np.where(up, mid, hi)
            lo = np.where(up, lo, mid)
        return 10**(0.5*(lo + hi))

    def nabla_ad(self, p, T):
        """dln T/dln p at constant entropy from de = (p/rho^2) d rho (finite differences
        of the dump in log rho, log T)."""
        rho = self.rho(p, T)[0]
        lr, lt = np.log10(rho), np.log10(T)
        h = 0.5*min(self.dx, self.dy)
        e = lambda a, b: 10**self._bil(self.le, a, b)
        pr = lambda a, b: 10**(a + self._bil(self.lpr, a, b))
        ln10 = np.log(10.0)
        e_r = (e(lr + h, lt) - e(lr - h, lt))/(2*h*ln10*rho)        # de/drho |_T
        e_t = (e(lr, lt + h) - e(lr, lt - h))/(2*h*ln10*T)          # de/dT |_rho
        p_r = (pr(lr + h, lt) - pr(lr - h, lt))/(2*h*ln10*rho)
        p_t = (pr(lr, lt + h) - pr(lr, lt - h))/(2*h*ln10*T)
        dT_drho = (p/rho**2 - e_r)/e_t                               # along the adiabat
        dp_drho = p_r + p_t*dT_drho
        return float(dT_drho/dp_drho*p/T)


def read_profile(fn):
    a = np.loadtxt(fn, comments='#')
    o = np.argsort(a[:, 0])
    return a[o, 0], a[o, 1]              # log10 p [barye] ascending, T [K]


class Column:
    """The pgen's column in potential form: Phi_col(p) = int_p^{250 bar} dp/rho (0 at
    250 bar, > 0 above, < 0 below), T held above the table top, the adiabat of the
    bottom row below the table bottom (problem/ic_adiabat_ext, dex decades)."""

    def __init__(self, profile, eos, ptop=1e-10, n=8001, dex=2.0, nad=801):
        lpt, Tt = read_profile(profile)
        if lpt[-1] < np.log10(P0):
            sys.exit('profile %s must reach 250 bar (as problem/ic_profile)' % profile)
        self.lp_bot, self.T_bot = lpt[-1], Tt[-1]
        # adiabat extension (RK4 in log10 p, as ic_adiabat_extension in the pgen)
        h = dex/(nad - 1)
        lpa = lpt[-1] + h*np.arange(nad)
        Ta = np.empty(nad)
        Ta[0] = Tt[-1]
        lt = np.log10(Tt[-1])
        sl = lambda lp, lt_: eos.nabla_ad(10**lp, 10**lt_)
        for k in range(nad - 1):
            k1 = sl(lpa[k], lt)
            k2 = sl(lpa[k] + 0.5*h, lt + 0.5*h*k1)
            k3 = sl(lpa[k] + 0.5*h, lt + 0.5*h*k2)
            k4 = sl(lpa[k] + h, lt + h*k3)
            lt += h*(k1 + 2*k2 + 2*k3 + k4)/6.0
            Ta[k + 1] = 10**lt
        self.nab_bot = sl(lpt[-1], np.log10(Tt[-1]))
        self.slope_bot = np.log10(Tt[-1]/Tt[-2])/(lpt[-1] - lpt[-2])
        lpx = np.concatenate([lpt, lpa[1:]])
        Tx = np.concatenate([Tt, Ta[1:]])
        self.lp = np.linspace(min(np.log10(ptop*BAR), lpt[0]), lpx[-1], n)
        self.T = np.interp(self.lp, lpx, Tx)          # held above the table top
        self.Tlin = np.interp(self.lp, lpt, Tt)       # (for reference: no extension)
        below = self.lp > lpt[-1]
        self.Tlin[below] = Tt[-1] + (Tt[-1] - Tt[-2])/(lpt[-1] - lpt[-2])*(self.lp[below]
                                                                           - lpt[-1])
        p = 10**self.lp
        self.rho = eos.rho(p, self.T)
        f = p/self.rho*np.log(10.0)                   # dPhi/dlog10 p
        cum = np.concatenate([[0], np.cumsum(0.5*(f[1:] + f[:-1])*np.diff(self.lp))])
        self.phi = np.interp(np.log10(P0), self.lp, cum) - cum

    def phi_at(self, pbar):
        return np.interp(np.log10(np.asarray(pbar)*BAR), self.lp, self.phi)

    def p_of_phi(self, ph):
        return 10**np.interp(-np.asarray(ph), -self.phi, self.lp)/BAR

    def T_of_p(self, pbar, lin=False):
        return np.interp(np.log10(np.asarray(pbar)*BAR), self.lp, self.Tlin if lin else self.T)


# ------------------------------------------------------------------------ potential
class Pot:
    def __init__(self, gm, om2, gms, a):
        self.gm, self.om2, self.gms, self.a = gm, om2, gms, a

    def tide(self, r, mu):
        if self.gms == 0.0:
            return 0.0*r
        d = np.sqrt(r*r + self.a**2 - 2*self.a*r*mu)
        return -self.gms/d + self.gms/self.a + self.gms*r*mu/self.a**2

    def phi(self, r, ap, s2, mu):
        """Phi_tot (the pgen's TotPotAt), s2 = sin^2 colatitude, mu = substellar cosine"""
        return self.gm*(1/ap - 1/r) - 0.5*self.om2*r*r*s2 + self.tide(r, mu)

    def surf(self, ap, ph, s2, mu, r0):
        """radius where Phi_tot = ph along (s2, mu), Newton"""
        r = np.full(np.broadcast(np.asarray(s2), np.asarray(mu)).shape, float(r0))
        for _ in range(60):
            hh = 1e-6*r
            F = self.phi(r, ap, s2, mu) - ph
            dF = (self.phi(r + hh, ap, s2, mu) - self.phi(r - hh, ap, s2, mu))/(2*hh)
            r = r - F/dF
        return r

    def l1(self):
        if self.gms == 0.0:
            return np.inf
        s = self.a*np.cbrt(self.gm/(3*(self.gm + self.gms)))
        for _ in range(100):
            F = self.gm/s**2 - self.gms/(self.a - s)**2 - self.om2*s + self.gms/self.a**2
            dF = -2*self.gm/s**3 - 2*self.gms/(self.a - s)**3 - self.om2
            s -= F/dF
        return s


PSI = np.linspace(-0.5*np.pi, 0.5*np.pi, 721)     # terminator ring angle (latitude on it)
S2_RING = np.cos(PSI)**2
MU_RING = np.zeros_like(PSI)


def r_eff(pot, col, ap, pref, rp):
    r = pot.surf(ap, col.phi_at(pref), S2_RING, MU_RING, rp)
    return np.sqrt(np.trapz(r**2, PSI)/np.pi)


def anchor(pot, col, rp, pref, ap0):
    ap = ap0
    for _ in range(60):
        R = r_eff(pot, col, ap, pref, rp)
        ap_new = 1.0/(1.0/ap + (1.0/rp - 1.0/R))   # d(1/R) ~ d(1/ap)
        if abs(ap_new/ap - 1) < 1e-14:
            ap = ap_new
            break
        ap = ap_new
    return ap


def wall_range(pot, col, ap, x1min, n=181):
    """Phi_tot on the inner-wall sphere r = x1min: min/max over the sphere -> p, T there"""
    th = np.linspace(0, np.pi, n)[:, None]
    ph = np.linspace(0, 2*np.pi, 2*n)[None, :]
    s2 = np.sin(th)**2*np.ones_like(ph)
    mu = np.sin(th)*np.cos(ph)               # substellar direction on the equator
    P = pot.phi(x1min, ap, s2, mu)
    return P.min(), P.max(), P


def main():
    A = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    A.add_argument('--mp', type=float, help='M_p [M_J]')
    A.add_argument('--gp', type=float, help='g at R_p [cm s^-2] (instead of --mp)')
    A.add_argument('--rp', type=float, required=True, help='R_p [R_J] at p_ref (terminator)')
    A.add_argument('--pref', type=float, required=True, help='transit pressure p_ref [bar]')
    A.add_argument('--a_au', type=float, required=True, help='semi-major axis a [AU]')
    A.add_argument('--porb', type=float, required=True, help='P_orb = P_rot [d]')
    A.add_argument('--teff', type=float, help='host T_eff [K] (for Teq)')
    A.add_argument('--rstar', type=float, help='R_* [R_sun] (for Teq)')
    A.add_argument('--mstar', type=float, help='literature M_* [M_sun] (Kepler check only)')
    A.add_argument('--profile', required=True, help='IC profile (log10 p [barye], T [K])')
    A.add_argument('--eos', required=True, help='general-EOS dump (dhjcs format)')
    A.add_argument('--x1max', type=float, help='deck x1max [cm] (Roche-lobe check)')
    A.add_argument('--adiabat_dex', type=float, default=2.0,
                   help='problem/ic_adiabat_ext_dex (default 2, as the pgen)')
    A.add_argument('--no_roche', action='store_true', help='rot_potential only')
    A.add_argument('--no_rot', action='store_true', help='no rot_potential (and no Roche)')
    a = A.parse_args()
    rp = a.rp*RJ
    if (a.mp is None) == (a.gp is None):
        sys.exit('give exactly one of --mp, --gp')
    gm = a.mp*GMJ if a.mp is not None else a.gp*rp**2
    om = 2*np.pi/(a.porb*DAY)
    om2 = 0.0 if a.no_rot else om**2
    aa = a.a_au*AU
    gms = 0.0 if (a.no_roche or a.no_rot) else om**2*aa**3 - gm
    eos = EOSDump(a.eos)
    col = Column(a.profile, eos, dex=a.adiabat_dex)
    pot = Pot(gm, om2, gms, aa)
    ap0 = 1.0/(1.0/rp + col.phi_at(a.pref)/gm)    # spherical anchor = start
    ap = anchor(pot, col, rp, a.pref, ap0)
    grav = gm/ap**2
    R = r_eff(pot, col, ap, a.pref, rp)
    ph = col.phi_at(a.pref)
    print('# G M_p = %.6e (M_p = %.4f M_J, g(R_p) = %.2f), R_p = %.4f R_J at p_ref = %.3e bar'
          % (gm, gm/GMJ, gm/rp**2, a.rp, a.pref))
    print('# omega = 2 pi/P = %.6e s^-1, a = %.6e cm' % (om, aa), end='')
    if gms > 0:
        print(', G M_* (Kepler, as the pgen) = %.6e -> M_* = %.4f M_sun%s'
              % (gms, gms/GMSUN, '' if a.mstar is None else
                 ' (literature %.4f: %+.2f %%)' % (a.mstar, 100*(gms/GMSUN/a.mstar - 1))))
    else:
        print(' (tide off)')
    if a.teff and a.rstar:
        print('# Teq = Teff sqrt(R*/2a) = %.1f K' % (a.teff*np.sqrt(a.rstar*RSUN/(2*aa))))
    print('# IC table bottom %.1f bar %.0f K: table slope dlnT/dlnp %.4f, nabla_ad %.4f'
          % (10**col.lp_bot/BAR, col.T_bot, col.slope_bot, col.nab_bot))
    mode = 'rot_potential only' if gms == 0 else 'Roche (rot_potential + tide)'
    if a.no_rot:
        mode = 'spherical (no rotation)'
    print('anchor [%s]: ap = x1min = %.6e cm, grav = %.4f cm/s^2, terminator R_eff(p_ref) '
          '= %.6f R_J' % (mode, ap, grav, R/RJ))
    dirs = dict(pole=(0.0, 0.0), terminator_eq=(1.0, 0.0), substellar=(1.0, 1.0),
                antistellar=(1.0, -1.0))
    print('   p_ref surface: ' + ', '.join('%s %.5f' % (k, pot.surf(ap, ph, v[0], v[1], rp)
                                                    / RJ) for k, v in dirs.items()) + ' R_J')
    pmin, pmax, _ = wall_range(pot, col, ap, ap)
    for tag, P in (('deepest (min Phi)', pmin), ('shallowest (max Phi)', pmax)):
        pb = col.p_of_phi(P)
        print('   inner wall r = x1min, %s: p = %.1f bar, T = %.0f K (adiabat ext.; the old '
              'linear extrapolation: %.0f K)%s' % (tag, pb, col.T_of_p(pb),
                                                   col.T_of_p(pb, lin=True),
                                                   '  [below the IC table]' if
                                                   np.log10(pb*BAR) > col.lp_bot else ''))
    for nm, (s2, mu) in dirs.items():
        pb = col.p_of_phi(pot.phi(ap, ap, s2, mu))
        print('      wall %-14s p = %8.1f bar  T = %6.0f K' % (nm, pb, col.T_of_p(pb)))
    if np.log10(col.p_of_phi(pmin)*BAR) > col.lp[-1] - 0.05:
        print('### WARNING: the wall reaches the end of the adiabat extension: raise '
              '--adiabat_dex and problem/ic_adiabat_ext_dex')
    l1 = pot.l1()
    if np.isfinite(l1):
        print('   L1 at %.5e cm = %.3f R_J; ap/r_L1 = %.3f%s' % (l1, l1/RJ, ap/l1, '' if a.x1max
              is None else ', x1max/r_L1 = %.3f' % (a.x1max/l1)))
        if a.x1max:
            ptop = max(pot.phi(a.x1max, ap, 0.0, 0.0), pot.phi(a.x1max, ap, 1.0, 1.0))
            print('   max Phi_tot on x1max %.4e vs Phi(L1) %.4e%s' % (
                ptop, pot.phi(l1, ap, 1.0, 1.0),
                '  (part of the top is outside the Roche lobe: the pgen WARNs, OK)'
                if ptop >= pot.phi(l1, ap, 1.0, 1.0) else ''))
    print('deck keys: <mesh> x1min = %.6e   <problem> ap = %.6e  grav = %.4f  omega = %.6e'
          '%s' % (ap, ap, grav, om, '' if gms == 0 else
                  '  roche = true  roche_a = %.6e' % aa))


if __name__ == '__main__':
    main()
