#!/usr/bin/env python3
"""T7 (radiative shocks): steady planar grey non-equilibrium radiative shock.

Design note section 8, (T7): Mach 2 (subcritical) and Mach 5 (supercritical)
Lowrie & Edwards (2008) shocks, coupled to hydro.  Pass: relative L1 below 2 %
in rho, T_gas and T_rad.

PROVENANCE.  Lowrie & Edwards, "Radiative shock solutions with grey
nonequilibrium diffusion", Shock Waves 18, 129 (2008), could not be read from
this machine.  The ODE system below is **derived here** from the steady grey
non-equilibrium-diffusion equations as they are written out in Skinner &
Ostriker 2013 (arXiv:1306.0010, their Eqs. 94-97) and Ferguson, Morel & Lowrie
2017 (arXiv:1612.06346, their Eqs. 9-11), which restate the Lowrie-Edwards
non-dimensionalisation.  The derivation is validated in `--selftest` and in the
README; no reference number is taken on faith.

NON-DIMENSIONALISATION (Lowrie-Edwards, as restated by Ferguson et al. 2017).
Lengths in L~, densities in rho0, temperatures in T0, velocities in the
upstream adiabatic sound speed a0, all evaluated in the far-upstream
equilibrium state, so that rho = T = 1 and v = M0 there.  Then p = rho T /
gamma, the sound speed is sqrt(T) and

    P0    = a_r T0^4 / (rho0 a0^2)          radiation/gas pressure ratio,
    C     = c / a0                          non-dimensional light speed,
    sigma = sigma~_t L~                     cross section per unit length,
    kappa = C / (3 sigma)                   radiative diffusivity,
    sigma_a                                 emission/absorption coupling, in
                                            the Lowrie-Edwards scaling
                                            sigma_a = sigma~_a L~ C.

With sigma_a = 1e6 and kappa = 1 (the standard case) one has C = sqrt(3
sigma_a) = 1732.05 and sigma~_t L~ = 577.35, which is the value Ferguson et al.
quote for the same problem.

STEADY EQUATIONS.  Writing Theta = T_rad^4, the three conservation integrals of
the steady planar system in the shock frame are

    m  = rho v                                                    (mass)
    Pm = rho v^2 + rho T / gamma + P0 Theta / 3                (momentum)
    Pe = m (v^2/2 + T/(gamma-1)) + (4/3) P0 v Theta + F_r        (energy)

with the diffusive radiative flux F_r = -P0 kappa dTheta/dx, and the gas
energy equation supplies the second ODE.  Eliminating rho and v gives, with
M = v/sqrt(T) the local Mach number and C_p = 1/(gamma-1),

    dTheta/dx = v [ 3 rho (v^2 - M0^2) + 6 C_p rho (T - 1)
                    + 8 P0 (Theta - rho) ] / (6 P0 kappa)
    dT/dx     = (gamma-1) P0 [ M0 dTheta/dx / 3
                    + rho sigma_a (Theta - T^4) (gamma M^2 - 1) ]
                / ( M0 rho (M^2 - 1) )

Both right-hand sides are singular at the adiabatic sonic point M = 1, so the
system is integrated in Mach space instead (as Lowrie & Edwards do): with
Theta(M, T) taken from the momentum integral, dM/dx follows by the chain rule
and dx/dM = (M^2-1)/W, dT/dM = U/W stay finite at M = 1, where U = (M^2-1)
dT/dx and W = (M^2-1) dM/dx.

ALGORITHM.  (1) The two far-field equilibrium states (T_rad = T_gas, F_r = 0)
follow from the three integrals by a 2-D root find.  (2) Both equilibria are
saddles; the Jacobian of (dT/dx, dTheta/dx) is formed numerically and the
precursor is launched from the upstream point along its unstable
eigendirection, the relaxation region from the downstream point along its
stable one.  (3) The embedded hydrodynamic shock conserves Theta and F_r (the
radiation field is continuous in the diffusion limit) as well as all three
integrals, so it maps the supersonic branch onto the subsonic one at equal
(Theta, F_r): the shock is located as the intersection of the two trajectories
in the (Theta, F_r) plane.  A subcritical shock (T just ahead of the shock
below the downstream T) and a supercritical one (equal to it) come out of the
same construction, as does the Zel'dovich spike behind the shock.

CAVEAT.  The Lowrie-Edwards solution uses Eddington closure, P_r = E_r/3.  An
M1 code in these optically thick shocks is in the diffusion regime and so
closes at 1/3 as well, which is why everyone gates on this solution, but the
comparison is a diffusion-limit comparison and not a test of the M1 closure
itself; Skinner & Ostriker attribute their own residual to exactly this.
"""

import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve, minimize_scalar

import common

KB_CGS = 1.380649e-16
MH_CGS = 1.67353e-24

# Skinner & Ostriker 2013 / Wibking & Krumholz 2022 cgs scaling of the same
# non-dimensional problem (QUOKKA RadhydroShockCGS): mu = m_p + m_e, a constant
# total cross section sigma~ = 577 cm^-1, rho0 = 5.69 g/cm^3, T0 = 2.18e6 K.
CGS_RHO0 = 5.69
CGS_T0 = 2.18e6
CGS_MU = MH_CGS
CGS_SIGMA = 577.0


class Params(object):
    """The non-dimensional problem."""

    def __init__(self, m0, p0=1.0e-4, gamma=5.0 / 3.0, sigma_a=1.0e6,
                 kappa=1.0):
        self.m0 = float(m0)
        self.p0 = float(p0)
        self.gamma = float(gamma)
        self.sigma_a = float(sigma_a)
        self.kappa = float(kappa)
        g, m0 = self.gamma, self.m0
        self.mass = m0
        self.mom = m0 ** 2 + 1.0 / g + self.p0 / 3.0
        self.ener = (m0 * (m0 ** 2 / 2.0 + 1.0 / (g - 1.0))
                     + (4.0 / 3.0) * self.p0 * m0)
        self.km = 3.0 * g * self.mom
        self.cspeed = np.sqrt(3.0 * self.sigma_a) * self.kappa

    def theta_of(self, mach, tgas):
        """Theta = T_rad^4 from the momentum integral at given (M, T)."""
        rho = self.m0 / (mach * np.sqrt(tgas))
        return ((self.km - 3.0 * rho * tgas - 3.0 * self.gamma * self.m0 ** 2
                 / rho) / (self.gamma * self.p0))

    def rho_branch(self, tgas, theta, sign):
        """rho from the momentum integral; sign=-1 supersonic, +1 subsonic."""
        b = self.km - self.gamma * self.p0 * theta
        disc = b * b - 36.0 * self.gamma * self.m0 ** 2 * tgas
        return (b + sign * np.sqrt(disc)) / (6.0 * tgas)

    def dtheta_dx(self, rho, vel, tgas, theta):
        return vel * (3.0 * rho * (vel ** 2 - self.m0 ** 2)
                      + 6.0 * rho * (tgas - 1.0) / (self.gamma - 1.0)
                      + 8.0 * self.p0 * (theta - rho)) \
            / (6.0 * self.p0 * self.kappa)

    def pieces(self, mach, tgas):
        """(v, rho, Theta, dTheta/dx, U, W) at a point of Mach space."""
        g = self.gamma
        vel = mach * np.sqrt(tgas)
        rho = self.m0 / vel
        theta = self.theta_of(mach, tgas)
        dth = self.dtheta_dx(rho, vel, tgas, theta)
        src = rho * self.sigma_a * (theta - tgas ** 4) * (g * mach ** 2 - 1.0)
        # U = (M^2-1) dT/dx, W = (M^2-1) dM/dx: both regular at M = 1.
        uu = (g - 1.0) * self.p0 * (self.m0 * dth / 3.0 + src) \
            / (self.m0 * rho)
        ww = mach ** 2 / (np.sqrt(tgas) * (g * mach ** 2 - 1.0)) \
            * (-g * self.p0 * dth * (mach ** 2 - 1.0) / (3.0 * self.m0)
               - (1.0 + g * mach ** 2) * uu / (2.0 * mach * np.sqrt(tgas)))
        return vel, rho, theta, dth, uu, ww

    def rhs_x(self, tgas, theta, sign):
        """(dT/dx, dTheta/dx) in x space, used only for the Jacobians."""
        rho = self.rho_branch(tgas, theta, sign)
        vel = self.m0 / rho
        mach = vel / np.sqrt(tgas)
        dth = self.dtheta_dx(rho, vel, tgas, theta)
        src = (rho * self.sigma_a * (theta - tgas ** 4)
               * (self.gamma * mach ** 2 - 1.0))
        dtg = (self.gamma - 1.0) * self.p0 * (self.m0 * dth / 3.0 + src) \
            / (self.m0 * rho * (mach ** 2 - 1.0))
        return np.array([dtg, dth])

    def equilibrium(self):
        """Downstream equilibrium (rho1, T1, v1) from the three integrals."""
        g, m0 = self.gamma, self.m0

        def residual(z):
            rho1, tg1 = z
            vel1 = m0 / rho1
            return [m0 * vel1 + rho1 * tg1 / g
                    + self.p0 * tg1 ** 4 / 3.0 - self.mom,
                    m0 * (vel1 ** 2 / 2.0 + tg1 / (g - 1.0))
                    + (4.0 / 3.0) * self.p0 * vel1 * tg1 ** 4 - self.ener]

        rg = (g + 1.0) * m0 ** 2 / ((g - 1.0) * m0 ** 2 + 2.0)
        tg = ((2.0 * g * m0 ** 2 - (g - 1.0))
              * ((g - 1.0) * m0 ** 2 + 2.0) / ((g + 1.0) ** 2 * m0 ** 2))
        sol = fsolve(residual, [rg, tg], full_output=True)
        z = sol[0]
        res = float(np.max(np.abs(residual(z))))
        return float(z[0]), float(z[1]), m0 / float(z[0]), res


def _eigendirection(par, tgas, theta, sign, growing):
    """dT/dTheta along the un/stable eigendirection at an equilibrium."""
    h = 1e-7
    jac = np.zeros((2, 2))
    for k, (dtg, dth) in enumerate([(h, 0.0), (0.0, h)]):
        jac[:, k] = (par.rhs_x(tgas + dtg, theta + dth, sign)
                     - par.rhs_x(tgas - dtg, theta - dth, sign)) / (2.0 * h)
    val, vec = np.linalg.eig(jac)
    want = [i for i in range(2) if (val[i].real > 0.0) == growing]
    if not want:
        raise RuntimeError("no %s eigendirection at the equilibrium "
                           "(eigenvalues %s)"
                           % ("growing" if growing else "decaying", val))
    i = want[0]
    if vec[1, i].real == 0.0:
        raise RuntimeError("degenerate eigendirection at the equilibrium")
    return vec[0, i].real / vec[1, i].real, float(val[i].real)


class Solution(object):
    """A solved steady radiative shock, with x = 0 at the embedded shock."""

    def __init__(self, par, eps=1e-6, rtol=1e-9, atol=1e-13, nsample=4000):
        self.par = par
        self.rho1, self.t1, self.v1, self.jump_res = par.equilibrium()
        theta1 = self.t1 ** 4

        sign1 = 1.0
        if abs(par.rho_branch(self.t1, theta1, 1.0) / self.rho1 - 1.0) > 1e-6:
            sign1 = -1.0
        self.sign_up, self.sign_dn = -1.0, sign1

        self.slope_up, self.rate_up = _eigendirection(par, 1.0, 1.0, -1.0,
                                                      True)
        self.slope_dn, self.rate_dn = _eigendirection(par, self.t1, theta1,
                                                      sign1, False)

        pre = self._integrate(1.0, 1.0, self.slope_up, eps, -1.0,
                              rtol, atol, nsample)
        rel = self._integrate(self.t1, theta1, self.slope_dn, -eps, sign1,
                              rtol, atol, nsample)
        self.pre, self.rel = pre, rel
        self._connect()

    def _integrate(self, tgas_eq, theta_eq, slope, eps, sign,
                   rtol, atol, nsample):
        par = self.par
        theta = theta_eq + eps
        tgas = tgas_eq + eps * slope
        rho = par.rho_branch(tgas, theta, sign)
        mach = par.m0 / rho / np.sqrt(tgas)

        def rhs(mm, yy):
            pcs = par.pieces(mm, yy[1])
            return [(mm * mm - 1.0) / pcs[5], pcs[4] / pcs[5]]

        sol = solve_ivp(rhs, (mach, 1.0), [0.0, tgas], method="LSODA",
                        rtol=rtol, atol=atol, dense_output=True)
        if sol.status != 0:
            raise RuntimeError("ODE integration failed: %s" % sol.message)
        # The approach to the equilibrium is exponential in x and therefore
        # linear in (M - M_eq): sample geometrically in that distance, or the
        # far tail (which sets the domain length) is not resolved at all.
        span = 1.0 - mach
        half = max(nsample // 2, 2)
        geo = np.geomspace(abs(span) * 1e-12, abs(span), half)
        mg = np.unique(np.concatenate([
            mach + np.sign(span) * geo,
            np.linspace(mach, 1.0, nsample - half)]))
        mg = np.clip(mg, min(mach, 1.0), max(mach, 1.0))
        if span < 0.0:
            mg = mg[::-1]
        nsample = mg.size
        yg = sol.sol(mg)
        out = {"mach": mg, "x": yg[0], "tgas": yg[1]}
        vel = np.empty(nsample)
        rho = np.empty(nsample)
        theta = np.empty(nsample)
        flux = np.empty(nsample)
        for i in range(nsample):
            pcs = par.pieces(mg[i], yg[1][i])
            vel[i], rho[i], theta[i] = pcs[0], pcs[1], pcs[2]
            flux[i] = -par.p0 * par.kappa * pcs[3]
        out.update(vel=vel, rho=rho, theta=theta, flux=flux)
        return out

    def _connect(self):
        """Embedded shock = equal (Theta, F_r) on the two branches."""
        pre, rel = self.pre, self.rel
        ia = np.argsort(pre["theta"])
        ib = np.argsort(rel["theta"])
        tha, fa = pre["theta"][ia], pre["flux"][ia]
        thb, fb = rel["theta"][ib], rel["flux"][ib]
        lo = max(tha[0], thb[0])
        hi = min(tha[-1], thb[-1])
        if not hi > lo:
            raise RuntimeError("precursor and relaxation branches do not "
                               "overlap in Theta; no shock can be located")

        def gap(th):
            return np.interp(th, tha, fa) - np.interp(th, thb, fb)

        grid = np.linspace(lo, hi, 2001)
        vals = gap(grid)
        cross = np.nonzero(np.diff(np.sign(vals)))[0]
        if cross.size == 0:
            raise RuntimeError("no (Theta, F_r) crossing: the shock is "
                               "continuous or the branches are wrong")
        k = cross[-1]
        self.theta_s = brentq(gap, grid[k], grid[k + 1], xtol=1e-14,
                              rtol=1e-14)
        self.x_pre_s = float(np.interp(self.theta_s, tha, pre["x"][ia]))
        self.x_rel_s = float(np.interp(self.theta_s, thb, rel["x"][ib]))
        self.tgas_pre_s = float(np.interp(self.theta_s, tha,
                                          pre["tgas"][ia]))
        self.tgas_post_s = float(np.interp(self.theta_s, thb,
                                           rel["tgas"][ib]))
        self.mach_pre_s = float(np.interp(self.theta_s, tha, pre["mach"][ia]))
        self.mach_post_s = float(np.interp(self.theta_s, thb, rel["mach"][ib]))
        self.trad_s = self.theta_s ** 0.25
        self.subcritical = self.tgas_pre_s < 0.999 * self.t1

        keep = pre["theta"] <= self.theta_s
        self.x_a = pre["x"][keep] - self.x_pre_s
        self.rho_a = pre["rho"][keep]
        self.vel_a = pre["vel"][keep]
        self.tgas_a = pre["tgas"][keep]
        self.theta_a = pre["theta"][keep]
        self.flux_a = pre["flux"][keep]

        keep = rel["theta"] >= self.theta_s
        self.x_b = rel["x"][keep] - self.x_rel_s
        self.rho_b = rel["rho"][keep]
        self.vel_b = rel["vel"][keep]
        self.tgas_b = rel["tgas"][keep]
        self.theta_b = rel["theta"][keep]
        self.flux_b = rel["flux"][keep]

        for nm in ("a", "b"):
            order = np.argsort(getattr(self, "x_" + nm))
            for q in ("x", "rho", "vel", "tgas", "theta", "flux"):
                setattr(self, "%s_%s" % (q, nm),
                        getattr(self, "%s_%s" % (q, nm))[order])
        self.tgas_spike = float(self.tgas_b.max())

    # ------------------------------------------------------------------
    def domain(self, dtfrac=1.0e-4):
        """Recommended (xmin, xmax): |T - T_eq| = dtfrac at both ends."""
        xmin = self.x_a[0]
        hit = np.nonzero(self.tgas_a - 1.0 >= dtfrac)[0]
        if hit.size:
            xmin = float(np.interp(1.0 + dtfrac, self.tgas_a, self.x_a))
        xmax = self.x_b[-1]
        dev = np.abs(self.tgas_b - self.t1)
        hit = np.nonzero(dev >= dtfrac)[0]
        if hit.size:
            xmax = float(self.x_b[hit[-1]])
        return xmin, xmax

    def sample(self, xgrid):
        """rho, v, T_gas, T_rad on an arbitrary x grid (shock at x = 0)."""
        xgrid = np.asarray(xgrid, dtype=float)
        out = {}
        left = xgrid < 0.0
        for name, qa, qb in (("rho", self.rho_a, self.rho_b),
                             ("vel", self.vel_a, self.vel_b),
                             ("tgas", self.tgas_a, self.tgas_b),
                             ("theta", self.theta_a, self.theta_b)):
            val = np.empty_like(xgrid)
            val[left] = np.interp(xgrid[left], self.x_a, qa)
            val[~left] = np.interp(xgrid[~left], self.x_b, qb)
            out[name] = val
        out["trad"] = out["theta"] ** 0.25
        return out

    def conservation(self):
        """max relative residual of the three integrals over the profile."""
        par = self.par
        rho = np.concatenate([self.rho_a, self.rho_b])
        vel = np.concatenate([self.vel_a, self.vel_b])
        tgas = np.concatenate([self.tgas_a, self.tgas_b])
        theta = np.concatenate([self.theta_a, self.theta_b])
        flux = np.concatenate([self.flux_a, self.flux_b])
        c1 = rho * vel / par.mass - 1.0
        c2 = (rho * vel ** 2 + rho * tgas / par.gamma
              + par.p0 * theta / 3.0) / par.mom - 1.0
        c3 = (par.mass * (vel ** 2 / 2.0 + tgas / (par.gamma - 1.0))
              + (4.0 / 3.0) * par.p0 * vel * theta + flux) / par.ener - 1.0
        return (float(np.abs(c1).max()), float(np.abs(c2).max()),
                float(np.abs(c3).max()))

    def endpoints(self):
        """Relative distance of the two profile ends from the equilibria."""
        up = max(abs(self.rho_a[0] - 1.0), abs(self.tgas_a[0] - 1.0),
                 abs(self.theta_a[0] - 1.0))
        dn = max(abs(self.rho_b[-1] / self.rho1 - 1.0),
                 abs(self.tgas_b[-1] / self.t1 - 1.0),
                 abs(self.theta_b[-1] / self.t1 ** 4 - 1.0))
        return float(up), float(dn)


def lowrie_edwards_rhs(par, mach, tgas):
    """(dx/dM, dT/dM) in the published Lowrie-Edwards closed form.

    Transcribed from the reference implementation of Lowrie & Edwards (2008)
    shipped with QUOKKA (``extern/LowrieEdwards/radshock.py``, function
    ``shock``), and used only to cross-check the derivation in this file's
    docstring, which was obtained independently.
    """
    g = par.gamma
    vel = mach * np.sqrt(tgas)
    rho = par.m0 / vel
    theta = par.theta_of(mach, tgas)
    dth = par.dtheta_dx(rho, vel, tgas, theta)
    rsrc = 3.0 * rho * par.sigma_a * (theta - tgas ** 4)
    zd = (par.m0 * dth
          + ((g - 1.0) / (g + 1.0)) * (g * mach ** 2 + 1.0) * rsrc)
    zn = par.m0 * dth + (g * mach ** 2 - 1.0) * rsrc
    dx_dm = (-6.0 * (par.m0 * rho * tgas) / ((g + 1.0) * par.p0 * mach)
             * (mach ** 2 - 1.0) / zd)
    dt_dm = -2.0 * ((g - 1.0) / (g + 1.0)) * (tgas / mach) * (zn / zd)
    return dx_dm, dt_dm


def cross_check(par, states=((2.5, 1.3), (1.4, 2.0), (0.6, 3.5))):
    """Max relative difference between our RHS and the published one."""
    worst = 0.0
    for mach, tgas in states:
        pcs = par.pieces(mach, tgas)
        mine = np.array([(mach ** 2 - 1.0) / pcs[5], pcs[4] / pcs[5]])
        theirs = np.array(lowrie_edwards_rhs(par, mach, tgas))
        worst = max(worst, float(np.max(np.abs(mine / theirs - 1.0))))
    return worst


class Units(object):
    """Scale factors from the non-dimensional solution to physical units."""

    def __init__(self, par, rho0=CGS_RHO0, t0=CGS_T0, mu=CGS_MU,
                 sigma_cgs=CGS_SIGMA, a_rad=common.AR_CGS, c_light=common.C_CGS,
                 nondim=False):
        self.par = par
        if nondim:
            self.rho0 = self.t0 = self.length = self.vel = 1.0
            # E_rad = P0 Theta and e_gas = rho T / (gamma (gamma-1)) in the
            # non-dimensional variables, so a_r -> P0 and c_v -> that.
            self.a_rad = par.p0
            self.cv = 1.0 / (par.gamma * (par.gamma - 1.0))
            self.p0_implied = par.p0
            self.c_light = par.cspeed
            self.sigma_cgs = par.cspeed / (3.0 * par.kappa)
            self.nondim = True
            return
        self.nondim = False
        self.rho0, self.t0, self.mu = rho0, t0, mu
        self.a_rad, self.c_light = a_rad, c_light
        self.vel = np.sqrt(par.gamma * KB_CGS * t0 / mu)
        self.cv = KB_CGS / (mu * (par.gamma - 1.0))
        self.p0_implied = a_rad * t0 ** 4 / (rho0 * self.vel ** 2)
        cnum = c_light / self.vel
        self.length = cnum / (3.0 * par.kappa * sigma_cgs)
        self.sigma_cgs = sigma_cgs


def solve(args):
    par = Params(args.m0, args.p0, args.gamma, args.sigma_a, args.kappa)
    return Solution(par, eps=args.eps, rtol=args.rtol, atol=args.atol,
                    nsample=args.nsample)


def describe(args, sol, units):
    par = sol.par
    common.report(args, "parameters: M0=%.6g P0=%.6g gamma=%.6g "
                        "sigma_a=%.6g kappa=%.6g  (C = %.6g)"
                  % (par.m0, par.p0, par.gamma, par.sigma_a, par.kappa,
                     par.cspeed))
    common.report(args, "upstream   : rho=1 v=%.8f T=%.8f" % (par.m0, 1.0))
    common.report(args, "downstream : rho=%.8f v=%.8f T=%.8f "
                        "(jump residual %.2e)"
                  % (sol.rho1, sol.v1, sol.t1, sol.jump_res))
    common.report(args, "embedded shock at x=0: Theta_s=%.8g T_rad,s=%.8f  "
                        "M %.5f -> %.5f" % (sol.theta_s, sol.trad_s,
                                            sol.mach_pre_s, sol.mach_post_s))
    common.report(args, "T_gas ahead=%.6f  behind=%.6f  spike max=%.6f  "
                        "(T1=%.6f) -> %s"
                  % (sol.tgas_pre_s, sol.tgas_post_s, sol.tgas_spike, sol.t1,
                     "SUBcritical" if sol.subcritical else "SUPERcritical"))
    c1, c2, c3 = sol.conservation()
    up, dn = sol.endpoints()
    common.report(args, "conservation residuals: mass %.2e momentum %.2e "
                        "energy %.2e" % (c1, c2, c3))
    common.report(args, "endpoint offsets from equilibrium: up %.2e dn %.2e"
                  % (up, dn))
    if not units.nondim:
        common.report(args, "cgs scaling: a0=%.6e cm/s  L=%.6e cm  "
                            "c_v=%.6e  P0(implied)=%.6e"
                      % (units.vel, units.length, units.cv, units.p0_implied))


def athinput_block(sol, units, dtfrac, ncells):
    par = sol.par
    xmin, xmax = sol.domain(dtfrac)
    lx = (xmax - xmin) * units.length
    lines = ["<problem>",
             "# Lowrie & Edwards steady radiative shock, M0 = %g, P0 = %g,"
             % (par.m0, par.p0),
             "# gamma = %g, sigma_a = %g, kappa = %g  (%s)"
             % (par.gamma, par.sigma_a, par.kappa,
                "subcritical" if sol.subcritical else "supercritical"),
             "m1_shock_rho_l = %.8e" % (1.0 * units.rho0),
             "m1_shock_v_l   = %.8e" % (par.m0 * units.vel),
             "m1_shock_t_l   = %.8e" % (1.0 * units.t0),
             "m1_shock_rho_r = %.8e" % (sol.rho1 * units.rho0),
             "m1_shock_v_r   = %.8e" % (sol.v1 * units.vel),
             "m1_shock_t_r   = %.8e" % (sol.t1 * units.t0),
             "m1_shock_xs    = %.8e" % (-xmin * units.length),
             "# <mesh>: x1min = 0  x1max = %.8e  nx1 = %d" % (lx, ncells),
             "#         shock sits at x1 = %.8e" % (-xmin * units.length),
             "#         inflow (Dirichlet) on both x1 faces",
             "# opacity: sigma = %.8e (per unit length, constant)"
             % (getattr(units, "sigma_cgs", 1.0 / units.length)),
             "</problem>"]
    return "\n".join(lines)


def write_reference(path, sol, units, xgrid):
    prof = sol.sample(xgrid)
    cols = np.column_stack([xgrid * units.length,
                            prof["rho"] * units.rho0,
                            prof["vel"] * units.vel,
                            prof["tgas"] * units.t0,
                            prof["trad"] * units.t0])
    par = sol.par
    head = ("Lowrie-Edwards steady radiative shock (t7_radshock.py)\n"
            "M0=%.10g P0=%.10g gamma=%.10g sigma_a=%.10g kappa=%.10g\n"
            "scales: rho0=%.10e T0=%.10e v0=%.10e L=%.10e\n"
            "shock at x=%.10e ; columns: x rho v T_gas T_rad"
            % (par.m0, par.p0, par.gamma, par.sigma_a, par.kappa,
               units.rho0, units.t0, units.vel, units.length,
               -sol.domain(1e-4)[0] * units.length))
    np.savetxt(path, cols, header=head, fmt="%.10e")


# ----------------------------------------------------------------------
def dump_profiles(dump, args, units):
    """(x, rho, T_gas, T_rad) from a dump, in the dump's own units."""
    x, rho = common.extract_1d(dump, "dens", axis=args.axis,
                               reduce=args.reduce)
    if dump.has("eint"):
        _, eint = common.extract_1d(dump, "eint", axis=args.axis,
                                    reduce=args.reduce)
    elif dump.has("press"):
        _, prs = common.extract_1d(dump, "press", axis=args.axis,
                                   reduce=args.reduce)
        eint = prs / (args.gamma - 1.0)
    else:
        raise KeyError("dump carries neither 'eint' nor 'press'")
    _, erad = common.extract_1d(dump, "m1_e", axis=args.axis,
                                reduce=args.reduce)
    cv = args.cv if args.cv is not None else units.cv
    arad = args.a_rad if args.a_rad is not None else units.a_rad
    tgas = eint / (rho * cv)
    trad = np.maximum(erad / arad, 0.0) ** 0.25
    return x, rho, tgas, trad


def compare(sol, units, x, rho, tgas, trad, args):
    """Align on x and return the three relative L1 errors and the shift."""
    scale = np.array([units.rho0, units.t0, units.t0])
    dat = np.vstack([rho, tgas, trad])

    def errors(shift):
        prof = sol.sample((x - shift) / units.length)
        ref = np.vstack([prof["rho"], prof["tgas"], prof["trad"]]) \
            * scale[:, None]
        return np.array([common.l1_rel(dat[i], ref[i]) for i in range(3)])

    dx = float(x[1] - x[0])
    guess = args.shock_x
    if guess is None:
        # steepest density gradient = the embedded hydrodynamic shock
        guess = float(0.5 * (x[1:] + x[:-1])[np.argmax(np.abs(np.diff(rho)))])
    span = args.max_shift * dx
    res = minimize_scalar(lambda s: float(errors(s).sum()),
                          bounds=(guess - span, guess + span),
                          method="bounded",
                          options={"xatol": 1e-4 * dx})
    shift = float(res.x)
    return errors(shift), shift, (shift - guess) / dx


# ----------------------------------------------------------------------
def selftest(sol, units, args):
    """Reference resampled onto a coarse grid, optionally corrupted."""
    xmin, xmax = sol.domain(args.dt_domain)
    n = args.selftest_nx
    xs = common.centres(0.0, (xmax - xmin) * units.length, n)
    xshift = -xmin * units.length
    prof = sol.sample((xs - xshift) / units.length)
    rho = prof["rho"] * units.rho0
    tgas = prof["tgas"] * units.t0
    trad = prof["trad"] * units.t0
    if args.selftest_fail:
        # 5 % error in the Zel'dovich spike / relaxation region and a
        # precursor stretched by 20 %.
        tgas = np.where(xs >= xshift, tgas * 1.05, tgas)
        pre = xs < xshift
        xstretch = xshift + (xs - xshift) * 1.2
        tg2 = np.interp(xstretch[pre], xs, tgas)
        tgas = tgas.copy()
        tgas[pre] = tg2
    cv = args.cv if args.cv is not None else units.cv
    arad = args.a_rad if args.a_rad is not None else units.a_rad
    data = {"dens": rho[None, None, :],
            "velx": (prof["vel"] * units.vel)[None, None, :],
            "eint": (rho * cv * tgas)[None, None, :],
            "m1_e": (arad * trad ** 4)[None, None, :]}
    return common.Dump(0.0, xs, np.array([0.5]), np.array([0.5]), data,
                       "<selftest>")


def main():
    p = common.base_parser(__doc__.splitlines()[0])
    p.add_argument("dump", nargs="?", default=None,
                   help="AthenaK bin dump of the steady shock")
    p.add_argument("--m0", type=float, default=2.0,
                   help="upstream Mach number (2 subcritical, 5 supercritical)")
    p.add_argument("--p0", type=float, default=1.0e-4,
                   help="P0 = a_r T0^4 / (rho0 a0^2)")
    p.add_argument("--gamma", type=float, default=5.0 / 3.0)
    p.add_argument("--sigma-a", type=float, default=1.0e6,
                   help="Lowrie-Edwards absorption coupling")
    p.add_argument("--kappa", type=float, default=1.0,
                   help="Lowrie-Edwards radiative diffusivity")
    p.add_argument("--nondim", action="store_true",
                   help="dump/reference are in the non-dimensional units")
    p.add_argument("--rho0", type=float, default=CGS_RHO0,
                   help="upstream density in the dump's units (cgs)")
    p.add_argument("--t0", type=float, default=CGS_T0,
                   help="upstream temperature in the dump's units (K)")
    p.add_argument("--mu", type=float, default=CGS_MU,
                   help="mean molecular mass in g")
    p.add_argument("--sigma-cgs", type=float, default=CGS_SIGMA,
                   help="total cross section per unit length (1/cm)")
    p.add_argument("--a-rad", type=float, default=None,
                   help="radiation constant in the dump's units")
    p.add_argument("--cv", type=float, default=None,
                   help="specific heat in the dump's units "
                        "(default: k_B / (mu (gamma-1)))")
    p.add_argument("--tol", type=float, default=0.02,
                   help="allowed relative L1 in rho, T_gas and T_rad")
    p.add_argument("--max-shift", type=float, default=20.0,
                   help="alignment search half-width, in cells")
    p.add_argument("--shock-x", type=float, default=None,
                   help="initial guess for the shock position in the dump")
    p.add_argument("--dt-domain", type=float, default=1.0e-4,
                   help="|T - T_eq| that defines the recommended domain")
    p.add_argument("--ncells", type=int, default=512,
                   help="recommended resolution printed in the IC block")
    p.add_argument("--write-ref", default=None,
                   help="write the reference profile to this text file")
    p.add_argument("--ref-nx", type=int, default=4096,
                   help="number of points in --write-ref")
    p.add_argument("--athinput", action="store_true",
                   help="print the <problem> block for the initial condition")
    p.add_argument("--eps", type=float, default=1.0e-6,
                   help="offset from the equilibria along the eigendirection")
    p.add_argument("--rtol", type=float, default=1.0e-9)
    p.add_argument("--atol", type=float, default=1.0e-13)
    p.add_argument("--nsample", type=int, default=4000,
                   help="points kept per branch of the ODE solution")
    p.add_argument("--cons-tol", type=float, default=1.0e-8,
                   help="allowed residual of the three conservation integrals")
    p.add_argument("--cross-tol", type=float, default=1.0e-10,
                   help="allowed difference from the published RHS")
    p.add_argument("--end-tol", type=float, default=1.0e-3,
                   help="allowed offset of the profile ends from equilibrium")
    p.add_argument("--selftest-nx", type=int, default=512)
    p.add_argument("--axis", type=int, default=1, choices=(1, 2, 3))
    p.add_argument("--reduce", default="mid", choices=("mid", "mean"))
    p.add_argument("--selftest-fail", action="store_true",
                   help="selftest variant that must FAIL")
    args = p.parse_args()

    sol = solve(args)
    units = Units(sol.par, args.rho0, args.t0, args.mu, args.sigma_cgs,
                  args.a_rad if args.a_rad is not None else common.AR_CGS,
                  common.C_CGS, nondim=args.nondim)
    describe(args, sol, units)

    c1, c2, c3 = sol.conservation()
    up, dn = sol.endpoints()
    xchk = cross_check(sol.par)
    common.report(args, "RHS vs the published Lowrie-Edwards closed form: "
                        "max relative difference %.2e" % xchk)
    internal = (max(c1, c2, c3) <= args.cons_tol
                and max(up, dn) <= args.end_tol
                and xchk <= args.cross_tol)

    xmin, xmax = sol.domain(args.dt_domain)
    lx = (xmax - xmin) * units.length
    dxc = lx / args.ncells
    common.report(args, "recommended domain: length %.8e (shock at %.8e), "
                        "%d cells" % (lx, -xmin * units.length, args.ncells))
    common.report(args, "  precursor %.4e (%.1f cells), relaxation %.4e "
                        "(%.1f cells) at that resolution"
                  % (-xmin * units.length, -xmin * units.length / dxc,
                     xmax * units.length, xmax * units.length / dxc))
    if args.write_ref:
        write_reference(args.write_ref, sol, units,
                        np.linspace(xmin, xmax, args.ref_nx))
        common.report(args, "wrote %s" % args.write_ref)
    if args.athinput:
        print(athinput_block(sol, units, args.dt_domain, args.ncells))

    if args.selftest:
        dump = selftest(sol, units, args)
    elif args.dump is None:
        common.verdict(internal,
                       "T7 reference M0=%g: conservation %.2e (<= %.3g)  "
                       "endpoints %.2e (<= %.3g)  [no dump given]"
                       % (args.m0, max(c1, c2, c3), args.cons_tol,
                          max(up, dn), args.end_tol))
        return
    else:
        dump = common.load_dump(args.dump, args.bin_convert_dir,
                                args.all_ranks)

    x, rho, tgas, trad = dump_profiles(dump, args, units)
    errs, shift, cells = compare(sol, units, x, rho, tgas, trad, args)
    ok = internal and bool((errs <= args.tol).all())
    common.verdict(ok, "T7 radshock M0=%g (%s): L1 rho=%.3e T_gas=%.3e "
                       "T_rad=%.3e (<= %.3g)  shift=%.2f cells  "
                       "conservation=%.2e"
                   % (args.m0, "sub" if sol.subcritical else "super",
                      errs[0], errs[1], errs[2], args.tol, cells,
                      max(c1, c2, c3)))


if __name__ == "__main__":
    sys.exit(main())
