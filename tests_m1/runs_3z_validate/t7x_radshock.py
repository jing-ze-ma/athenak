#!/usr/bin/env python3
"""t7_radshock.py with the two branches integrated in x instead of Mach space.

For strong (M0 = 27, 50) Lowrie-Edwards shocks the Mach-space integration of
t7_radshock.py stalls (LSODA step underflow on the precursor branch).  Each
branch is integrated here directly in x from its saddle (the same eigen-
directions, the same (Theta, F_r) connection and everything downstream of the
Solution object), and stopped just before M = 1 or where the momentum
quadratic loses its root.  Usage is that of t7_radshock.py (same options).
"""
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/tests_m1')
import t7_radshock as t7  # noqa: E402

XSPAN = 50.0   # max |x| of a branch in L~ units
DM = 1.0e-4    # stop at |M - 1| = DM


def _integrate_x(self, tgas_eq, theta_eq, slope, eps, sign, rtol, atol,
                 nsample):
    par = self.par
    # state (M, T_gas): regular everywhere except M = 1 (the (T, Theta) pair
    # folds at the isothermal sonic point M^2 = 1/gamma)
    # eps RELATIVE to Theta_eq (Theta1 = T1^4 is 1e7 at M0 = 50)
    # the eigenvector (dT, dTheta) = (slope, 1) scaled to a RELATIVE size eps
    scl = eps / max(abs(slope) / tgas_eq, 1.0 / theta_eq)
    theta = theta_eq + scl
    tgas = tgas_eq + scl * slope
    rho = par.rho_branch(tgas, theta, sign)
    y0 = [par.m0 / rho / np.sqrt(tgas), tgas]
    # precursor (sign -1, supersonic): the growing mode runs forward in x;
    # relaxation (subsonic): the decaying mode, integrated backward in x
    xend = XSPAN if sign < 0 else -XSPAN

    def rhs(x, y):
        pcs = par.pieces(y[0], y[1])
        d = y[0] * y[0] - 1.0
        return [pcs[5] / d, pcs[4] / d]

    def ev_sonic(x, y):
        return abs(y[0] - 1.0) - DM
    ev_sonic.terminal = True

    sol = solve_ivp(rhs, (0.0, xend), y0, method="Radau", rtol=rtol,
                    atol=atol, events=(ev_sonic,), dense_output=True)
    if sol.status < 0:
        raise RuntimeError("x-space ODE integration failed: %s" % sol.message)
    xe = sol.t[-1]
    # geometric sampling from the saddle (exponential approach) + linear
    half = max(nsample // 2, 2)
    geo = np.geomspace(abs(xe) * 1e-12, abs(xe), half)
    xg = np.unique(np.concatenate([np.sign(xe) * geo,
                                   np.linspace(0.0, xe, nsample - half)]))
    if xe < 0.0:
        xg = xg[::-1]
    yg = sol.sol(xg)
    n = xg.size
    out = {"x": xg, "tgas": yg[1]}
    for q in ("mach", "vel", "rho", "theta", "flux"):
        out[q] = np.empty(n)
    for i in range(n):
        pcs = par.pieces(yg[0][i], yg[1][i])
        out["mach"][i] = yg[0][i]
        out["vel"][i], out["rho"][i], out["theta"][i] = pcs[0], pcs[1], pcs[2]
        out["flux"][i] = -par.p0 * par.kappa * pcs[3]
    print("branch sign %+d: x_end = %.6g, M_end = %.6f, status %d (%s)"
          % (sign, xe, out["mach"][-1], sol.status, sol.message),
          file=sys.stderr)
    return out


def _rho_branch_clamped(self, tgas, theta, sign):
    """t7's rho_branch with the discriminant clamped at 0: identical where the
    branch exists; keeps the implicit solver's trial points finite past the
    ev_disc event."""
    b = self.km - self.gamma * self.p0 * theta
    disc = np.maximum(b * b - 36.0 * self.gamma * self.m0 ** 2 * tgas, 0.0)
    return (b + sign * np.sqrt(disc)) / (6.0 * tgas)


def _eigendirection_rel(par, tgas, theta, sign, growing):
    """t7._eigendirection with RELATIVE difference steps (h = 1e-7 absolute is
    round-off at Theta = 1e7)."""
    ht, hh = 1e-7 * tgas, 1e-7 * theta
    jac = np.zeros((2, 2))
    for k, (dtg, dth) in enumerate([(ht, 0.0), (0.0, hh)]):
        jac[:, k] = (par.rhs_x(tgas + dtg, theta + dth, sign)
                     - par.rhs_x(tgas - dtg, theta - dth, sign)) \
            / (2.0 * (dtg + dth))
    val, vec = np.linalg.eig(jac)
    want = [i for i in range(2) if (val[i].real > 0.0) == growing]
    if not want:
        raise RuntimeError("no eigendirection (eigenvalues %s)" % val)
    i = want[0]
    return vec[0, i].real / vec[1, i].real, float(val[i].real)


_connect_jump = t7.Solution._connect


def _connect_any(self):
    """The embedded-shock connection of t7_radshock.py, or, when the branches
    have no (Theta, F_r) crossing but meet at the sonic point with equal
    (T, Theta, F_r), the CONTINUOUS (radiation-dominated) shock: the two
    branches joined at M = 1, which is put at x = 0."""
    try:
        _connect_jump(self)
        self.continuous = False
        return
    except RuntimeError:
        pre, rel = self.pre, self.rel
        mis = max(abs(pre["tgas"][-1] / rel["tgas"][-1] - 1.0),
                  abs(pre["theta"][-1] / rel["theta"][-1] - 1.0),
                  abs(pre["flux"][-1] / rel["flux"][-1] - 1.0))
        print("no embedded shock; sonic-point mismatch of the branches %.2e"
              % mis, file=sys.stderr)
        if mis > 2.0e-3:
            raise
    self.continuous = True
    self.theta_s = 0.5 * (pre["theta"][-1] + rel["theta"][-1])
    self.x_pre_s, self.x_rel_s = pre["x"][-1], rel["x"][-1]
    self.tgas_pre_s, self.tgas_post_s = pre["tgas"][-1], rel["tgas"][-1]
    self.mach_pre_s, self.mach_post_s = pre["mach"][-1], rel["mach"][-1]
    self.trad_s = self.theta_s ** 0.25
    self.subcritical = False
    for nm, br, xs in (("a", pre, self.x_pre_s), ("b", rel, self.x_rel_s)):
        order = np.argsort(br["x"] - xs)
        setattr(self, "x_" + nm, (br["x"] - xs)[order])
        for q in ("rho", "vel", "tgas", "theta", "flux"):
            setattr(self, "%s_%s" % (q, nm), br[q][order])
    self.tgas_spike = float(self.tgas_b.max())


t7.Solution._connect = _connect_any
t7._eigendirection = _eigendirection_rel
t7.Solution._integrate = _integrate_x
t7.Params.rho_branch = _rho_branch_clamped

if __name__ == "__main__":
    sys.exit(t7.main())
