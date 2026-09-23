#!/usr/bin/env python3
"""DISCRETE 1-D model of the implicit rad_m1 step with the gas velocity in the enthalpy
flux (4/3) v E made implicit (option 1 of tests_m1/runs_3u_time2impl/README.md).

Space = the code's staggered layout (rad_m1_implicit.cpp): E, T, rho, v at cells, F0 on
x1 faces.  Fourier symbol per mode theta = k dx on a periodic grid:
  face flux   F0'_f  from (E_{i+1} - E_i)/dx                 (Eddington, w = 1/3)
  E row       -(F0_{i+1/2} - F0_{i-1/2})/dx - (A_{i+1/2} - A_{i-1/2})/dx + coupling
  enthalpy    A_f = (4/3) E0 v_up : 'cen' = face mean of the two cells,
              'upw' = left cell (the upwind choice frozen for v_f > 0)
  force       dv/dt += (kappa/c) * mean of the two face F0 (the write-back's dm1)
  hydro rows  centred i sin(theta)/dx + Rusanov dissipation at the adiabatic gas speed
The linearisation is about v0 = 0, where a(v')E' = (4/3) E0 v' exactly, so the Newton
form a(v^k)E' + (4/3)E^k (v' - v^k) of the plan IS the fully implicit 'imp' here.

Time schemes (time2_model conventions, FSAL carried exactly):
  imp     H-ESDIRK2, stage solve with v' in the linear system
  rhs     H-ESDIRK2, v lagged at the stage's old gas state (the code's M1_IW_V1)
  code    present scheme: Heun gas-only, then backward Euler with v lagged
  codeimp present scheme but v' implicit in the BE solve (the phase-B first gate)
Order: omega_ref = the SEMI-DISCRETE acoustic root of the same spatial operator (so p is
the temporal order), N = 64 cells per wavelength, 12 radwave cases.
Stability: max over theta = 2 pi m/64, m = 1..32, of rho(step map) - 1; dt = cfl dx/c_s.
Part 3 (python3 vimp_model.py iters): right-preconditioned BiCGStab iterations of the
face-eliminated E system on a 2-D grid, x1 line preconditioner tri vs penta."""
import math
import sys
import numpy as np
sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_3t_time2")
sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_3r_radwave")
import time2_model as tm  # noqa: E402
import radwave_disp as rd  # noqa: E402

G = tm.GAM
B2 = G / (2 * (1 - G))
I5 = np.eye(5)


class DOps(object):
    """discrete symbol for q = (rho, v, T, E, F0_face) at theta = k dx"""

    def __init__(self, bg, th, dx, enth="cen", diss=True):
        r0, t0, g, c, kap = bg.rho0, bg.t0, bg.gamma, bg.c, bg.kappa
        e0, a = bg.e0, bg.arad
        b = 4.0 * a * t0 ** 3
        s = 2j * math.sin(th / 2) / dx          # staggered difference
        av = math.cos(th / 2)                   # face <-> cell mean
        dc = 1j * math.sin(th) / dx             # centred cell difference
        dis = 0.0
        if diss:
            dis = -2.0 * math.sqrt(g * t0) / dx * math.sin(th / 2) ** 2
        Lh = np.zeros((5, 5), complex)
        Lh[0, 1] = -dc * r0
        Lh[1, 0] = -dc * t0 / r0
        Lh[1, 2] = -dc
        Lh[2, 1] = -dc * (g - 1.0) * t0
        for q in range(3):
            Lh[q, q] += dis
        Lr = np.zeros((5, 5), complex)
        Lr[1, 4] = kap / c * av
        Lr[2, 3] = (g - 1.0) * c * kap
        Lr[2, 2] = -(g - 1.0) * c * kap * b
        Lr[3, 4] = -s
        Lr[3, 3] = -c * r0 * kap
        Lr[3, 2] = c * r0 * kap * b
        Lr[4, 3] = -s * c * c / 3.0
        Lr[4, 4] = -c * r0 * kap
        Le = np.zeros((5, 5), complex)
        up = av if enth == "cen" else np.exp(-0.5j * th)
        Le[3, 1] = -(4.0 / 3.0) * e0 * s * up
        self.Lh, self.Lr, self.Le = Lh, Lr, Le
        self.L = Lh + Lr + Le


def heun(o, q, dt):
    y1 = q + dt * (o.Lh @ q)
    return 0.5 * q + 0.5 * (y1 + dt * (o.Lh @ y1))


def stage(o, mode, rhs, dt):
    g = G
    if mode == "imp":
        y = np.linalg.solve(I5 - g * dt * (o.Lr + o.Le), rhs)
    else:                                   # rhs: v lagged at the stage's old state
        y = np.linalg.solve(I5 - g * dt * o.Lr, rhs + g * dt * (o.Le @ rhs))
    return y, (y - rhs) / (g * dt)


def step(o, mode, q, k1, dt):
    if mode == "code":
        qh = heun(o, q, dt)
        return np.linalg.solve(I5 - dt * o.Lr, qh + dt * (o.Le @ qh)), k1
    if mode == "codeimp":
        qh = heun(o, q, dt)
        return np.linalg.solve(I5 - dt * (o.Lr + o.Le), qh), k1
    lh = o.Lh
    s1 = q + dt * (lh @ q)
    y2, k2 = stage(o, mode, s1 + (1 - G) * dt * k1, dt)
    s2 = 0.5 * q + 0.5 * (y2 + dt * (lh @ y2))
    y3, k3 = stage(o, mode, s2 + dt * (0.5 * G * k1 + (B2 - 0.5 * G) * k2), dt)
    return y3, k3


def smap(o, mode, dt):
    M = np.zeros((10, 10), complex)
    for j in range(10):
        e = np.zeros(10, complex)
        e[j] = 1
        a, b = step(o, mode, e[:5], e[5:], dt)
        M[:5, j], M[5:, j] = a, b
    return M


def semidisc_root(o, wcont):
    om = 1j * np.linalg.eigvals(o.L)
    return om[np.argmin(np.abs(om - wcont))]


def err(o, mode, dt, wref):
    lam = np.linalg.eigvals(smap(o, mode, dt))
    lm = lam[np.argmin(np.abs(lam - np.exp(-1j * wref * dt)))]
    return abs(1j * np.log(lm) / dt - wref) / abs(wref)


MODES = ("imp", "rhs", "code", "codeimp")


def part1():
    print("Part 1: temporal order vs the semi-discrete root, N = 64 cells/wavelength")
    print("  e128 / e1024 / p(1024->2048) per (P, tau_lambda)")
    nx = 64
    for enth in ("cen", "upw"):
        for mode in MODES:
            print("-- enthalpy %s, scheme %s" % (enth, mode))
            for pr, tau in tm.PT:
                bg = rd.Bg(pr, tau)
                wc, _ = rd.moment_root(bg, 0.0, "f0")
                o = DOps(bg, 2 * math.pi / nx, 1.0 / nx, enth)
                wref = semidisc_root(o, wc)
                per = 2 * np.pi / wref.real
                e = [err(o, mode, per / nt, wref) for nt in (128, 1024, 2048)]
                print("  P=%6g tau=%7g  %.2e %.2e %.2f" % (pr, tau, e[0], e[1],
                                                           math.log2(e[1] / e[2])))


def part2():
    print("\nPart 2: max over theta = 2 pi m/64 (m = 1..32), tau in 1,1e2,1e4,1e6 "
          "(per unit length, dx = 1/64) of rho(G) - 1")
    print("dt = cfl dx / c_s,gas.  'full' = with the Rusanov hydro rows; 'Lh=0' as C3")
    nx = 64
    dx = 1.0 / nx
    for enth in ("cen", "upw"):
        for lh in ("full", "Lh=0"):
            for cfl in (0.15, 0.3, 0.8):
                line = "%s %-4s cfl %.2f " % (enth, lh, cfl)
                for pr in (1.0, 10.0, 100.0):
                    w = {md: -1.0 for md in MODES}
                    for tau in (1.0, 1e2, 1e4, 1e6):
                        bg = rd.Bg(pr, 1.0)
                        bg.kappa = tau / bg.rho0
                        dt = cfl * dx / math.sqrt(bg.gamma * bg.t0)
                        for m in range(1, nx // 2 + 1):
                            o = DOps(bg, 2 * math.pi * m / nx, dx, enth)
                            if lh == "Lh=0":
                                o.Lh = np.zeros_like(o.Lh)
                            for md in MODES:
                                r = max(abs(np.linalg.eigvals(smap(o, md, dt)))) - 1
                                w[md] = max(w[md], r)
                    line += "| P=%g " % pr + " ".join("%s %.1e" % (md, w[md])
                                                      for md in MODES)
                print(line)


# ------------------------------------------------------------------ part 3: iterations
def build2d(bg, a, dx, n1, n2, enth, kvar=None):
    """E-only operator after eliminating T, F0 (faces), v on an n1 x n2 periodic grid,
    for a stage solve (I - a L).  Returns (M, M_old) as scipy sparse: M_old is the
    present operator without the v' enthalpy terms (the existing line preconditioner's
    source).  kvar: optional per-cell kappa multiplier (n2, n1)."""
    import scipy.sparse as sp
    c, e0, g, t0 = bg.c, bg.e0, bg.gamma, bg.t0
    b = 4.0 * bg.arad * t0 ** 3
    n = n1 * n2
    kap = bg.kappa * (np.ones((n2, n1)) if kvar is None else kvar)

    def idx(i, j):
        return (j % n2) * n1 + (i % n1)

    def shift(d):   # periodic shift operator along axis d by +1: (S x)_c = x_{c+e_d}
        r, cc = [], []
        for j in range(n2):
            for i in range(n1):
                r.append(idx(i, j))
                cc.append(idx(i + 1, j) if d == 0 else idx(i, j + 1))
        return sp.csr_matrix((np.ones(n), (r, cc)), shape=(n, n))
    kc = kap.ravel()
    Id = sp.identity(n, format="csr")
    # diagonal: 1 + a c kappa (1 - b beta), beta from the T elimination
    beta = a * (g - 1.0) * c * kc / (1.0 + a * (g - 1.0) * c * kc * b)
    diag = 1.0 + a * c * kc * (1.0 - b * beta)
    M0 = sp.diags(diag)
    Menth = sp.csr_matrix((n, n))
    for d in (0, 1):
        S = shift(d)
        St = S.T.tocsr()
        # face f = i+1/2 stored at cell i: face kappa = mean of the two cells
        kf = 0.5 * (kc + S @ kc)
        th = 1.0 / (1.0 + a * c * kf)
        Gr = (S - Id) / dx                        # cell -> face gradient
        Dv = (Id - St) / dx                       # face -> cell divergence
        Fop = sp.diags(-th * a * c * c / 3.0) @ Gr          # F0' = Fop E' (+ rhs)
        M0 = M0 + a * (Dv @ Fop)
        Avg = 0.5 * (Id + St)                     # face -> cell mean (v' from F0')
        Vop = (a / c) * Avg @ sp.diags(kf) @ Fop  # v' = Vop E' (+ rhs), the dm1 rule
        if enth == "cen":
            Up = 0.5 * (Id + S)                   # cell -> face mean
        else:
            Up = Id                               # left cell
        Menth = Menth + a * (4.0 / 3.0) * e0 * (Dv @ Up @ Vop)
    return (M0 + Menth).tocsr(), M0.tocsr()


def line_prec(M, n1, n2, width):
    """block-diagonal x1-line part of M, bandwidth `width` (1 = tri, 2 = penta)"""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl
    Mc = M.tocoo()
    i1r, jr = Mc.row % n1, Mc.row // n1
    i1c, jc = Mc.col % n1, Mc.col // n1
    dd = (i1c - i1r) % n1
    dd = np.where(dd > n1 // 2, dd - n1, dd)
    keep = (jr == jc) & (np.abs(dd) <= width)
    P = sp.csc_matrix((Mc.data[keep], (Mc.row[keep], Mc.col[keep])), shape=M.shape)
    return spl.splu(P)


def bicgstab(M, b, prec, tol=1e-10, maxit=500):
    x = np.zeros_like(b)
    r = b.copy()
    rh = r.copy()
    rho = al = om = 1.0
    v = p = np.zeros_like(b)
    nb = np.max(np.abs(b))
    for it in range(1, maxit + 1):
        rn = rh @ r
        be = (rn / rho) * (al / om)
        rho = rn
        p = r + be * (p - om * v)
        ph = prec.solve(p)
        v = M @ ph
        al = rho / (rh @ v)
        s = r - al * v
        if np.max(np.abs(s)) <= tol * nb:
            return it - 0.5
        sh = prec.solve(s)
        t = M @ sh
        om = (t @ s) / (t @ t)
        x = x + al * ph + om * sh
        r = s - om * t
        if np.max(np.abs(r)) <= tol * nb:
            return it
    return float("nan")


def part3():
    print("\nPart 3: BiCGStab iterations (right line preconditioner along x1, "
          "max-norm tol 1e-10)")
    print("2-D 64 x 64 periodic, dx = 1/64, Eddington, c = 1e3; stage factor a = g dt "
          "(H-ESDIRK2) and a = dt (BE)")
    print("prec: old = present tridiagonal (no v' terms), tri = tridiagonal part of the "
          "new operator, penta = +-2 x1 terms added; 'none' = operator without v' and "
          "old prec (today's count)")
    n1 = n2 = 64
    dx = 1.0 / n1
    rng = np.random.default_rng(1)
    bvec = rng.standard_normal(n1 * n2)
    kv = np.exp(rng.uniform(-1.15, 1.15, (n2, n1)))   # x10 random kappa spread
    for enth in ("cen", "upw"):
        for gname, gf in (("hesdirk2", G), ("BE", 1.0)):
            for cfl in (0.15, 0.3):
                for pr in (1.0, 10.0, 100.0):
                    out = []
                    for taucell in (1e-2, 1.0, 1e2, 1e4):
                        bg = rd.Bg(pr, 1.0)
                        bg.kappa = taucell / dx
                        dt = cfl * dx / math.sqrt(bg.gamma * bg.t0)
                        M, M0 = build2d(bg, gf * dt, dx, n1, n2, enth, kvar=kv)
                        p1 = line_prec(M0, n1, n2, 1)
                        n_none = bicgstab(M0, bvec, p1)
                        n_old = bicgstab(M, bvec, p1)
                        n_tri = bicgstab(M, bvec, line_prec(M, n1, n2, 1))
                        n_pen = bicgstab(M, bvec, line_prec(M, n1, n2, 2))
                        out.append("tc=%g %g/%g/%g/%g"
                                   % (taucell, n_none, n_old, n_tri, n_pen))
                    print("%s %-8s cfl %.2f P=%-5g none/old/tri/penta: %s"
                          % (enth, gname, cfl, pr, "  ".join(out)))


if __name__ == "__main__":
    what = sys.argv[1:] or ["order", "stab", "iters"]
    if "order" in what:
        part1()
    if "stab" in what:
        part2()
    if "iters" in what:
        part3()
