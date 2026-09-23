#!/usr/bin/env python3
"""H-ESDIRK2 with FSAL, treatments of the enthalpy term -(4/3) E0 ik v in a stage solve
Y = rhs + g dt (Lr Y + Le v*):
  imp   v* = v(Y)                     (fully implicit: needs v inside the linear solve)
  rhs   v* = v(rhs)                   (the code's V1: lagged at the solve's old gas state)
  pred  v* = v(rhs) + g dt Kv_prev    (rhs velocity + the force slope of the previous
                                       stage solve: K1 (FSAL) for Y2, K2 for Y3)
  exp   Le explicit with the Heun weights (design option 1)
and the present code (Heun gas-only, then BE, v* = v after the hydro).
FSAL is modelled exactly: K1 = (Y3 - rhs3)/(g dt) is carried as state.
Part 1: order (test A grid).  Part 2: max spectral radius - 1 of the full step map
(Lh included; and Lh = 0 as C3) over k = 2 pi m, m = 1..32 (N = 64), tau 1e2..1e6."""
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
MODES = ("imp", "rhs", "pred", "exp", "code")


def stage(o, mode, rhs, dt, kprev):
    g = G
    if mode == "imp":
        y = np.linalg.solve(I5 - g * dt * (o.Lr + o.Le), rhs)
    elif mode in ("rhs", "pred"):
        vs = rhs.copy()
        if mode == "pred":
            vs = vs + g * dt * kprev
        y = np.linalg.solve(I5 - g * dt * o.Lr, rhs + g * dt * (o.Le @ vs))
    else:
        y = np.linalg.solve(I5 - g * dt * o.Lr, rhs)
    return y, (y - rhs) / (g * dt)


def step(o, mode, q, k1, dt):
    if mode == "code":
        qh = tm.heun(o, q, dt, False, gas_only=True)
        rhs = qh + dt * (o.Le @ qh)
        return np.linalg.solve(I5 - dt * o.Lr, rhs), k1
    lh = o.Lh if mode != "exp" else o.Lh + o.Le
    s1 = q + dt * (lh @ q)
    y2, k2 = stage(o, mode, s1 + (1 - G) * dt * k1, dt, k1)
    s2 = 0.5 * q + 0.5 * (y2 + dt * (lh @ y2))
    y3, k3 = stage(o, mode, s2 + dt * (0.5 * G * k1 + (B2 - 0.5 * G) * k2), dt, k2)
    return y3, k3


def smap(o, mode, dt):
    M = np.zeros((10, 10), complex)
    for j in range(10):
        e = np.zeros(10, complex)
        e[j] = 1
        a, b = step(o, mode, e[:5], e[5:], dt)
        M[:5, j], M[5:, j] = a, b
    return M


def err(o, mode, dt, wref):
    lam = np.linalg.eigvals(smap(o, mode, dt))
    lm = lam[np.argmin(np.abs(lam - np.exp(-1j * wref * dt)))]
    return abs(1j * np.log(lm) / dt - wref) / abs(wref)


print("Part 1: order.  e128 / e1024 / p(1024->2048) per (P, tau)")
for mode in MODES:
    print("-- enthalpy", mode)
    for pr, tau in tm.PT:
        bg = rd.Bg(pr, tau)
        wref, _ = rd.moment_root(bg, 0.0, "f0")
        per = 2 * np.pi / wref.real
        o = tm.Ops(bg)
        e = [err(o, mode, per / nt, wref) for nt in (128, 1024, 2048)]
        print("  P=%6g tau=%7g  %.2e %.2e %.2f" % (pr, tau, e[0], e[1],
                                                   math.log2(e[1] / e[2])))

print("\nPart 2: max over k (m = 1..32, N = 64) and tau in 1e2,1e4,1e6 of rho(G) - 1")
print("dt = cfl dx / c_s,gas.  'full' = with the gas pressure; 'Lh=0' as C3")
for lh in ("full", "Lh=0"):
    for cfl in (0.15, 0.3, 0.8):
        line = "%-5s cfl %.2f " % (lh, cfl)
        for pr in (1.0, 10.0, 100.0):
            w = {md: -1.0 for md in MODES}
            for tau in (1e2, 1e4, 1e6):
                for m in (1, 2, 4, 8, 16, 32):
                    bg = rd.Bg(pr, tau, lam=1.0 / m)
                    bg.kappa = tau / bg.rho0
                    o = tm.Ops(bg)
                    if lh == "Lh=0":
                        o.Lh = np.zeros_like(o.Lh)
                    dt = cfl * (1.0 / 64) / math.sqrt(bg.gamma * bg.t0)
                    for md in MODES:
                        r = max(abs(np.linalg.eigvals(smap(o, md, dt)))) - 1
                        w[md] = max(w[md], r)
            line += " | P=%g " % pr + " ".join("%s %.1e" % (md, w[md]) for md in MODES)
        print(line)
