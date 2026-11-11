#!/usr/bin/env python3
"""rw_p333-1112: acoustic root of the P 3.33 tau 10 cell under model variants of the reference."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import radwave_disp as rd  # noqa: E402

bg = rd.Bg(10.0 / 3.0, 10.0, lam=1.0, c=10.0)
quad = rd.quad_code(4, 8, (1, 0, 0))


def root_L(L):
    lam, vec = np.linalg.eig(L)
    om = 1j * lam
    w = rd.pick_acoustic(bg, om)
    return w


def L_vet(qT=0.0, qv=0.0, f_qs=False, scale=1.0):
    """lab moment system with P1 = E1/3 + qT T1 + qv v1."""
    L = rd.matrix(bg, scale * qT)
    k, c = bg.k, bg.c
    L[4, 1] += -1j * k * c * c * qv
    return L


chi = bg.rho0 * bg.kappa
b = 4.0 * bg.arad
j0, j1, j2, j3 = rd.jfun_quad(chi + 0j, bg.k, *quad)
qT = chi * b * (j2 - j0 / 3.0)
# static formal solution WITH the O(v/c) Doppler terms: I1 = chi[b T1 + 4 E0 mu v1/c]/(chi+ik mu)
# (source S(1+3 mu beta) - I(1 - mu beta), I0 = E0: +3 mu beta E0 + mu beta E0 = 4 mu beta E0)
qv = chi * 4.0 * bg.e0 * (j3 - j1 / 3.0) / bg.c
print("vetqs_d             ", root_L(L_vet(qT)))
print("vetqs_d + Doppler   ", root_L(L_vet(qT, qv)))
for s in (0.99, 1.01):
    print("qT x %.2f          " % s, root_L(L_vet(qT, scale=s)))
print("eddington           ", root_L(L_vet(0.0)))
print("qT", qT, "qv", qv)
