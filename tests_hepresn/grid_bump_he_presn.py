#!/usr/bin/env python3
"""Radial grid with LOCAL refinement at the FeCZ edges: the 8-coefficient poly stretch
plus two tanh bumps (mesh/f_stretch_r_b<n>_amp/_x/_w, coordinates/grid_stretch.hpp).
usage: grid_bump_he_presn.py BASE_NPY NXBASE RTOP DEPTH[,DEPTH2] WIDTH_R OUT_NPY [NX]
Target cell width = the base grid's dx(r) (BASE_NPY on NXBASE cells: it already meets the
cells-per-scale-height rules) times 1 - DEPTH sech^2((r - r_e)/WIDTH_R R) at the FeCZ
edges r_e = 0.636 R and 0.933 R
(WIDTH_R > 1: half width in base cells).  The relative width misfit is minimised: the poly
coefficients linearly, the bump parameters by least_squares.  OUT_NPY holds 14 numbers:
c1..c8, then (amp, x, w) of bump 1 and 2."""
import sys

import numpy as np
from scipy.optimize import least_squares

R = 2.3717e11
EDGES = (0.636, 0.933)


def mapping(p, xi, deriv=False):
    """u(xi) (or du/dxi) of StretchRPoly with p = c1..c8, (a, xb, w) x 2."""
    c, u = p[:8], (np.zeros_like(xi) if deriv else xi.copy())
    if deriv:
        u += 1.0
    for k in range(1, 9):
        u += c[k-1] * ((k * xi**(k-1) - (k+1) * xi**k) if deriv else xi**k * (1 - xi))
    for b in range(2):
        a, xb, w = p[8+3*b:11+3*b]
        if a == 0.0:
            continue
        if deriv:
            u += a * (1 / np.cosh((xi - xb) / w)**2
                      - w * (np.tanh((1 - xb) / w) + np.tanh(xb / w)))
        else:
            u += a * w * (np.tanh((xi - xb) / w) - (1 - xi) * np.tanh(-xb / w)
                          - xi * np.tanh((1 - xb) / w))
    return u


def basis_d(xi):
    return np.stack([k * xi**(k-1) - (k+1) * xi**k for k in range(1, 9)], 1)


def main():
    cb = np.load(sys.argv[1])
    nb = int(sys.argv[2])
    rtop = float(sys.argv[3])
    deps = [float(v) for v in sys.argv[4].split(',')]
    deps = deps * 2 if len(deps) == 1 else deps
    wid, out = float(sys.argv[5]), sys.argv[6]
    r0 = 0.5 * R
    pb = np.concatenate([cb[:8], np.zeros(6)]) if len(cb) == 8 else cb
    eb = r0 + (rtop - r0) * mapping(pb, np.arange(nb + 1) / nb)
    rcb, dxb = 0.5 * (eb[1:] + eb[:-1]), np.diff(eb)
    rr = np.linspace(r0, rtop, 40001)
    dxt = np.exp(np.interp(rr, rcb, np.log(dxb)))
    wids = []
    for e, dep in zip(EDGES, deps):
        # WIDTH_R > 1: the half width in BASE cells at that edge (keeps the stretch
        # per cell the same at both edges)
        wids.append(wid * np.interp(e * R, rcb, dxb) / R if wid > 1 else wid)
        dxt *= 1 - dep / np.cosh((rr / R - e) / wids[-1])**2
    cnt = np.concatenate([[0], np.cumsum(np.diff(rr) / (0.5 * (dxt[1:] + dxt[:-1])))])
    nx = int(sys.argv[7]) if len(sys.argv) > 7 else int(round(cnt[-1]))
    xs = np.linspace(0, 1, 4001)
    ut = np.interp(xs, cnt / cnt[-1], (rr - r0) / (rtop - r0))
    dut = np.gradient(ut, xs)                      # target du/dxi
    Bd = basis_d(xs)

    def solve_c(q):
        p = np.concatenate([np.zeros(8), q])
        rhs = (dut - mapping(p, xs, deriv=True)) / dut
        c, *_ = np.linalg.lstsq(Bd / dut[:, None], rhs, rcond=None)
        return np.concatenate([c, q])

    def resid(q):
        p = solve_c(q)
        return mapping(p, xs, deriv=True) / dut - 1

    q0 = []
    for e, wd, dep in zip(EDGES, wids, deps):
        xb = np.interp(e * R, rr, cnt / cnt[-1])
        w = wd * R / (rtop - r0) / np.interp(xb, xs, dut)
        q0 += [-dep * np.interp(xb, xs, dut), xb, w]
    sol = least_squares(resid, np.array(q0),
                        bounds=([-5, 0, 1e-3, -5, 0, 1e-3], [0, 1, 0.5, 0, 1, 0.5]))
    p = solve_c(sol.x)
    print('target cells %.1f -> nx1 %d; max relative width misfit %.3f'
          % (cnt[-1], nx, np.abs(sol.fun).max()))
    assert (mapping(p, np.linspace(0, 1, 4097), deriv=True) > 0).all(), 'fold-over'
    print('c1..c8 =', ' '.join('%.10e' % v for v in p[:8]))
    for b in range(2):
        print('b%d amp %.10e x %.10e w %.10e' % ((b + 1,) + tuple(p[8+3*b:11+3*b])))
    np.save(out, p)
    print(nx)


if __name__ == '__main__':
    main()
