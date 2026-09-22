#!/usr/bin/env python3
"""Discrete Gauss identity  sum_f sigma_f A_f nhat_f = 0  for a cubed-sphere cell,
as AthenaK builds A_f (exact areas) and nhat_f (face-CENTRE normals).
The identity is exact for the true closed surface, so what is printed is the code's
quadrature error: a cell for which it does not vanish leaks momentum."""
import numpy as np

A = np.array([0., 1., 0.])
B = np.array([0., 0., 1.])
N = np.array([1., 0., 0.])


def rhat(xi, eta):
    x, y = np.tan(xi), np.tan(eta)
    d = np.sqrt(1.0 + x*x + y*y)
    return (A*x + B*y + N)/d


def tangents(xi, eta):
    x, y = np.tan(xi), np.tan(eta)
    d = np.sqrt(1.0 + x*x + y*y)
    r = (A*x + B*y + N)/d
    e1 = (A*d - r*x)/np.sqrt(1.0 + y*y)
    e2 = (B*d - r*y)/np.sqrt(1.0 + x*x)
    return e1, e2


def nxi(xi, eta):
    e1, e2 = tangents(xi, eta)
    c = e1@e2
    return (e1 - c*e2)/np.sqrt(1.0 - c*c)


def neta(xi, eta):
    e1, e2 = tangents(xi, eta)
    c = e1@e2
    return (e2 - c*e1)/np.sqrt(1.0 - c*c)


def solid_angle(xl, xr, yl, yr):
    """exact solid angle of the gnomonic patch (tan coords), Gauss-Legendre in xi,eta"""
    gx, gw = np.polynomial.legendre.leggauss(60)
    xis = 0.5*(xr-xl)*gx + 0.5*(xr+xl)
    ets = 0.5*(yr-yl)*gx + 0.5*(yr+yl)
    tot = 0.0
    for a, wa in zip(xis, gw):
        for b, wb in zip(ets, gw):
            x, y = np.tan(a), np.tan(b)
            d2 = 1.0 + x*x + y*y
            jac = (1+x*x)*(1+y*y)/d2**1.5
            tot += wa*wb*jac
    return tot*0.25*(xr-xl)*(yr-yl)


def arc(xa, xb, yy):   # angle between (1,xa,yy) and (1,xb,yy)
    return np.arccos((1+xa*xb+yy*yy)/np.sqrt(1+xa*xa+yy*yy)/np.sqrt(1+xb*xb+yy*yy))


def cell(xil, xir, etl, etr, rl, rr):
    xic, etc = 0.5*(xil+xir), 0.5*(etl+etr)
    dom = solid_angle(xil, xir, etl, etr)
    S = np.zeros(3)                      # code's quadrature
    # radial faces: normal = rhat(cell centre angles), area = r^2 dOmega
    S += rr*rr*dom*rhat(xic, etc)
    S -= rl*rl*dom*rhat(xic, etc)
    # xi faces: area = 0.5(rr^2-rl^2)*arc_eta at that xi face
    hr = 0.5*(rr*rr - rl*rl)
    S += hr*arc(np.tan(etl), np.tan(etr), np.tan(xir))*nxi(xir, etc)
    S -= hr*arc(np.tan(etl), np.tan(etr), np.tan(xil))*nxi(xil, etc)
    # eta faces
    S += hr*arc(np.tan(xil), np.tan(xir), np.tan(etr))*neta(xic, etr)
    S -= hr*arc(np.tan(xil), np.tan(xir), np.tan(etl))*neta(xic, etl)

    return S, dom, hr


for n in (16, 32, 64):
    d = 0.5*np.pi/2/n     # half cell width in xi (panel spans [-pi/4,pi/4])
    dd = np.pi/2/n
    print(f"--- n = {n}, dxi = {dd:.5f} ---")
    for name, (jc, kc) in (("INTERIOR(centre)", (0, 0)),
                           ("SEAM   (xi edge)", (n-1, 0)),
                           ("VERTEX (corner) ", (n-1, n-1))):
        xil = -np.pi/4 + jc*dd if jc else -0.5*dd
        etl = -np.pi/4 + kc*dd if kc else -0.5*dd
        xir, etr = xil+dd, etl+dd
        S, dom, hr = cell(xil, xir, etl, etr, 1.0, 1.0+dd)
        scale = dom*1.0 + 4*hr*dd   # rough total surface area
        print(f"  {name} |sum A n| code = {np.linalg.norm(S):.6e}"
              f"  rel(area) = {np.linalg.norm(S)/scale:.6e}")
