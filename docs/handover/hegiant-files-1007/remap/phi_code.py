#!/usr/bin/env python3
"""The he_star_m1 he_gm_column potential, replicated: the IC column (he_ic_file) resampled in log onto the pgen's
uniform fine grid (rlo = x1f(0) - 2 dx_lo, rhi = x1f(top ghost face) + 2 dx_hi, nf = he_nfine), G m(r) = he_gm +
G int 4 pi r^2 rho dr from r_in (trapezoid), Phi = int G m/r^2 dr from r_in (trapezoid), then linear interpolation
at the cell volume centroids x1v (HsLinInterp).  phi_at(faces, ...) -> Phi at every cell incl. ghosts."""
import numpy as np

G_CODE = 6.67408e-8


def phi_at(faces, ic_file, he_gm, rin, nf):
    """faces: all x1 faces incl. ghosts (n1 + 2 ng + 1); returns (x1v, Phi(x1v))"""
    d = np.loadtxt(ic_file, usecols=(0, 1))
    fr, fd = d[:, 0], d[:, 1]
    dxlo = faces[1] - faces[0]
    dxhi = faces[-1] - faces[-2]
    # the pgen uses x1f(0) and x1f(n1m1+1) (the outermost ghost faces) and the FIRST/LAST active widths;
    # with ng ghosts: dxlo = x1f(is+1) - x1f(is), dxhi = x1f(ie+1) - x1f(ie)
    return faces, fr, fd, dxlo, dxhi


def phi_cells(faces, ng, ic_file, he_gm, nf):
    d = np.loadtxt(ic_file, usecols=(0, 1))
    fr, fd = d[:, 0], d[:, 1]
    n1 = faces.size - 1 - 2 * ng
    rin = faces[ng]
    dxlo = faces[ng + 1] - faces[ng]
    dxhi = faces[ng + n1] - faces[ng + n1 - 1]
    rlo = faces[0] - 2.0 * dxlo
    rhi = faces[-1] + 2.0 * dxhi
    dr = (rhi - rlo) / (nf - 1)
    hr = rlo + dr * np.arange(nf)
    k = np.clip(np.searchsorted(fr, hr, side='left') - 1, 0, fr.size - 2)
    w = (hr - fr[k]) / (fr[k + 1] - fr[k])
    hd = np.exp((1.0 - w) * np.log(fd[k]) + w * np.log(fd[k + 1]))
    cm = np.concatenate([[0.0], np.cumsum(0.5 * dr * 4.0 * np.pi * (hr[:-1]**2 * hd[:-1]
                                                                      + hr[1:]**2 * hd[1:]))])
    x = (rin - rlo) / dr
    i0 = int(min(max(np.floor(x), 0), nf - 2))
    w0 = x - i0
    hgm = he_gm + G_CODE * (cm - ((1 - w0) * cm[i0] + w0 * cm[i0 + 1]))
    cp = np.concatenate([[0.0], np.cumsum(0.5 * dr * (hgm[:-1] / hr[:-1]**2 + hgm[1:] / hr[1:]**2))])
    phig = cp - ((1 - w0) * cp[i0] + w0 * cp[i0 + 1])
    f0, f1 = faces[:-1], faces[1:]
    x1v = 0.75 * (f1**4 - f0**4) / (f1**3 - f0**3)
    xx = (x1v - rlo) / dr
    ii = np.clip(np.floor(xx).astype(int), 0, nf - 2)
    ww = np.clip(xx - ii, 0.0, 1.0)
    return x1v, (1.0 - ww) * phig[ii] + ww * phig[ii + 1]
