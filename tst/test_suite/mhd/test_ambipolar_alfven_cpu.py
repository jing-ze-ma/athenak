"""
Ambipolar-diffusion regression test (ambipolar-1101): a linear Alfven wave damped by AD.

1-D linear_wave, wave_flag = 5, background B0 = (1, sqrt2, 0.5), <mhd>/ambipolar =
constant with eta_A = eta_ad |B|^2.  For an Alfven wave J is perpendicular to k and to
dB and lies in the (k, B0) plane, so J.B0 != 0 and the AD projection
J - (J.B)B/B^2 matters: the wave damps with the diffusivity
    eta_eff = eta_A cos^2(k, B0) = eta_ad B_par^2      (= 0.02 here),
NOT eta_A = 0.0650 that an Ohmic-like E = eta_A J would give.  The exact linear system
for the Fourier amplitudes (v_e, B_e) along the polarisation e = k x B0 is
    dv/dt = i k B_par B/rho,  dB/dt = i k B_par v - eta_eff k^2 B,
whose dominant mode decays at Gamma = eta_eff k^2/2.  The test measures Gamma from the
first and last tab dump (projecting the state onto that mode) at nx1 = 32 and 64.

Measured on ambipolar-1101 (all-dump fit over 2 periods, RESULTS.md): relative error of
Gamma 2.8e-2 (32), 3.0e-3 (64), 2.8e-4 (128); this two-dump, one-period version
measures 2.76e-2 (32), 2.96e-3 (64), 2.67e-4 (128).  An Ohmic-like E = eta_A J (no
projection) would damp 3.25x faster and fails by two orders of magnitude.
"""

import numpy as np
import pytest
import test_suite.testutils as testutils

_input = "inputs/ambipolar_alfven.athinput"
_eta_eff = 0.02
_bx0, _by0, _bz0, _d0 = 1.0, 1.4142136, 0.5, 1.0
# |Gamma - Gamma_exact|/Gamma_exact upper bounds (measured value + margin)
_max_err = {32: 4.0e-2, 64: 5.0e-3}
_min_ratio = 4.0      # err(32)/err(64): second order or better


def read_tab(fname):
    with open(fname) as f:
        f.readline()
        names = f.readline().lstrip("#").split()
    data = np.loadtxt(fname)
    if data.shape[1] < len(names):
        # 1-D tables print no j/x2v/k/x3v columns, but the header still names them
        names = [n for n in names if n not in ("j", "x2v", "k", "x3v")]
    return {n: data[:, c] for c, n in enumerate(names)}


def amplitudes(fname, k):
    d = read_tab(fname)
    e = np.array([0.0, -_bz0, _by0])
    e /= np.linalg.norm(e)
    ph = np.exp(-1j*k*d["x1v"])
    ve = e[1]*d["vely"] + e[2]*d["velz"]
    be = e[1]*d["bcc2"] + e[2]*d["bcc3"]
    n = len(ph)
    return np.array([2.0*np.sum(ve*ph)/n, 2.0*np.sum(be*ph)/n])


def gamma(nx1):
    testutils.cleanup()
    testutils.run(_input, ["mesh/nx1=" + str(nx1), "meshblock/nx1=" + str(nx1)])
    k = 2.0*np.pi
    x0 = amplitudes("tab/AdAlfven.mhd_w_bcc.00000.tab", k)
    x1 = amplitudes("tab/AdAlfven.mhd_w_bcc.00001.tab", k)
    with open("tab/AdAlfven.mhd_w_bcc.00001.tab") as f:
        t1 = float(f.readline().split("time=")[1].split()[0])
    m = np.array([[0.0, 1j*k*_bx0/_d0], [1j*k*_bx0, -_eta_eff*k*k]])
    lam, r = np.linalg.eig(m)
    a0 = np.linalg.solve(r, x0)
    a1 = np.linalg.solve(r, x1)
    i = int(np.argmax(np.abs(a0)))
    g = -np.log(np.abs(a1[i])/np.abs(a0[i]))/t1
    return g, -lam[i].real


def test_ambipolar_alfven_damping():
    try:
        err = {}
        for n in (32, 64):
            g, gex = gamma(n)
            assert np.isfinite(g), "non-finite damping rate at nx1={:}".format(n)
            err[n] = abs(g - gex)/gex
            if err[n] > _max_err[n]:
                pytest.fail("AD Alfven damping at nx1={:}: Gamma={:.6e} exact={:.6e} "
                            "relerr={:.3e} > {:.1e}".format(n, g, gex, err[n],
                                                            _max_err[n]))
        if err[32]/err[64] < _min_ratio:
            pytest.fail("AD Alfven damping converges too slowly: err32/err64 = "
                        "{:.2f} < {:}".format(err[32]/err[64], _min_ratio))
    finally:
        testutils.cleanup()
