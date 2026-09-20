#!/usr/bin/env python3
"""tests_r12/relax2.py: thermally pre-relaxed IC (2026-09-19).

A FROZEN probe (problem/vdamp_all_time=0.05, rho fixed) lets the run pull T to its own
discrete radiative equilibrium: dT/T = +2e-4 (0.91 R) ... -1.2e-3 (0.96 R) within ~50 s.
With p ~ T^4 that is the 0.5-3 % of rho g dipole of tests_r12/fbud.py.  Here:
  1. dlnT(r) = T_probe(t_end)/T_probe(0) - 1 on the cells, interpolated to the IC grid,
     zero below R_LO (and optionally above R_HI, --rhi);
  2. hydrostatic re-solve of rho above R_LO as a DIFFERENCE against the old IC, so the
     offline force law's own error cancels:
       pg_new[i+1]-pg_new[i] = pg_old[i+1]-pg_old[i] + int (rhs_new - rhs_old) dr
  3. write ic_rx<tag>.txt (r rho eint).
Usage: relax2.py <probe_arm> <out_tag> [--start ic] [--rlo 0.62] [--rhi 1.1] [--gain 1]
"""
import sys

import numpy as np

import relax_ic as R


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv \
        else default


def rd(fn):
    import struct
    out = []
    f = open(fn, 'rb')
    while True:
        h = f.read(16)
        if len(h) < 16:
            break
        t, n1, nv = struct.unpack('<dii', h)
        b = f.read(8*(n1 + nv*n1))
        if len(b) < 8*(n1 + nv*n1):
            break
        a = np.frombuffer(b, '<f8')
        out.append((t, a[:n1].copy(), a[n1:].reshape(nv, n1).copy()))
    return out


def main():
    arm, tag = sys.argv[1], sys.argv[2]
    ic = arg('--start', R.IC0)
    rlo, rhi = arg('--rlo', 0.62)*R.RSTAR, arg('--rhi', 1.1)*R.RSTAR
    gain = arg('--gain', 1.0)
    P = rd(arm + '/rt_profile.bin')
    rc = P[0][1]
    kd = arg("--idump", -1)
    dl = P[kd][2][5]/P[0][2][5] - 1.0
    dprev = P[-3][2][5]/P[0][2][5] - 1.0
    print('probe %s: t = %.0f s ; max |dlnT| %.2e ; still moving (last 2 dumps) %.2e'
          % (arm, P[kd][0], np.abs(dl).max(), np.abs(dl - dprev).max()))
    print('  max drho/rho in the probe %.2e' % np.abs(P[-1][2][0]/P[0][2][0] - 1).max())

    d = np.loadtxt(ic, comments='#')
    r, rho, eint = d[:, 0], d[:, 1], d[:, 2]
    T = R.temp_of(rho, eint)
    win = R.smooth((r - rlo)/(0.03*R.RSTAR))*(1.0 - R.smooth((r - rhi)/(0.01*R.RSTAR)))
    Tn = T*(1.0 + gain*win*np.interp(r, rc, dl))

    def rhs(rh, Tt, dpr, rr):
        w = R.taper_w(rh, Tt)
        return (-rh*R.GM/rr**2 - w*dpr
                + (1.0 - w)*rh*R.kappa(Tt, rh)*R.LSTAR/(4*np.pi*rr**2)/R.CLIGHT)
    dpo = np.gradient(R.A_RAD*T**4/3.0, r)
    dpn = np.gradient(R.A_RAD*Tn**4/3.0, r)
    pgo = R.pgas(rho, T)
    rho_o = rhs(rho, T, dpo, r)
    rn = rho.copy()
    i0 = int(np.searchsorted(r, rlo))
    pg = pgo[i0]
    one = np.ones(1)
    for i in range(i0, len(r) - 1):
        h = r[i+1] - r[i]
        k1 = rhs(rn[i]*one, Tn[i]*one, dpn[i]*one, r[i]*one)[0] - rho_o[i]
        base = pgo[i+1] - pgo[i]
        pe = max(pg + base + h*k1, 1e-3*pg)
        re = float(R.rho_of(pe, Tn[i+1], rn[i]))
        k2 = rhs(re*one, Tn[i+1]*one, dpn[i+1]*one, r[i+1]*one)[0] - rho_o[i+1]
        pg = max(pg + base + 0.5*h*(k1 + k2), 1e-3*pg)
        rn[i+1] = float(R.rho_of(pg, Tn[i+1], re))
    for s in (0.70, 0.85, 0.90, 0.925, 0.95, 0.97, 0.99, 1.00, 1.01):
        k = int(np.argmin(np.abs(r/R.RSTAR - s)))
        print('  %.3f R: T_new/T %.5f  rho_new/rho %.5f' % (s, Tn[k]/T[k], rn[k]/rho[k]))
    print('  top rho %.3e (old %.3e)' % (rn[-1], rho[-1]))
    en = R.egas(rn, Tn) + R.taper_w(rn, Tn)*R.A_RAD*Tn**4
    # below R_LO the file is untouched; above, write the DIFFERENCE of the offline eint
    # onto the old one so the offline EOS error cancels here too
    eo = R.egas(rho, T) + R.taper_w(rho, T)*R.A_RAD*T**4
    en = eint + (en - eo)
    out = R.HERE + '/ic_rx%s.txt' % tag
    with open(out, 'w') as fh:
        fh.write('# tests_r12/relax2.py %s from %s: r[cm] rho[g/cm^3] eint[erg/cm^3]\n'
                 % (tag, arm))
        for a, b_, c_ in zip(r, rn, en):
            fh.write('%.10e %.10e %.10e\n' % (a, b_, c_))
    print('  wrote', out)


if __name__ == '__main__':
    main()
