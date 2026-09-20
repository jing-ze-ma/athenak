#!/usr/bin/env python3
"""tests_r12: put the He4 ic_profile in equilibrium WITH THE CODE'S OWN two-stream.

Why (tests_r11, 2026-09-18/19): ic_he4_presn_sph.txt keeps the column's rho and total
pressure, so its table T is no longer the grey solution; the run's t = 0 two-stream flux
is 1.04 / 1.21 / 1.28 of L/(4 pi r^2) at 0.99 / 1.00 / 1.01 R.  With rt_rad_force the
momentum equation carries (1-w) rho kappa F_2s/c there, so the top ~3 % of the radius
starts with up to ~20 % of g of net outward force and loses 0.5-3 % of L: the secular
expansion no damping removes (damp3).

Each pass:
  1. PROBE: a one-cycle run of the 1-D column writes column_<tag>.txt (the run's own T,
     rho, kappa, tau on the fine grid) and mltfaces_<tag>.txt (F_req, F_2s on the faces).
  2. UNSOLD flux correction of T above R_A:  d(sigma T^4) = 3/4 int_0^tau dF dtau' +
     1/2 dF(0),  dF = F_req - F_2s, under-relaxed, windowed to zero below R_A.
  3. HYDROSTATIC re-solve of rho above R_A with the run's force law,
       d p_gas/dr = -rho g - w dPrad/dr + (1-w) rho kappa F_req/c,  Prad = a T^4/3,
     (the Prad grad w pieces cancel between the EOS and rt_rad_force), p_gas and e_gas
     from the EOS table dump, kappa from the Rosseland table, w = rad_taper WeightGated.
  4. write ic_pass<N>.txt (r rho eint), eint = e_gas + w a T^4.

Usage: relax_ic.py <npass> [--start ic_file] [--first N]
"""
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B = '/viper/u2/jinma/ATHENAK/bench/wt_he4'
IC0 = '/viper/u2/jinma/ATHENAK/bench/hestar_presn/ic_he4_presn_sph.txt'
DUMP = ('/viper/u2/jinma/ATHENAK/bench/analysis_0916/iso3d_w5/hestar/'
        'eos_table_box_w5.txt')
OPAC = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/rosseland_he_x0.0_z0.02.txt'
RSTAR, LSTAR, GM = 2.3717e11, 2.3066e38, 6.674e-8*6.2651e33
A_RAD, CLIGHT, SIGSB = 7.5657332503e-15, 2.99792458e10, 5.670374419e-5
RHO_HI, RHO_LO, T_HI, T_LO = 1.93e-9, 4.97e-10, 6.3281e4, 4.5213e4
R_A, R_W = 0.935*RSTAR, 0.90*RSTAR     # correct above R_A, window to zero at R_W
RELAX = 0.35


def read_dump(fn):
    with open(fn) as fh:
        lines = [fh.readline() for _ in range(6)]
    nx, ny, xmin, dx, ymin, dy = [float(v) for v in lines[2].split()[1:]]
    nx, ny = int(nx), int(ny)
    d = np.loadtxt(fn, comments='#')
    return (xmin + dx*np.arange(nx), ymin + dy*np.arange(ny),
            d[:, 0].reshape(ny, nx), d[:, 1].reshape(ny, nx))


def bilin(x, y, f, xq, yq):
    ix = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    iy = np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2)
    tx = (xq - x[ix])/(x[ix+1] - x[ix])
    ty = (yq - y[iy])/(y[iy+1] - y[iy])
    return ((1-tx)*(1-ty)*f[iy, ix] + tx*(1-ty)*f[iy, ix+1]
            + (1-tx)*ty*f[iy+1, ix] + tx*ty*f[iy+1, ix+1])


def load_opac(fn):
    hdr, vals = [], []
    for ln in open(fn):
        (hdr if ln.startswith('#') else vals).append(ln)
    g = [h for h in hdr if h.strip().startswith('# 2')]
    nT, nD, lTmin, dlT, lDmin, dlD = [float(v) for v in g[-1].strip('# \n').split()]
    k = np.array([float(v) for v in vals]).reshape(int(nT), int(nD))
    return k, lTmin, dlT, int(nT), lDmin, dlD, int(nD)


EX, EY, ELE, ELP = read_dump(DUMP)
OP = load_opac(OPAC)


def kappa(T, rho):
    k, lTmin, dlT, nT, lDmin, dlD, nD = OP
    x = np.clip((np.log10(T) - lTmin)/dlT, 0, nT - 1.000001)
    y = np.clip((np.log10(rho) - lDmin)/dlD, 0, nD - 1.000001)
    i0, j0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - i0, y - j0
    return 10.0**((1-fx)*(1-fy)*k[i0, j0] + fx*(1-fy)*k[i0+1, j0]
                  + (1-fx)*fy*k[i0, j0+1] + fx*fy*k[i0+1, j0+1])


def smooth(s):
    s = np.clip(s, 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def taper_w(rho, T):
    wr = smooth((np.log10(rho) - np.log10(RHO_LO))/(np.log10(RHO_HI) - np.log10(RHO_LO)))
    wt = smooth((np.log10(T) - np.log10(T_LO))/(np.log10(T_HI) - np.log10(T_LO)))
    return np.maximum(wr, wt)


def pgas(rho, T):
    return rho*10.0**bilin(EX, EY, ELP, np.log10(rho), np.log10(T))


def egas(rho, T):
    return rho*10.0**bilin(EX, EY, ELE, np.log10(rho), np.log10(T))


def rho_of(pg, T, guess):
    lo, hi = np.log10(guess) - 1.5, np.log10(guess) + 1.5
    for _ in range(60):
        mid = 0.5*(lo + hi)
        if pgas(10.0**mid, T) < pg:
            lo = mid
        else:
            hi = mid
    return 10.0**(0.5*(lo + hi))


def temp_of(rho, eint):
    """T from (rho, eint = e_gas + w a T^4), the run's tapered EOS, by bisection."""
    lo = np.full_like(rho, EY[0] + 1e-9)
    hi = np.full_like(rho, EY[-1] - 1e-9)
    for _ in range(80):
        mid = 0.5*(lo + hi)
        T = 10.0**mid
        f = egas(rho, T) + taper_w(rho, T)*A_RAD*T**4 - eint
        lo = np.where(f < 0.0, mid, lo)
        hi = np.where(f < 0.0, hi, mid)
    return 10.0**(0.5*(lo + hi))


def probe(ic, tag):
    arm = os.path.join(HERE, tag)
    log = os.path.join(HERE, tag + '.log')
    out = subprocess.run(['sbatch', '--parsable', '-o', log,
                          os.path.join(HERE, 'r12_arm.sh'), tag, 'time/nlim=1',
                          'time/cfl_number=0.15', 'problem/rt_top_vacuum=true',
                          'problem/ic_profile=' + ic],
                         capture_output=True, text=True, cwd=HERE)
    jid = out.stdout.strip().split(';')[0]
    while subprocess.run(['squeue', '-h', '-j', jid], capture_output=True,
                         text=True).stdout.strip():
        time.sleep(20)
    mf = np.loadtxt(os.path.join(arm, 'mltfaces_%s.txt' % tag), comments='#')
    col = np.loadtxt(os.path.join(arm, 'column_%s.txt' % tag), comments='#')
    return mf, col


def main():
    npass = int(sys.argv[1])
    ic = sys.argv[sys.argv.index('--start') + 1] if '--start' in sys.argv else IC0
    first = int(sys.argv[sys.argv.index('--first') + 1]) if '--first' in sys.argv else 1
    for n in range(first, first + npass):
        tag = 'p%02d' % n
        mf, col = probe(ic, tag)
        rfc, freq, f2s, fmlt = mf[:, 1], mf[:, 10], mf[:, 17], mf[:, 9]
        q = (f2s + fmlt)/freq
        top = rfc > 0.93*RSTAR
        print('PASS %d  ic = %s' % (n, os.path.basename(ic)))
        print('  (F_2s+F_mlt)/F_req on faces above 0.93 R: max |q-1| = %.4f ; at'
              ' 0.96/0.98/0.99/1.00/1.01 R: %s'
              % (np.abs(q[top] - 1).max(),
                 ' '.join('%.4f' % q[int(np.argmin(np.abs(rfc/RSTAR - s)))]
                          for s in (0.96, 0.98, 0.99, 1.00, 1.01))))
        sys.stdout.flush()

        d = np.loadtxt(ic, comments='#')
        r, rho, eint = d[:, 0], d[:, 1], d[:, 2]
        T = temp_of(rho, eint)
        # the run's own tau on the fine grid (column dump: r p T rho kappa grad tau)
        o = np.argsort(col[:, 0])
        tau = np.interp(r, col[o, 0], col[o, 6])
        Tc = np.interp(r, col[o, 0], col[o, 2])
        kc = np.interp(r, col[o, 0], col[o, 4])
        m = r > R_W
        print('  offline T vs the run\'s dump T above 0.90 R: max |dT/T| = %.2e ; kappa'
              ' table vs dump: max |dk/k| = %.2e'
              % (np.abs(T[m]/Tc[m] - 1).max(),
                 np.abs(kappa(T[m], rho[m])/kc[m] - 1).max()))

        # ---- Unsold correction, integrated DOWN from the top in the run's tau
        freq_r = LSTAR/(4*np.pi*r**2)
        dF = freq_r*(1.0 - np.interp(r, rfc, q))
        dF[r < R_A] = 0.0
        idx = np.argsort(tau)                      # tau ascending = top first
        integ = np.zeros_like(r)
        integ[idx] = np.concatenate(([0.0], np.cumsum(
            0.5*(dF[idx][1:] + dF[idx][:-1])*np.diff(tau[idx]))))
        dB = 0.75*integ + 0.5*dF[idx][0]
        win = smooth((r - R_W)/(R_A - R_W))
        T4 = SIGSB*T**4 + RELAX*win*dB
        Tn = (np.maximum(T4, 0.05*SIGSB*T**4)/SIGSB)**0.25
        print('  T correction: top %.4f, at 0.99 R %.4f, at 0.96 R %.4f (T_new/T_old)'
              % (Tn[-1]/T[-1], (Tn/T)[int(np.argmin(np.abs(r/RSTAR - 0.99)))],
                 (Tn/T)[int(np.argmin(np.abs(r/RSTAR - 0.96)))]))

        # ---- hydrostatic re-solve of rho above R_W (Heun, in ln p_gas)
        rn = rho.copy()
        i0 = int(np.searchsorted(r, R_W))
        prad = A_RAD*Tn**4/3.0
        dprad = np.gradient(prad, r)

        def rhs(i, rh):
            w = taper_w(np.array([rh]), np.array([Tn[i]]))[0]
            kap = kappa(np.array([Tn[i]]), np.array([rh]))[0]
            return (-rh*GM/r[i]**2 - w*dprad[i]
                    + (1.0 - w)*rh*kap*LSTAR/(4*np.pi*r[i]**2)/CLIGHT)
        pg = pgas(np.array([rn[i0]]), np.array([Tn[i0]]))[0]
        for i in range(i0, len(r) - 1):
            h = r[i+1] - r[i]
            k1 = rhs(i, rn[i])
            pe = max(pg + h*k1, 1e-3*pg)
            re = float(rho_of(pe, Tn[i+1], rn[i]))
            k2 = rhs(i+1, re)
            pg = max(pg + 0.5*h*(k1 + k2), 1e-3*pg)
            rn[i+1] = float(rho_of(pg, Tn[i+1], re))
        print('  rho_new/rho_old: at 0.96 R %.4f, 0.99 R %.4f, 1.00 R %.4f, top %.4f ;'
              ' top rho %.3e'
              % tuple([(rn/rho)[int(np.argmin(np.abs(r/RSTAR - s)))]
                       for s in (0.96, 0.99, 1.00)] + [rn[-1]/rho[-1], rn[-1]]))
        en = egas(rn, Tn) + taper_w(rn, Tn)*A_RAD*Tn**4
        en[:i0] = eint[:i0]
        rn[:i0] = rho[:i0]
        ic = os.path.join(HERE, 'ic_pass%02d.txt' % n)
        with open(ic, 'w') as fh:
            fh.write('# tests_r12/relax_ic.py pass %d: r[cm] rho[g/cm^3] eint[erg/cm^3]\n'
                     % n)
            for a, b_, c_ in zip(r, rn, en):
                fh.write('%.10e %.10e %.10e\n' % (a, b_, c_))
        print('  wrote', ic)
        sys.stdout.flush()


if __name__ == '__main__':
    main()
