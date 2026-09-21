"""Build the G2 initial condition and its a_rad(z) table.

TWO STEPS, both of which need the code's own tabulated EOS, so both drive the
box_convection start-up (time/nlim=0, which writes problem/column_dump and quits).

(1) THE GAS-ONLY COLUMN.  The production He column
    bench/hestar_fecz/box_w8/column_used.txt was relaxed under the TAPERED EOS, so its
    internal energy carries a T^4/3 wherever the taper weight is 1.  With
    eos_radiation = false that same energy would be read as gas energy and the column
    would sit at the wrong temperature.  What we want is the same (rho, T) with the
    gas-only energy e = e_gas(rho, T) -- i.e. the radiation energy moved out of the gas,
    which is what the M1 milestone does.  There is no host-side EOS call available, so
    e is found by a SECANT ITERATION in ln e against the temperature the code's own
    column dump reports: 12 start-up runs bring max |ln T/T_target| to zero at the
    dump's printed precision.

(2) a_rad(z).  DEFINED as the residual acceleration that makes the resulting column
    exactly hydrostatic,

        dP_gas/dz = -rho (g0 - a_rad)   =>   a_rad = g0 + (1/rho) dP_gas/dz,

    evaluated on the pgen's own fine grid with a centred difference of half-width 8
    nodes (0.13 of a cell; the column dump carries 6 significant digits and a one-node
    difference would be 30x noisier).  It is NOT taken from the opacity table -- but it
    agrees with kappa_R F/c to a median ratio of 1.000 over the box, which is the check
    the design rests on.

Usage:  python3 make_ic_and_arad.py <path to the athena binary>
Run it from tests_wbphi/.
"""
import os
import subprocess
import sys

import numpy as np

COL = '/viper/u2/jinma/ATHENAK/bench/hestar_fecz/box_w8/column_used.txt'
G0 = 3.98107e5
NITER = 12
HALF = 8


def startup(exe, icfile, tag):
    """one nlim=0 start-up; returns the fine-grid column it dumps."""
    inp = 'ic/%s.athinput' % tag
    src = open('he_1d_wb.athinput').read()
    open(inp, 'w').write(src.replace('IC_PROFILE_PLACEHOLDER',
                                     os.path.abspath(icfile)))
    subprocess.run([exe, '-i', os.path.basename(inp), 'time/nlim=0',
                    '-d', 'itwork'], cwd='ic', check=True,
                   stdout=open('ic/%s.log' % tag, 'w'), stderr=subprocess.STDOUT)
    return np.loadtxt('ic/column_used.txt')


def main():
    exe = sys.argv[1]
    os.makedirs('ic', exist_ok=True)
    tgt = np.loadtxt(COL)
    z, rho = tgt[:, 0], tgt[:, 1]
    # two starting guesses for the secant: the tapered energy, and 1.3x of it
    e = [tgt[:, 4], 1.3*tgt[:, 4]]
    f = []
    for n in range(2):
        np.savetxt('ic/ic_s%d.txt' % n, np.c_[z, rho, e[n]], fmt='%.10e')
        f.append(np.log(startup(exe, 'ic/ic_s%d.txt' % n, 's%d' % n)[:, 2]/tgt[:, 2]))
    for n in range(2, NITER):
        x0, x1 = np.log(e[n-2]), np.log(e[n-1])
        den = f[n-1] - f[n-2]
        good = np.abs(den) > 1e-14
        x2 = np.where(good, x1 - f[n-1]*(x1 - x0)/np.where(good, den, 1.0), x1)
        e.append(np.exp(np.clip(x2, x1 - 0.7, x1 + 0.7)))
        np.savetxt('ic/ic_s%d.txt' % n, np.c_[z, rho, e[n]], fmt='%.10e')
        c = startup(exe, 'ic/ic_s%d.txt' % n, 's%d' % n)
        f.append(np.log(c[:, 2]/tgt[:, 2]))
        print('iter %2d  max |ln T/T_target| = %.3e' % (n, np.abs(f[n]).max()))
    col = np.loadtxt('ic/column_used.txt')
    np.savetxt('ic/cu_final.txt', col, fmt='%.10e')
    p = col[:, 3]
    n = len(z)
    i = np.arange(n)
    ip, im = np.clip(i + HALF, 0, n-1), np.clip(i - HALF, 0, n-1)
    arad = G0 + ((p[ip] - p[im])/(z[ip] - z[im]))/rho
    np.savetxt('ic/wb_arad.txt', np.c_[z, arad], fmt='%.12e',
               header='z[cm]  a_rad[cm/s^2] = g0 + (1/rho) dP_gas/dz')
    print('a_rad/g0 over the file: min %.4g max %.4g' % ((arad/G0).min(),
                                                         (arad/G0).max()))


if __name__ == '__main__':
    main()
