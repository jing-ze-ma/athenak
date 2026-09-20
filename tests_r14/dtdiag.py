"""tests_r14: WHERE dt collapses in the tall 1-D He4 columns.  Read-only."""
import struct
import sys

import numpy as np

RS = 2.3717e11
TURN = 4705.0
ARAD = 7.5657e-15
KB = 1.380649e-16
MU = 4.0026/3.0      # fully ionized He4 ~ 4/3
MH = 1.6605390e-24
CFL = 0.15


def rd(fn):
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


def faces(r):
    rf = np.empty(len(r)+1)
    rf[1:-1] = 0.5*(r[1:] + r[:-1])
    rf[0] = r[0] - (rf[1]-r[0])
    rf[-1] = r[-1] + (r[-1]-rf[-2])
    return rf


def weight(rho, lo=1e-11, hi=1e-10):
    x = np.log10(np.maximum(rho, 1e-300))
    s = np.clip((x-np.log10(lo))/(np.log10(hi)-np.log10(lo)), 0, 1)
    return s*s*(3-2*s)


def diag(arm, times):
    R = rd(arm + '/rt_profile.bin')
    r = R[0][1]
    rf = faces(r)
    dx = np.diff(rf)
    vol = 4*np.pi/3*(rf[1:]**3 - rf[:-1]**3)
    for tt in times:
        k = int(np.argmin([abs(x[0]-tt*TURN) for x in R]))
        t, _, q = R[k]
        rho, v1, T, ei = q[0], q[1], q[5], q[6]
        w = weight(rho)
        # gas + tapered radiation pressure, sound speed estimate
        pg = rho*KB*T/(MU*MH)
        prad = w*ARAD*T**4/3.0
        p = pg + prad
        # Gamma1 ~ (5/3 gas, 4/3 rad) mixture; crude: cs^2 = Gam1 p/rho
        beta = np.where(p > 0, pg/np.maximum(p, 1e-300), 1.0)
        gam1 = beta*(5./3.) + (1-beta)*(4./3.)
        cs = np.sqrt(gam1*p/np.maximum(rho, 1e-300))
        sig = CFL*dx/(np.abs(v1)+cs)
        j = int(np.argmin(sig))
        print("\n== %s  t=%.1f s = %.3f turn   dt_min(est)=%.3g s at i=%d r/R=%.4f"
              % (arm, t, t/TURN, sig[j], j, r[j]/RS))
        print("   i  r/R      dx        rho        T         eint       v1"
              "        w    cs        p_rad/p   dt_i")
        sel = sorted(set(list(np.argsort(sig)[:8]) +
                         [i for i in range(len(r)) if 0.80 < r[i]/RS < 1.30][::3]))
        for i in sel:
            print("  %3d %7.4f %9.3e %10.3e %9.3e %10.3e %9.2e %5.2f %9.3e %8.3f %9.3e"
                  % (i, r[i]/RS, dx[i], rho[i], T[i], ei[i], v1[i], w[i], cs[i],
                     (1-beta[i]), sig[i]))
    return R, r, rf, vol


if __name__ == '__main__':
    arm = sys.argv[1]
    diag(arm, [float(x) for x in sys.argv[2:]])
