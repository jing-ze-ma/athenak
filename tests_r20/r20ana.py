"""tests_r20: mlt_relax_down_time / mlt_sync_passes arm report.  Read-only.

usage: python3 r20ana.py ARM [ARM ...]
"""
import glob
import sys

import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r14')
sys.path.insert(0, '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r17')
from dtdiag import rd                      # noqa: E402
from mltts import rd as rdmlt, zig         # noqa: E402

TURN = 4705.0
RS = 2.3717e11
REF = {0.25: (1.374e30, 1.066e35), 0.5: (3.801e30, 3.437e35),
       1.0: (9.596e30, 20.27e35)}


def band(a, lo, hi, nlo, nhi):
    """fraction of the variance of a[lo:hi+1] in wavelengths nlo..nhi faces"""
    x = np.asarray(a[lo:hi+1], dtype=float)
    x = x - x.mean()
    n = len(x)
    f = np.abs(np.fft.rfft(x))**2
    tot = f[1:].sum()
    if tot <= 0.0:
        return 0.0
    k = np.arange(len(f))
    lam = np.full(len(f), np.inf)
    lam[1:] = n/k[1:]
    s = (lam >= nlo) & (lam <= nhi)
    s[0] = False
    return f[s].sum()/tot


def main(arms):
    print("== hst: |momentum| (col 4) and KEr (col 8); ref = tests_r15/S0.075")
    for arm in arms:
        d = np.loadtxt(arm + '/he4.hydro.hst')
        out = []
        for fr in (0.25, 0.5, 1.0, 1.5):
            k = int(np.argmin(np.abs(d[:, 0] - fr*TURN)))
            if abs(d[k, 0] - fr*TURN) > 0.05*TURN:
                out.append("%.2f -- " % fr)
                continue
            m, ke = abs(d[k, 3]), d[k, 7]
            if fr in REF:
                out.append("%.2f mom %.3e (%+.1f%%) KEr %.3e (%+.1f%%)"
                           % (fr, m, 100*(m/REF[fr][0] - 1.0),
                              ke, 100*(ke/REF[fr][1] - 1.0)))
            else:
                out.append("%.2f mom %.3e KEr %.3e" % (fr, m, ke))
        print("  %-6s %s" % (arm, "\n         ".join(out)))

    print("\n== zigzag(rho) over i = 15..40 from rt_profile.bin")
    for arm in arms:
        try:
            R = rd(arm + '/rt_profile.bin')
        except IOError:
            print("  %-6s (no rt_profile.bin)" % arm)
            continue
        row = []
        for fr in (0.25, 0.5, 1.0, 1.5):
            k = int(np.argmin([abs(x[0] - fr*TURN) for x in R]))
            if abs(R[k][0] - fr*TURN) > 0.05*TURN:
                row.append("%.2f --" % fr)
                continue
            rho = R[k][2][0]
            row.append("%.2f %.3e" % (fr, zig(rho)[15:41].max()))
        print("  %-6s %s" % (arm, " | ".join(row)))

    print("\n== mlt face table at the last record: variance fraction of F_2s and")
    print("   F_used over i = 15..40 in the 2-face and the 3-6-face bands")
    for arm in arms:
        g = glob.glob(arm + '/mltfaces_*.txt.ts')
        if not g:
            print("  %-6s (no .ts)" % arm)
            continue
        R = rdmlt(g[0])
        t, ii, q = R[-1]
        lo = int(np.where(ii == 15)[0][0])
        hi = int(np.where(ii == 40)[0][0])
        f2, fu = q[16], q[8]
        print("  %-6s t/turn %.3f  F_2s: 2f %.3f 3-6f %.3f | "
              "F_used: 2f %.3f 3-6f %.3f | zz(F_2s) %.2e zz(F_used) %.2e"
              % (arm, t/TURN, band(f2, lo, hi, 1.9, 2.1),
                 band(f2, lo, hi, 2.5, 6.5), band(fu, lo, hi, 1.9, 2.1),
                 band(fu, lo, hi, 2.5, 6.5), zig(f2)[lo:hi+1].max(),
                 zig(fu)[lo:hi+1].max()))


if __name__ == '__main__':
    main(sys.argv[1:])
