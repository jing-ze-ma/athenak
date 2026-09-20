"""tests_r17: read the periodic shell-mean MLT face table (problem/mlt_dump_dt).

Columns per row, after the face index i:
  0 r 1 T_f 2 p_f 3 grad 4 grad_ad 5 x 6 F_mlt 7 F_rad 8 F_used 9 F_req 10 F_raddiff
  11 w_blend 12 F_rad_col0 13 F_conv_res 14 grad_rad 15 D 16 F_2s
  17 F_deposited (fconv after the per-column cap)  18 F_cap
Read-only.
"""
import sys

import numpy as np

RS = 2.3717e11
TURN = 4705.0


def rd(fn):
    """-> list of (t, i_array, q[nq, nf])"""
    out = []
    t = None
    rows = []
    idx = []
    for ln in open(fn):
        if ln.startswith('# t '):
            if t is not None:
                out.append((t, np.array(idx), np.array(rows).T))
            t = float(ln.split()[2])
            rows, idx = [], []
        elif ln.startswith('#'):
            continue
        else:
            f = ln.split()
            idx.append(int(f[0]))
            rows.append([float(x) for x in f[1:]])
    if t is not None:
        out.append((t, np.array(idx), np.array(rows).T))
    return out


def zig(a):
    """alternation: |a_i - (a_{i-1}+a_{i+1})/2| / |a_i|, padded"""
    z = np.zeros_like(a)
    z[1:-1] = np.abs(a[1:-1] - 0.5*(a[2:] + a[:-2]))/np.maximum(np.abs(a[1:-1]), 1e-300)
    return z


if __name__ == '__main__':
    arm = sys.argv[1]
    fn = sys.argv[2] if len(sys.argv) > 2 else None
    if fn is None:
        import glob
        fn = glob.glob(arm + '/mltfaces_*.txt.ts')[0]
    R = rd(fn)
    print("%s : %d records, t = %.0f .. %.0f s" % (fn, len(R), R[0][0], R[-1][0]))
    lo, hi = 15, 40
    print(" t/turn  zz(F_used) zz(F_2s) zz(F_req) zz(F_used+F_2s) zz(T) "
          "  max|F_dep-F_used|/F_used   <F_used>/<F_req>")
    for t, ii, q in R:
        s = (ii >= lo) & (ii <= hi)
        fu, f2, fq, td = q[8], q[16], q[9], q[1]
        dep = q[17]
        tot = fu + f2
        cap = np.max(np.abs(dep[s] - fu[s])/np.maximum(np.abs(fu[s]), 1e-300))
        print("%7.3f  %9.2e %8.2e %9.2e %13.2e %9.2e %10.2e %12.4f"
              % (t/TURN, zig(fu)[s].max(), zig(f2)[s].max(), zig(fq)[s].max(),
                 zig(tot)[s].max(), zig(td)[s].max(), cap,
                 np.mean(fu[s])/np.mean(fq[s])))
    t, ii, q = R[-1]
    print("\n  i    r/R     F_req      F_2s       F_used     F_dep      "
          "F_req-F_2s-F_used   T")
    for n, i in enumerate(ii):
        if not (lo-2 <= i <= hi+2):
            continue
        print("%4d %7.4f %10.3e %10.3e %10.3e %10.3e %12.3e %12.5e"
              % (i, q[0][n]/RS, q[9][n], q[16][n], q[8][n], q[17][n],
                 q[9][n]-q[16][n]-q[8][n], q[1][n]))
