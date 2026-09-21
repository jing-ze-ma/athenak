"""Gate (5): T(p) of a relaxed column, ck_spherical off vs on.

The correlated-k column dump carries the EOS temperature (column 3) and the pressure in
bar (column 2), so this is a real T(p) in kelvin -- the earlier version of this gate used
e/rho from the binary output because it did not have the relaxed column dump.
"""
import sys
import numpy as np


def load(path):
    d = np.loadtxt(path)
    return d[:, 2], d[:, 3], d[:, 1]           # p[bar], T[K], r_face


if __name__ == "__main__":
    tag = sys.argv[1]
    p0, T0, r = load(f"h_{tag}_off/col.txt")
    p1, T1, _ = load(f"h_{tag}_on/col.txt")
    print(f"=== {tag}: T(p), off vs on")
    print("     r[cm]     p_off[bar]  p_on[bar]    T_off[K]   T_on[K]    dT[K]   dT/T")
    n = p0.size
    for i in range(0, n, max(1, n//22)):
        print(f"  {r[i]:.4e} {p0[i]:.4e} {p1[i]:.4e} {T0[i]:9.2f} {T1[i]:9.2f} "
              f"{T1[i]-T0[i]:+9.2f} {T1[i]/T0[i]-1:+7.3f}")
    dT = T1 - T0
    for lab, m in [("p > 10 bar (below the ck cut)", p0 > 10.0),
                   ("1 .. 10 bar", (p0 > 1.0) & (p0 <= 10.0)),
                   ("1e-3 .. 1 bar (photosphere)", (p0 > 1.0e-3) & (p0 <= 1.0)),
                   ("< 1e-3 bar (upper atmosphere)", p0 <= 1.0e-3)]:
        if m.any():
            print(f"  {lab:32s} dT: min {dT[m].min():+8.1f} K  max "
                  f"{dT[m].max():+8.1f} K  mean {dT[m].mean():+8.1f} K")
