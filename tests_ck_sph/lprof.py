"""Gates (3) and (4): L(r) = A(r) F_net(r) normalised by the injected internal luminosity.

On this path the internal flux enters as sigma T_int^4 PER UNIT AREA through the face
i = is, i.e. AT x1min (conduction.cpp:1533-1537, `if (i == is && fin != 0.0)
flx1(m,IEN,k,j,i) += fin`), and the RK update multiplies it by area1(m,k,j,is).  So the
injected luminosity is L_int = A(x1min) sigma T_int^4 and the steady-state statement is

    A(r) F_net(r) = L_int                      (night column, no beam)
    A_top F_top,thermal = L_int + P_beam       (day column)

The two-stream's face flux is defined only for i >= icut; below the cut the flux is
carried by the radiative-conduction operator, which does not expose a per-face flux to
this dump.  So L(r) is reported over the faces the two-stream owns, and the blend weight
w_diff (column 9) is printed so the handover band is visible.
"""
import sys
import numpy as np

SIG = 5.6704e-5


def load(path):
    hdr = {}
    with open(path) as f:
        for ln in f:
            if not ln.startswith("#"):
                break
            for k in ["mu0", "icut", "T_int", "T_irr"]:
                if k + " =" in ln:
                    hdr[k] = float(ln.split(k + " =")[1].split(",")[0]
                                   .split("(")[0].replace("K", "").strip())
    return hdr, np.loadtxt(path)


def prof(path, label):
    hdr, d = load(path)
    i = d[:, 0].astype(int)
    rf, p, F, Q, Ac, V, w = (d[:, 1], d[:, 2], d[:, 4], d[:, 5],
                             d[:, 11], d[:, 12], d[:, 9])
    icut = int(hdr["icut"])
    lo = int(np.where(i == icut)[0][0])
    hi = len(i) - 1
    dom = 3.0*V[lo]/(rf[lo+1]**3 - rf[lo]**3)
    A = dom*rf**2
    Lint = A[0]*SIG*hdr["T_int"]**4          # A(x1min) sigma T_int^4
    L = A*F/Lint
    mu0 = hdr["mu0"]
    dz = V/np.where(Ac > 0, Ac, 1.0)
    Pb = np.sum(V[lo:hi]*Q[lo:hi])/Lint
    print(f"--- {label}  mu0 = {mu0:+.4f}  icut = {icut}  "
          f"L_int = A(x1min) sigma T_int^4 = {Lint:.5e}")
    print(f"    beam power absorbed / L_int = {Pb:.4f}")
    pure = np.where((w[lo:hi+1] == 0.0))[0] + lo
    if pure.size:
        print(f"    faces the two-stream owns alone (w = 0): "
              f"{pure.size}, r {rf[pure[0]]:.4e}..{rf[pure[-1]]:.4e}, "
              f"p {p[pure[0]]:.3e}..{p[pure[-1]]:.3e} bar")
        print(f"    L(r)/L_int there: min {L[pure].min():.4f} "
              f"max {L[pure].max():.4f}  max/min {L[pure].max()/L[pure].min():.4f}")
    band = np.where((w[lo:hi+1] > 0.0) & (w[lo:hi+1] < 1.0))[0] + lo
    if band.size:
        print(f"    handover band 0<w<1: {band.size} faces, "
              f"r {rf[band[0]]:.4e}..{rf[band[-1]]:.4e}")
    print(f"    L(top)/L_int = {L[hi]:.4f}   "
          f"(day: target 1 + P_beam = {1.0 + max(Pb, 0.0):.4f})")
    print("      i      r[cm]     p[bar]      w     F_net       L/L_int")
    for k in range(lo, hi+1, max(1, (hi-lo)//14)):
        print(f"    {i[k]:4d} {rf[k]:.4e} {p[k]:.3e} {w[k]:5.2f} "
              f"{F[k]:+.4e} {L[k]:+.4f}")
    print(f"    {i[hi]:4d} {rf[hi]:.4e} {p[hi]:.3e} {w[hi]:5.2f} "
          f"{F[hi]:+.4e} {L[hi]:+.4f}")
    return L, rf, hi, Pb


if __name__ == "__main__":
    for a in sys.argv[1:]:
        prof(a, a)
