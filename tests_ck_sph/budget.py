"""Energy budget and dilution check for a correlated-k column dump.

Columns of the dump (see two_stream_rt.hpp, rt_dumpcol):
  0 i   1 r_face  2 p[bar]  3 T[K]  4 F_lw_net  5 Q_sw  6 Gamma_1  7 grad_ad
  8 tau_R(face)  9 w_diff(face)  10 Src_lw  11 A_cell  12 V_cell

The face areas are rebuilt from the cell geometry the solver itself reported: for any
mesh whose cells are radial shells of a fixed solid angle, V_i = dOm (r_{i+1}^3-r_i^3)/3
and A(f) = dOm r_f^2, so dOm comes out of V_i and the two r_face values.  A_cell is
V_i/dx_i, which is what the scheme uses as the constant frame area inside a cell.

Reported:
  budget_rel  = [sum_i V_i Src_i + (A_top F_top - A_bot F_bot)] / max|A F|
                -- zero to round-off iff the deposit telescopes (the spherical form).
  ppar_rel    = the same with every area set to 1 (what the old form telescopes to).
  AF spread   = max/min of A(f) F(f) over the optically thin faces, the
                transparent-shell statement L = 4 pi r^2 F = const.
"""
import sys
import numpy as np


def load(path):
    with open(path) as f:
        head = [ln for ln in f if ln.startswith("#")]
    icut = int([w for ln in head for w in [ln] if "icut =" in w][0]
               .split("icut =")[1].split("(")[0])
    d = np.loadtxt(path)
    return icut, d


def report(path, label):
    icut, d = load(path)
    i = d[:, 0].astype(int)
    rf, F, Src, Vc = d[:, 1], d[:, 4], d[:, 10], d[:, 12]
    ie = i[-1] - 1                     # the last row is the top FACE, ie+1
    lo = int(np.where(i == icut)[0][0])
    hi = int(np.where(i == ie + 1)[0][0])
    # solid angle from the first active cell, then the face areas
    dom = 3.0*Vc[lo]/(rf[lo+1]**3 - rf[lo]**3)
    A = dom*rf**2
    dep = np.sum(Vc[lo:hi]*Src[lo:hi])
    flux = A[hi]*F[hi] - A[lo]*F[lo]
    scale = max(abs(A[hi]*F[hi]), abs(A[lo]*F[lo]))
    dep_pp = np.sum(Vc[lo:hi]*Src[lo:hi])       # same deposit
    flux_pp = (F[hi] - F[lo])*np.mean(A[lo:hi+1])
    print(f"--- {label}  ({path})")
    print(f"    icut={icut} ie={ie}  r_in={rf[lo]:.4e} r_out={rf[hi]:.4e} "
          f"r_out/r_in={rf[hi]/rf[lo]:.4f}  A_out/A_in={A[hi]/A[lo]:.4f}")
    print(f"    sum V Src      = {dep:+.8e}")
    print(f"    -(A_t F_t - A_b F_b) = {-flux:+.8e}")
    print(f"    BUDGET_REL (spherical) = {(dep + flux)/scale:+.3e}")
    print(f"    budget_rel (plane-parallel bookkeeping) = "
          f"{(dep_pp + flux_pp)/max(abs(flux_pp), 1e-300):+.3e}")
    # transparent statement: A F over the thin faces
    thin = d[:, 8] if np.any(d[:, 8] > 0) else None
    sel = slice(hi - 12, hi + 1)
    AF = A[sel]*F[sel]
    print(f"    A F over the top 12 faces: min {AF.min():.6e} max {AF.max():.6e} "
          f"max/min {AF.max()/AF.min():.5f}")
    print(f"    F   over the top 12 faces: min {F[sel].min():.6e} "
          f"max {F[sel].max():.6e} max/min {F[sel].max()/F[sel].min():.5f}")
    return dict(A=A, F=F, Src=Src, V=Vc, rf=rf, lo=lo, hi=hi,
                p=d[:, 2], T=d[:, 3], Q=d[:, 5])


if __name__ == "__main__":
    for a in sys.argv[1:]:
        report(a, a)
