"""Per-cell / column budgets for a correlated-k column dump, thermal and beam separately.

Dump columns (two_stream_rt.hpp, rt_dumpcol):
  0 i  1 r_face  2 p[bar]  3 T[K]  4 F_lw_net  5 Q_sw  6 Gamma_1  7 grad_ad
  8 tau_R(face)  9 w_diff(face)  10 Src_lw  11 A_cell  12 V_cell

Face areas come from the solver's own geometry: for radial shells of fixed solid angle,
V_i = dOm (r_{i+1}^3 - r_i^3)/3 and A(f) = dOm r_f^2, so dOm falls out of V_i and the two
r_face values.  dz_i = V_i / A_cell_i.

THERMAL:  Phi(f) = A(f) F_lw_net(f); the identity is V_i Src_i = Phi(f_lo) - Phi(f_hi).
BEAM:     Q_sw is a volumetric absorption rate.  The beam's flux PER UNIT AREA at a face
          is rebuilt downward from the top, F_b(f) = F_b(top) - sum_{i>=f} Q_i dz_i, which
          is exactly how the kernel forms it (one running transmission per chain), and
          F_b(top) = (1-albedo) F_star mu0 T_ghost.  The absorbed power is sum_i V_i Q_i.
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
                    hdr[k] = float(ln.split(k + " =")[1].split(",")[0].split("(")[0]
                                    .replace("K", "").strip())
    d = np.loadtxt(path)
    return hdr, d


def report(path):
    hdr, d = load(path)
    i = d[:, 0].astype(int)
    rf, F, Q, S, Ac, V = d[:, 1], d[:, 4], d[:, 5], d[:, 10], d[:, 11], d[:, 12]
    icut = int(hdr["icut"])
    lo = int(np.where(i == icut)[0][0])
    hi = len(i) - 1
    dom = 3.0*V[lo]/(rf[lo+1]**3 - rf[lo]**3)
    A = dom*rf**2
    dz = V/np.where(Ac > 0, Ac, 1.0)
    Phi = A*F
    sc = np.max(np.abs(Phi))
    res = V[lo:hi]*S[lo:hi] - (Phi[lo:hi] - Phi[lo+1:hi+1])
    print(f"--- {path}   mu0 = {hdr['mu0']:.4f}, icut = {icut}, "
          f"r {rf[lo]:.4e}..{rf[hi]:.4e}, A_out/A_in = {A[hi]/A[lo]:.3f}")
    print(f"  THERMAL  max per-cell |V dep - dPhi|/max|Phi| = "
          f"{np.max(np.abs(res))/sc:.3e}")
    print(f"  THERMAL  column budget  [sum V Src + (A_t F_t - A_b F_b)]/max|A F| = "
          f"{(np.sum(V[lo:hi]*S[lo:hi]) + (Phi[hi]-Phi[lo]))/sc:.3e}")
    # ---- beam
    Fstar = SIG*(hdr["T_irr"]**4)
    mu0 = hdr["mu0"]
    if mu0 <= 0.0:
        print("  BEAM     night column, no beam")
        return
    absorbed = np.sum(V[lo:hi]*Q[lo:hi])
    # beam flux per unit area, rebuilt downward from the top face
    Fb = np.zeros(hi+1)
    for f in range(hi, lo-1, -1):
        Fb[f] = Fb[f+1] + Q[f]*dz[f] if f < hi else 0.0
    # Fb currently holds the CUMULATIVE absorption below the top; the top value is
    # whatever entered, which is absorbed_total + what survives to icut.  Use the
    # deposit-free statement instead: F_b(f) = F_b(top) - cumulative(from top to f)
    cum = np.zeros(hi+1)
    for f in range(hi-1, lo-1, -1):
        cum[f] = cum[f+1] + Q[f]*dz[f]
    Ftop_beam = cum[lo] + 0.0            # everything absorbed, if the column is opaque
    print(f"  BEAM     absorbed power sum V Q      = {absorbed:.5e}")
    print(f"           / (A_top mu0 F*)            = "
          f"{absorbed/(A[hi]*mu0*Fstar):.5f}")
    # radius where the slant optical depth reaches 1, from the rebuilt beam flux
    prof = Ftop_beam - cum                # beam flux per unit area at each face
    tgt = Ftop_beam/np.e
    f1 = hi
    for f in range(hi, lo-1, -1):
        if prof[f] <= tgt:
            f1 = f
            break
    print(f"           tau_slant = 1 at r = {rf[f1]:.5e} "
          f"(i = {i[f1]}, p = {d[f1,2]:.3e} bar)")
    print(f"           / (A(tau=1) mu0 F*)         = "
          f"{absorbed/(A[f1]*mu0*Fstar):.5f}")
    print(f"           beam flux entering the top face / ((1-A_b) F* mu0): "
          f"{prof[hi]/(Fstar*mu0):.5f}  [albedo not removed]")


if __name__ == "__main__":
    for a in sys.argv[1:]:
        report(a)
