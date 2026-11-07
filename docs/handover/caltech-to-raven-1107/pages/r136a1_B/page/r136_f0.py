"""r136_f0.py (r136kf-1106): the solver-consistent comoving radial flux F0 for the R136a1 page.

WHY.  On a backward-Euler-path step (time2_auto_be / time2_reject redo, about 1/3 of the cycles since 3.5 d) the stored
cell flux is F = <F0 faces> + (E + P_rr) v*, v* = the velocity BEFORE that solve's radiative kick (rad_m1_implicit.cpp
write-back: vfx = t2st && impl_vimp && t2_fvnew), while the dumped v already holds the kick
dv = dt kappa (F0 - F_wb) / c.  F - (E + P_rr) v then is F0 - G (F0 - F_wb), G = (4/3) kappa E dt / c: G = 50 at the wall.
Stage steps re-form F at the kicked velocity, so a dump after a stage step has F - (E + P_rr) v = <F0 faces> to round-off.

WHAT.  Where G can matter (deep, thick) take F0 from the face equation the solve converges to in the thick limit:
F0_f = -(c / 3 chi_f) dE/dr with chi_f the harmonic face mean (implicit_face_opac = harmonic), F0 = mean of the two faces;
wall face = the deck's fixed flux.  Elsewhere keep F - (E + P_rr) v, where G is negligible.
Weight: w = 1 for G(dt_ref) >= 3e-3, 0 for <= 1e-3 (linear in log G), times a guard on the cell optical depth (0 at
dtau <= 3, 1 at >= 10).  dt_ref = 8 s bounds the run's dt.
"""
import numpy as np

C = 2.99792458e10
FWALL = 1.131943656797e15    # <rad_m1>/implicit_flux_x1min of decks/r136a1_128_B.athinput
DT_REF = 8.0


def f0_solver(f0_lab, E, kr, rho, x1f, r, fwall=FWALL, dt_ref=DT_REF):
    """f0_lab = F - (E + P_rr) v_r [k,j,i]; returns (F0, w), same shape and dtype as f0_lab."""
    chi = kr * rho
    dr = np.diff(x1f).astype(chi.dtype)
    G = (4. / 3.) * kr * E * dt_ref / C
    w = np.clip(np.log10(np.maximum(G, 1e-30) / 1e-3) / np.log10(3.), 0., 1.) * np.clip((chi * dr - 3.) / 7., 0., 1.)
    del G
    iw = np.where(w.max(axis=(0, 1)) > 0.)[0]
    F0 = f0_lab.copy()
    if len(iw) == 0:
        return F0, w
    n = int(iw[-1]) + 1                     # cells 0..n-1 can take the thick form
    nf = min(n + 1, E.shape[2])
    ch = 2. / (1. / chi[..., :nf - 1] + 1. / chi[..., 1:nf])
    dx = np.diff(r[:nf]).astype(chi.dtype)
    ff = -(C / (3. * ch)) * (E[..., 1:nf] - E[..., :nf - 1]) / dx      # interior faces 1 .. nf-1
    del ch
    fl = np.empty_like(E[..., :n]); fr = np.empty_like(E[..., :n])
    fl[..., 0] = fwall; fl[..., 1:] = ff[..., :n - 1]
    fr[..., :nf - 1] = ff[..., :nf - 1]
    if nf - 1 < n:
        fr[..., n - 1] = f0_lab[..., n - 1]  # the outermost cell of the mesh (w = 0 there in practice)
    F0[..., :n] = w[..., :n] * 0.5 * (fl + fr) + (1. - w[..., :n]) * f0_lab[..., :n]
    return F0, w
