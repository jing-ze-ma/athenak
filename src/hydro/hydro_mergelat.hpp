#ifndef HYDRO_HYDRO_MERGELAT_HPP_
#define HYDRO_HYDRO_MERGELAT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_mergelat.hpp
//! \brief mrglat-1019: MERGED LATERAL CELLS near the inner radius (H-AMR style), opt-in.
//!
//! The global dt of the spherical (sp wedge, cubed sphere) stellar setups is set by the
//! LATERAL hydro CFL at r_in (r dtheta, r sin(theta) dphi are 4x smaller than dr there).
//! For r < r_c every 1 x f x f (x1, x2, x3) group of cells is treated as ONE hydro cell:
//!  (1) after every hydro stage update (RKUpdate + source terms + the conduction
//!      operators, before SendU) and after every write of the gas by the implicit
//!      radiation solve (RadiationM1::FluidRefresh), the conserved state of each group is
//!      replaced by its group value: rho, the total energy and the scalars by their
//!      VOLUME average; the momentum as a Cartesian vector field m_i = P/V + Omega x d_i
//!      (d_i = cell centre - group volume centroid) that keeps the group's linear momentum
//!      P = sum V_i m_i AND its angular momentum L = sum V_i x_i x m_i exactly (both in
//!      Cartesian components, x_i = x1v rhat at the cell centre, the arm of the code's L_z
//!      = sum V R m_phi), <hydro>/merge_lat_mom = rigid (default) | cart (Omega = 0); the
//!      energy E_i = e_bar + |m_i|^2/(2 rho_bar), e_bar so that sum V_i E_i is unchanged
//!      (falls back to Omega = 0 where the rotation part would take more than half of the
//!      thermal energy).  The vector is moved between bases exactly: sp orthonormal
//!      (r, theta, phi); cs radial + COVARIANT on the unit panel tangents (raised with
//!      cos_cell).  So every stage ends with group-uniform rho, e, P and a group-rigid
//!      velocity for r < r_c; the fluxes between the cells of a group cancel in the group
//!      sums, i.e. the group value is updated by its OUTER face fluxes only.
//!  (2) the x2 / x3 face states of the merged cells come from a GROUP-level limited
//!      linear reconstruction (group values of the two neighbouring groups, van Leer
//!      slope on the coarse spacing, evaluated at the fine faces), so the core is second
//!      order on the coarse lateral grid; x1 is untouched (fine, incl. wb_x1).
//!  (3) Hydro::NewTimeStep uses the GROUP widths f dx2, f dx3 in the merged region.
//!  (4) <hydro>/merge_lat_rad = true also replaces the radiation E (volume average) and F
//!      (Cartesian volume average, same basis transforms) of each group after the implicit
//!      solve; default false (radiation fine).
//! FOFC runs unchanged on the merged states (fine trial update, fine flags; a flagged cell
//! makes its faces first order; the group average after the update keeps it conservative).
//! Groups never straddle a MeshBlock: each block's nx2, nx3 must be divisible by every f
//! (startup FATAL), groups start at js / ks.
//!
//! Keys (read only when named; without merge_lat_rc nothing is constructed and every
//! kernel is the default one -- bitwise):
//!   <hydro>/merge_lat_rc   = r_1[,r_2,...]   shell outer radii (increasing, code units)
//!   <hydro>/merge_lat_f    = f_1[,f_2,...]   lateral merge factor per shell (x1v < r_1:
//!                                            f_1, r_1 <= x1v < r_2: f_2, ...)
//!   <hydro>/merge_lat_mom  = rigid | cart
//!   <hydro>/merge_lat_rad  = false | true
//!   <hydro>/merge_lat_check = true  (test) every average re-sums its groups: max relative
//!                                change of M, E, |P|, |L| printed as MLATDIAG chk
//!   <hydro>/merge_lat_diag = N   every N cycles rank 0 prints MLATDIAG: the global gas
//!                                M, Etot(IEN), P (Cartesian), L (Cartesian) sums (may be
//!                                named without merge_lat_rc: diagnostic of a fine run)

#include <string>
#include "athena.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"

class MeshBlockPack;
class ParameterInput;

namespace hydro {

constexpr int kMlatMaxShell = 8;

//! the geometry the merge kernels capture (all Views / values: no host pointer)
struct MergeLatGeo {
  DvceArray4D<Real> vol;
  DvceArray2D<Real> x1v, x2v, x3v;
  DvceArray3D<Real> ccell;
  DualArray1D<RegionSize> size;
  DualArray1D<int> panel;
  int js, ks, nx2, nx3;
  bool cs;
};

//! cell-centre direction rhat, the two tangential unit vectors (sp e_theta, e_phi; cs
//! e_xi, e_eta) and their cosine cc (sp 0)
KOKKOS_INLINE_FUNCTION
void MlatBasis(const MergeLatGeo &g, const int m, const int k, const int j,
               Real *rh, Real *b1, Real *b2, Real &cc) {
  if (g.cs) {
    const Real x2l = g.size.d_view(m).x2min, x2u = g.size.d_view(m).x2max;
    const Real x3l = g.size.d_view(m).x3min, x3u = g.size.d_view(m).x3max;
    const Real xc = 0.25*M_PI*CellCenterX(j-g.js, g.nx2, x2l, x2u);
    const Real ec = 0.25*M_PI*CellCenterX(k-g.ks, g.nx3, x3l, x3u);
    const int pnl = g.panel.d_view(m);
    Real qc[3];
    cubed_sphere::PanelToCart(pnl, xc, ec, qc);
    cubed_sphere::PanelTangents(pnl, xc, ec, b1, b2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    for (int q = 0; q < 3; ++q) {rh[q] = qc[q]*qn;}
    cc = g.ccell(m,k,j);
  } else {
    const Real t = g.x2v(m,j), p = g.x3v(m,k);
    const Real st = sin(t), ct = cos(t), sp = sin(p), cp = cos(p);
    rh[0] = st*cp; rh[1] = st*sp; rh[2] = ct;
    b1[0] = ct*cp; b1[1] = ct*sp; b1[2] = -st;
    b2[0] = -sp;   b2[1] = cp;    b2[2] = 0.0;
    cc = 0.0;
  }
}

//! stored components (radial, tangential pair) -> Cartesian vector
KOKKOS_INLINE_FUNCTION
void MlatToCart(const Real *rh, const Real *b1, const Real *b2, const Real cc,
                const Real a1, const Real a2, const Real a3, Real *v) {
  const Real s2 = 1.0 - cc*cc;
  const Real u2 = (a2 - cc*a3)/s2, u3 = (a3 - cc*a2)/s2;
  for (int q = 0; q < 3; ++q) {v[q] = a1*rh[q] + u2*b1[q] + u3*b2[q];}
}

//! Cartesian vector -> stored components (radial, covariant pair; sp orthonormal)
KOKKOS_INLINE_FUNCTION
void MlatFromCart(const Real *rh, const Real *b1, const Real *b2, const Real *v,
                  Real &a1, Real &a2, Real &a3) {
  a1 = v[0]*rh[0] + v[1]*rh[1] + v[2]*rh[2];
  a2 = v[0]*b1[0] + v[1]*b1[1] + v[2]*b1[2];
  a3 = v[0]*b2[0] + v[1]*b2[1] + v[2]*b2[2];
}

//! group-level limited linear reconstruction at the faces of fine cell s (0..f-1) of a
//! group of f cells: qm, qc, qp the group values of the lower neighbour group, this group
//! (the cell's own value) and the upper neighbour group (coarse spacing f).  Index units:
//! the group centre sits at (f-1)/2, the cell's faces at s -+ 1/2.  qlo = state at the
//! cell's lower face (the RIGHT state of that face), qhi = at its upper face.
KOKKOS_INLINE_FUNCTION
void MlatCoarsePLM(const Real qm, const Real qc, const Real qp, const int s, const int f,
                   Real &qlo, Real &qhi) {
  const Real dql = qc - qm, dqr = qp - qc;
  const Real dq2 = dql*dqr;
  // van Leer slope per COARSE cell (= 2 dqm of the fine PLM)
  const Real sl = (dq2 > 0.0) ? 2.0*dq2/(dql + dqr) : 0.0;
  const Real rf = 1.0/static_cast<Real>(f);
  const Real c0 = 0.5*static_cast<Real>(f - 1);
  qlo = qc + sl*(static_cast<Real>(s) - 0.5 - c0)*rf;
  qhi = qc + sl*(static_cast<Real>(s) + 0.5 - c0)*rf;
}

class MergeLat {
 public:
  MergeLat(MeshBlockPack *ppack, ParameterInput *pin);
  static bool Named(ParameterInput *pin);

  bool on = false;          // merging active (merge_lat_rc named)
  int nsh = 0;
  Real rc[kMlatMaxShell];
  int fs[kMlatMaxShell];
  int fmax = 1;
  bool rigid = true;
  bool rad = false;
  int diag_every = 0;
  bool check = false;       // <hydro>/merge_lat_check: group-sum errors of every average
  double chkerr[4] = {0.0, 0.0, 0.0, 0.0};   // max rel. errors M, E, P, L (MLATDIAG chk)
  // per (m, i) merge factor (1 = not merged), incl. the x1 ghosts
  DvceArray2D<int> fi;
  int imax_loc = -1;        // largest merged interior i on this rank (-1 none)
  long ncall = 0, nfallback = 0;

  MergeLatGeo Geo();
  void Build();
  //! kind 0: hydro u (nvar = nhydro + nscalars, IDN IM1..3 IEN, scalars);
  //! kind 1: radiation (E, F1..3)
  //! phi/grav: <hydro>/etotgrav (IEN carries rho Phi; Phi per cell)
  void Average(DvceArray5D<Real> &u, const int nhyd, const int nvar, const int kind,
               const DvceArray4D<Real> &phi, const bool grav);
  //! number of groups whose density (kind 0) / E (kind 1) is not bitwise uniform
  int CountNonUniform(const DvceArray5D<Real> &u, const int n0);
  void Diag(const DvceArray5D<Real> &u, const char *tag);

 private:
  MeshBlockPack *pmy_pack;
  bool built = false;
};

}  // namespace hydro
#endif  // HYDRO_HYDRO_MERGELAT_HPP_
