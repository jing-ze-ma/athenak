#ifndef RECONSTRUCT_PLM_C1_HPP_
#define RECONSTRUCT_PLM_C1_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file plm_c1.hpp
//! \brief lhllc2nd-1102: PLM with a C1 slope limiter, <hydro>/plm_c1 = true (and
//! <mhd>/plm_c1; EOS_Data::plm_c1, read only when named, default false = the van Leer
//! PLM of plm.hpp, bit for bit).
//!
//! WHY: the van Leer (harmonic-mean) slope dl dr/(dl + dr), set to 0 for dl dr <= 0, has
//! a KINK at every local extremum of every reconstructed variable (slope 2 dl on one side
//! of dl = 0, slope 0 on the other).  An extremum that does not move (stationary vortex,
//! a slowly evolving low-Mach flow, the deep layers of a hot Jupiter) keeps the state ON
//! that kink, step after step, and the RK2 step loses second order in dt there (see
//! hydro/rsolvers/lhllc_c1_hyd.hpp for the argument).  The low-dissipation low-Mach
//! solvers (LHLLC/LHLLD) show it most, but it is a property of the reconstruction.
//!
//! WHAT: the van Leer slope times theta (2 - theta), theta = 2 dl dr/(dl^2 + dr^2) (> 0;
//! else 0).  For smooth data theta = 1 - O(dx^2), so the slope keeps second order; at
//! an extremum theta -> 0 linearly and the van Leer slope -> 0 linearly, so the product
//! -> 0 QUADRATICALLY and the limiter is C1 (except at dl = dr = 0, where it is still
//! Lipschitz).  It never exceeds the van Leer slope (theta (2 - theta) <= 1), so it is as
//! TVD as van Leer.  The non-uniform (Mignone 2014) slope of Hydro/MHD::PLM_nonuniform
//! gets the same factor.

#include <math.h>
#include "athena.hpp"

KOKKOS_INLINE_FUNCTION
void PLMC1(const Real &q_im1, const Real &q_i, const Real &q_ip1,
           Real &ql_ip1, Real &qr_i) {
  const Real dql = (q_i - q_im1);
  const Real dqr = (q_ip1 - q_i);
  const Real dq2 = dql*dqr;
  Real dqm = 0.0;
  if (dq2 > 0.0) {
    const Real th = 2.0*dq2/(dql*dql + dqr*dqr);
    dqm = dq2/(dql + dqr)*th*(2.0 - th);
  }
  ql_ip1 = q_i + dqm;
  qr_i   = q_i - dqm;
  return;
}

KOKKOS_INLINE_FUNCTION
void PiecewiseLinearC1X1(TeamMember_t const &member, const int m, const int k,
     const int j, const int il, const int iu, const DvceArray5D<Real> &q,
     ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
  int nvar = q.extent_int(1);
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      PLMC1(q(m,n,k,j,i-1), q(m,n,k,j,i), q(m,n,k,j,i+1), ql(n,i+1), qr(n,i));
    });
  }
  return;
}

KOKKOS_INLINE_FUNCTION
void PiecewiseLinearC1X2(TeamMember_t const &member, const int m, const int k,
     const int j, const int il, const int iu, const DvceArray5D<Real> &q,
     ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
  int nvar = q.extent_int(1);
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      PLMC1(q(m,n,k,j-1,i), q(m,n,k,j,i), q(m,n,k,j+1,i), ql_jp1(n,i), qr_j(n,i));
    });
  }
  return;
}

KOKKOS_INLINE_FUNCTION
void PiecewiseLinearC1X3(TeamMember_t const &member, const int m, const int k,
     const int j, const int il, const int iu, const DvceArray5D<Real> &q,
     ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
  int nvar = q.extent_int(1);
  for (int n=0; n<nvar; ++n) {
    par_for_inner(member, il, iu, [&](const int i) {
      PLMC1(q(m,n,k-1,j,i), q(m,n,k,j,i), q(m,n,k+1,j,i), ql_kp1(n,i), qr_k(n,i));
    });
  }
  return;
}
#endif  // RECONSTRUCT_PLM_C1_HPP_
