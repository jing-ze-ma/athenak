#ifndef DIFFUSION_CS_EDGE_JB_HPP_
#define DIFFUSION_CS_EDGE_JB_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_edge_jb.hpp
//! \brief seamexact-1103: the current J and field B at a cubed-sphere cell EDGE, the
//! inputs of the resistive edge EMF, in one place.  Used by the ambipolar edge EMF
//! (ambipolar.cpp, ad_cs_emf) and by the single-valued seam current of
//! <mhd>/cs_seam_res = mean_j (MHD::CsSeamResist), so both evaluate J and B at an edge
//! with the same expressions in the same order.
//!
//! Components (see AddEMFAmbipolarGnomonic): contravariant (jr, ju1, ju2), (br, bu1, bu2)
//! on {rhat, e_xi, e_eta}; covariant bd1 = B.e_xi, bd2 = B.e_eta; c = e_xi.e_eta and
//! sn = sqrt(1 - c^2) at the edge; eta_a the edge's ambipolar coefficient.

#include "athena.hpp"

//----------------------------------------------------------------------------------------
// B.e_xi at the centre of an x2 face / B.e_eta at the centre of an x3 face, from the
// face-normal components (the partner component from the four faces that share the
// face's corners).  THE one copy: the resistive curl (resistivity_gnomonic.cpp, pass 1)
// and the cubed-sphere seam current (mhd_cs_seam.cpp, cs_seam_res = own | stokes) both
// use it (stokes-1103: moved here from resistivity_gnomonic.cpp, the same expressions).

KOKKOS_INLINE_FUNCTION
Real CsBcovXi(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
              const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
              const int i) {
  const Real b3a = 0.25*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i)
                       + b.x3f(m,k,j-1,i) + b.x3f(m,k+1,j-1,i));
  return (b.x2f(m,k,j,i) + cosf(m,k,j)*b3a)/sinf(m,k,j);
}

KOKKOS_INLINE_FUNCTION
Real CsBcovEta(const DvceFaceFld4D<Real> &b, const DvceArray3D<Real> &cosf,
               const DvceArray3D<Real> &sinf, const int m, const int k, const int j,
               const int i) {
  const Real b2a = 0.25*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i)
                       + b.x2f(m,k-1,j,i) + b.x2f(m,k-1,j+1,i));
  return (b.x3f(m,k,j,i) + cosf(m,k,j)*b2a)/sinf(m,k,j);
}

struct CsEdgeJB {
  Real c, sn, eta_a;
  Real jr, ju1, ju2;
  Real br, bu1, bu2, bd1, bd2;
};

struct CsEdgeOps {
  DvceFaceFld4D<Real> b;
  DvceArray5D<Real> bc;
  DvceEdgeFld4D<Real> jn;
  DvceArray3D<Real> cxi, sxi, cet, set, ccl, scl;
  DvceArray4D<Real> etc;   // the cell eta_A (ad_eta_c); unused when has_ad = false
  bool has_ad = true;
  CsEdgeOps(const DvceFaceFld4D<Real> &b_, const DvceArray5D<Real> &bc_,
            const DvceEdgeFld4D<Real> &jn_) : b(b_), bc(bc_), jn(jn_) {}

  KOKKOS_INLINE_FUNCTION
  Real Eta(const int m, const int k, const int j, const int i) const {
    return has_ad ? etc(m,k,j,i) : 0.0;
  }
  // x1 (radial) edge at (k,j,i): r centre, xi face j, eta face k
  KOKKOS_INLINE_FUNCTION
  void X1(const int m, const int k, const int j, const int i, CsEdgeJB &s) const {
    s.c = 0.5*(cxi(m,k,j) + cxi(m,k-1,j));
    s.sn = sqrt(1.0 - s.c*s.c);
    s.eta_a = 0.25*(Eta(m,k,j,i) + Eta(m,k,j-1,i) + Eta(m,k-1,j,i) + Eta(m,k-1,j-1,i));
    s.br = 0.25*(bc(m,IBX,k,j,i) + bc(m,IBX,k,j-1,i)
               + bc(m,IBX,k-1,j,i) + bc(m,IBX,k-1,j-1,i));
    s.bu1 = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k-1,j,i))/s.sn;   // B^xi
    s.bu2 = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k,j-1,i))/s.sn;   // B^eta
    s.bd1 = s.bu1 + s.c*s.bu2; s.bd2 = s.c*s.bu1 + s.bu2;
    s.jr = jn.x1e(m,k,j,i);
    s.ju1 = 0.25*(jn.x2e(m,k,j,i) + jn.x2e(m,k,j,i+1)
                + jn.x2e(m,k,j-1,i) + jn.x2e(m,k,j-1,i+1))/s.sn;
    s.ju2 = 0.25*(jn.x3e(m,k,j,i) + jn.x3e(m,k,j,i+1)
                + jn.x3e(m,k-1,j,i) + jn.x3e(m,k-1,j,i+1))/s.sn;
  }
  // x2 edge at (k,j,i): r face i, xi centre j, eta face k; ia/ib the radial-edge columns
  KOKKOS_INLINE_FUNCTION
  void X2(const int m, const int k, const int j, const int i, const int ia, const int ib,
          CsEdgeJB &s) const {
    s.c = cet(m,k,j); s.sn = set(m,k,j);
    s.eta_a = 0.25*(Eta(m,k,j,i) + Eta(m,k,j,i-1) + Eta(m,k-1,j,i) + Eta(m,k-1,j,i-1));
    s.br = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k-1,j,i));
    s.bu2 = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k,j,i-1))/s.sn;   // B^eta
    s.bd1 = 0.25*(bc(m,IBY,k,j,i) + bc(m,IBY,k,j,i-1)       // B_xi
                + bc(m,IBY,k-1,j,i) + bc(m,IBY,k-1,j,i-1));
    s.bu1 = s.bd1 - s.c*s.bu2; s.bd2 = s.c*s.bu1 + s.bu2;
    s.jr = 0.25*(jn.x1e(m,k,j,ib) + jn.x1e(m,k,j,ia)
               + jn.x1e(m,k,j+1,ib) + jn.x1e(m,k,j+1,ia));
    s.ju1 = jn.x2e(m,k,j,i)/s.sn;
    s.ju2 = 0.25*(jn.x3e(m,k,j,i) + jn.x3e(m,k,j+1,i)
                + jn.x3e(m,k-1,j,i) + jn.x3e(m,k-1,j+1,i))/s.sn;
  }
  // x3 edge at (k,j,i): r face i, xi face j, eta centre k
  KOKKOS_INLINE_FUNCTION
  void X3(const int m, const int k, const int j, const int i, const int ia, const int ib,
          CsEdgeJB &s) const {
    s.c = cxi(m,k,j); s.sn = sxi(m,k,j);
    s.eta_a = 0.25*(Eta(m,k,j,i) + Eta(m,k,j-1,i) + Eta(m,k,j,i-1) + Eta(m,k,j-1,i-1));
    s.br = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j-1,i));
    s.bu1 = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j,i-1))/s.sn;   // B^xi
    // B.e_eta per cell = c_cell B_xi + s_cell N_eta, then averaged
    auto bde = [&](int jj, int ii) {
      return ccl(m,k,jj)*bc(m,IBY,k,jj,ii) + scl(m,k,jj)*bc(m,IBZ,k,jj,ii);
    };
    s.bd2 = 0.25*(bde(j,i) + bde(j,i-1) + bde(j-1,i) + bde(j-1,i-1));
    s.bu2 = s.bd2 - s.c*s.bu1; s.bd1 = s.bu1 + s.c*s.bu2;
    s.jr = 0.25*(jn.x1e(m,k,j,ib) + jn.x1e(m,k,j,ia)
               + jn.x1e(m,k+1,j,ib) + jn.x1e(m,k+1,j,ia));
    s.ju2 = jn.x3e(m,k,j,i)/s.sn;
    s.ju1 = 0.25*(jn.x2e(m,k,j,i) + jn.x2e(m,k+1,j,i)
                + jn.x2e(m,k,j-1,i) + jn.x2e(m,k+1,j-1,i))/s.sn;
  }
};

// the ambipolar EMF increment eta_A (J_a - (J.B) B_a/B^2) along the edge (dir 1, 2, 3),
// exactly the expressions of the former ad_cs_emf kernel
KOKKOS_INLINE_FUNCTION
Real CsEdgeAdEmf(const int dir, const CsEdgeJB &s) {
  const Real bsq = s.br*s.br + s.bu1*s.bd1 + s.bu2*s.bd2 + 1.0e-300;
  const Real jdb = s.jr*s.br + s.ju1*s.bd1 + s.ju2*s.bd2;
  if (dir == 1) return s.eta_a*(s.jr - jdb*s.br/bsq);
  if (dir == 2) {
    const Real jd1 = s.ju1 + s.c*s.ju2;
    return s.eta_a*(jd1 - jdb*s.bd1/bsq);
  }
  const Real jd2 = s.c*s.ju1 + s.ju2;
  return s.eta_a*(jd2 - jdb*s.bd2/bsq);
}

// <mhd>/cs_seam_res = mean_j: the cell-centred field the ambipolar seam current uses
// (the stage bcc of the caller; nullptr = MHD::bcc0).  Set by AddEMFAmbipolarGnomonic
// before MHD::CsSeamResist (mhd_cs_seam.cpp).
void CsSeamResJBcc(const DvceArray5D<Real> *bc);

#endif  // DIFFUSION_CS_EDGE_JB_HPP_
