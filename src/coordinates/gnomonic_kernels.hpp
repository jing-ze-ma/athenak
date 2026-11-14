#ifndef COORDINATES_GNOMONIC_KERNELS_HPP_
#define COORDINATES_GNOMONIC_KERNELS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_kernels.hpp
//! \brief Device functions that rotate reconstructed states, magnetic fields, EMFs and
//! fluxes between the non-orthogonal gnomonic equiangle basis and the orthonormal frame
//! of each sweep of the cubed-sphere Riemann solver. They read only the trig arrays held
//! in a GnomonicTrig, which kernels capture by value.

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \struct GnomonicTrig
//! \brief One (sin, cos) pair of the angle between the xi and eta axes, indexed (m,k,j),
//! for the cell centres (x1 sweep), the xi faces (x2 sweep) or the eta faces (x3 sweep).
//! Built on the host by Coordinates::GnomonicTrig{Cell,FaceXi,FaceEta}().

struct GnomonicTrig {
  DvceArray3D<Real> sn, cs;
};

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquianglePrimFaceX1
//! \brief Rotates the tangential velocities of the x1-face states into the orthonormal
//! frame {rhat, e_xi, (e_eta - c e_xi)/s}.

KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql,
     ScrArray2D<Real> &qr) {
    const Real sin_theta = gt.sn(m,k,j);
    const Real cos_theta = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql(IVY,i+1) += cos_theta * ql(IVZ,i+1);
      ql(IVZ,i+1) *= sin_theta;
      qr(IVY,i)   += cos_theta * qr(IVZ,i);
      qr(IVZ,i)   *= sin_theta;
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquianglePrimFaceX2
//! \brief Rotates the velocities of the xi-face states (normal IVY, tangential IVZ) into
//! the orthonormal frame of the x2 sweep.

KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql_jp1,
     ScrArray2D<Real> &qr_j) {
    const Real sin_jp1 = gt.sn(m,k,j+1);
    const Real cos_jp1 = gt.cs(m,k,j+1);
    const Real sin_j = gt.sn(m,k,j);
    const Real cos_j = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql_jp1(IVZ,i) += cos_jp1 * ql_jp1(IVY,i);
      ql_jp1(IVY,i) *= sin_jp1;
      qr_j(IVZ,i)   += cos_j * qr_j(IVY,i);
      qr_j(IVY,i)   *= sin_j;
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquianglePrimFaceX3
//! \brief Rotates the velocities of the eta-face states (normal IVZ, tangential IVY) into
//! the orthonormal frame of the x3 sweep.

KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX3(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql_kp1,
     ScrArray2D<Real> &qr_k) {
    const Real sin_kp1 = gt.sn(m,k+1,j);
    const Real cos_kp1 = gt.cs(m,k+1,j);
    const Real sin_k = gt.sn(m,k,j);
    const Real cos_k = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql_kp1(IVY,i) += cos_kp1 * ql_kp1(IVZ,i);
      ql_kp1(IVZ,i) *= sin_kp1;
      qr_k(IVY,i)   += cos_k * qr_k(IVZ,i);
      qr_k(IVZ,i)   *= sin_k;
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquiangleFaceBX2
//! \brief Rotates the face-parallel magnetic field of the xi-face states into the x2
//! sweep frame {nhat_xi, e_eta, rhat}. bcc holds the orthonormal triple
//! (B.rhat, B.e_xi, B.f2) with f2 = (e_eta - c e_xi)/s, so only slot 2 changes, to
//! B.e_eta = c*(B.e_xi) + s*(B.f2). The x1 and x3 sweep frames need no rotation.

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFaceBX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &bl_jp1,
     ScrArray2D<Real> &br_j) {
    const Real sin_jp1 = gt.sn(m,k,j+1);
    const Real cos_jp1 = gt.cs(m,k,j+1);
    const Real sin_j = gt.sn(m,k,j);
    const Real cos_j = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      bl_jp1(2,i) = cos_jp1*bl_jp1(1,i) + sin_jp1*bl_jp1(2,i);
      br_j(2,i)   = cos_j*br_j(1,i)     + sin_j*br_j(2,i);
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquiangleEmfX1
//! \brief Projects the x1-sweep EMF onto the eta edge, E.e_eta = c*e21 + s*e31, as needed
//! by CT. The x2 and x3 sweep EMFs are already along edge directions.

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleEmfX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     DvceArray4D<Real> e31,
     DvceArray4D<Real> e21) {
    const Real sin_theta = gt.sn(m,k,j);
    const Real cos_theta = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      e31(m,k,j,i) = cos_theta*e21(m,k,j,i) + sin_theta*e31(m,k,j,i);
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquiangleFluxX1
//! \brief Converts the tangential momentum fluxes of the x1 sweep from the orthonormal
//! frame back to covariant components.

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sn(m,k,j);
    const Real cos_theta = gt.cs(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      Real fb = flx(m,IM3,k,j,i)/sin_theta;
      Real fa = flx(m,IM2,k,j,i) - fb*cos_theta;
      flx(m,IM2,k,j,i) = fa + fb*cos_theta;
      flx(m,IM3,k,j,i) = fb + fa*cos_theta;
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquiangleFluxX2
//! \brief Converts the tangential momentum fluxes of the x2 sweep from the orthonormal
//! frame back to covariant components.

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sn(m,k,j);
    const Real cos_theta = gt.cs(m,k,j);
    const Real T22 = 1.0/sin_theta;
    const Real T32 = -cos_theta/sin_theta;
    par_for_inner(member, il, iu, [&](const int i) {
      Real fb = flx(m,IM3,k,j,i) + T32 * flx(m,IM2,k,j,i);
      Real fa = T22 * flx(m,IM2,k,j,i);
      flx(m,IM2,k,j,i) = fa + fb*cos_theta;
      flx(m,IM3,k,j,i) = fb + fa*cos_theta;
    });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void GnomonicEquiangleFluxX3
//! \brief Converts the tangential momentum fluxes of the x3 sweep from the orthonormal
//! frame back to covariant components.

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX3(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sn(m,k,j);
    const Real cos_theta = gt.cs(m,k,j);
    const Real T23 = -cos_theta/sin_theta;
    const Real T33 = 1.0/sin_theta;
    par_for_inner(member, il, iu, [&](const int i) {
      Real fa = flx(m,IM2,k,j,i) + T23 * flx(m,IM3,k,j,i);
      Real fb = T33 * flx(m,IM3,k,j,i);
      flx(m,IM2,k,j,i) = fa + fb*cos_theta;
      flx(m,IM3,k,j,i) = fb + fa*cos_theta;
    });
  return;
}
#endif // COORDINATES_GNOMONIC_KERNELS_HPP_
