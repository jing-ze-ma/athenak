#ifndef COORDINATES_GNOMONIC_KERNELS_HPP_
#define COORDINATES_GNOMONIC_KERNELS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gnomonic_kernels.hpp
//! \brief the cubed-sphere rotations a flux kernel applies, as FREE functions.
//!
//! These were members of Coordinates, reached from inside the flux kernels as
//! `pmy_pack->pcoord->GnomonicEquiangle...`. That made each of those lambdas capture
//! `this` and dereference two host pointers on the device -- hipcc says so
//! (-Wgpu-maybe-wrong-side: "capture host side class data member by this pointer in
//! device or host device lambda function may result in invalid memory access if this
//! pointer is not accessible on device side"), and it is an illegal access on a GPU not
//! share memory with the host. It survived only because the machine this was developed on
//! has unified-memory APUs.
//!
//! Between them the eight helpers read just SIX Views, so they need no Coordinates
//! object at all -- only `GnomonicTrig`, a small POD of those Views a kernel captures
//! BY VALUE the way it captures any other View. Build one on the host with
//! `Coordinates::GnomonicTrigData()` before the parallel dispatch and pass it in.

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \struct GnomonicTrig
//! \brief the gnomonic-equiangle trig a device kernel needs, capturable by value.
//!
//! Indexed (m,k,j): these depend on the two PANEL-TANGENTIAL angles only -- xi = x2 and
//! eta = x3 -- so they have no radial extent. The _xi arrays live on x2 faces (staggered
//! in j), the _eta arrays on x3 faces (staggered in k).

struct GnomonicTrig {
  DvceArray3D<Real> sin_cell, cos_cell;
  DvceArray3D<Real> sin_face_xi, cos_face_xi;
  DvceArray3D<Real> sin_face_eta, cos_face_eta;
};

KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql,
     ScrArray2D<Real> &qr) {
    const Real sin_theta = gt.sin_cell(m,k,j);
    const Real cos_theta = gt.cos_cell(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql(IVY,i+1) += cos_theta * ql(IVZ,i+1);
      ql(IVZ,i+1) *= sin_theta;
      qr(IVY,i)   += cos_theta * qr(IVZ,i);
      qr(IVZ,i)   *= sin_theta;
    });
  return;
}

// --- x2 sweep: XI face, normal component is IVY, tangential is IVZ
KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql_jp1,
     ScrArray2D<Real> &qr_j) {
    const Real sin_jp1 = gt.sin_face_xi(m,k,j+1);
    const Real cos_jp1 = gt.cos_face_xi(m,k,j+1);
    const Real sin_j = gt.sin_face_xi(m,k,j);
    const Real cos_j = gt.cos_face_xi(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql_jp1(IVZ,i) += cos_jp1 * ql_jp1(IVY,i);
      ql_jp1(IVY,i) *= sin_jp1;
      qr_j(IVZ,i)   += cos_j * qr_j(IVY,i);
      qr_j(IVY,i)   *= sin_j;
    });
  return;
}

// --- x3 sweep: ETA face, normal component is IVZ, tangential is IVY
KOKKOS_INLINE_FUNCTION
void GnomonicEquianglePrimFaceX3(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &ql_kp1,
     ScrArray2D<Real> &qr_k) {
    const Real sin_kp1 = gt.sin_face_eta(m,k+1,j);
    const Real cos_kp1 = gt.cos_face_eta(m,k+1,j);
    const Real sin_k = gt.sin_face_eta(m,k,j);
    const Real cos_k = gt.cos_face_eta(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      ql_kp1(IVY,i) += cos_kp1 * ql_kp1(IVZ,i);
      ql_kp1(IVZ,i) *= sin_kp1;
      qr_k(IVY,i)   += cos_k * qr_k(IVZ,i);
      qr_k(IVZ,i)   *= sin_k;
    });
  return;
}

// --- MAGNETIC FIELD into the sweep's orthonormal frame ----------------------------
//
// bcc holds the ORTHONORMAL triple (B.rhat, B.e_xi, B.f2), f2 = (e_eta - c e_xi)/s
// -- see Coordinates::GnomonicEquiangleRaiseVelMHD for why the field is stored that
// way and not as the raw face-normal average. Slot ordering matches wl/wr.
//
// Each sweep's frame is {nhat, then its two face-parallel axes}, and in every sweep
// the NORMAL component is taken from b0.x*f rather than from bcc, so only the two
// face-parallel slots matter:
//   x1, frame {rhat, e_xi, f2}      -> slots 1,2 are already exactly right;
//   x2, frame {nhat_xi, e_eta, rhat} -> slot 0 is right, slot 2 must become B.e_eta;
//   x3, frame {nhat_eta, rhat, e_xi} -> slots 0,1 are already exactly right.
// So x1 and x3 need nothing at all, and only this one routine exists.
//
// B.e_eta = c*(B.e_xi) + s*(B.f2), the same projection the x1 EMF needs below.
KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFaceBX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     ScrArray2D<Real> &bl_jp1,
     ScrArray2D<Real> &br_j) {
    const Real sin_jp1 = gt.sin_face_xi(m,k,j+1);
    const Real cos_jp1 = gt.cos_face_xi(m,k,j+1);
    const Real sin_j = gt.sin_face_xi(m,k,j);
    const Real cos_j = gt.cos_face_xi(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      bl_jp1(2,i) = cos_jp1*bl_jp1(1,i) + sin_jp1*bl_jp1(2,i);
      br_j(2,i)   = cos_j*br_j(1,i)     + sin_j*br_j(2,i);
    });
  return;
}

// --- EMF returned by the x1 sweep -------------------------------------------------
// CT integrates E around a face, so e2/e3 must be E projected on the EDGE directions
// e_xi and e_eta, whose lengths dxedge.x2e/x3e it multiplies them by. Of the three
// sweeps only x1 returns something else: its orthonormal frame is
// {rhat, e_xi, (e_eta - c e_xi)/s}, whose third axis is NOT the eta edge, so
//   E.e_eta = c*(E.e_xi) + s*(E.f2)  =  c*e21 + s*e31.
// The x2 and x3 frames are {nhat, e_eta, rhat} and {nhat, rhat, e_xi}: both of their
// face-parallel axes ARE edge directions, so those EMFs need no rotation.
KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleEmfX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
     DvceArray4D<Real> e31,
     DvceArray4D<Real> e21) {
    const Real sin_theta = gt.sin_cell(m,k,j);
    const Real cos_theta = gt.cos_cell(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      e31(m,k,j,i) = cos_theta*e21(m,k,j,i) + sin_theta*e31(m,k,j,i);
    });
  return;
}

// --- x1 flux: radial. IM2 comes back unchanged, which is correct: the orthonormal
// frame used above is {e_xi, n-hat} with n-hat perpendicular to e_xi, so the first
// component already IS the covariant m_2.
KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX1(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sin_cell(m,k,j);
    const Real cos_theta = gt.cos_cell(m,k,j);
    par_for_inner(member, il, iu, [&](const int i) {
      Real fb = flx(m,IM3,k,j,i)/sin_theta;
      Real fa = flx(m,IM2,k,j,i) - fb*cos_theta;
      flx(m,IM2,k,j,i) = fa + fb*cos_theta;
      flx(m,IM3,k,j,i) = fb + fa*cos_theta;
    });
  return;
}

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX2(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sin_face_xi(m,k,j);
    const Real cos_theta = gt.cos_face_xi(m,k,j);
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

KOKKOS_INLINE_FUNCTION
void GnomonicEquiangleFluxX3(const GnomonicTrig &gt, TeamMember_t const &member,
     const int m, const int k, const int j, const int il, const int iu,
                             DvceArray5D<Real> flx) {
    const Real sin_theta = gt.sin_face_eta(m,k,j);
    const Real cos_theta = gt.cos_face_eta(m,k,j);
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
