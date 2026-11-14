#ifndef COORDINATES_CUBED_SPHERE_HPP_
#define COORDINATES_CUBED_SPHERE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cubed_sphere.hpp
//! \brief Global geometry of the gnomonic equiangular cubed sphere: the orientation of
//! each panel on the cube and its tangent basis, plus the transformations of vectors
//! between the charts of two neighbouring panels used at panel seams.
//! On panel p a point is the unit vector (a*x + b*y + n)/delta, with x = tan(xi),
//! y = tan(eta), delta = sqrt(1 + x^2 + y^2), (a,b,n) the right-handed panel frame, and
//! xi = x2, eta = x3 in [-pi/4, pi/4]. The map is also used beyond that range, in the
//! ghost zones of a panel. The frames are consistent with Mesh::panel_neighbors.

#include <math.h>

#include "athena.hpp"

namespace cubed_sphere {
//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::PanelFrame
//! \brief Right-handed frame (a,b,n) of panel p. Panels 0,1,2,4 lie on the equator,
//! panels 3 and 5 are the +z and -z poles (see the face diagram in mesh.hpp).

KOKKOS_INLINE_FUNCTION
void PanelFrame(const int p, Real a[3], Real b[3], Real n[3]) {
  switch (p) {
    case 0:
      a[0]= 0.0; a[1]= 1.0; a[2]= 0.0;
      b[0]= 0.0; b[1]= 0.0; b[2]= 1.0;
      n[0]= 1.0; n[1]= 0.0; n[2]= 0.0;
      break;
    case 1:
      a[0]=-1.0; a[1]= 0.0; a[2]= 0.0;
      b[0]= 0.0; b[1]= 0.0; b[2]= 1.0;
      n[0]= 0.0; n[1]= 1.0; n[2]= 0.0;
      break;
    case 2:
      a[0]= 0.0; a[1]=-1.0; a[2]= 0.0;
      b[0]= 0.0; b[1]= 0.0; b[2]= 1.0;
      n[0]=-1.0; n[1]= 0.0; n[2]= 0.0;
      break;
    case 3:
      a[0]= 0.0; a[1]= 1.0; a[2]= 0.0;
      b[0]=-1.0; b[1]= 0.0; b[2]= 0.0;
      n[0]= 0.0; n[1]= 0.0; n[2]= 1.0;
      break;
    case 4:
      a[0]= 1.0; a[1]= 0.0; a[2]= 0.0;
      b[0]= 0.0; b[1]= 0.0; b[2]= 1.0;
      n[0]= 0.0; n[1]=-1.0; n[2]= 0.0;
      break;
    default:
      a[0]= 0.0; a[1]= 1.0; a[2]= 0.0;
      b[0]= 1.0; b[1]= 0.0; b[2]= 0.0;
      n[0]= 0.0; n[1]= 0.0; n[2]=-1.0;
      break;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::PanelToCart
//! \brief Converts panel coordinates (xi,eta) on panel p to a Cartesian unit vector.

KOKKOS_INLINE_FUNCTION
void PanelToCart(const int p, const Real xi, const Real eta, Real q[3]) {
  Real a[3], b[3], n[3];
  PanelFrame(p, a, b, n);
  const Real x = tan(xi);
  const Real y = tan(eta);
  const Real delta = sqrt(1.0 + x*x + y*y);
  for (int c=0; c<3; ++c) {
    q[c] = (a[c]*x + b[c]*y + n[c])/delta;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::CartToPanel
//! \brief Converts a Cartesian direction q to (xi,eta) in the chart of panel p (valid
//! for q.n > 0).

KOKKOS_INLINE_FUNCTION
void CartToPanel(const int p, const Real q[3], Real &xi, Real &eta) {
  Real a[3], b[3], n[3];
  PanelFrame(p, a, b, n);
  Real qa = 0.0, qb = 0.0, qn = 0.0;
  for (int c=0; c<3; ++c) {
    qa += q[c]*a[c];
    qb += q[c]*b[c];
    qn += q[c]*n[c];
  }
  xi  = atan(qa/qn);
  eta = atan(qb/qn);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int cubed_sphere::FindPanel
//! \brief Returns the panel whose outward normal is most aligned with direction q.

KOKKOS_INLINE_FUNCTION
int FindPanel(const Real q[3]) {
  int best = 0;
  Real bdot = -2.0;
  for (int p=0; p<6; ++p) {
    Real a[3], b[3], n[3];
    PanelFrame(p, a, b, n);
    const Real d = q[0]*n[0] + q[1]*n[1] + q[2]*n[2];
    if (d > bdot) { bdot = d; best = p; }
  }
  return best;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::PanelTangents
//! \brief Unit tangent vectors e1 = d(rhat)/d(xi), e2 = d(rhat)/d(eta) (normalized) of
//! panel p at (xi,eta) in Cartesian components. This is the basis of the stored velocity
//! components; it is not orthogonal (e1.e2 = -sin(xi) sin(eta)).

KOKKOS_INLINE_FUNCTION
void PanelTangents(const int p, const Real xi, const Real eta,
                   Real e1[3], Real e2[3]) {
  Real a[3], b[3], n[3];
  PanelFrame(p, a, b, n);
  const Real x = tan(xi);
  const Real y = tan(eta);
  const Real delta = sqrt(1.0 + x*x + y*y);
  const Real inv1 = 1.0/sqrt(1.0 + y*y);
  const Real inv2 = 1.0/sqrt(1.0 + x*x);
  for (int c=0; c<3; ++c) {
    const Real rhat = (a[c]*x + b[c]*y + n[c])/delta;
    e1[c] = (a[c]*delta - rhat*x)*inv1;
    e2[c] = (b[c]*delta - rhat*y)*inv2;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \struct cubed_sphere::SeamXform
//! \brief Transforms the two tangential covariant momentum components of a cell on panel
//! psrc at (xi,eta) into the tangent basis of panel pdst at the same physical point.
//! The components are raised with the source metric, rebuilt as a Cartesian vector and
//! projected onto the destination basis. The radial component is common to both charts.
//! The geometry (SeamXformAt) depends only on (psrc,pdst,xi,eta) and is computed once per
//! (k,j) by the seam packers; ApplyMomentumXform does the per-cell arithmetic.

struct SeamXform {
  Real s1[3], s2[3], d1[3], d2[3], cs, det;
};

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::SeamXformAt
//! \brief Geometry of TransformMomentum for one panel pair and (xi,eta).

KOKKOS_INLINE_FUNCTION
void SeamXformAt(const int psrc, const int pdst, const Real xi, const Real eta,
                 SeamXform &t) {
  PanelTangents(psrc, xi, eta, t.s1, t.s2);

  Real q[3];
  PanelToCart(psrc, xi, eta, q);
  Real xid, etad;
  CartToPanel(pdst, q, xid, etad);
  PanelTangents(pdst, xid, etad, t.d1, t.d2);

  // raise on the source panel: g = [[1,c],[c,1]]
  t.cs = t.s1[0]*t.s2[0] + t.s1[1]*t.s2[1] + t.s1[2]*t.s2[2];
  t.det = 1.0 - t.cs*t.cs;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::ApplyMomentumXform
//! \brief Per-cell part of TransformMomentum: raise, rebuild V, project on the
//! destination basis.

KOKKOS_INLINE_FUNCTION
void ApplyMomentumXform(const SeamXform &t, const Real m2, const Real m3,
                        Real &m2_out, Real &m3_out) {
  const Real p2 = (m2 - t.cs*m3)/t.det;
  const Real p3 = (m3 - t.cs*m2)/t.det;

  // V = p2*s1 + p3*s2, then lower on the destination panel
  m2_out = 0.0;
  m3_out = 0.0;
  for (int c=0; c<3; ++c) {
    const Real vc = p2*t.s1[c] + p3*t.s2[c];
    m2_out += vc*t.d1[c];
    m3_out += vc*t.d2[c];
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::TransformMomentum
//! \brief Full momentum transformation: SeamXformAt followed by ApplyMomentumXform.

KOKKOS_INLINE_FUNCTION
void TransformMomentum(const int psrc, const int pdst, const Real xi, const Real eta,
                       const Real m2, const Real m3, Real &m2_out, Real &m3_out) {
  SeamXform t;
  SeamXformAt(psrc, pdst, xi, eta, t);
  ApplyMomentumXform(t, m2, m3, m2_out, m3_out);
  return;
}

//----------------------------------------------------------------------------------------
//! \struct cubed_sphere::SeamFieldXform
//! \brief Transforms the tangential magnetic field given as its projections on the face
//! normals of panel psrc at (xi,eta), (B.nhat_xi, B.nhat_eta) as stored in b0.x2f/x3f,
//! into its projections on the face normals of panel pdst at the same point. The route is
//! normals -> contravariant components -> Cartesian -> destination normals. The radial
//! component is common to both charts and is not handled here. As for TransformMomentum,
//! the geometry (SeamFieldXformAt) is split from the per-cell work (ApplyFieldXform).

struct SeamFieldXform {
  Real s1[3], s2[3], d1[3], d2[3], sns, csd, snd;
};

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::SeamFieldXformAt
//! \brief Geometry of TransformFieldToDstNormals for one panel pair and (xi,eta).

KOKKOS_INLINE_FUNCTION
void SeamFieldXformAt(const int psrc, const int pdst, const Real xi, const Real eta,
                      SeamFieldXform &t) {
  PanelTangents(psrc, xi, eta, t.s1, t.s2);
  const Real css = t.s1[0]*t.s2[0] + t.s1[1]*t.s2[1] + t.s1[2]*t.s2[2];
  t.sns = sqrt(1.0 - css*css);

  Real q[3];
  PanelToCart(psrc, xi, eta, q);
  Real xid, etad;
  CartToPanel(pdst, q, xid, etad);
  PanelTangents(pdst, xid, etad, t.d1, t.d2);
  t.csd = t.d1[0]*t.d2[0] + t.d1[1]*t.d2[1] + t.d1[2]*t.d2[2];
  t.snd = sqrt(1.0 - t.csd*t.csd);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::ApplyFieldXform
//! \brief Per-cell part of TransformFieldToDstNormals.

KOKKOS_INLINE_FUNCTION
void ApplyFieldXform(const SeamFieldXform &t, const Real bn_xi, const Real bn_eta,
                     Real &bn_xi_out, Real &bn_eta_out) {
  // B.nhat_xi = sn * B^xi, so the contravariant pair on the unit basis is just this
  const Real p2 = bn_xi/t.sns;
  const Real p3 = bn_eta/t.sns;

  // V = p2*s1 + p3*s2, then project on the destination normals
  // nhat_xi' = (d1 - csd d2)/snd,  nhat_eta' = (d2 - csd d1)/snd
  bn_xi_out = 0.0;
  bn_eta_out = 0.0;
  for (int c=0; c<3; ++c) {
    const Real vc = p2*t.s1[c] + p3*t.s2[c];
    bn_xi_out  += vc*(t.d1[c] - t.csd*t.d2[c]);
    bn_eta_out += vc*(t.d2[c] - t.csd*t.d1[c]);
  }
  bn_xi_out /= t.snd;
  bn_eta_out /= t.snd;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void cubed_sphere::TransformFieldToDstNormals
//! \brief Full field transformation: SeamFieldXformAt followed by ApplyFieldXform.

KOKKOS_INLINE_FUNCTION
void TransformFieldToDstNormals(const int psrc, const int pdst, const Real xi,
                                const Real eta, const Real bn_xi, const Real bn_eta,
                                Real &bn_xi_out, Real &bn_eta_out) {
  SeamFieldXform t;
  SeamFieldXformAt(psrc, pdst, xi, eta, t);
  ApplyFieldXform(t, bn_xi, bn_eta, bn_xi_out, bn_eta_out);
  return;
}

} // namespace cubed_sphere
#endif // COORDINATES_CUBED_SPHERE_HPP_
