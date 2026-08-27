#ifndef COORDINATES_CUBED_SPHERE_HPP_
#define COORDINATES_CUBED_SPHERE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cubed_sphere.hpp
//! \brief Global geometry of the gnomonic equiangular cubed sphere: where each panel
//! sits on the cube, and the tangent basis it carries.
//!
//! Coordinates::CoordGnomonicEquiangle only ever works INSIDE one panel, and every
//! panel's metric is the same function of (xi,eta) -- nothing there knows how a panel is
//! oriented in space. The orientation lives in exactly two places: Mesh::panel_neighbors
//! (which panel you reach across which face, and how the two tangential indices line up)
//! and the frames below. They must agree, and they do: propagating panel_neighbors around
//! the cube starting from panel 0's frame reproduces all six frames uniquely, and the
//! resulting 24 directed edges reproduce Mesh::GetPanelBoundary entry for entry.
//!
//! On panel p a point is the unit vector (a*x + b*y + n)/delta, with x = tan(xi),
//! y = tan(eta), delta = sqrt(1 + x^2 + y^2), and (a,b,n) the right-handed frame below;
//! xi = x2 and eta = x3 both run over [-pi/4, pi/4], and n is the outward face normal.
//! The map is used OUTSIDE that range too: a panel's ghost zone is its own gnomonic map
//! extended past the seam, which is what makes the halo geometry well defined.
//!
//! WHY THIS IS SHARED. Two panel frames written out twice drift apart: an independent
//! copy in a problem generator had panels 3 and 4 interchanged, which is self-consistent
//! within each panel and therefore invisible to any single-panel test, but silently makes
//! a global initial condition a different field on two of the six panels. There must be
//! one copy, and this is it.

#include <math.h>

#include "athena.hpp"

namespace cubed_sphere {

//----------------------------------------------------------------------------------------
//! \brief The right-handed frame (a,b,n) of panel p. Panels 4,0,1,2 walk in longitude
//! (normals +x, +y, -x, -y) and 3/5 are the two poles (+z, -z), matching the face diagram
//! in mesh.hpp.

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
//! \brief Panel coordinates (xi,eta) on panel p to a Cartesian unit vector. Valid for
//! |xi|,|eta| beyond pi/4, i.e. in the ghost zones.

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
//! \brief The inverse: the (xi,eta) of the direction q in panel p's chart. Well defined
//! whenever q.n > 0, which includes a neighbouring panel's near ghost zones.

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
//! \brief Which panel the direction q belongs to: the one whose outward normal it is
//! most aligned with.

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
//! \brief The two UNIT tangent vectors of panel p at (xi,eta), in Cartesian components.
//! These are d(rhat)/d(xi) and d(rhat)/d(eta) normalised, which is the basis the code
//! stores velocity components on -- see the note above GnomonicEquianglePrimFaceX1.
//! They are NOT orthogonal: e1.e2 = -sin(xi) sin(eta), which is Coordinates::cos_cell.

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
//! \fn cubed_sphere::TransformMomentum
//! \brief Re-express the two tangential COVARIANT momentum components of a cell on panel
//! `psrc` at (xi,eta) in the tangent basis of panel `pdst` at the SAME physical point.
//!
//! Across a panel seam the two charts do NOT share a tangent basis. Their eta axes agree
//! along the shared edge, but the xi axes differ by a shear that grows away from the edge
//! midpoint (for the 0/1 seam, e_xi = e_xi' - sqrt(2) sin(eta) e_eta'), reaching O(1) at
//! the panel corners. Copying the components across with a signed axis permutation -- the
//! permutation is exact for the INDICES, and exact for the components only on the seam
//! midline -- therefore leaves an O(1) error in the halo velocity, which makes the scheme
//! inconsistent at every seam. The radial component is untouched: rhat is common to both
//! charts.
//!
//! m_i = V.e_i, so the round trip is: raise with the source Gram matrix to get the
//! contravariant pair, rebuild V in Cartesian, then dot into the destination basis.

KOKKOS_INLINE_FUNCTION
void TransformMomentum(const int psrc, const int pdst, const Real xi, const Real eta,
                       const Real m2, const Real m3, Real &m2_out, Real &m3_out) {
  Real s1[3], s2[3];
  PanelTangents(psrc, xi, eta, s1, s2);

  Real q[3];
  PanelToCart(psrc, xi, eta, q);
  Real xid, etad;
  CartToPanel(pdst, q, xid, etad);
  Real d1[3], d2[3];
  PanelTangents(pdst, xid, etad, d1, d2);

  // raise on the source panel: g = [[1,c],[c,1]]
  const Real cs = s1[0]*s2[0] + s1[1]*s2[1] + s1[2]*s2[2];
  const Real det = 1.0 - cs*cs;
  const Real p2 = (m2 - cs*m3)/det;
  const Real p3 = (m3 - cs*m2)/det;

  // V = p2*s1 + p3*s2, then lower on the destination panel
  m2_out = 0.0;
  m3_out = 0.0;
  for (int c=0; c<3; ++c) {
    const Real vc = p2*s1[c] + p3*s2[c];
    m2_out += vc*d1[c];
    m3_out += vc*d2[c];
  }
  return;
}

} // namespace cubed_sphere
#endif // COORDINATES_CUBED_SPHERE_HPP_
