//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_geometry.cpp
//! \brief Option checks for implicit VET on curvilinear meshes (spherical-polar wedge or
//! full sphere, cubed sphere) and the cubed-sphere seam set-up (CubedS1Init).
//!
//! The geometry lives in the implicit kernels: the E row takes dt A_f/V_i per face and
//! each face equation the centre-to-centre distance dxface (dr, r dtheta, r sin dphi).
//! F is stored in physical orthonormal components (r, theta, phi), as the momentum.
//!
//!  called from: rad_vet.cpp, rad_vet_init.cpp, rad_vet_steady.cpp,
//!    rays/sc_cart_setup.cpp
//!  calls into: rad_vet_bvals.cpp

#include <float.h>
#include <math.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_faces.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanCsGeometry
//! \brief cubed sphere: the cell geometry of the finite-volume divergence in Cartesian y.
//! Lateral face: in the great-circle plane of its edge points P, P' on the cube surface,
//! normal along P x P', area (r_r^2 - r_l^2) th/2 (th the angle P, P'), first moment
//! (r_r^3 - r_l^3)/3 2 sin(th/2) along the bisector.  Radial face: r_f^2 S (S the vector
//! solid angle) minus half the sum of outward th n of the lateral faces, so a cell's six
//! vector areas sum to zero to round-off; centroid (r_r^4 - r_l^4) S/(4 V).
//! Cube points come from the global panel index (exactly +-1 on edges), so both panels at
//! a seam build a face from the same numbers, bitwise.  The frame of y is the hydro's:
//! r_hat and the panel tangents e_xi, e_eta.  Centroids, r_hat, the owner's x2/x3
//! normals, volume and gid (CLG_F, for CleanCsSlot) take the halo exchange.

void RadiationVet::CleanCsGeometry() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, ks = indcs.ks, nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*indcs.ng, n2 = indcs.nx2 + 2*indcs.ng;
  const int n3 = indcs.nx3 + 2*indcs.ng, np = pmy_pack->pmesh->mesh_indcs.nx2;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  for (int d = 0; d < 3; ++d) {
    Kokkos::realloc(cln_ga[d], nmb1 + 1, n3, n2, n1);
    Kokkos::realloc(cln_gd[d], nmb1 + 1, n3, n2, n1);
  }
  // all CLG_NV in gn first (exchanges, distances), then into cln_gn's two arrays
  DvceArray5D<Real> gn("cln_csgn", nmb1 + 1, CLG_NV, n3, n2, n1);
  auto ga0 = cln_ga[0], ga1 = cln_ga[1], ga2 = cln_ga[2];
  auto gd0 = cln_gd[0], gd1 = cln_gd[1], gd2 = cln_gd[2];
  auto size = pmy_pack->pmb->mb_size.d_view;
  auto pnl = pmy_pack->pmb->mb_panel.d_view;
  auto rf = pmy_pack->pcoord->xx1f;
  auto vol = pmy_pack->pcoord->volume;
  auto gid = pmy_pack->pmb->mb_gid.d_view;
  par_for("cln_csgeo", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real fa[3], fb[3], fn[3];
    cubed_sphere::PanelFrame(pnl(m), fa, fb, fn);
    // tangent of the panel angle at face q along x2 (d = 0) or x3 (d = 1)
    auto tq = [&](const int d, const int q) {
      const Real x0 = (d == 0) ? size(m).x2min : size(m).x3min;
      const int s = 2*(q - ((d == 0) ? js : ks) +
                       static_cast<int>(floor(0.5*(x0 + 1.0)*np + 0.5))) - np;
      return (s == np || s == -np) ? static_cast<Real>(s/np)
                                   : tan(0.25*M_PI*static_cast<Real>(s)/np);
    };
    // angle and unit normal (+index) of face q along x2 (d = 0) or x3 (d = 1) across
    // cell o of the other direction: P(eta_r) x P(eta_l), P(xi_l) x P(xi_r)
    auto lat = [&](const int d, const int q, const int o, Real &th, Real nv[3],
                   Real bv[3]) {
      const Real t = tq(d, q), u0 = tq(1 - d, o), u1 = tq(1 - d, o + 1);
      Real p[3], pp[3];
      for (int c = 0; c < 3; ++c) {
        p[c] = (d == 0) ? fa[c]*t + fb[c]*u1 + fn[c] : fa[c]*u0 + fb[c]*t + fn[c];
        pp[c] = (d == 0) ? fa[c]*t + fb[c]*u0 + fn[c] : fa[c]*u1 + fb[c]*t + fn[c];
      }
      nv[0] = p[1]*pp[2] - p[2]*pp[1];
      nv[1] = p[2]*pp[0] - p[0]*pp[2];
      nv[2] = p[0]*pp[1] - p[1]*pp[0];
      const Real s = std::sqrt(nv[0]*nv[0] + nv[1]*nv[1] + nv[2]*nv[2]);
      th = atan2(s, p[0]*pp[0] + p[1]*pp[1] + p[2]*pp[2]);
      const Real ip = 1.0/std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
      const Real ipp = 1.0/std::sqrt(pp[0]*pp[0] + pp[1]*pp[1] + pp[2]*pp[2]);
      for (int c = 0; c < 3; ++c) {
        nv[c] /= s;
        bv[c] = p[c]*ip + pp[c]*ipp;
      }
      const Real ib = 2.0*sin(0.5*th)/std::sqrt(bv[0]*bv[0] + bv[1]*bv[1] + bv[2]*bv[2]);
      for (int c = 0; c < 3; ++c) bv[c] *= ib;
    };
    Real th[4], nv[4][3], bv[4][3], sv[3];
    lat(0, j, k, th[0], nv[0], bv[0]);
    lat(0, j + 1, k, th[1], nv[1], bv[1]);
    lat(1, k, j, th[2], nv[2], bv[2]);
    lat(1, k + 1, j, th[3], nv[3], bv[3]);
    for (int c = 0; c < 3; ++c) {
      sv[c] = -0.5*(th[1]*nv[1][c] - th[0]*nv[0][c] + th[3]*nv[3][c] - th[2]*nv[2][c]);
    }
    const Real ss = std::sqrt(sv[0]*sv[0] + sv[1]*sv[1] + sv[2]*sv[2]);
    const Real rl = rf(m,i), rr = rf(m,i+1), al = 0.5*(rr*rr - rl*rl);
    ga0(m,k,j,i) = rl*rl*ss;
    ga1(m,k,j,i) = al*th[0];
    ga2(m,k,j,i) = al*th[2];
    for (int c = 0; c < 3; ++c) {
      gn(m,CLG_N+c,k,j,i) = sv[c]/ss;
      gn(m,CLG_N+3+c,k,j,i) = nv[0][c];
      gn(m,CLG_N+6+c,k,j,i) = nv[2][c];
      gn(m,CLG_X+c,k,j,i) = 0.25*(rr*rr*rr*rr - rl*rl*rl*rl)*sv[c]/vol(m,k,j,i);
      gn(m,CLG_M+c,k,j,i) = (rr*rr*rr - rl*rl*rl)/3.0*bv[0][c];
      gn(m,CLG_M+3+c,k,j,i) = (rr*rr*rr - rl*rl*rl)/3.0*bv[2][c];
      gn(m,CLG_F+c,k,j,i) = nv[0][c];
      gn(m,CLG_F+4+c,k,j,i) = nv[2][c];
    }
    gn(m,CLG_F+3,k,j,i) = vol(m,k,j,i);
    gn(m,CLG_F+7,k,j,i) = gid(m);
    gn(m,CLG_A,k,j,i) = 0.5*(rl*rl + rr*rr)*ss;
    gn(m,CLG_A+1,k,j,i) = 0.5*al*(th[0] + th[1]);
    gn(m,CLG_A+2,k,j,i) = 0.5*al*(th[2] + th[3]);
    const Real xc = 0.25*M_PI*CellCenterX(j - js, nx2, size(m).x2min, size(m).x2max);
    const Real ec = 0.25*M_PI*CellCenterX(k - ks, nx3, size(m).x3min, size(m).x3max);
    Real qc[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(pnl(m), xc, ec, qc);
    cubed_sphere::PanelTangents(pnl(m), xc, ec, e1, e2);
    const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
    for (int c = 0; c < 3; ++c) {
      gn(m,CLG_R+c,k,j,i) = qc[c]*qn;
      gn(m,CLG_E+c,k,j,i) = e1[c];
      gn(m,CLG_E+3+c,k,j,i) = e2[c];
    }
  });
  CleanGhostFill(gn, CLG_X, 6, false);
  CleanGhostFill(gn, CLG_F, 11, false);
  par_for("cln_csq", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int b = 0;
    for (int d = 1; d < 3; ++d) {
      for (int lo = 0; lo < 2 && (lo == 0 || (d == 1 ? j : k) > 0); ++lo) {
        Real s;
        const int q = CleanCsSlotQ(gn, d, m, k, j, i, k - lo*(d == 2), j - lo*(d == 1),
                                   i, s);
        b += ((q - 1) + 2*(s > 0.0)) << (4*(d - 1) + 2*lo);
      }
    }
    gn(m,CLG_Q,k,j,i) = b;
  });
  par_for("cln_csdist", DevExeSpace(), 0, nmb1, 1, n3-1, 1, n2-1, 1, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int d = 0; d < 3; ++d) {
      const int kl = k - (d == 2), jl = j - (d == 1), il = i - (d == 0);
      Real s = 0.0;
      for (int c = 0; c < 3; ++c) {
        s += SQR(gn(m,CLG_X+c,k,j,i) - gn(m,CLG_X+c,kl,jl,il));
      }
      auto gd = (d == 0) ? gd0 : ((d == 1) ? gd1 : gd2);
      gd(m,k,j,i) = std::sqrt(s);
    }
  });
  // the directions of CleanGeo::a are bitwise the same along each x1 range of a line
  const int is = indcs.is, ie = indcs.ie;
  Kokkos::realloc(cln_gn.r, nmb1 + 1, CLG_NV - CLG_NA, n3, n2, n1);
  Kokkos::realloc(cln_gn.a, nmb1 + 1, CLG_NA, n3, n2, 3);
  cln_gn.is = is;
  cln_gn.ie = ie;
  auto gr = cln_gn.r;
  auto gt = cln_gn.a;
  int nbad = 0;
  Kokkos::parallel_reduce("cln_csgeo2",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, DevExeSpace>(
      {0, 0, 0, 0}, {nmb1 + 1, n3, n2, n1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nb) {
    for (int n = 0; n < CLG_NV; ++n) {
      const int s = CleanGeo::Slot(n);
      const int c = (i < is) ? 0 : ((i > ie) ? 2 : 1), i0 = (c == 0) ? 0 : ((c == 1) ? is
                                                                        : ie + 1);
      const Real v = gn(m,n,k,j,i), v0 = gn(m,n,k,j,i0);
      if (s < 0) {
        gr(m,-1-s,k,j,i) = v;
      } else if (v != v0 || Kokkos::signbit(v) != Kokkos::signbit(v0)) {
        nb += 1;
      } else if (i == i0) {
        gt(m,s,k,j,c) = v;
      }
    }
  }, nbad);
  if (nbad > 0) {
    CleanFatal("<rad_vet> clean cubed sphere: " + std::to_string(nbad) +
               " geometry directions change along x1 (CleanGeo)");
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::MeshBlock0Dx
//! \brief cell width along x(d+1) of the mesh's MeshBlock 0: the same on every rank
//! (a rank's own first MeshBlock differs by round-off)

Real RadiationVet::MeshBlock0Dx(const int d) const {
  auto *pm = pmy_pack->pmesh;
  auto &ms = pm->mesh_size;
  auto &id = pm->mb_indcs;
  const int lv = pm->lloc_eachmb[0].level - pm->root_level;
  const int nr = ((d == 0) ? pm->nmb_rootx1 :
                  ((d == 1) ? pm->nmb_rootx2 : pm->nmb_rootx3)) << lv;
  const Real a = (d == 0) ? ms.x1min : ((d == 1) ? ms.x2min : ms.x3min);
  const Real b = (d == 0) ? ms.x1max : ((d == 1) ? ms.x2max : ms.x3max);
  const int n = (d == 0) ? id.nx1 : ((d == 1) ? id.nx2 : id.nx3);
  return (((nr > 1) ? LeftEdgeX(1, nr, a, b) : b) - a)/static_cast<Real>(n);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanSeamImage
//! \brief vet_gd_theta_seam_sym: the geometry the clean scheme reads, with the ghost rows
//! of a theta-periodic seam (MeshBlocks at x2min / x2max) taking that of the cells they
//! hold, their periodic images.  On a theta wedge centred on the equator the image of the
//! ghost at depth g beyond a wall is the equatorial mirror of the row at depth g inside
//! it: vol, the x1 / x3 face lengths and areas of that row, the x2 faces mirrored, the
//! centroid reflected about the wall (as BuildX2vRecon).  The seam face's two sides then
//! have mirror-paired stencils, as they have across any interior face; without it the
//! one copy of the face computed on the x2min side sees a ghost nearer the pole than the
//! cell it holds (a net F_theta that grows with the wedge)

void RadiationVet::CleanSeamImage() {
  auto *pc = pmy_pack->pcoord;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je, ng = indcs.ng, nmb = pmy_pack->nmb_thispack;
  if (indcs.nx2 < ng) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "vet_gd_theta_seam_sym needs nx2 >= nghost per MeshBlock" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto cp = [](const DvceArray4D<Real> &a, const char *nm) {
    DvceArray4D<Real> b(std::string(nm), a.extent(0), a.extent(1), a.extent(2),
                        a.extent(3));
    Kokkos::deep_copy(b, a);
    return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);
  };
  auto v = cp(pc->volume, "cln_sgv"), a0 = cp(pc->area.x1f, "cln_sga0");
  auto a1 = cp(pc->area.x2f, "cln_sga1"), d1 = cp(pc->dxface.x2f, "cln_sgd1");
  auto d2 = cp(pc->dxface.x3f, "cln_sgd2");
  cln_sgx = DvceArray2D<Real>("cln_sgx", pc->x2v.extent(0), pc->x2v.extent(1));
  cln_sgtg = DvceArray3D<Real>("cln_sgtg", cln_tg.extent(0), cln_tg.extent(1),
                               cln_tg.extent(2));
  Kokkos::deep_copy(cln_sgx, pc->x2v);
  Kokkos::deep_copy(cln_sgtg, cln_tg);
  auto x2v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_sgx);
  auto tg = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cln_sgtg);
  auto x2f = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->xx2f);
  auto x1v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pc->x1v);
  auto bc = pmy_pack->pmb->mb_bcs.h_view;
  const BoundaryFlag pf = BoundaryFlag::periodic;
  for (int m = 0; m < nmb; ++m) {
    for (int sd = 0; sd < 2; ++sd) {
      if (bc(m, sd ? BoundaryFace::outer_x2 : BoundaryFace::inner_x2) != pf) {continue;}
      const int w = sd ? je + 1 : js;   // the seam face
      for (int g = 0; g < ng; ++g) {
        const int jg = sd ? je + 1 + g : js - 1 - g, jm = sd ? je - g : js + g;
        // (cells; the x2 face beyond the ghost, depth g + 1 from the wall)
        const int fg = sd ? je + 2 + g : js - 1 - g, fm = sd ? je - g : js + 1 + g;
        x2v(m,jg) = 2.0*x2f(m,w) - x2v(m,jm);
        tg(m,0,jg) = tg(m,0,jm);   // (tan, sin of the cell and of its outer x2 face)
        tg(m,2,jg) = tg(m,2,jm);
        tg(m,1,fg) = tg(m,1,fm);
        tg(m,3,fg) = tg(m,3,fm);
        for (int k = 0; k < static_cast<int>(v.extent(1)); ++k) {
          for (int i = 0; i < static_cast<int>(a0.extent(3)); ++i) {
            if (i < static_cast<int>(v.extent(3))) {
              v(m,k,jg,i) = v(m,k,jm,i);
              a1(m,k,fg,i) = a1(m,k,fm,i);
              d2(m,k,jg,i) = d2(m,k,jm,i);
            }
            a0(m,k,jg,i) = a0(m,k,jm,i);
          }
        }
        for (int i = 0; i < static_cast<int>(v.extent(3)); ++i) {   // (d2: k = nk too)
          d2(m,d2.extent(1)-1,jg,i) = d2(m,d2.extent(1)-1,jm,i);
        }
      }
      // the centre distances across the seam face and the ghosts' own x2 faces
      for (int g = 0; g <= ng - 1; ++g) {
        const int fg = sd ? w + g : w - g;
        for (int k = 0; k < static_cast<int>(d1.extent(1)); ++k) {
          for (int i = 0; i < static_cast<int>(d1.extent(3)); ++i) {
            d1(m,k,fg,i) = x1v(m,i)*(x2v(m,fg) - x2v(m,fg-1));
          }
        }
      }
    }
  }
  auto dv = [](auto &h, DvceArray4D<Real> &d, const char *nm) {
    d = DvceArray4D<Real>(std::string(nm), h.extent(0), h.extent(1), h.extent(2),
                          h.extent(3));
    Kokkos::deep_copy(d, h);
  };
  dv(v, cln_sgv, "cln_sgv");
  dv(a0, cln_sga0, "cln_sga0");
  dv(a1, cln_sga1, "cln_sga1");
  dv(d1, cln_sgd1, "cln_sgd1");
  dv(d2, cln_sgd2, "cln_sgd2");
  Kokkos::deep_copy(cln_sgx, x2v);
  Kokkos::deep_copy(cln_sgtg, tg);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SphericalS1Check
//! \brief fatal on any option the sp / cs implicit operator does not support

void RadiationVet::SphericalS1Check(ParameterInput *pin) {
  std::string why;
  auto *pm = pmy_pack->pmesh;
  if (pm->multilevel) {why += " SMR/AMR;";}
  if (pm->mesh_indcs.nx1 != pm->mb_indcs.nx1) {
    why += " more than one MeshBlock along x1 (meshblock/nx1 must equal mesh/nx1);";
  }
  // cubed sphere
  if (cs_geom) {
    std::string wcs;
    // seam geometry needs the equiangular panel [-1, 1]^2 with nx2 = nx3, and the
    // one-sided tangential derivatives at a seam need >= 3 cells per block
    if (pm->mesh_indcs.nx2 != pm->mesh_indcs.nx3 ||
        pm->mesh_size.x2min != -1.0 || pm->mesh_size.x2max != 1.0 ||
        pm->mesh_size.x3min != -1.0 || pm->mesh_size.x3max != 1.0) {
      wcs += " a panel other than x2, x3 in [-1, 1] with mesh nx2 = nx3;";
    }
    if (pm->mb_indcs.nx2 < 3 || pm->mb_indcs.nx3 < 3) {
      wcs += " meshblock nx2 or nx3 < 3;";
    }
    if (!wcs.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet> on the cubed sphere supports only "
        << "time_scheme = be | hesdirk2 (and the spherical-polar list below); this "
        << "input has:"
        << wcs << why << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // full sphere: a polar face has zero area and must be reflecting (F_theta = 0)
  if (!cs_geom && pm->use_polar_boundary) {
    const bool p2lo = (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar);
    const bool p2hi = (pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar);
    if ((p2lo && ibc_x2min != VET_IBC_REFLECT) || (p2hi && ibc_x2max != VET_IBC_REFLECT))
        {
      why += " implicit_bc_x2min|max != auto | reflect at a polar face;";
    }
  }
  if (!why.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> on a spherical-polar mesh supports "
      << "only time_scheme = be | hesdirk2, one MeshBlock along x1, no SMR; this "
      << "input has:" << why
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (global_variable::my_rank == 0 && cs_geom) {
    std::cout << "<rad_vet>: cubed sphere (covariant cell F and gas coupling; skewed "
              << "transverse faces with the lagged cross term, mirror-pair seam faces "
              << "without the resample; implicit with areas/volumes/face distances from "
              << "Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  ", stretched radial grid)" : ")") << std::endl;
  } else if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: spherical-polar "
              << (pm->use_polar_boundary ? "FULL SPHERE [polar x2 faces reflecting at "
                                           "zero area] " : "wedge ")
              << (pm->use_grid_stretch_theta ? "[stretched theta] " : "")
              << "("
              << "implicit with areas/volumes/face distances from Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  " (stretched radial grid)" : "") << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CubedS1Init
//! \brief cubed sphere, after ImplicitInit: open the panel seams of the implicit
//! operator, tabulate the seam geometry (VetCsSeamGeom, incl. ghost rows), effective
//! transverse distances/areas and the face and cell-centre normals

void RadiationVet::CubedS1Init() {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ng = indcs.ng;
  const int n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  const int npan = pm->mesh_indcs.nx2;   // cells per panel edge (= nx3, checked)
  auto &bc = pmy_pack->pmb->mb_bcs;
  auto &msz = pmy_pack->pmb->mb_size;
  vetbcs = DualArray2D<BoundaryFlag>("vetbcs", nmb, 6);
  cs_seam = DualArray2D<int>("vetcsseam", nmb, 4);
  csg2 = DvceArray4D<Real>("vetcsg2", nmb, 4, n3, 2);
  csg3 = DvceArray4D<Real>("vetcsg3", nmb, 4, n2, 2);
  auto g2h = Kokkos::create_mirror_view(csg2);
  auto g3h = Kokkos::create_mirror_view(csg3);
  Kokkos::deep_copy(g2h, 0.0);
  Kokkos::deep_copy(g3h, 0.0);
  for (int m = 0; m < nmb; ++m) {
    for (int f = 0; f < 6; ++f) {
      BoundaryFlag b = bc.h_view(m,f);
      vetbcs.h_view(m,f) = (b == BoundaryFlag::panel) ? BoundaryFlag::block : b;
    }
    for (int f = 0; f < 4; ++f) {
      cs_seam.h_view(m,f) = (bc.h_view(m,2+f) == BoundaryFlag::panel) ? 1 : 0;
    }
    // the first active cell's index along the panel edge
    const int j0 = static_cast<int>(std::lround(0.5*(msz.h_view(m).x2min + 1.0)*npan));
    const int k0 = static_cast<int>(std::lround(0.5*(msz.h_view(m).x3min + 1.0)*npan));
    for (int sd = 0; sd < 2; ++sd) {
      for (int k = 0; k < n3; ++k) {
        const int q = k - indcs.ks + k0;
        if (q < 0 || q >= npan) continue;
        Real a, b, c, d;
        VetCsSeamGeom(q, npan, a, b, c, d);
        g2h(m,0,k,sd) = a; g2h(m,1,k,sd) = b; g2h(m,2,k,sd) = c; g2h(m,3,k,sd) = d;
      }
      for (int j = 0; j < n2; ++j) {
        const int q = j - indcs.js + j0;
        if (q < 0 || q >= npan) continue;
        Real a, b, c, d;
        VetCsSeamGeom(q, npan, a, b, c, d);
        g3h(m,0,j,sd) = a; g3h(m,1,j,sd) = b; g3h(m,2,j,sd) = c; g3h(m,3,j,sd) = d;
      }
    }
  }
  vetbcs.modify_host();
  vetbcs.sync_device();
  cs_seam.modify_host();
  cs_seam.sync_device();
  Kokkos::deep_copy(csg2, g2h);
  Kokkos::deep_copy(csg3, g3h);
  // effective transverse distances and areas (rad_vet.hpp)
  {
    auto *pc = pmy_pack->pcoord;
    const int n1 = indcs.nx1 + 2*ng;
    cs_dxf_eff = DvceFaceFld4D<Real>("vetcsdxf", nmb, n3, n2, n1);
    cs_area_eff = DvceFaceFld4D<Real>("vetcsarea", nmb, n3, n2, n1);
    Kokkos::deep_copy(cs_dxf_eff.x1f, pc->dxface.x1f);
    Kokkos::deep_copy(cs_dxf_eff.x2f, pc->dxface.x2f);
    Kokkos::deep_copy(cs_dxf_eff.x3f, pc->dxface.x3f);
    Kokkos::deep_copy(cs_area_eff.x1f, pc->area.x1f);
    Kokkos::deep_copy(cs_area_eff.x2f, pc->area.x2f);
    Kokkos::deep_copy(cs_area_eff.x3f, pc->area.x3f);
    auto dx2 = cs_dxf_eff.x2f, dx3 = cs_dxf_eff.x3f;
    auto a2 = cs_area_eff.x2f, a3 = cs_area_eff.x3f;
    auto sn2 = pc->sin_face_xi, sn3 = pc->sin_face_eta;
    auto x1v = pc->x1v, x1f = pc->xx1f;
    auto g2 = csg2, g3 = csg3;
    auto cse = cs_seam.d_view;
    const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
    const int is = indcs.is, ie = indcs.ie;
    par_for("vet_cs_eff2", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool sm = (j == js && cse(m,0) != 0) || (j == je+1 && cse(m,1) != 0);
      const int sd = (j == js) ? 0 : 1;
      if (sm) {
        dx2(m,k,j,i) = x1v(m,i)*g2(m,1,k,sd);
        a2(m,k,j,i) = 0.5*(x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i))*g2(m,0,k,sd);
      } else {
        dx2(m,k,j,i) = dx2(m,k,j,i)*sn2(m,k,j);
      }
    });
    par_for("vet_cs_eff3", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool sm = (k == ks && cse(m,2) != 0) || (k == ke+1 && cse(m,3) != 0);
      const int sd = (k == ks) ? 0 : 1;
      if (sm) {
        dx3(m,k,j,i) = x1v(m,i)*g3(m,1,j,sd);
        a3(m,k,j,i) = 0.5*(x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i))*g3(m,0,j,sd);
      } else {
        dx3(m,k,j,i) = dx3(m,k,j,i)*sn3(m,k,j);
      }
    });
  }
  // planar face normals nhat_xi = (a - tan(xi_f) q)/sqrt(1 + tan^2), likewise eta; the
  // face plane holds the origin, so the normal is constant along the column
  {
    Kokkos::realloc(cs_fn2, nmb, n2 + 1, 3);
    Kokkos::realloc(cs_fn3, nmb, n3 + 1, 3);
    auto f2h = Kokkos::create_mirror_view(cs_fn2);
    auto f3h = Kokkos::create_mirror_view(cs_fn3);
    auto &pan = pmy_pack->pmb->mb_panel;
    for (int m = 0; m < nmb; ++m) {
      Real fa[3], fb[3], fq[3];
      cubed_sphere::PanelFrame(pan.h_view(m), fa, fb, fq);
      for (int j = 0; j <= n2; ++j) {
        const Real x = std::tan(0.25*M_PI*LeftEdgeX(j - indcs.js, indcs.nx2,
                                                   msz.h_view(m).x2min,
                                                   msz.h_view(m).x2max));
        const Real nn = 1.0/std::sqrt(1.0 + x*x);
        for (int c = 0; c < 3; ++c) {f2h(m,j,c) = (fa[c] - x*fq[c])*nn;}
      }
      for (int k = 0; k <= n3; ++k) {
        const Real y = std::tan(0.25*M_PI*LeftEdgeX(k - indcs.ks, indcs.nx3,
                                                   msz.h_view(m).x3min,
                                                   msz.h_view(m).x3max));
        const Real nn = 1.0/std::sqrt(1.0 + y*y);
        for (int c = 0; c < 3; ++c) {f3h(m,k,c) = (fb[c] - y*fq[c])*nn;}
      }
    }
    Kokkos::deep_copy(cs_fn2, f2h);
    Kokkos::deep_copy(cs_fn3, f3h);
    // cell-centre normals (same formula at xi_c, eta_c)
    Kokkos::realloc(cs_cn2, nmb, n2, 3);
    Kokkos::realloc(cs_cn3, nmb, n3, 3);
    auto c2h = Kokkos::create_mirror_view(cs_cn2);
    auto c3h = Kokkos::create_mirror_view(cs_cn3);
    for (int m = 0; m < nmb; ++m) {
      Real fa[3], fb[3], fq[3];
      cubed_sphere::PanelFrame(pan.h_view(m), fa, fb, fq);
      for (int j = 0; j < n2; ++j) {
        const Real x = std::tan(0.25*M_PI*CellCenterX(j - indcs.js, indcs.nx2,
                                                     msz.h_view(m).x2min,
                                                     msz.h_view(m).x2max));
        const Real nn = 1.0/std::sqrt(1.0 + x*x);
        for (int c = 0; c < 3; ++c) {c2h(m,j,c) = (fa[c] - x*fq[c])*nn;}
      }
      for (int k = 0; k < n3; ++k) {
        const Real y = std::tan(0.25*M_PI*CellCenterX(k - indcs.ks, indcs.nx3,
                                                     msz.h_view(m).x3min,
                                                     msz.h_view(m).x3max));
        const Real nn = 1.0/std::sqrt(1.0 + y*y);
        for (int c = 0; c < 3; ++c) {c3h(m,k,c) = (fb[c] - y*fq[c])*nn;}
      }
    }
    Kokkos::deep_copy(cs_cn2, c2h);
    Kokkos::deep_copy(cs_cn3, c3h);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CheckGeometry
//! \brief mesh checks (Cartesian, spherical polar, cubed sphere), removed options

void RadiationVet::CheckGeometry(ParameterInput *pin) {
  MeshBlockPack *ppack = pmy_pack;
  // geometry guard: uniform Cartesian, a spherical-polar wedge or full sphere, or the
  // cubed sphere (checks in SphericalS1Check, VetRayInit); no relativistic coordinates
  sph_geom = false;
  cs_geom = false;
  {
    Mesh *pm_ = ppack->pmesh;
    std::string why;
    const bool sp_ = pm_->use_spherical_polar && !pm_->use_cubed_sphere;
    const bool cs_ = pm_->use_cubed_sphere && !pm_->use_spherical_polar;
    if (pm_->use_spherical_polar && !sp_) {why += " mesh/use_spherical_polar";}
    // full sphere (use_polar_boundary) and theta stretch allowed on sp: a polar face is a
    // physical face of area 0 (reflect is exact there)
    const bool full_ok = sp_;
    if (pm_->use_polar_boundary) {
      if (!full_ok) {
        why += " mesh/use_polar_boundary (the full sphere needs transport = implicit on "
               "spherical polar)";
      }
    } else if (sp_ && (pm_->mesh_size.x2min <= 1.0e-10 ||
                       pm_->mesh_size.x2max >= M_PI - 1.0e-10)) {
      why += " a theta range that reaches a pole without mesh/use_polar_boundary (use "
             "ix2_bc = ox2_bc = polar with use_polar_boundary, or a wedge clear of the "
             "poles)";
    }
    if (pm_->use_grid_stretch_r && !sp_ && !cs_) {why += " mesh/use_grid_stretch_r";}
    if (pm_->use_grid_stretch_r_poly && !sp_ && !cs_) {
      why += " mesh/use_grid_stretch_r_poly";
    }
    if (pm_->use_grid_stretch_theta && !full_ok) {
      why += " mesh/use_grid_stretch_theta (needs transport = implicit on spherical "
             "polar)";
    }
    if (ppack->pcoord != nullptr &&
        (ppack->pcoord->is_special_relativistic ||
         ppack->pcoord->is_general_relativistic ||
         ppack->pcoord->is_dynamical_relativistic)) {
      why += " relativistic coordinates";
    }
    if (!why.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet> (explicit and implicit transport) supports a "
        << "uniform Cartesian mesh, a spherical-polar wedge or the full sphere and the "
        << "cubed sphere; this input sets:"
        << why << "." << std::endl
        << "The VET kernels would run with Cartesian uniform-dx arithmetic on it."
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    sph_geom = sp_;
    if (cs_) {
      sph_geom = true;
      cs_geom = true;
    }
  }
  nstage = VET_NSTAGE;
  impl_cfl = -1.0;
  marshak_q = 0.5;
  ibc_x1min = VET_IBC_MARSHAK;
  ibc_x1max = VET_IBC_MARSHAK;
  iflux_x1min = 0.0;
  iflux_x1max = 0.0;
  iebath_x1min = 0.0;
  iebath_x1max = 0.0;
  impl_nstep = 0.0;
  impl_itsum = 0.0;
  impl_itmax = 0.0;
  trans_on = false;
  trans_x3 = false;
  impl_prec_fb = -1;
  bcg_capped = false;
  bcg_r0rel = 0.0;
  ew_fprev = 0.0;
  ew_etaprev = 0.0;
  pred_ok = false;
  pred_dt = 0.0;
  vimp_now = false;
  time_scheme = 0;
  t2_ok = false;
  t2_solve = 0;
  t2_fail = false;
  t2_dbg_fail = -1;
  t2_nstep = 0.0;
  t2_nbe = 0.0;
  t2_nfall = 0.0;
  t2_dtprev = 0.0;
  t2_vprev = false;
  t2_nclip = 0.0;
  pred2_ok = false;
  pred2_dt = 0.0;
  t2_lin_save = 0.0;
  // Anderson acceleration
  aa_nh = 0;
  aa_head = 0;
  aa_hasp = false;
  aa_fnp = -1.0;
  st_edges = true;
  impl_prec = 0;
  e_floor = pin->GetOrAddReal("rad_vet","e_floor",(FLT_MIN));

  // removed explicit face-flux options: only their old defaults are accepted
  if (pin->DoesParameterExist("rad_vet","scaled_prefactor") &&
      pin->GetReal("rad_vet","scaled_prefactor") != 20.0) {
    std::cout << "### FATAL ERROR: <rad_vet>/scaled_prefactor = "
              << pin->GetString("rad_vet","scaled_prefactor") << " was removed (only 20)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

}  // namespace radvet
