//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_sph2.inc
//! \brief accsp-1102 (USER 11-02 "try sp with stretched theta with gd and wb tide"):
//! vet_gd on a STRETCHED theta (mesh/theta_plateau, mesh/use_grid_stretch_theta) and on
//! the FULL sphere (mesh/use_polar_boundary).  Included once at the end of
//! rad_m1_vetgd.cpp (member definitions; no new translation unit).
//!
//! (a) STRETCHED THETA.  The sweep took the band rows at the uniform positions x2min +
//!     (j - js + 0.5) dx2 and the upwind row as j + (thu - th)/dx2.  With vgd_sth the
//!     rows are the face MIDPOINTS of the real faces (the node table vgd_thn, computed
//!     from the mesh's own theta map), and the fractional row is linear in theta between
//!     the two bracketing nodes (VgdThFj).  On a uniform grid this is the old formula up
//!     to round-off; vgd_sth is false there (bitwise), env VGD_THTAB=1 forces the table
//!     path for the round-off gate.
//! (b) FULL SPHERE.  The band rows beyond a pole continue the theta of the EXTENDED
//!     chart (-theta north, 2 pi - theta south, i.e. the point (|theta|, phi + pi)): the
//!     node table is mirrored there, the halo slot beyond the pole takes the pole partner
//!     (lx2, lx3 + dk + nb3/2) with its rows reversed (VetGdHaloInit / VetGdHalo; an
//!     on-rank partner is read in-kernel with VgdPoleRow), and an upwind point behind the
//!     axis is located in that chart (VgdPoleChart).  The directions are global Cartesian,
//!     so no direction map is needed across the pole.  The band depth is computed from
//!     sin(vet_gd_band_colat) instead of sin(0); the longer near-tangent reads of the
//!     polar rings are clamped to the shell's stored band and counted (as every clamp).

namespace radm1 {

namespace {
// the mesh's theta face of logical index n (0 .. nx2 the mesh faces), the map of
// Coordinates (LeftEdgeX + ApplyThetaStretch)
double VgdThetaFace(const Mesh *pm, const int n) {
  const int nx = pm->mesh_indcs.nx2;
  const double x0 = pm->mesh_size.x2min, x1 = pm->mesh_size.x2max;
  const double rx = static_cast<double>(n)/static_cast<double>(nx);
  Real t = static_cast<Real>(rx*x1 - rx*x0 + x0);
  if (pm->use_grid_stretch_theta) {
    ApplyThetaStretch(pm->fStretchTheta, pm->fStretchThetaPlat, t);
  }
  return static_cast<double>(t);
}
} // namespace

int RadiationM1::VgdPoleSlot(const int o, const int lx2) const {
  if (!vgd_pol) {return 0;}
  const int oo = (o < 4) ? o : (o + 1);
  const int dj = oo%3 - 1;
  if (dj < 0 && lx2 == 0) {return 1;}
  if (dj > 0 && lx2 == vgd_nb2 - 1) {return 2;}
  return 0;
}

//----------------------------------------------------------------------------------------
//! the narrowest lateral spacing (radians of arc on the unit sphere) of the band-depth
//! formula: min over the mesh's theta cells of dtheta, and sin(theta_lim) dphi with
//! theta_lim the theta limit of a wedge (as before) or, on the full sphere,
//! vet_gd_band_colat (the polar rings poleward of it read beyond the band: clamped)

double RadiationM1::VetGdThetaDmin() const {
  Mesh *pm = pmy_pack->pmesh;
  const int nx = pm->mesh_indcs.nx2;
  double dth = 1.0e300;
  for (int n = 0; n < nx; ++n) {
    dth = std::min(dth, VgdThetaFace(pm, n + 1) - VgdThetaFace(pm, n));
  }
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  const double smin = vgd_pol ? std::sin(vgd_bcolat*M_PI/180.0)
                              : std::min(std::sin(t0), std::sin(t1));
  auto &mbs = pmy_pack->pmb->mb_size;
  return std::min(dth, smin*mbs.h_view(0).dx3);
}

//----------------------------------------------------------------------------------------
//! the node table vgd_thn(m, band row) (after vgd_w is known)

void RadiationM1::VetGdThetaInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int w = vgd_w, nx2b = indcs.nx2, c2w = nx2b + 2*w;
  const int nx = pm->mesh_indcs.nx2;
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  std::vector<double> node(nx);
  for (int n = 0; n < nx; ++n) {
    node[n] = 0.5*(VgdThetaFace(pm, n) + VgdThetaFace(pm, n + 1));
  }
  // a global row J (any integer within one period/mirror of the mesh) -> its theta
  auto nodeg = [&](const int J) -> double {
    if (J >= 0 && J < nx) {return node[J];}
    if (vgd_pol) {   // the extended chart beyond a pole: the mirror row
      if (J < 0) {return -node[-1 - J];}
      return 2.0*M_PI - node[2*nx - 1 - J];
    }
    // a periodic wedge: the image of the other side's row (the band holds its copy)
    if (J < 0) {return node[J + nx] - (t1 - t0);}
    return node[J - nx] + (t1 - t0);
  };
  vgd_thn_h.assign(static_cast<size_t>(nmb)*c2w, 0.0);
  Kokkos::realloc(vgd_thn, nmb, c2w);
  auto tn_h = Kokkos::create_mirror_view(vgd_thn);
  const int g0 = pmy_pack->gids;
  for (int m = 0; m < nmb; ++m) {
    const int lx2 = pm->lloc_eachmb[g0 + m].lx2;
    for (int jb = 0; jb < c2w; ++jb) {
      const int J = lx2*nx2b + (jb - w);
      const double v = nodeg(J);
      vgd_thn_h[static_cast<size_t>(m)*c2w + jb] = v;
      tn_h(m,jb) = v;
    }
    for (int jb = 1; jb < c2w; ++jb) {
      if (!(tn_h(m,jb) > tn_h(m,jb-1))) {
        VgdFatal("accsp-1102: the theta node table is not increasing (block " +
                 std::to_string(m) + ", band row " + std::to_string(jb) +
                 "): the band is deeper than one period/mirror of the mesh");
      }
    }
  }
  Kokkos::deep_copy(vgd_thn, tn_h);
  // (the pole slot table vgd_hpol is made by VetGdHaloInit, which knows nb2)
  {
    const int c1 = indcs.nx1 + 2*indcs.ng;
    Kokkos::realloc(vgd_wsd, 2, c1);
    auto wd_h = Kokkos::create_mirror_view(vgd_wsd);
    for (int i = 0; i < c1; ++i) {
      wd_h(0,i) = (i < static_cast<int>(vgd_wsi.size())) ? vgd_wsi[i] : 0;
      wd_h(1,i) = (i < static_cast<int>(vgd_wso.size())) ? vgd_wso[i] : 0;
    }
    Kokkos::deep_copy(vgd_wsd, wd_h);
  }
  if (global_variable::my_rank == 0) {
    double dmin = 1.0e300, dmax = 0.0;
    for (int n = 0; n < nx; ++n) {
      const double d = VgdThetaFace(pm, n + 1) - VgdThetaFace(pm, n);
      dmin = std::min(dmin, d);
      dmax = std::max(dmax, d);
    }
    std::cout << "<rad_m1> vet_gd (accsp-1102): theta node table ("
              << (pm->use_grid_stretch_theta ? "STRETCHED theta" : "uniform theta")
              << ", rows " << dmin*180.0/M_PI << " .. " << dmax*180.0/M_PI << " deg)"
              << (vgd_pol ? "; FULL SPHERE: pole slots take the partner (lx2, lx3 + "
                            "nb3/2) with reversed rows, extended chart behind the axis, "
                            "band depth from vet_gd_band_colat = " : "")
              << (vgd_pol ? std::to_string(vgd_bcolat) + " deg (reads beyond the shell's"
                            " band are clamped and counted)" : std::string(""))
              << std::endl;
  }
}

} // namespace radm1
