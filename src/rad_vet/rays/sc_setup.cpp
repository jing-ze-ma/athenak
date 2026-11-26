//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_setup.cpp
//! \brief curved-mesh rays: keys, arrays, halo partners, GPU IPC handles.
//!
//!  called from: rays/sc_build.cpp, rays/sc_cs_sweep.cpp, rays/sc_source.cpp,
//!    rays/sc_sph_globalband.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rad_vet_init.cpp, rays/sc_cs_setup.cpp, rays/sc_directions.cpp,
//!    rays/sc_sph_globalband.cpp, rays/sc_sph_walls.cpp, rays/sc_tables.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {

// moment errors of a direction set (rows x y z w): max |sum w - 1|, |sum w n|,
// |sum w n_i n_j - delta_ij/3|, |sum w n_i n_j n_k n_l - (dd + dd + dd)/15|
void ScMomentCheck(const std::vector<double> &d, double e[4]) {
  const int n = static_cast<int>(d.size()/4);
  double m0 = 0.0, m1[3] = {0, 0, 0}, m2[3][3] = {{0}}, m4[3][3][3][3] = {{{{0}}}};
  for (int q = 0; q < n; ++q) {
    const double *x = &d[4*q], w = d[4*q+3];
    m0 += w;
    for (int i = 0; i < 3; ++i) {
      m1[i] += w*x[i];
      for (int j = 0; j < 3; ++j) {
        m2[i][j] += w*x[i]*x[j];
        for (int k = 0; k < 3; ++k) {
          for (int l = 0; l < 3; ++l) {m4[i][j][k][l] += w*x[i]*x[j]*x[k]*x[l];}
        }
      }
    }
  }
  e[0] = std::fabs(m0 - 1.0);
  e[1] = e[2] = e[3] = 0.0;
  for (int i = 0; i < 3; ++i) {
    e[1] = std::max(e[1], std::fabs(m1[i]));
    for (int j = 0; j < 3; ++j) {
      e[2] = std::max(e[2], std::fabs(m2[i][j] - ((i == j) ? 1.0/3.0 : 0.0)));
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          const double iso = ((i == j && k == l) + (i == k && j == l) +
                              (i == l && j == k))/15.0;
          e[3] = std::max(e[3], std::fabs(m4[i][j][k][l] - iso));
        }
      }
    }
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScInit
//! \brief direction table, intensity array, exchange object, wall maps (from VetRayInit)

void RadiationVet::ScInit() {
  Mesh *pm = pmy_pack->pmesh;
  // nx3 = 1: axisymmetric; the phi band columns are the cell's own column with the
  // directions rotated by their phi offset (must map the set onto itself)
  if (!pm->multi_d) {ScFatal("needs a 3-D wedge or a 2-D axisymmetric mesh");}
  if (!pm->three_d && vgd_slab != 1 && (impl_beam_hr || impl_beam_fs)) {
    ScFatal("2-D axisymmetric: no implicit_flux_beam");
  }
  if (ibc_x1max != VET_IBC_MARSHAK &&
      !(ibc_x1max == VET_IBC_REFLECT || ibc_x1max == VET_IBC_OUTFLOW) &&
      !(vgd_slab == 1 && vet_x1per)) {
    // the sweep has a vacuum top, or a mirror (ScTopMirror) for the clean scheme
    ScFatal("outer x1 needs implicit_bc_x1max = marshak or "
            "reflect with scheme = clean");
  }
  if (vgd_slab != 1 && (ibc_x2min == VET_IBC_MARSHAK || ibc_x2max == VET_IBC_MARSHAK ||
                         ibc_x3min == VET_IBC_MARSHAK || ibc_x3max == VET_IBC_MARSHAK)) {
    // the curved sweep has no vacuum transverse face (the Cartesian slab has)
    ScFatal("implicit_bc_x2*/x3* = marshak is not wired into the vet_gd sweep");
  }
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  if (vgd_slab == 1) {ScSlabGeom();}   // (the Cartesian mesh keeps no x1 arrays)
  auto cx1v = (vgd_slab == 1) ? vgd_sx1v : pmy_pack->pcoord->x1v;
  auto cx1f = (vgd_slab == 1) ? vgd_sx1f : pmy_pack->pcoord->xx1f;
  // cubed sphere (rays/sc_cs_sweep.cpp): wedge-only options are refused; the cs band
  // is filled by its own exact copies (compact lists, see ScCsSweep); no walls
  vgd_csg = pm->use_cubed_sphere;
  // <mesh>/layout_invariant: the sweep's (theta, phi) from the global index
  // (rays/sc_band.hpp); a uniform single-level sp mesh
  vgd_lg.on = false;
  if (!vgd_csg && pm->layout_invariant) {
    vgd_lg.on = true;
    vgd_lg.x2a = pm->mesh_size.x2min;
    vgd_lg.dx2 = (pm->mesh_size.x2max - pm->mesh_size.x2min)/pm->mesh_indcs.nx2;
    vgd_lg.x3a = pm->mesh_size.x3min;
    vgd_lg.dx3 = (pm->mesh_size.x3max - pm->mesh_size.x3min)/pm->mesh_indcs.nx3;
    Kokkos::realloc(vgd_lg.o, pmy_pack->nmb_thispack, 2);
    vgd_lg.oh = Kokkos::create_mirror_view(vgd_lg.o);
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      const auto &ll = pm->lloc_eachmb[pmy_pack->gids + m];
      if (ll.level != pm->root_level) {
        ScFatal("layout_invariant: the global vet_gd geometry needs a single-level "
                "mesh");
      }
      vgd_lg.oh(m,0) = static_cast<int>(ll.lx2)*pm->mb_indcs.nx2;
      vgd_lg.oh(m,1) = static_cast<int>(ll.lx3)*pm->mb_indcs.nx3;
    }
    Kokkos::deep_copy(vgd_lg.o, vgd_lg.oh);
  }
  if (vgd_csg) {
    std::string bad;
    if (!bad.empty()) {ScFatal("not available on the cubed sphere:" + bad);}
    vgd_wint = false;
    vgd_ipc_on = false;
    vgd_hx_on = false;
    // the unfused twin always uses this build's cut (as the fused twin); the previous
    // build's cut would make it decomposition-dependent
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd on the CUBED SPHERE: seam "
                << "strip + vertex triangle upwind stencil, exact per-shell band copies "
                << "(the sp halo lists / IPC / overlap keys do not apply)"
                << std::endl;
    }
  }
  // stretched theta / full sphere (rays/sc_sph_globalband.cpp): options that assume the
  // uniform index formula are refused; the compact halo falls back to the dense one
  // (same values)
  if (vgd_sth) {
    std::string bad;
    if (vgd_thmir || vgd_thvac > 0) {bad += " vet_gd_theta_wall = mirror|vacuum|extrap;";}
    if (!bad.empty()) {
      ScFatal("not available on a stretched theta / the full sphere:" +
              bad);
    }
    // with vet_gd_band_global the compact halo is the keep-tested global-band exchange
    // (ScGbHaloC); otherwise the dense halo (the slot keep test assumes uniform theta)
    if (vgd_hcomp > 0 && !vgd_gband) {
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd: vet_gd_halo_compact = " << vgd_hcomp
                  << " -> 0 (the dense per-shell slot halo; the compact keep test is "
                  << "ported only to vet_gd_band_global; speed only)" << std::endl;
      }
      vgd_hcomp = 0;
    }
    if (vgd_hcomp > 0) {
      // the slot-path speed machinery does not apply to the global band
      vgd_hx_on = false;
    }
    vgd_ipc_on = false;
  }
  // several MeshBlocks along x1: levels swept in turn, the boundary shell handed to the
  // next (ScX1Hand); flags 1 innermost, 2 outermost
  vgd_nb1 = 1;
  for (int g = 0; g < pm->nmb_total; ++g) {
    vgd_nb1 = std::max(vgd_nb1, pm->lloc_eachmb[g].lx1 + 1);
  }
  vgd_x1nb.assign(2*nmb, -1);
  vgd_x1lv.assign(nmb, 0);
  Kokkos::realloc(vgd_x1f, nmb);
  auto xf_h = Kokkos::create_mirror_view(vgd_x1f);
  for (int m = 0; m < nmb; ++m) {
    const LogicalLocation &l = pm->lloc_eachmb[pmy_pack->gids + m];
    vgd_x1lv[m] = l.lx1;
    xf_h(m) = (l.lx1 == 0 ? 1 : 0) + (l.lx1 == vgd_nb1 - 1 ? 2 : 0);
    if (l.lx1 == vgd_nb1 - 1 && ibc_x1max == VET_IBC_REFLECT &&
        (vgd_slab != 1 || vet_bc_mirror_hi)) {xf_h(m) += 4;}
    xf_h(m) += 8*l.lx1;
    for (int g = 0; g < pm->nmb_total; ++g) {
      const LogicalLocation &o = pm->lloc_eachmb[g];
      if (o.level == l.level && o.panel == l.panel && o.lx2 == l.lx2 && o.lx3 == l.lx3 &&
          std::abs(o.lx1 - l.lx1) == 1) {
        vgd_x1nb[2*m + (o.lx1 > l.lx1 ? 1 : 0)] = g;
      }
    }
  }
  Kokkos::deep_copy(vgd_x1f, xf_h);
  Kokkos::realloc(vgd_x1mo, nmb);
  auto mo_h = Kokkos::create_mirror_view(vgd_x1mo);
  for (int m = 0; m < nmb; ++m) {mo_h(m) = m;}
  std::stable_sort(mo_h.data(), mo_h.data() + nmb,
                   [&](const int a, const int b) {return vgd_x1lv[a] < vgd_x1lv[b];});
  Kokkos::deep_copy(vgd_x1mo, mo_h);
  if (vgd_nb1 > 1) {
    if (vgd_csg || vgd_pol || indcs.nx1 < 2) {
      ScFatal("several MeshBlocks along x1 need a wedge (no full sphere) and "
              "MeshBlock nx1 >= 2");
    }
    vgd_hcomp = 0;   // the dense per-shell halo (the compact lists assume one level)
  }
  // vgd_raw: the raw band where the exact-list band would clamp (same values where
  // nothing clamps)
  const bool raw_ok = !vgd_csg && !vgd_pol && !vgd_sth;
  vgd_raw = false;
  // nx3 = 1: the dense halo (the compact / exact lists mark 3-D band reads)
  if (!pm->three_d) {vgd_hcomp = 0;}
  auto raw_halo = [&]() {
    if (vgd_raw || !pm->three_d) {   // the compact halo: ScGbHaloC
      vgd_hx_on = false;
      vgd_ipc_on = vgd_ipc_on && vgd_raw && pm->three_d;   // ScGbHaloC's own IPC
    }
  };
  raw_halo();
  // vet_gd_wall_sym / vet_gd_theta_seam_sym do nothing without a phi wall or a
  // periodic theta seam
  const bool sp3 = !vgd_csg && !vgd_pol && pm->three_d;
  vgd_wsym = vgd_wsym && sp3 &&
             pm->mesh_size.x3max - pm->mesh_size.x3min < 2.0*M_PI - 1.0e-12;
  vgd_wst = vgd_wst && sp3 &&
            pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic;
  // the direction set: any table of unit vectors + weights with sum w = 1,
  // sum w nn = I/3
  if (vgd_gl_nmu > 0 && vgd_slab != 1) {
    // product set: the azimuth step must divide the phi period (exact walls);
    // nphi >= 4 on the cubed sphere
    const double per = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const double m = per*vgd_gl_nphi/(2.0*M_PI);
    if (vgd_gl_nphi < 4 ||
        (!vgd_csg && !vgd_wsym && std::fabs(m - std::round(m)) > 1.0e-9)) {
      ScFatal("vet_gd_gl_nphi must be >= 4 with 2 pi/nphi dividing the phi period of "
              "the wedge (exact wall maps), or vet_gd_wall_sym = true");
    }
  }
  if (vgd_slab == 1) {
    // Cartesian slab: the vet_sc set (mu1, mu2, mu3, w), x1 in the role of r
    auto ah = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vet_ang);
    vgd_base.resize(4*vet_nray);
    for (int q = 0; q < 4*vet_nray; ++q) {vgd_base[q] = ah(q/4, q%4);}
  } else {
    ScBuildSet(vgd_base);
  }
  // the set's image under (x, y, z) -> (c x - t s y, s x + t c y, tz z): a turn about z
  // (t = 1), a meridian mirror (t = -1) or the equatorial mirror (tz = -1); -1 where
  // the image is not in the set
  const size_t nb = vgd_base.size()/4;
  auto img = [&](const double c, const double s, const double t, const double tz) {
    std::vector<int> p(nb, -1);
    for (size_t q = 0; q < nb; ++q) {
      const double *u = &vgd_base[4*q];
      const double mx = c*u[0] - t*s*u[1], my = s*u[0] + t*c*u[1];
      for (size_t r = 0; r < nb && p[q] < 0; ++r) {
        const double *v = &vgd_base[4*r];
        if (mx*v[0] + my*v[1] + tz*u[2]*v[2] > 1.0 - 1.0e-12 &&
            std::fabs(v[3] - u[3]) <= 1.0e-12*u[3]) {p[q] = static_cast<int>(r);}
      }
    }
    return p;
  };
  auto dev = [&](const std::vector<int> &p, const char *nm) {
    if (*std::min_element(p.begin(), p.end()) < 0) {
      ScFatal(std::string(nm) + ": the set is not mirror-symmetric about the wedge "
              "mid-plane");
    }
    DvceArray1D<int> v(std::string(nm), nb);
    auto h = Kokkos::create_mirror_view(v);
    for (size_t q = 0; q < nb; ++q) {h(q) = p[q];}
    Kokkos::deep_copy(v, h);
    return v;
  };
  if (vgd_wsym && ScRotOn()) {
    // a set turning about z on a wedge of whole azimuthal periods: the turn by the period
    // W maps every turned set onto itself (the two turns commute), so its phi wall maps
    // are exact and there is nothing to symmetrise (the key does nothing)
    const double w = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const std::vector<int> pw = img(std::cos(w), std::sin(w), 1.0, 1.0);
    if (!vgd_wint || vgd_rot3 || *std::min_element(pw.begin(), pw.end()) < 0) {
      ScFatal("vet_gd_wall_sym with a rotating set needs vet_gd_wall_interp = true, a "
              "turn about z and a wedge of whole azimuthal periods of the set");
    }
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_wall_sym: the turning set's phi wall maps are exact "
                << "(no effect)" << std::endl;
    }
    vgd_wsym = false;
  }
  if (vgd_wsym) {
    // vet_gd_wall_sym: turn the set about z by the least angle that puts one of its
    // mirror planes (GL k pi/nphi, the octahedral sets k pi/4) on the wedge mid-plane
    if (!vgd_wint) {
      ScFatal("vet_gd_wall_sym needs vet_gd_wall_interp = true");
    }
    const double pc = 0.5*(pm->mesh_size.x3min + pm->mesh_size.x3max);
    const double pp = (vgd_gl_nmu > 0) ? M_PI/vgd_gl_nphi : 0.25*M_PI;
    const double da = pc - pp*std::round(pc/pp);
    const double ca = std::cos(da), sa = std::sin(da);
    if (std::fabs(da) > 1.0e-12) {
      for (size_t q = 0; q < nb; ++q) {
        const double x = vgd_base[4*q], y = vgd_base[4*q+1];
        vgd_base[4*q] = ca*x - sa*y;
        vgd_base[4*q+1] = sa*x + ca*y;
      }
    }
    vgd_mir = dev(img(std::cos(2.0*pc), std::sin(2.0*pc), -1.0, 1.0), "vet_gd_wall_sym");
    // inexact phi wall map (the turn by the period is not a set map): the clean
    // scheme's seam faces in both frames (RayFaceSums)
    vgd_wsw = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const std::vector<int> pw = img(std::cos(vgd_wsw), std::sin(vgd_wsw), 1.0, 1.0);
    vgd_wsx = *std::min_element(pw.begin(), pw.end()) < 0;
  }
  if (vgd_wst) {
    // vet_gd_theta_seam_sym: the clean scheme's theta seam faces in both frames, with the
    // equatorial mirror (a theta wedge centred on the equator).  A turn about z keeps the
    // set's equatorial mirror (the permutation of every turned set); a 3-D one does not
    if (!vgd_wint || (ScRotOn() && vgd_rot3) ||
        std::fabs(pm->mesh_size.x2min + pm->mesh_size.x2max - M_PI) > 1.0e-12) {
      ScFatal("vet_gd_theta_seam_sym needs a theta wedge centred on the equator, "
              "vet_gd_wall_interp = true and no 3-D set rotation");
    }
    vgd_mirz = dev(img(1.0, 0.0, 1.0, -1.0), "vet_gd_theta_seam_sym");
  }
  // the moments of the set, as built and after a fixed random 3-D rotation
  if (global_variable::my_rank == 0) {
    double e[4], er[4];
    ScMomentCheck(vgd_base, e);
    const double ax[3] = {0.3713906763541037, -0.5570860145311556, 0.7427813527082074};
    const double an = 0.987654321, ca = std::cos(an), sa = std::sin(an);
    std::vector<double> rb(vgd_base);
    for (size_t q = 0; q < rb.size()/4; ++q) {   // Rodrigues
      const double *x = &vgd_base[4*q];
      const double kx = ax[1]*x[2] - ax[2]*x[1], ky = ax[2]*x[0] - ax[0]*x[2],
                   kz = ax[0]*x[1] - ax[1]*x[0];
      const double kd = ax[0]*x[0] + ax[1]*x[1] + ax[2]*x[2];
      rb[4*q] = x[0]*ca + kx*sa + ax[0]*kd*(1.0 - ca);
      rb[4*q+1] = x[1]*ca + ky*sa + ax[1]*kd*(1.0 - ca);
      rb[4*q+2] = x[2]*ca + kz*sa + ax[2]*kd*(1.0 - ca);
    }
    ScMomentCheck(rb, er);
    std::cout << "<rad_vet> vet_gd direction set: " << vgd_base.size()/4
              << " directions; moment errors |sum w - 1| " << e[0] << ", |sum w n| "
              << e[1]
              << ", |sum w nn - I/3| " << e[2] << " (rotated " << er[2]
              << "), 4th-moment anisotropy " << e[3] << " (rotated " << er[3] << ")"
              << std::endl;
  }
  vgd_n = static_cast<int>(vgd_base.size()/4);
  // the set's z-rotation period: product 2 pi/nphi, geodesic pi, HEALPix / LQ_N pi/2
  vgd_zper = (vgd_gl_nmu > 0) ? 2.0*M_PI/vgd_gl_nphi : ((vgd_ls < 0) ? M_PI : 0.5*M_PI);
  if (!pm->three_d && vgd_slab != 1) {
    const double per = pm->mesh_size.x3max - pm->mesh_size.x3min;
    const double c = std::cos(per), s = std::sin(per);
    for (int q = 0; q < vgd_n; ++q) {
      const double *u = &vgd_base[4*q];
      double bd = -2.0;
      for (int r = 0; r < vgd_n; ++r) {
        const double *v = &vgd_base[4*r];
        bd = std::max(bd, (c*u[0] - s*u[1])*v[0] + (s*u[0] + c*u[1])*v[1] + u[2]*v[2]);
      }
      if (bd < 1.0 - 1.0e-12) {
        ScFatal("2-D axisymmetric: the phi extent of the mesh must rotate the direction "
                "set onto itself (e.g. 2 pi/vet_gd_gl_nphi, pi/2 for HEALPix)");
      }
    }
  }
  if (vgd_n > VET_VGD_NMAX) {ScFatal("too many directions");}
  const int n = vgd_n;
  Kokkos::realloc(vgd_dir, n, 4);
  // (one shell's n intensities, or the 6 face-sum columns of RayFaceSums: 4 sums, the
  // volume below and the volume of the column's levels)
  // (also the slab's 12 face-sum slots, ScSlabFaceX1)
  Kokkos::realloc(vgd_x1b, nmb, std::max(n, 12)*indcs.nx3*indcs.nx2);
  Kokkos::realloc(vgd_x1r, nmb, std::max(n, 12)*indcs.nx3*indcs.nx2);
  // lateral band: the deepest reach of an upwind point in cells over every shell (cell
  // sizes r dtheta, r sin(theta) dphi); on the cubed sphere ScCsInit
  if (vgd_csg) {
    ScCsInit();
  } else {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cx1v);
    auto &mbs = pmy_pack->pmb->mb_size;
    const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
    const double smin = std::min(std::sin(t0), std::sin(t1));
    // stretched theta / full sphere: the narrowest row and the band latitude limit
    // (ScThetaDmin)
    const double dmin = vgd_sth ? ScThetaDmin()
                             : std::min(vgd_lg.Dx2H(mbs, 0),
                                        smin*vgd_lg.Dx3H(mbs, 0));
    // shell s is read (a) by s-1 inward, chord <= sqrt(r_s^2 - r_{s-1}^2), (b) by s+1
    // outward, chord <= sqrt(r_{s+1}^2 - r_s^2), (c) by s itself (turn / inner face),
    // chord <= 2 sqrt(r_s^2 - r_lo^2); angular offset <= chord / r_lo; + 2 cells
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cx1f);
    // Cartesian slab: lateral cells per unit x1 chord, max over the set (constant feet)
    const bool slb = (vgd_slab == 1);
    double tsl = 0.0;
    // (a set turning about x1: the bound of every turn, |n_perp| over the finer side)
    const bool srot = slb && vsc_rot && ScRotOn();
    const double dxm = std::min(vgd_lg.Dx2H(mbs, 0),
                                pm->three_d ? vgd_lg.Dx3H(mbs, 0) : 1.0e300);
    // (a turning vet_ring changes n_1: both ends of its turn range, VetRotBeta, where
    // the rays come closest to the x1 layers)
    const bool rring = srot && vet_ring_n > 0;
    for (int e = 0; e < (rring ? 2 : 1); ++e) {
      std::vector<double> sb(vgd_base);
      if (rring) {
        VetRotSet(((e == 0) ? -0.4 : 0.4)*2.0*M_PI/vet_ring_n);
        auto ah = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vet_ang);
        for (int q = 0; q < 4*vet_nray; ++q) {sb[q] = ah(q/4, q%4);}
      }
      for (size_t q = 0; slb && q < sb.size()/4; ++q) {
        const double *u = &sb[4*q];
        const double t = (srot && !rring) ? std::sqrt(u[1]*u[1] + u[2]*u[2])/dxm
                         : std::max(std::fabs(u[1])/vgd_lg.Dx2H(mbs, 0), pm->three_d ?
                                    std::fabs(u[2])/vgd_lg.Dx3H(mbs, 0) : 0.0);
        tsl = std::max(tsl, t/std::fabs(u[0]));
      }
    }
    if (rring) {VetRotSet(vsc_beta);}
    // pass 1: the raw band, where pass 0's exact lists would clamp
    for (int pass = 0; pass < 2; ++pass) {
      int nmax = vgd_raw ? std::numeric_limits<int>::max()
                         : (pm->three_d ? std::min(indcs.nx2, indcs.nx3) : indcs.nx2);
      // vet_gd_band_global: the band is filled from any block, so its depth is capped
      // by the sphere (half a ring in phi, pole-to-pole in theta): layout independent
      if (vgd_gband) {
        if (!vgd_pol) {
          ScFatal("vet_gd_band_global needs the full sphere (mesh/use_polar_boundary)");
        }
        nmax = std::min(pm->mesh_indcs.nx3/2, pm->mesh_indcs.nx2);
      }
      vgd_wsh.assign(indcs.nx1 + 2*indcs.ng, 0);
      vgd_wsi.assign(indcs.nx1 + 2*indcs.ng, 0);
      vgd_wso.assign(indcs.nx1 + 2*indcs.ng, 0);
      int wneed = 0;
      for (int m = 0; m < ((vgd_nb1 > 1) ? nmb : 1); ++m) {
        // shells next to another level: radii from the x1 ghosts
        const bool lo = (xf_h(m) & 1), hi = (xf_h(m) & 2);
        for (int l = indcs.is; l <= indcs.ie; ++l) {
          const double r = x1v_h(m,l);
          const double rlo = (l > indcs.is || !lo) ? x1v_h(m,l-1) : x1f_h(m,l);
          const bool up = (l < indcs.ie || !hi);
          // inward-pass halo of shell s serves (a) and (c), the outward-pass halo (b);
          // slab: the x1 chords to the shells below and above
          const double ain = slb ? (r - rlo)
                                 : 2.0*std::sqrt(std::max(r*r - rlo*rlo, 0.0));
          const double aout = !up ? 0.0 : slb ? (x1v_h(m,l+1) - r)
                                      : std::sqrt(x1v_h(m,l+1)*x1v_h(m,l+1) - r*r);
          auto cells = [&](const double a) {
            return static_cast<int>(std::ceil(slb ? a*tsl : a/(rlo*dmin))) + 2 + vgd_xr;
          };
          const int wl = cells(std::max(ain, aout));
          wneed = std::max(wneed, wl);
          vgd_wsh[l] = std::min(wl, nmax);
          vgd_wsi[l] = std::min(cells(ain), nmax);
          vgd_wso[l] = std::min(cells(aout), nmax);
        }
      }
  #if MPI_PARALLEL_ENABLED
      if (vgd_nb1 > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &wneed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      }
  #endif
      vgd_w = std::min(wneed, nmax);
      if (pass == 0 && raw_ok && wneed > nmax) {
        vgd_raw = true;
        raw_halo();
        if (global_variable::my_rank == 0) {
          std::cout << "<rad_vet> vet_gd: the exact-list band needs "
                    << wneed << " cells, MeshBlocks of " << indcs.nx2 << " x "
                    << indcs.nx3 << " (theta x phi) hold " << nmax << ": the raw "
                    << "global band (no clamped reads)" << std::endl;
        }
        continue;
      }
      // several MeshBlocks along x1, x1 periodic slab (the exit shells' bands are read
      // by the next pass): every shell at the full band, the x1 ghosts too
      if (vgd_nb1 > 1 || (slb && vet_x1per)) {
        vgd_wsh.assign(c1, vgd_w);
        vgd_wsi = vgd_wsh;
        vgd_wso = vgd_wsh;
      }
      vgd_capped = (vgd_w < wneed);
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_vet> vet_gd: exact per-shell lateral halo, band " << vgd_w
                  << " cells (needed " << wneed << "), per-shell depth mean "
                  << std::accumulate(vgd_wsh.begin(), vgd_wsh.end(), 0.0)/indcs.nx1
                  << (vgd_capped ? ": MeshBlocks too small, the deepest near-tangent "
                                   "reads are clamped (counted)" : "") << std::endl;
      }
      break;
    }
  }
  const int c2w = indcs.nx2 + 2*vgd_w, c3w = indcs.nx3 + 2*vgd_w;
  if (vgd_sth) {ScThetaInit();}   // rays/sc_sph_globalband.cpp
  Kokkos::realloc(vgd_cs, nmb, 2, c3w, c2w, c1);
  Kokkos::realloc(vgd_wall, nmb, c3w, c2w);
  Kokkos::realloc(vgd_map, nmb, c3w, c2w, n);
  Kokkos::realloc(vgd_mr, nmb, c3w, c2w, n);
  if (vgd_wint) {
    Kokkos::realloc(vgd_m3, nmb, c3w, c2w, n, 3);
    Kokkos::realloc(vgd_w3, nmb, c3w, c2w, n, 3);
  }
  if (!vgd_csg) {ScHaloInit();}   // cs lists made by ScCsInit
  if (vgd_gband || vgd_raw) {ScGlobalBandInit();}   // sph2.inc: after the slots
  if (vgd_thmir) {
    // vet_gd_theta_wall = mirror: the theta neighbour across a wall is the block itself
    // (one MeshBlock in theta): corner slots are phi-edge slots
    if (vgd_csg) {
      ScFatal("vet_gd_theta_wall = mirror: no theta walls on the cubed sphere");
    }
    if (vgd_nb2 != 1) {
      ScFatal("vet_gd_theta_wall = mirror needs ONE MeshBlock along x2 (theta)");
    }
    if (vgd_ipc_on) {ScFatal("vet_gd_theta_wall = mirror: vet_gd_halo_ipc not ported");}
  }
  if (vgd_thvac > 0) {
    // vet_gd_theta_wall = vacuum / extrap: only the theta ghosts beyond a physical wall
    // change (ScTabCompute flags 3 / 4); vacuum allows any MeshBlock count in theta
    if (vgd_csg) {ScFatal("vet_gd_theta_wall = vacuum/extrap: not on the cubed sphere");}
    if (vgd_ipc_on) {
      ScFatal("vet_gd_theta_wall = vacuum/extrap: vet_gd_halo_ipc not ported");
    }
    if (vgd_thvac == 2 && vgd_nb2 != 1) {
      ScFatal("vet_gd_theta_wall = extrap needs ONE MeshBlock along x2 (theta)");
    }
  }
  {
    // vet_gd_band_ring: only where every band value the sweep reads comes from this
    // sweep's exact-list (cs: compact; vgd_raw: own) exchange and some block has a
    // remote side; without the cs compact lists the clean sc rays keep the ragged band
    if (vgd_csg && vgd_hcomp == 0) {
      ScFatal("vet_gd_halo_compact = 0 on the cubed sphere with the clean sc rays (no "
              "compact lists, no band ring: the full band); name vet_gd_halo_compact "
              "= 2");
    }
    const bool hxok = vgd_raw ? (vgd_hmpi && vgd_nb1 == 1)
                    : vgd_csg ? (vgd_cscomp && vgd_hmpi)
                              : (vgd_hx_on && vgd_hmpi && (vgd_hcomp > 0));
    bool anyside = false;
    auto hl_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_hloc);
    for (size_t q = 0; q < hl_h.extent(0); ++q) {anyside = anyside || (hl_h(q) < 0);}
    vgd_ring = vgd_ring_par && hxok && anyside;
    // the fused two-set sweep keeps the ragged band on sp; on the cubed sphere the ring
    // records both sets' compact lists (ScCsSweep's ringrec)
  }
  // combinations the fused sweep never covers (ScXffOk) keep two sweeps and no
  // second array (the clean ring, RingMoments, sweeps one set at a time)
  if (vgd_xff && ((cln_ring && !(cln_rfuse && cln_ring == 1 &&
                                          ibc_x1max != VET_IBC_REFLECT)) ||
                  (vgd_gband && !vgd_raw) || vgd_nb1 > 1)) {
    vgd_xff = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: not covered here ("
                << "xfade_m, ...): two sweeps per build" << std::endl;
    }
  }
  // the shells below the ray cut are filled just before a set's closure reads them
  // (ClosureFromSc): the fused sweep's two sets share their slabs.  From the cut of
  // the restart file, with a margin; a lower cut takes two sweeps (ScXffOk)
  vgd_rag_cut = (vgd_xff && vgd_nb1 == 1 && marena_on && !marena_sealed) ?
                std::max(cln_rcut - indcs.nx1/32, 0) : 0;
  ScRagAlloc(vgd_i);   // after the halo init: band sides follow vgd_hloc
  if (vgd_xff) {
    // the other set's intensities, same layout and band sides
    ScRagAlloc(vgd_i_xb);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: second-set intensities "
                << static_cast<double>(vgd_i_xb.d.extent(0))*sizeof(Real)/1.0e9
                << " GB on rank 0" << std::endl;
    }
  }
  ScRagThick();
  // vet_gd_halo_exact marks the reads in the shell-list sweep
  if (vgd_hx_on && (!vgd_hmpi || vgd_hcomp <= 0)) {
    vgd_hx_on = false;         // does not apply (or nothing to exchange)
  }
  vgd_alpha = -1.0;
  // table prefetch (own device instance for the upload)
  if (ScRotOn()) {
    vgd_tex = Kokkos::Experimental::partition_space(DevExeSpace(), 1)[0];
  }
  ScTables((vgd_slab == 1) ? ScSlabKey(RotBlock(pm->ncycle)) : RotAngle(pm->ncycle));
}

namespace {

// second, high-priority stream of the cubed sphere's early send (ScEs2Make)
#if defined(KOKKOS_ENABLE_HIP)
DevExeSpace ScPrioSpace() {
  int lo = 0, hi = 0;
  (void)hipDeviceGetStreamPriorityRange(&lo, &hi);
  hipStream_t s;
  (void)hipStreamCreateWithPriority(&s, hipStreamNonBlocking, hi);
  return DevExeSpace(s, Kokkos::Impl::ManageStream::yes);
}
#elif defined(KOKKOS_ENABLE_CUDA)
DevExeSpace ScPrioSpace() {
  int lo = 0, hi = 0;
  (void)cudaDeviceGetStreamPriorityRange(&lo, &hi);
  cudaStream_t s;
  (void)cudaStreamCreateWithPriority(&s, cudaStreamNonBlocking, hi);
  return DevExeSpace(s, Kokkos::Impl::ManageStream::yes);
}
#else
DevExeSpace ScPrioSpace() {return DevExeSpace();}
#endif

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScEs2Make
//! \brief the early send's second, high-priority stream vgd_es2 (made once)

void RadiationVet::ScEs2Make() {
  if (!vgd_es2_made) {
    vgd_es2 = ScPrioSpace();
    vgd_es2_made = true;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScHaloInit
//! \brief the 8 lateral neighbour slots of every local MeshBlock (same x1 level; theta
//! and phi periodic by wrapping the logical location) and the message buffers

void RadiationVet::ScHaloInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  int nb2 = 1, nb3 = 1;
  for (int g = 0; g < pm->nmb_total; ++g) {
    nb2 = std::max(nb2, pm->lloc_eachmb[g].lx2 + 1);
    nb3 = std::max(nb3, pm->lloc_eachmb[g].lx3 + 1);
  }
  std::vector<int> gof(vgd_nb1*nb2*nb3, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx1*nb3 + l.lx3)*nb2 + l.lx2] = g;
  }
  for (int t = 0; t < vgd_nb1*nb2*nb3; ++t) {
    if (gof[t] < 0) {ScFatal("incomplete MeshBlock grid (SMR/AMR not supported)");}
  }
  Kokkos::realloc(vgd_hloc, 8*nmb);
  auto hl_h = Kokkos::create_mirror_view(vgd_hloc);
  vgd_hrank.assign(8*nmb, me);
  vgd_hlid.assign(8*nmb, -1);
  vgd_hmpi = false;
  vgd_nb2 = nb2;
  vgd_nb3 = nb3;
  Kokkos::realloc(vgd_lxy, nmb, 2);
  {
    auto lx_h = Kokkos::create_mirror_view(vgd_lxy);
    for (int m = 0; m < nmb; ++m) {
      lx_h(m,0) = pm->lloc_eachmb[g0 + m].lx2;
      lx_h(m,1) = pm->lloc_eachmb[g0 + m].lx3;
    }
    Kokkos::deep_copy(vgd_lxy, lx_h);
  }
  // full sphere: a slot beyond a pole takes the pole partner (lx2, lx3 + dk + nb3/2),
  // rows reversed (vgd_hpol 1 north, 2 south, else 0)
  if (vgd_pol && (nb3 % 2 != 0)) {
    ScFatal("the full sphere needs an EVEN number of MeshBlocks in phi (pole partners)");
  }
  Kokkos::realloc(vgd_hpol, 8*nmb);
  auto hp_h = Kokkos::create_mirror_view(vgd_hpol);
  for (int m = 0; m < nmb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    for (int o = 0; o < 8; ++o) {
      const int oo = (o < 4) ? o : (o + 1);       // skip the centre of the 3x3
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      const int pol = ScPoleSlot(o, l.lx2);
      hp_h(8*m + o) = pol;
      const int l2 = (pol != 0) ? l.lx2 : (l.lx2 + dj + nb2) % nb2;
      const int l3 = (pol != 0) ? (l.lx3 + dk + nb3/2 + nb3) % nb3
                                : (l.lx3 + dk + nb3) % nb3;
      const int g = gof[(l.lx1*nb3 + l3)*nb2 + l2];
      const int rk = pm->rank_eachmb[g];
      if (rk == me) {
        hl_h(8*m + o) = g - g0;
      } else {
        hl_h(8*m + o) = -1;
        vgd_hrank[8*m + o] = rk;
        vgd_hlid[8*m + o] = g - pm->gids_eachrank[rk];
        vgd_hmpi = true;
      }
    }
  }
  Kokkos::deep_copy(vgd_hloc, hl_h);
  Kokkos::deep_copy(vgd_hpol, hp_h);
  // the largest region (edge slot: nx x W) times the most variables per message
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nreg = std::max(indcs.nx2, indcs.nx3)*vgd_w;
  // one shell of all directions (the multi-shell ln chi, ln S exchange is chunked to
  // fit, ScHaloCs)
  (void) c1;
  vgd_maxcnt = nreg*std::max(vgd_n, 2);
  if (vgd_hmpi) {
    Kokkos::realloc(vgd_sbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    Kokkos::realloc(vgd_rbuf, static_cast<size_t>(nmb)*8*vgd_maxcnt);
    // message plan: pieces sorted by (rank, receiver block, receiver slot) on both
    // sides; piece size kn*jn with kn, jn in {nx, ws}
    std::vector<std::array<int, 4>> ps, pr;
    for (int m = 0; m < nmb; ++m) {
      for (int o = 0; o < 8; ++o) {
        if (hl_h(8*m + o) >= 0) {continue;}
        // the receiver's slot pointing back here: the opposite one, or across a pole
        // the one with the same theta side and the opposite phi side
        int orc = 7 - o;
        if (hp_h(8*m + o) != 0) {
          const int oo = (o < 4) ? o : (o + 1);
          const int dk = oo/3 - 1, dj = oo%3 - 1;
          const int ob = 3*(-dk + 1) + (dj + 1);
          orc = (ob < 4) ? ob : (ob - 1);
        }
        ps.push_back({vgd_hrank[8*m + o], vgd_hlid[8*m + o], orc, 8*m + o});
        pr.push_back({vgd_hrank[8*m + o], m, o, 8*m + o});
      }
    }
    std::sort(ps.begin(), ps.end());
    std::sort(pr.begin(), pr.end());
    vgd_prk.clear();
    for (const auto &x : ps) {
      if (vgd_prk.empty() || vgd_prk.back() != x[0]) {vgd_prk.push_back(x[0]);}
    }
    const int np = static_cast<int>(vgd_prk.size());
    Kokkos::realloc(vgd_soff, vgd_w + 1, 8*nmb);
    Kokkos::realloc(vgd_roff, vgd_w + 1, 8*nmb);
    auto so_h = Kokkos::create_mirror_view(vgd_soff);
    auto ro_h = Kokkos::create_mirror_view(vgd_roff);
    vgd_pdsp.assign(vgd_w + 1, std::vector<int>(np, 0));
    vgd_pscnt = vgd_pdsp;
    vgd_prdsp = vgd_pdsp;
    vgd_prcnt = vgd_pdsp;
    auto sz = [&](const int o, const int ws) {
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      return ((dk == 0) ? indcs.nx3 : ws)*((dj == 0) ? indcs.nx2 : ws);
    };
    for (int ws = 0; ws <= vgd_w; ++ws) {
      for (int t = 0; t < 8*nmb; ++t) {so_h(ws,t) = -1; ro_h(ws,t) = -1;}
      for (int side = 0; side < 2; ++side) {
        const auto &L = (side == 0) ? ps : pr;
        int off = 0, p = -1, prev = -1;
        for (const auto &x : L) {
          if (x[0] != prev) {
            ++p;
            prev = x[0];
            ((side == 0) ? vgd_pdsp : vgd_prdsp)[ws][p] = off;
          }
          const int o = x[3] % 8;
          ((side == 0) ? so_h : ro_h)(ws, x[3]) = off;
          // piece size: the receiver's slot 7 - o has the shape of my slot o
          off += sz(o, ws);
          ((side == 0) ? vgd_pscnt : vgd_prcnt)[ws][p] += sz(o, ws);
        }
      }
    }
    Kokkos::deep_copy(vgd_soff, so_h);
    Kokkos::deep_copy(vgd_roff, ro_h);
    if (vgd_hcomp > 0) {
      // partner boundaries of both dense layouts per depth ws (in units of nv*ni)
      const int npp = static_cast<int>(vgd_prk.size());
      Kokkos::realloc(vgd_pbd, vgd_w + 1, 2*(npp + 1));
      auto pb_h = Kokkos::create_mirror_view(vgd_pbd);
      for (int ws = 0; ws <= vgd_w; ++ws) {
        for (int p = 0; p < npp; ++p) {
          pb_h(ws,p) = vgd_pdsp[ws][p];
          pb_h(ws,npp + 1 + p) = vgd_prdsp[ws][p];
        }
        pb_h(ws,npp) = (npp > 0) ? (vgd_pdsp[ws][npp-1] + vgd_pscnt[ws][npp-1]) : 0;
        pb_h(ws,2*npp + 1) = (npp > 0) ? (vgd_prdsp[ws][npp-1] + vgd_prcnt[ws][npp-1])
                                       : 0;
      }
      Kokkos::deep_copy(vgd_pbd, pb_h);
      for (int sl = 0; sl < 2; ++sl) {
        Kokkos::realloc(vgd_hpb[sl], 2*(npp + 1));
      }
      auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       pmy_pack->pcoord->x1v);
      auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                       pmy_pack->pcoord->xx1f);
      vgd_r1v.assign(x1v_h.extent(1), 0.0);
      vgd_r1f.assign(x1f_h.extent(1), 0.0);
      for (size_t q = 0; q < vgd_r1v.size(); ++q) {vgd_r1v[q] = x1v_h(0,q);}
      for (size_t q = 0; q < vgd_r1f.size(); ++q) {vgd_r1f[q] = x1f_h(0,q);}
      // exact per-shell sizes (ni = 1, nv = vgd_n) at the deepest band: sums of the
      // remote pieces (only the dense sbuf / rbuf carry the multi-shell exchange)
      const size_t ns = static_cast<size_t>(pb_h(vgd_w,npp))*vgd_n;
      const size_t nr = static_cast<size_t>(pb_h(vgd_w,2*npp + 1))*vgd_n;
      Kokkos::realloc(vgd_csb, std::max<size_t>(ns, 1));
      Kokkos::realloc(vgd_crb, std::max<size_t>(nr, 1));
      // mask slots: this shell's and the next one's
      for (int sl = 0; sl < 2; ++sl) {
        Kokkos::realloc(vgd_hfs[sl], ns + 1);
        Kokkos::realloc(vgd_hfr[sl], nr + 1);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScIpcInit
//! \brief vet_gd_halo_ipc: node-local CUDA IPC halo of the exact lists: raw receive
//! buffer (2 parities x half, cs the partner stride) and flags, handles swapped with
//! every partner and opened. Off unless every partner is on this node. Collective.

void RadiationVet::ScIpcInit(const std::vector<int> &prk, size_t cs, size_t half) {
  vgd_ipc_state = -1;
#if MPI_PARALLEL_ENABLED && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  const int me = global_variable::my_rank;
  const int np = static_cast<int>(prk.size());
  int ok = (np <= VGD_IPC_NS) ? 1 : 0;
  MPI_Comm node;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, me, MPI_INFO_NULL, &node);
  int nn;
  MPI_Comm_size(node, &nn);
  std::vector<int> nodeg(nn);
  MPI_Allgather(&me, 1, MPI_INT, nodeg.data(), 1, MPI_INT, node);
  MPI_Comm_free(&node);
  // vgd_comm is MPI_COMM_WORLD or a duplicate (vet_gd_async): the same ranks
  for (int p = 0; p < np; ++p) {
    if (std::find(nodeg.begin(), nodeg.end(), prk[p]) == nodeg.end()) {ok = 0;}
  }
  vgd_ipc_cs = cs;
  vgd_ipc_half = half;
  if (ok) {
    if (VETG(Malloc)(reinterpret_cast<void **>(&vgd_ipc_rbuf),
                     2*vgd_ipc_half*sizeof(Real))
        != VETG(Success)) {ok = 0;}
#if defined(KOKKOS_ENABLE_HIP)
    if (ok && hipExtMallocWithFlags(reinterpret_cast<void **>(&vgd_ipc_sig),
                                    2*VGD_IPC_NS*sizeof(long long),
                                    hipDeviceMallocUncached) != hipSuccess) {ok = 0;}
#else
    if (ok && cudaMalloc(reinterpret_cast<void **>(&vgd_ipc_sig),
                         2*VGD_IPC_NS*sizeof(long long)) != cudaSuccess) {ok = 0;}
#endif
    if (ok && VETG(Malloc)(reinterpret_cast<void **>(&vgd_ipc_err), sizeof(int))
        != VETG(Success)) {ok = 0;}
    if (ok) {
      (void)VETG(Memset)(vgd_ipc_sig, 0, 2*VGD_IPC_NS*sizeof(long long));
      (void)VETG(Memset)(vgd_ipc_err, 0, sizeof(int));
      (void)VETG(Memset)(vgd_ipc_rbuf, 0, 2*vgd_ipc_half*sizeof(Real));
      (void)VETG(DeviceSynchronize)();
    }
  }
  {
    int g = ok;
    MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    ok = g;
  }
  if (ok) {
    VETG(IpcMemHandle_t) hb, hs;
    (void)VETG(IpcGetMemHandle)(&hb, vgd_ipc_rbuf);
    (void)VETG(IpcGetMemHandle)(&hs, vgd_ipc_sig);
    struct Info {VETG(IpcMemHandle_t) b, s; long long half, cs; int slot;};
    std::vector<Info> mine(np), theirs(np);
    std::vector<MPI_Request> rq(2*np);
    for (int p = 0; p < np; ++p) {
      mine[p] = Info{hb, hs, static_cast<long long>(vgd_ipc_half),
                     static_cast<long long>(vgd_ipc_cs), p};
      MPI_Irecv(&theirs[p], sizeof(Info), MPI_BYTE, prk[p], 7010, vgd_comm, &rq[p]);
    }
    for (int p = 0; p < np; ++p) {
      MPI_Isend(&mine[p], sizeof(Info), MPI_BYTE, prk[p], 7010, vgd_comm,
                &rq[np + p]);
    }
    MPI_Waitall(2*np, rq.data(), MPI_STATUSES_IGNORE);
    vgd_ipc_prb = DualArray1D<unsigned long long>("vet_vgd_ipc_prb", std::max(np, 1));
    vgd_ipc_psig = DualArray1D<unsigned long long>("vet_vgd_ipc_psig", std::max(np, 1));
    vgd_ipc_phalf = DualArray1D<long long>("vet_vgd_ipc_phalf", std::max(np, 1));
    vgd_ipc_pcs = DualArray1D<long long>("vet_vgd_ipc_pcs", std::max(np, 1));
    vgd_ipc_pslot = DualArray1D<int>("vet_vgd_ipc_pslot", std::max(np, 1));
    vgd_ipc_open.assign(2*np, nullptr);
    for (int p = 0; p < np; ++p) {
      void *pb = nullptr, *ps = nullptr;
      if (VETG(IpcOpenMemHandle)(&pb, theirs[p].b, VETG(IpcMemLazyEnablePeerAccess))
          != VETG(Success)) {ok = 0;}
      if (VETG(IpcOpenMemHandle)(&ps, theirs[p].s, VETG(IpcMemLazyEnablePeerAccess))
          != VETG(Success)) {ok = 0;}
      vgd_ipc_open[2*p] = pb;
      vgd_ipc_open[2*p + 1] = ps;
      vgd_ipc_prb.h_view(p) = reinterpret_cast<unsigned long long>(pb);
      vgd_ipc_psig.h_view(p) = reinterpret_cast<unsigned long long>(ps);
      vgd_ipc_phalf.h_view(p) = theirs[p].half;
      vgd_ipc_pcs.h_view(p) = theirs[p].cs;
      vgd_ipc_pslot.h_view(p) = theirs[p].slot;
    }
    (void)VETG(GetLastError)();
    vgd_ipc_prb.modify_host(); vgd_ipc_prb.sync_device();
    vgd_ipc_psig.modify_host(); vgd_ipc_psig.sync_device();
    vgd_ipc_phalf.modify_host(); vgd_ipc_phalf.sync_device();
    vgd_ipc_pcs.modify_host(); vgd_ipc_pcs.sync_device();
    vgd_ipc_pslot.modify_host(); vgd_ipc_pslot.sync_device();
    int g = ok;
    MPI_Allreduce(&ok, &g, 1, MPI_INT, MPI_MIN, vgd_comm);
    ok = g;
    MPI_Barrier(vgd_comm);   // every peer's flags are zero before any exchange
  }
  vgd_ipc_state = ok ? 1 : -1;
  vgd_ipc_seq = 0;
  if (me == 0) {
    std::cout << "<rad_vet> vet_gd_halo_ipc: "
              << (ok ? "ON (node-local GPU IPC, GPU-driven per-shell halo of the exact "
                       "lists, no host wait)"
                     : "not possible (a partner off-node, or IPC failed): MPI path")
              << std::endl;
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScIpcCheck
//! \brief vet_gd_halo_ipc: after a (fenced) sweep, fatal if a spin timed out

void RadiationVet::ScIpcCheck() {
#if MPI_PARALLEL_ENABLED && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  if (vgd_ipc_state != 1) {return;}
  int h = 0;
  (void)VETG(Memcpy)(&h, vgd_ipc_err, sizeof(int), VETG(MemcpyDeviceToHost));
  if (h != 0) {
    std::cout << "### FATAL ERROR in rad_vet vet_gd_halo_ipc: a halo spin timed out "
                 "(rank "
              << global_variable::my_rank << ", exchange " << vgd_ipc_seq << ")"
              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScRagAlloc
//! \brief the ragged intensity array (ScRag, rad_vet.hpp): shell i gets band depth
//! vgd_wsh[i] (its data's deepest read), ghost shells none; zero-filled

void RadiationVet::ScRagAlloc(ScRag &a) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  a.nmb = pmy_pack->nmb_thispack;
  a.n = vgd_n;
  a.nx2 = indcs.nx2;
  a.nx3 = indcs.nx3;
  a.w = vgd_w;
  a.c1 = c1;
  Kokkos::realloc(a.off, a.nmb, c1);
  Kokkos::realloc(a.wi, c1);
  Kokkos::realloc(a.sb, a.nmb, 4);
  Kokkos::realloc(a.oos, 1);
  auto off_h = Kokkos::create_mirror_view(a.off);
  auto wi_h = Kokkos::create_mirror_view(a.wi);
  auto sb_h = Kokkos::create_mirror_view(a.sb);
  // a side carries the band when any of its three neighbours (edge, two corners) is
  // not on this rank (vgd_hloc < 0)
  auto hl_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_hloc);
  for (int m = 0; m < a.nmb; ++m) {
    int side[4] = {0, 0, 0, 0};
    for (int o = 0; o < 8; ++o) {
      if (hl_h(8*m + o) >= 0) {continue;}
      const int oo = (o < 4) ? o : (o + 1);
      const int dk = oo/3 - 1, dj = oo%3 - 1;
      if (dk < 0) {side[0] = 1;}
      if (dk > 0) {side[1] = 1;}
      if (dj < 0) {side[2] = 1;}
      if (dj > 0) {side[3] = 1;}
    }
    for (int q = 0; q < 4; ++q) {sb_h(m,q) = side[q];}
  }
  if (vgd_raw) {
    // raw band (every slot remote): a side carries it only where a band position has
    // no on-rank owner (vgd_own < 0); on-rank owners are read directly
    auto ow_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_own);
    const int w = vgd_w;
    for (int m = 0; m < a.nmb; ++m) {
      for (int q = 0; q < 4; ++q) {sb_h(m,q) = 0;}
      for (int kk = 0; kk < a.nx3 + 2*w; ++kk) {
        for (int jj = 0; jj < a.nx2 + 2*w; ++jj) {
          if (ow_h(m,kk,jj,0) >= 0) {continue;}
          sb_h(m,0) |= (kk < w);
          sb_h(m,1) |= (kk >= w + a.nx3);
          sb_h(m,2) |= (jj < w);
          sb_h(m,3) |= (jj >= w + a.nx2);
        }
      }
    }
  }
  for (int i = 0; i < c1; ++i) {
    int wl = 0;
    if ((i >= indcs.is && i <= indcs.ie) || vgd_nb1 > 1) {
      wl = std::max(vgd_wsh[i], std::max(vgd_wsi[i], vgd_wso[i]));
      wl = std::min(wl, vgd_w);
    }
    wi_h(i) = wl;
  }
  int64_t tot = 0, tth = 0;
  const int cg = (vgd_rag_cut > 0) ? indcs.is + vgd_rag_cut : 0;
  if (vgd_ring) {
    // vet_gd_band_ring: every shell compact (depth 0) at [0, ...), then vgd_rg_k wide
    // slots of one shell each at the deepest band; shell i's full depth in
    // vgd_rg_wfull (same for both arrays)
    vgd_rg_wfull.assign(c1, 0);
    for (int i = 0; i < c1; ++i) {vgd_rg_wfull[i] = wi_h(i);}
    const int64_t cap = static_cast<int64_t>(a.nx3)*a.nx2*a.n;
    for (int m = 0; m < a.nmb; ++m) {
      for (int i = 0; i < c1; ++i) {
        off_h(m,i) = (i < cg) ? -1 - tth : tot;   // (below cg: ScRagThick)
        ((i < cg) ? tth : tot) += cap;
      }
    }
    const int64_t ncomp = tot;
    vgd_rg_k = VGD_RING_K;
    vgd_rg_soff.assign(static_cast<size_t>(vgd_rg_k)*a.nmb, 0);
    for (int s = 0; s < vgd_rg_k; ++s) {
      for (int m = 0; m < a.nmb; ++m) {
        int64_t mx = 0;
        for (int i = 0; i < c1; ++i) {
          const int wl = wi_h(i);
          mx = std::max(mx, static_cast<int64_t>(a.n)*
                            (a.nx3 + (sb_h(m,0) + sb_h(m,1))*wl)*
                            (a.nx2 + (sb_h(m,2) + sb_h(m,3))*wl));
        }
        vgd_rg_soff[static_cast<size_t>(s)*a.nmb + m] = tot;
        tot += mx;
      }
    }
    for (int i = 0; i < c1; ++i) {wi_h(i) = 0;}
    vgd_rg_state.assign(c1, -1);
    vgd_rg_sown.assign(vgd_rg_k, -1);
    vgd_rg_so.assign(c1, 0);
    vgd_rg_sn.assign(c1, 0);
    vgd_rg_sk.assign(c1, 0);
    vgd_rg_se.assign(c1, -1);
    vgd_rg_pend.clear();
    Kokkos::realloc(a.offc, a.nmb, c1);
    Kokkos::realloc(vgd_rg_dsoff, vgd_rg_k, a.nmb);
    auto so_h = Kokkos::create_mirror_view(vgd_rg_dsoff);
    for (int m = 0; m < a.nmb; ++m) {
      for (int s = 0; s < vgd_rg_k; ++s) {
        so_h(s,m) = vgd_rg_soff[static_cast<size_t>(s)*a.nmb + m];
      }
    }
    Kokkos::deep_copy(a.offc, off_h);
    Kokkos::deep_copy(vgd_rg_dsoff, so_h);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd_band_ring: compact intensities "
                << ncomp*sizeof(Real)/1.0e9 << " GB + "
                << vgd_rg_k << " band slots "
                << (tot - ncomp)*sizeof(Real)/1.0e9 << " GB on rank 0" << std::endl;
    }
  } else {
  for (int m = 0; m < a.nmb; ++m) {
    for (int i = 0; i < c1; ++i) {
      const int wl = wi_h(i);
      off_h(m,i) = (i < cg) ? -1 - tth : tot;
      ((i < cg) ? tth : tot) += static_cast<int64_t>(a.n)*
             (a.nx3 + (sb_h(m,0) + sb_h(m,1))*wl)*(a.nx2 + (sb_h(m,2) + sb_h(m,3))*wl);
    }
  }
  }
  a.dummy = tot;
  a.nth = tth;
  Kokkos::deep_copy(a.off, off_h);
  Kokkos::deep_copy(a.wi, wi_h);
  Kokkos::deep_copy(a.sb, sb_h);
  const int64_t nd = tot + 1;
  if (!MemArenaRays(a, nd)) {Kokkos::realloc(a.d, nd);}
  Kokkos::deep_copy(a.d, 0.0);
  Kokkos::deep_copy(a.oos, 0.0);
  if (global_variable::my_rank == 0) {
    const double dense = static_cast<double>(a.nmb)*a.n*(a.nx3 + 2*vgd_w)*
                         (a.nx2 + 2*vgd_w)*c1;
    std::cout << "<rad_vet> vet_gd: ragged per-shell band intensity array "
              << nd*sizeof(Real)/1.0e9 << " GB on rank 0 (dense band "
              << dense*sizeof(Real)/1.0e9 << " GB)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScRagThick
//! \brief the slabs of the shells below vgd_rag_cut, one set of them for both fused
//! arrays, behind both in the arena: the arrays' offsets there (their own bases)

void RadiationVet::ScRagThick() {
  if (vgd_i.nth == 0) return;
  MemArenaGrow(marena_sw + vgd_i.nth, false);   // (both arrays hold zeros)
  Real *t = marena.data() + marena_sw;
  marena_sw += vgd_i.nth + ((64 - (vgd_i.nth % 64)) % 64);
  for (ScRag *a : {&vgd_i, &vgd_i_xb}) {
    const int64_t b = t - a->d.data();
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), a->off);
    for (size_t m = 0; m < h.extent(0); ++m) {
      for (size_t i = 0; i < h.extent(1); ++i) {
        if (h(m,i) < 0) h(m,i) = b + (-1 - h(m,i));
      }
    }
    Kokkos::deep_copy(a->off, h);
    if (vgd_ring) Kokkos::deep_copy(a->offc, h);
    a->d = DvceArray1D<Real>(a->d.data(), b + a->nth);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse: the shells below " << vgd_rag_cut
              << " (ray cut " << cln_rcut << ") share one slab, "
              << vgd_i.nth*sizeof(Real)/1.0e9 << " GB on rank 0" << std::endl;
  }
}

}  // namespace radvet
