//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_walls.cpp
//! \brief curved-mesh rays: x1 hand-over, mirrors, theta walls and poles.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_faces.cpp,
//!    rad_vet_closure_ring.cpp, rays/sc_setup.cpp, rays/sc_sph_globalband.cpp,
//!    rays/sc_sph_halo.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rays/sc_sph_halo.cpp

#include <math.h>

#include <algorithm>
#include <cmath>
#include <iostream>
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
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdThVacCs
//! \brief vet_gd_theta_wall = vacuum / extrap: ln chi, ln S of the theta ghost rows
//! beyond a physical wall (flags 3, 4) = the edge row's (zero gradient)

void RadiationVet::VetGdThVacCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int wb = vgd_w, nx2 = indcs.nx2;
  const int c3 = indcs.nx3 + 2*vgd_w;
  const int c2 = nx2 + 2*wb;
  const Real t0 = pmy_pack->pmesh->mesh_size.x2min;
  const Real t1 = pmy_pack->pmesh->mesh_size.x2max;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const VgdLG lg_ = vgd_lg;
  auto a_ = a;
  auto wl_ = vgd_wall;
  par_for("vet_vgd_thvac_cs", vgd_cur, 0, nmb1, 0, 1, 0, c3 - 1, 0, 2*wb - 1, i0, i1,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int s, const int i) {
    const int g = (s < wb) ? s : (s - wb);
    const int j = (s < wb) ? (wb - 1 - g) : (wb + nx2 + g);
    if (wl_(m,k,j) < 3) {return;}
    // the band row next to the wall (the block's edge row with one theta MeshBlock)
    int je = (s < wb) ? static_cast<int>(lg_.Jw(mbsize, m, wb, t0))
                      : (static_cast<int>(lg_.Jw(mbsize, m, wb, t1)) - 1);
    je = (je < 0) ? 0 : ((je > c2 - 1) ? (c2 - 1) : je);
    a_(m,c,k,j,i) = a_(m,c,k,je,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdThMirCs
//! \brief vet_gd_theta_wall = mirror: ln chi, ln S of the theta ghost rows = their
//! mirror rows (after the exchange filled the phi band)

void RadiationVet::VetGdThMirCs(DvceArray5D<Real> &a, const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int wb = vgd_w, nx2 = indcs.nx2;
  const int c3 = indcs.nx3 + 2*vgd_w;
  auto a_ = a;
  par_for("vet_vgd_thmir_cs", vgd_cur, 0, nmb1, 0, 1, 0, c3 - 1, 0, 2*wb - 1, i0, i1,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int s, const int i) {
    const int g = (s < wb) ? s : (s - wb);
    const int j = (s < wb) ? (wb - 1 - g) : (wb + nx2 + g);
    const int jm = VgdMirRow(g, s >= wb, wb, nx2);
    a_(m,c,k,j,i) = a_(m,c,k,jm,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdX1Hand
//! \brief several MeshBlocks along x1: the boundary shell of every level to the x1
//! ghost shell of the level above (up) or below, its band exchanged as a sweep shell

void RadiationVet::VetGdX1Hand(const bool up) {
  if (vgd_scut != 0) {VgdFatal("several MeshBlocks along x1 need vray_icut = 0 (the "
                               "sweep from the inner face)");}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w;
  const int nx2 = indcs.nx2, nk = indcs.nx3*nx2;
  const int isrc = up ? indcs.ie : indcs.is, idst = up ? (indcs.is - 1) : (indcs.ie + 1);
  auto vi_ = vgd_i;
  auto b_ = vgd_x1b;
  auto r_ = vgd_x1r;
  par_for("vet_vgd_x1pack", vgd_cur, 0, nmb1, 0, n - 1, 0, nk - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int t) {
    b_(m,d*nk + t) = vi_(m,d,w + t/nx2,w + t % nx2,isrc);
  });
  VetGdX1Pass(up);
  par_for("vet_vgd_x1unpack", vgd_cur, 0, nmb1, 0, n - 1, 0, nk - 1,
  KOKKOS_LAMBDA(const int m, const int d, const int t) {
    vi_(m,d,w + t/nx2,w + t % nx2,idst) = r_(m,d*nk + t);
  });
  VetGdHalo(vgd_i, n, idst, idst, w, true);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTopMirror
//! \brief reflecting outer x1 (clean scheme): the inflow of every inward direction at
//! the top face = the top cell's last-sweep intensity in the specular image
//! n - 2 (n . r) r (nearest set direction), scaled so the inflow returns the outflow's
//! radial flux; the top cell's isotropic E when there is no outflow yet

void RadiationVet::VetGdTopMirror() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, ie = indcs.ie;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, ks = indcs.ks, js = indcs.js;
  auto dir_ = vgd_dir;
  auto tg_ = vgd_trig;
  auto mp_ = vgd_rmap;
  if (mp_.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(vgd_rmap, nmb1 + 1, nx3, nx2, n);
    mp_ = vgd_rmap;
  }
  // every sweep: a rotating set changes the directions
  par_for("vet_vgd_rmap", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1, 0, n - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
    Real best = -2.0;
    for (int q = 0; q < n; ++q) {
      Real c = 0.0;
      for (int a = 0; a < 3; ++a) {c += dir_(q,a)*(dir_(d,a) - 2.0*mr*r[a]);}
      if (c > best + 1.0e-12) {best = c; mp_(m,k,j,d) = q;}
    }
  });
  auto vi_ = vgd_i;
  auto iw_ = iw;
  auto xf_ = vgd_x1f;
  // clean sc rays: the last build's top shell from its copy cln_top (kept in restart
  // files; clean_rot_ring: the set's slot), else the top cell's E. cln_top holds either
  // every direction per cell (older files) or n1 per packed slot
  auto fi_ = cln_top;
  const int ns = cln_ring ? 2 : 1, n1 = fi_.extent_int(4);
  const bool all = (fi_.extent_int(1) == ns*n);
  const int off = (cln_ring ? cln_rqm : 0)*(all ? n : (n + n1 - 1)/std::max(n1, 1));
  const bool cp = cln_scmem && fi_.extent_int(0) == nmb1 + 1 && fi_.extent_int(1) > 0;
  const bool e1 = cln_scmem && !cp;
  par_for("vet_vgd_mirror", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (!(xf_(m) & 4)) {return;}
    auto top = [&](const int q) {
      if (!cp) {return vi_(m,q,k+w,j+w,ie);}
      return all ? fi_(m,off+q,k+ks,j+js,ie) : fi_(m,off+q/n1,k+ks,j+js,q % n1);
    };
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    Real hp = 0.0, hm = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
      if (mr >= 0.0) {
        hp += dir_(d,3)*mr*top(d);
      } else {
        hm -= dir_(d,3)*mr*top(mp_(m,k,j,d));
      }
    }
    const Real e0 = iw_(m,VET_IW_EN,k+ks,j+js,ie);
    for (int d = 0; d < n; ++d) {   // outward ones are not read
      vi_(m,d,k+w,j+w,ie+1) = (hp > 0.0 && hm > 0.0 && !e1) ? hp/hm*top(mp_(m,k,j,d))
                                                            : e0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdTopCopy
//! \brief cln_scmem: the top shell of this build's sweep into its slot of cln_top
//! (blk mod 2 on clean_rot_ring), for the next build's VetGdTopMirror

void RadiationVet::VetGdTopCopy() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, ie = indcs.ie;
  const int ks = indcs.ks, js = indcs.js, n1 = cln_s.extent_int(4), nt = (n + n1 - 1)/n1;
  const int ns = cln_ring ? 2 : 1, off = (cln_ring ? cln_rqm : 0)*nt;
  if (cln_top.extent_int(0) != nmb1 + 1 || cln_top.extent_int(1) != ns*nt) {
    Kokkos::realloc(cln_top, nmb1 + 1, ns*nt, cln_s.extent_int(2), cln_s.extent_int(3),
                    n1);
  }
  auto t_ = cln_top;
  auto vi_ = vgd_i;
  par_for("vet_vgd_top", DevExeSpace(), 0, nmb1, 0, indcs.nx3 - 1, 0, indcs.nx2 - 1, 0,
          n - 1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    t_(m,off+d/n1,k+ks,j+js,d % n1) = vi_(m,d,k+w,j+w,ie);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdInnerMirror
//! \brief reflecting inner x1 (clean scheme), after the inward pass: per column of the
//! first shell, inward radial flux / outward flux of its specular images (nearest set
//! directions, robust to round-off): the scale of the mirrored inner-face inflow
//! (1 where either vanishes)

void RadiationVet::VetGdInnerMirror() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb1 = pmy_pack->nmb_thispack - 1, n = vgd_n, w = vgd_w, is = indcs.is;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  if (vgd_isc.extent_int(0) != nmb1 + 1) {Kokkos::realloc(vgd_isc, nmb1 + 1, nx3, nx2);}
  auto dir_ = vgd_dir;
  auto tg_ = vgd_trig;
  auto vi_ = vgd_i;
  auto xf_ = vgd_x1f;
  auto sc_ = vgd_isc;
  par_for("vet_vgd_imirror", DevExeSpace(), 0, nmb1, 0, nx3 - 1, 0, nx2 - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    sc_(m,k,j) = 1.0;
    if (!(xf_(m) & 1)) {return;}
    const Real r[3] = {tg_(m,k,j,0)*tg_(m,k,j,3), tg_(m,k,j,0)*tg_(m,k,j,2),
                       tg_(m,k,j,1)};
    Real hi = 0.0, ho = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real mr = dir_(d,0)*r[0] + dir_(d,1)*r[1] + dir_(d,2)*r[2];
      if (mr < 0.0) {
        hi -= dir_(d,3)*mr*vi_(m,d,k+w,j+w,is);
        continue;
      }
      int dm = 0;
      Real best = -2.0;
      for (int q = 0; q < n; ++q) {
        // the nearest inward set direction: outward ones are swept only after this
        Real c = 0.0;
        for (int a = 0; a < 3; ++a) {c += dir_(q,a)*(dir_(d,a) - 2.0*mr*r[a]);}
        const Real mq = dir_(q,0)*r[0] + dir_(q,1)*r[1] + dir_(q,2)*r[2];
        if (mq < 0.0 && c > best + 1.0e-12) {best = c; dm = q;}
      }
      ho += dir_(d,3)*mr*vi_(m,dm,k+w,j+w,is);
    }
    if (hi > 0.0 && ho > 0.0) {sc_(m,k,j) = hi/ho;}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdX1Pass
//! \brief vgd_x1b(m, :) of every MeshBlock to vgd_x1r of its x1 neighbour above (up)
//! or below (on-rank copy or one message per block)

void RadiationVet::VetGdX1Pass(const bool up) {
  Mesh *pm = pmy_pack->pmesh;
  const int nmb = pmy_pack->nmb_thispack, g0 = pmy_pack->gids;
  const int cnt = vgd_x1b.extent_int(1);
  Kokkos::fence();
  auto bh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_x1b);
  auto rh = Kokkos::create_mirror_view(vgd_x1r);
#if MPI_PARALLEL_ENABLED
  std::vector<MPI_Request> req;
#endif
  for (int m = 0; m < nmb; ++m) {
    const int gr = vgd_x1nb[2*m + (up ? 1 : 0)], gs = vgd_x1nb[2*m + (up ? 0 : 1)];
    if (gr >= 0 && pm->rank_eachmb[gr] == global_variable::my_rank) {
      for (int t = 0; t < cnt; ++t) {rh(gr - g0,t) = bh(m,t);}
    }
#if MPI_PARALLEL_ENABLED
    if (gr >= 0 && pm->rank_eachmb[gr] != global_variable::my_rank) {
      req.emplace_back();
      MPI_Isend(&bh(m,0), cnt, MPI_ATHENA_REAL, pm->rank_eachmb[gr], gr, vgd_comm,
                &req.back());
    }
    if (gs >= 0 && pm->rank_eachmb[gs] != global_variable::my_rank) {
      req.emplace_back();
      MPI_Irecv(&rh(m,0), cnt, MPI_ATHENA_REAL, pm->rank_eachmb[gs], g0 + m, vgd_comm,
                &req.back());
    }
#endif
  }
#if MPI_PARALLEL_ENABLED
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
#endif
  Kokkos::deep_copy(vgd_x1r, rh);
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::VgdPoleSlot
//! \brief pole slot of halo neighbour o of a block in theta row lx2: 1 N, 2 S, else 0

int RadiationVet::VgdPoleSlot(const int o, const int lx2) const {
  if (!vgd_pol) {return 0;}
  const int oo = (o < 4) ? o : (o + 1);
  const int dj = oo%3 - 1;
  if (dj < 0 && lx2 == 0) {return 1;}
  if (dj > 0 && lx2 == vgd_nb2 - 1) {return 2;}
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn double RadiationVet::VetGdThetaDmin
//! \brief narrowest lateral spacing (arc on the unit sphere) for the band depth:
//! min of dtheta and sin(theta_lim) dphi, theta_lim = wedge limit or vet_gd_band_colat

double RadiationVet::VetGdThetaDmin() const {
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
  return std::min(dth, smin*vgd_lg.Dx3H(mbs, 0));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdThetaInit
//! \brief node table vgd_thn(m, band row), after vgd_w is known

void RadiationVet::VetGdThetaInit() {
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
  // global row J (within one period / mirror of the mesh) -> theta
  auto nodeg = [&](const int J) -> double {
    if (J >= 0 && J < nx) {return node[J];}
    if (vgd_pol) {   // extended chart beyond a pole: the mirror row
      if (J < 0) {return -node[-1 - J];}
      return 2.0*M_PI - node[2*nx - 1 - J];
    }
    // periodic wedge: image of the other side's row
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
        VgdFatal("the theta node table is not increasing (block " +
                 std::to_string(m) + ", band row " + std::to_string(jb) +
                 "): the band is deeper than one period/mirror of the mesh");
      }
    }
  }
  Kokkos::deep_copy(vgd_thn, tn_h);
  // the pole slot table vgd_hpol is made by VetGdHaloInit (knows nb2)
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
    std::cout << "<rad_vet> vet_gd: theta node table ("
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

}  // namespace radvet
