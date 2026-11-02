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


//----------------------------------------------------------------------------------------
//! the polar RING buffers (full sphere): every phi block's polar rows (node colatitude <
//! vet_gd_band_colat, plus one row), gathered so that an upwind stencil inside the polar
//! zone reads its true phi on the whole ring (the band cannot reach there: the phi cells
//! of a polar ring are sin(theta) narrow).  Items: (block, pole) for every block at a
//! pole, in gid order; each rank packs its own items, MPI_Allgatherv, every rank unpacks.

void RadiationM1::VetGdRingInit() {
  vgd_rbon = false;
  if (!vgd_pol) {return;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nx = pm->mesh_indcs.nx2, nx3t = pm->mesh_indcs.nx3;
  const int nx2b = indcs.nx2, nx3b = indcs.nx3;
  std::vector<double> node(nx);
  for (int n = 0; n < nx; ++n) {
    node[n] = 0.5*(VgdThetaFace(pm, n) + VgdThetaFace(pm, n + 1));
  }
  const double tlim = vgd_bcolat*M_PI/180.0;
  int nz = 0;
  while (nz < nx/2 && node[nz] < tlim) {++nz;}
  nz = std::min(nz + 1, std::min(nx2b, nx/2));
  vgd_nzr = nz;
  vgd_nx3t = nx3t;
  vgd_nt2 = nx;
  Kokkos::realloc(vgd_rbth, 2, nz);
  auto th_h = Kokkos::create_mirror_view(vgd_rbth);
  for (int q = 0; q < nz; ++q) {th_h(0,q) = node[q]; th_h(1,q) = node[nx - 1 - q];}
  Kokkos::deep_copy(vgd_rbth, th_h);
  // the items in gid order (ranks own contiguous gid ranges)
  const int me = global_variable::my_rank, nr = global_variable::nranks;
  const int g0 = pmy_pack->gids;
  std::vector<std::array<int, 4>> it;
  vgd_rbcnt.assign(nr, 0);
  vgd_rbdsp.assign(nr, 0);
  const int U = nz*nx3b;
  for (int g = 0; g < pm->nmb_total; ++g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    const int rk = pm->rank_eachmb[g];
    for (int p = 0; p < 2; ++p) {
      if ((p == 0 && l.lx2 != 0) || (p == 1 && l.lx2 != vgd_nb2 - 1)) {continue;}
      const int off = static_cast<int>(it.size())*U;
      if (vgd_rbcnt[rk] == 0) {vgd_rbdsp[rk] = off;}
      vgd_rbcnt[rk] += U;
      it.push_back({(rk == me) ? (g - g0) : -1, p, l.lx3, off});
    }
  }
  vgd_rbnit = static_cast<int>(it.size());
  vgd_rbnloc = 0;
  Kokkos::realloc(vgd_rbit, std::max(vgd_rbnit, 1), 4);
  auto e_h = Kokkos::create_mirror_view(vgd_rbit);
  int first = -1;
  for (int e = 0; e < vgd_rbnit; ++e) {
    for (int c = 0; c < 4; ++c) {e_h(e,c) = it[e][c];}
    if (it[e][0] >= 0) {if (first < 0) {first = e;} ++vgd_rbnloc;}
  }
  Kokkos::deep_copy(vgd_rbit, e_h);
  vgd_rbf0 = std::max(first, 0);
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int n = vgd_n;
  const int nsh = c1;
  const size_t ni = static_cast<size_t>(c1)*2*nz*nx3t*n;
  const size_t nc = static_cast<size_t>(c1)*2*2*nz*nx3t;
  Kokkos::realloc(vgd_rbi, ni);
  Kokkos::realloc(vgd_rbc, nc);
  Kokkos::deep_copy(vgd_rbi, 0.0);
  Kokkos::deep_copy(vgd_rbc, 0.0);
  const size_t per = static_cast<size_t>(U)*std::max(n, 2*nsh);
  Kokkos::realloc(vgd_rbs, std::max<size_t>(per*vgd_rbnloc, 1));
  Kokkos::realloc(vgd_rbr, std::max<size_t>(per*vgd_rbnit, 1));
  vgd_rbon = true;
  if (me == 0) {
    std::cout << "<rad_m1> vet_gd (accsp-1102): polar RING buffers: " << nz
              << " rows per pole (node colatitude < " << vgd_bcolat << " deg + 1 row, to "
              << node[nz - 1]*180.0/M_PI << " deg), " << vgd_rbnit << " (block, pole) items; "
              << (ni + nc)*sizeof(Real)/1.0e9 << " GB per rank; stencils inside the polar "
              << "zone read the whole ring (no band, no clamp)" << std::endl;
  }
}

void RadiationM1::VetGdRingGatherI(const int i) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nz = vgd_nzr, nx3b = indcs.nx3, nx2b = indcs.nx2, w = vgd_w, n = vgd_n;
  const int nx3t = vgd_nx3t, U = nz*nx3b;
  auto e_ = vgd_rbit;
  auto sb_ = vgd_rbs;
  auto vi_ = vgd_i;
  const int nloc = vgd_rbnloc, f0 = vgd_rbf0;   // the local items are contiguous
  if (nloc > 0) {
    Kokkos::parallel_for("m1_vgd_rb_packi", Kokkos::RangePolicy<>(vgd_cur, 0, nloc*U*n),
    KOKKOS_LAMBDA(const int t) {
      const int d = t % n;
      const int r1 = t/n;
      const int k = r1 % nx3b;
      const int r2 = r1/nx3b;
      const int q = r2 % nz;
      const int el = r2/nz;
      const int m = e_(f0 + el,0), p = e_(f0 + el,1);
      const int row = (p == 0) ? (w + q) : (w + nx2b - 1 - q);
      sb_(t) = vi_(m,d,w + k,row,i);
    });
  }
  auto rb_ = sb_;
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    vgd_cur.fence();
    std::vector<int> cnt(vgd_rbcnt.size()), dsp(vgd_rbdsp.size());
    for (size_t r = 0; r < cnt.size(); ++r) {cnt[r] = vgd_rbcnt[r]*n; dsp[r] = vgd_rbdsp[r]*n;}
    MPI_Allgatherv(vgd_rbs.data(), nloc*U*n, MPI_ATHENA_REAL, vgd_rbr.data(), cnt.data(),
                   dsp.data(), MPI_ATHENA_REAL, MPI_COMM_WORLD);
    rb_ = vgd_rbr;
  }
#endif
  auto ri_ = vgd_rbi;
  const bool loc = !(global_variable::nranks > 1);
  Kokkos::parallel_for("m1_vgd_rb_unpi",
                       Kokkos::RangePolicy<>(vgd_cur, 0, vgd_rbnit*U*n),
  KOKKOS_LAMBDA(const int t) {
    const int d = t % n;
    const int r1 = t/n;
    const int k = r1 % nx3b;
    const int r2 = r1/nx3b;
    const int q = r2 % nz;
    const int e = r2/nz;
    const int p = e_(e,1), lx3 = e_(e,2);
    if (loc && e_(e,0) < 0) {return;}
    const int src = loc ? (((e - f0)*U + q*nx3b + k)*n + d) : ((e_(e,3) + q*nx3b + k)*n + d);
    const size_t at = ((static_cast<size_t>(i)*2*nz + p*nz + q)*nx3t + lx3*nx3b + k)*n + d;
    ri_(at) = rb_(src);
  });
}

void RadiationM1::VetGdRingGatherC(const int i0, const int i1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nz = vgd_nzr, nx3b = indcs.nx3, nx2b = indcs.nx2, w = vgd_w;
  const int nx3t = vgd_nx3t, U = nz*nx3b, ns = i1 - i0 + 1;
  auto e_ = vgd_rbit;
  auto sb_ = vgd_rbs;
  auto cs_ = vgd_cs;
  const int nloc = vgd_rbnloc, f0 = vgd_rbf0;
  const int B = U*2*ns;   // values per item
  if (nloc > 0) {
    Kokkos::parallel_for("m1_vgd_rb_packc", Kokkos::RangePolicy<>(vgd_cur, 0, nloc*B),
    KOKKOS_LAMBDA(const int t) {
      const int k = t % nx3b;
      int r = t/nx3b;
      const int q = r % nz;
      r /= nz;
      const int c = r % 2;
      r /= 2;
      const int is = r % ns;
      const int el = r/ns;
      const int m = e_(f0 + el,0), p = e_(f0 + el,1);
      const int row = (p == 0) ? (w + q) : (w + nx2b - 1 - q);
      sb_(t) = cs_(m,c,w + k,row,i0 + is);
    });
  }
  auto rb_ = sb_;
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    vgd_cur.fence();
    std::vector<int> cnt(vgd_rbcnt.size()), dsp(vgd_rbdsp.size());
    for (size_t r = 0; r < cnt.size(); ++r) {
      cnt[r] = vgd_rbcnt[r]/U*B;
      dsp[r] = vgd_rbdsp[r]/U*B;
    }
    MPI_Allgatherv(vgd_rbs.data(), nloc*B, MPI_ATHENA_REAL, vgd_rbr.data(), cnt.data(),
                   dsp.data(), MPI_ATHENA_REAL, MPI_COMM_WORLD);
    rb_ = vgd_rbr;
  }
#endif
  auto rc_ = vgd_rbc;
  const bool loc = !(global_variable::nranks > 1);
  Kokkos::parallel_for("m1_vgd_rb_unpc", Kokkos::RangePolicy<>(vgd_cur, 0, vgd_rbnit*B),
  KOKKOS_LAMBDA(const int t) {
    const int k = t % nx3b;
    int r = t/nx3b;
    const int q = r % nz;
    r /= nz;
    const int c = r % 2;
    r /= 2;
    const int is = r % ns;
    const int e = r/ns;
    if (loc && e_(e,0) < 0) {return;}
    const int p = e_(e,1), lx3 = e_(e,2);
    const int src = loc ? ((e - f0)*B + ((is*2 + c)*nz + q)*nx3b + k)
                        : (e_(e,3)/U*B + ((is*2 + c)*nz + q)*nx3b + k);
    const size_t at = ((static_cast<size_t>(i0 + is)*2 + c)*2*nz + p*nz + q)*nx3t
                      + lx3*nx3b + k;
    rc_(at) = rb_(src);
  });
}

} // namespace radm1
