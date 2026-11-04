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
  if (!vgd_pol || !vgd_rbkey) {return;}
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nx = pm->mesh_indcs.nx2, nx3t = pm->mesh_indcs.nx3;
  const int nx2b = indcs.nx2, nx3b = indcs.nx3;
  std::vector<double> node(nx);
  for (int n = 0; n < nx; ++n) {
    node[n] = 0.5*(VgdThetaFace(pm, n) + VgdThetaFace(pm, n + 1));
  }
  const double tlim = ((vgd_rbcolat > 0.0) ? vgd_rbcolat : vgd_bcolat)*M_PI/180.0;
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
              << " rows per pole (node colatitude < " << tlim*180.0/M_PI << " deg + 1 row, to "
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

//----------------------------------------------------------------------------------------
//! accpolar-1103 <rad_m1>/vet_gd_band_global: the band lists keyed on the GLOBAL cell
//! (the cubed sphere's exact band copies, VetGdCsInit, on the full sphere's index space).
//! Band position (kk, jj) of block (lx2, lx3) is the global cell J = lx2 nx2 + jj - w, K =
//! lx3 nx3 + kk - w, continued beyond a pole in the extended chart (J -> -1 - J north,
//! 2 nt2 - 1 - J south, K -> K + nx3t/2: the pole partner's reversed rows, as the slot
//! halo), K periodic.  Every position at depth d >= 1 (Chebyshev distance from the
//! interior) takes its owner's interior value: on-rank copies and one message per
//! partner rank, both sorted by depth so that the per-shell halo of depth ws is a prefix
//! (VetGdCsHalo with ws).  Every slot is marked remote (vgd_hloc = -1): the sweep reads
//! every band position from the band (no direct neighbour reads, which reach one block).

void RadiationM1::VetGdGlobalBandInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  const int ngb = pm->nmb_total;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, w = vgd_w;
  const int c2w = nx2 + 2*w, c3w = nx3 + 2*w;
  const int nt2 = pm->mesh_indcs.nx2, nt3 = pm->mesh_indcs.nx3;
  const int nb2 = vgd_nb2, nb3 = vgd_nb3;
  if (pm->multilevel) {VgdFatal("vet_gd_band_global: SMR/AMR not supported");}
  if (nb2*nx2 != nt2 || nb3*nx3 != nt3 || (nt3 % 2) != 0) {
    VgdFatal("vet_gd_band_global: the MeshBlocks must tile the mesh, nx3 even");
  }
  std::vector<int> gof(nb2*nb3, -1);
  for (int g = 0; g < ngb; ++g) {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    gof[l.lx3*nb2 + l.lx2] = g;
  }
  // the owner (gid, band k, band j) of band position (kk, jj) of block g; false: none
  auto owner = [&](const int g, const int kk, const int jj, int &go, int &ko,
                   int &jo) -> bool {
    const LogicalLocation &l = pm->lloc_eachmb[g];
    int J = l.lx2*nx2 + jj - w, K = l.lx3*nx3 + kk - w;
    if (J < 0) {
      J = -1 - J;
      K += nt3/2;
    } else if (J >= nt2) {
      J = 2*nt2 - 1 - J;
      K += nt3/2;
    }
    if (J < 0 || J >= nt2) {return false;}
    K = ((K % nt3) + nt3) % nt3;
    const int b2 = J/nx2, b3 = K/nx3;
    go = gof[b3*nb2 + b2];
    ko = K - b3*nx3 + w;
    jo = J - b2*nx2 + w;
    return true;
  };
  auto depth = [&](const int kk, const int jj) -> int {
    const int dk = std::max(std::max(w - kk, kk - (w + nx3 - 1)), 0);
    const int dj = std::max(std::max(w - jj, jj - (w + nx2 - 1)), 0);
    return std::max(dk, dj);
  };
  // buckets by depth (1..w): on-rank (6 ints), per partner rank send / recv (3 ints)
  std::vector<std::vector<int>> lb(w + 1);
  std::map<int, std::vector<std::vector<int>>> sb, rb;
  for (int gr = 0; gr < ngb; ++gr) {
    const int rr = pm->rank_eachmb[gr];
    for (int kk = 0; kk < c3w; ++kk) {
      for (int jj = 0; jj < c2w; ++jj) {
        const int d = depth(kk, jj);
        if (d == 0) {continue;}
        int go, ko, jo;
        if (!owner(gr, kk, jj, go, ko, jo)) {
          VgdFatal("vet_gd_band_global: a band position beyond the sphere (band deeper "
                   "than the theta rows)");
        }
        const int ro = pm->rank_eachmb[go];
        if (rr != me && ro != me) {continue;}
        if (rr == me && ro == me) {
          lb[d].insert(lb[d].end(), {gr - g0, kk, jj, go - g0, ko, jo, d});
        } else if (rr == me) {
          auto &v = rb[ro];
          if (v.empty()) {v.resize(w + 1);}
          v[d].insert(v[d].end(), {gr - g0, kk, jj});
        } else {
          auto &v = sb[rr];
          if (v.empty()) {v.resize(w + 1);}
          v[d].insert(v[d].end(), {go - pm->gids_eachrank[me], ko, jo});
        }
      }
    }
  }
  // the on-rank list in depth order, its prefix counts
  std::vector<int> loc;
  vgd_gdnl.assign(w + 1, 0);
  for (int d = 1; d <= w; ++d) {
    for (size_t q = 0; q < lb[d].size(); q += 7) {
      loc.insert(loc.end(), lb[d].begin() + q, lb[d].begin() + q + 6);
    }
    vgd_gdnl[d] = static_cast<int>(loc.size()/6);
  }
  vgd_csnloc = static_cast<int>(loc.size()/6);
  Kokkos::realloc(vgd_csloc, std::max(vgd_csnloc, 1), 6);
  {
    auto l_h = Kokkos::create_mirror_view(vgd_csloc);
    for (int e = 0; e < vgd_csnloc; ++e) {
      for (int c = 0; c < 6; ++c) {l_h(e,c) = loc[6*e + c];}
    }
    Kokkos::deep_copy(vgd_csloc, l_h);
  }
  // the partners (sorted), their segments in depth order
  vgd_csprk.clear();
  for (const auto &x : sb) {vgd_csprk.push_back(x.first);}
  for (const auto &x : rb) {
    if (sb.find(x.first) == sb.end()) {vgd_csprk.push_back(x.first);}
  }
  std::sort(vgd_csprk.begin(), vgd_csprk.end());
  const int np = static_cast<int>(vgd_csprk.size());
  vgd_csps.assign(np, 0);
  vgd_cspn.assign(np, 0);
  vgd_cspr.assign(np, 0);
  vgd_cspm.assign(np, 0);
  vgd_gdsn.assign(w + 1, std::vector<int>(np, 0));
  vgd_gdrn.assign(w + 1, std::vector<int>(np, 0));
  std::vector<int> sl, rl, sd, rd;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    vgd_csps[p] = static_cast<int>(sl.size()/3);
    vgd_cspr[p] = static_cast<int>(rl.size()/3);
    for (int d = 1; d <= w; ++d) {
      if (sb.count(rk)) {
        const auto &v = sb[rk][d];
        sl.insert(sl.end(), v.begin(), v.end());
        sd.insert(sd.end(), v.size()/3, d);
      }
      if (rb.count(rk)) {
        const auto &v = rb[rk][d];
        rl.insert(rl.end(), v.begin(), v.end());
        rd.insert(rd.end(), v.size()/3, d);
      }
      vgd_gdsn[d][p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
      vgd_gdrn[d][p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
    }
    vgd_cspn[p] = static_cast<int>(sl.size()/3) - vgd_csps[p];
    vgd_cspm[p] = static_cast<int>(rl.size()/3) - vgd_cspr[p];
  }
  vgd_csnsnd = static_cast<int>(sl.size()/3);
  vgd_csnrcv = static_cast<int>(rl.size()/3);
  Kokkos::realloc(vgd_cssl, std::max(vgd_csnsnd, 1), 3);
  Kokkos::realloc(vgd_csrl, std::max(vgd_csnrcv, 1), 3);
  Kokkos::realloc(vgd_gdps, std::max(vgd_csnsnd, 1));
  Kokkos::realloc(vgd_gdpr, std::max(vgd_csnrcv, 1));
  {
    auto s_h = Kokkos::create_mirror_view(vgd_cssl);
    auto r_h = Kokkos::create_mirror_view(vgd_csrl);
    auto ds_h = Kokkos::create_mirror_view(vgd_gdps);
    auto dr_h = Kokkos::create_mirror_view(vgd_gdpr);
    for (int e = 0; e < vgd_csnsnd; ++e) {
      for (int c = 0; c < 3; ++c) {s_h(e,c) = sl[3*e + c];}
      ds_h(e) = sd[e];
    }
    for (int e = 0; e < vgd_csnrcv; ++e) {
      for (int c = 0; c < 3; ++c) {r_h(e,c) = rl[3*e + c];}
      dr_h(e) = rd[e];
    }
    Kokkos::deep_copy(vgd_cssl, s_h);
    Kokkos::deep_copy(vgd_csrl, r_h);
    Kokkos::deep_copy(vgd_gdps, ds_h);
    Kokkos::deep_copy(vgd_gdpr, dr_h);
  }
  // message buffers: one shell of every direction, or ln chi and ln S of some shells
  const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
  const int nvb = std::max(vgd_n, 2*8);
  // device buffers (CUDA-aware MPI, as the sp slot halo; VetGdCsHalo with ws >= 0);
  // the pinned host buffers of the cubed sphere get one element
  Kokkos::realloc(vgd_sbuf, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_rbuf, static_cast<size_t>(nent)*nvb);
  Kokkos::realloc(vgd_cssb, 1);
  Kokkos::realloc(vgd_csrb, 1);
  vgd_hmpi = (np > 0);
  // the band columns' phi from the GLOBAL index (P0 + (K + 0.5) D, K unwrapped): the
  // per-block x3min + (k + 0.5) dx3 differs in the last bit between layouts, and a
  // tangent cell's branch test n . r_hat >= 0 may flip with it
  {
    const double p0 = pm->mesh_size.x3min, dp = (pm->mesh_size.x3max - p0)/nt3;
    vgd_phc_h.assign(static_cast<size_t>(nmb)*c3w, 0.0);
    Kokkos::realloc(vgd_phc, nmb, c3w);
    auto ph_h = Kokkos::create_mirror_view(vgd_phc);
    for (int m = 0; m < nmb; ++m) {
      const int lx3 = pm->lloc_eachmb[g0 + m].lx3;
      for (int kb = 0; kb < c3w; ++kb) {
        const int K = lx3*nx3 + kb - w;
        const double v = p0 + (K + 0.5)*dp;
        vgd_phc_h[static_cast<size_t>(m)*c3w + kb] = v;
        ph_h(m,kb) = v;
      }
    }
    Kokkos::deep_copy(vgd_phc, ph_h);
  }
  // every slot reads the band (no direct on-rank neighbour reads)
  Kokkos::deep_copy(vgd_hloc, -1);
  for (int t = 0; t < 8*nmb; ++t) {vgd_hrank[t] = me;}
  if (me == 0) {
    std::cout << "<rad_m1> vet_gd_band_global (accpolar-1103): band " << w
              << " cells keyed on the global cell (exact copies from any block; rank 0: "
              << "on-rank " << vgd_csnloc << ", sent " << vgd_csnsnd << ", received "
              << vgd_csnrcv << " positions per shell at full depth, partners " << np
              << "); stencil clamped to the cell's own window (layout independent)"
              << std::endl;
  }
}

} // namespace radm1
