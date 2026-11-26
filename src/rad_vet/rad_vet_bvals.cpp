//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_bvals.cpp
//! \brief ghost zones of the solve's arrays (deferred and batched exchanges), global
//! sums.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_faces.cpp,
//!    rad_vet_closure_rays.cpp, rad_vet_closure_ring.cpp, rad_vet_fct.cpp,
//!    rad_vet_flux.cpp, rad_vet_flux_bounds.cpp, rad_vet_flux_schur.cpp,
//!    rad_vet_flux_seam.cpp, rad_vet_gas.cpp, rad_vet_geometry.cpp, rad_vet_ledger.cpp,
//!    rad_vet_newton.cpp, rad_vet_oldvector.cpp, rad_vet_rhs.cpp, rad_vet_solver.cpp,
//!    rad_vet_steady.cpp, rad_vet_tasks.cpp, rad_vet_update.cpp

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_args.hpp"
#include "rad_vet/rad_vet_reduce.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanGhostFill
//! \brief fill the ghosts of components [c0, c0+nc) by the cell-centred exchange
//! (cubed-sphere seams without resample), then (phys) copy the edge cell on physical
//! boundaries.  part 1 starts, part 2 finishes (inner cells may run between).  Deferred
//! fills (CleanGhostDefer) join the next exchange where they fit.  D12, D13, D23 are odd
//! across a reflecting wall normal to one of their axes and D_r,theta, D_r,phi across
//! the pole, so d_b(D E) puts no shear traction on the wall

void RadiationVet::CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                                  const bool phys, const int part, const int dep,
                                  const std::array<uint32_t, 4> fm) {
  CleanGhostFillT(a, c0, nc, phys, part, dep, fm);
}

void RadiationVet::CleanGhostDefer(DvceArray5D<Real> a, const int c0, const int nc,
                                   const int dep) {
  int np = nc;
  for (auto &s : cln_hpend) np += s.nc;
  if (np > cln_hmax) CleanGhostFlush();
  cln_hpend.push_back({a, c0, nc, true, dep});
}

// (a range lies in one of the two arrays: CLW_NS)
void RadiationVet::CleanGhostFill(const CleanW &a, const int c0, const int nc,
                                  const int part, const int dep,
                                  const std::array<uint32_t, 4> fm) {
  CleanGhostFillT((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc, true, part,
                  dep, fm);
}

void RadiationVet::CleanGhostDefer(const CleanW &a, const int c0, const int nc,
                                   const int dep) {
  CleanGhostDefer((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc, dep);
}

void RadiationVet::CleanGhostFlush() {
  if (cln_hpend.empty()) return;
  std::vector<CleanGSeg<Real>> sg;
  sg.swap(cln_hpend);
  CleanGhostExchange(sg, 0);
}

template <typename T>
void RadiationVet::CleanGhostFillT(DvceArray5D<T> a, const int c0, const int nc,
                                   const bool phys, const int part, const int dep,
                                   const std::array<uint32_t, 4> fm) {
  if (nc > cln_hmax) {   // (the exchange moves at most cln_hmax variables at once)
    const int h = cln_hmax;
    // (a full mask stays full for every chunk)
    auto sh = [&](const uint32_t m) {return (m == ~0u || h >= 32) ? m : m >> h;};
    CleanGhostFillT(a, c0, h, phys, 0, dep, fm);
    CleanGhostFillT(a, c0 + h, nc - h, phys, 0, dep, {sh(fm[0]), sh(fm[1]), sh(fm[2]),
                                                       sh(fm[3])});
    return;
  }
  std::vector<CleanGSeg<T>> sg;
  if constexpr (std::is_same<T, Real>::value) {
    int np = nc;
    for (auto &s : cln_hpend) np += s.nc;
    // (part 2: the segments part 1 left in cln_hpend, nothing deferred in between)
    if (part == 2 || np <= cln_hmax) sg.swap(cln_hpend);
    if (part == 2) return CleanGhostExchange(sg, 2);
  }
  CleanGhostFlush();   // (those that do not fit: on their own first)
  sg.push_back({a, c0, nc, phys, dep, fm});
  CleanGhostExchange(sg, part);
  if constexpr (std::is_same<T, Real>::value) {
    if (part == 1) cln_hpend.swap(sg);
  }
}

template <typename T>
void RadiationVet::CleanGhostExchange(const std::vector<CleanGSeg<T>> &sg,
                                      const int part) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = sg[0].a.extent_int(4), n2 = sg[0].a.extent_int(3);
  const int n3 = sg[0].a.extent_int(2);
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  if (!cln_hpb) {
    cln_hpb.reset(new MeshBoundaryValuesCC(pmy_pack, nullptr, false));
    cln_hpb->InitializeBuffers(cln_hmax);
    cln_hpb->SetVectorPairs(cln_hmax, {});
    cln_hpb->cs_noresample = true;
    // (16 teams per buffer instead of one per MeshBlock face; bitwise)
    cln_hpb->cc_split = 16;
  }
  if (cln_hdl.state == 0) CleanHaloInit(n3, n2, n1);
  const bool dl = (cln_hdl.state == 1);
  if (dl && (n1 != cln_hdl.n1 || n2 != cln_hdl.n2 || n3 != cln_hdl.n3)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CleanGhostExchange: an array of another shape than the "
              << "lists of clean_halo_direct" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (dl && part != 2) {   // (the depth: the deepest its segments need)
    int d = 0;
    for (auto &s : sg) d = (s.dep <= 0) ? ng : std::max(d, s.dep);
    cln_hdl.dcur = std::min(std::max(d, 1), ng);
    CleanHaloPost(sg);
  }
  if (dl && part == 1) return;
  if (dl) CleanHaloFinish(sg);
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  // one MeshBlock and one Real array: exchange on the array itself, no copies
  bool own = false;
  if constexpr (std::is_same<T, Real>::value) own = (nmb1 == 0 && sg.size() == 1);
  // (clean_halo_direct, doubles: the cube-vertex corners on the arrays themselves)
  const bool dc = dl && cs && std::is_same<T, Real>::value;
  // (exchanges run one at a time: their scratch shares one allocation, grown on demand)
  if (!own && (dl ? cs && !dc : true) && cln_hbuf.extent_int(1) < nc) {
    Kokkos::realloc(cln_hbuf, nmb1 + 1, nc, n3, n2, n1);
    cln_halo.clear();
  }
  CleanHalo &hl = cln_halo[nc];
  if (hl.ca.extent_int(0) == 0) {
    hl.a = DvceArray5D<Real>(cln_hbuf.data(), nmb1 + 1, nc, n3, n2, n1);
    Kokkos::realloc(hl.ca, nmb1 + 1, nc, 1, 1, 1);
  }
  auto h = hl.a;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  MeshBoundaryValuesCC *pb = cln_hpb.get();
  DvceArray5D<Real> ha = hl.a;
  if constexpr (std::is_same<T, Real>::value) {
    const size_t off = static_cast<size_t>(sg[0].c0)*n3*n2*n1;
    if (own) ha = DvceArray5D<Real>(sg[0].a.data() + off, 1, nc, n3, n2, n1);
  }
  // (only the shell outside the inner box is copied; clean_halo_direct: only on the
  // cubed sphere, for the cube-vertex corners)
  if (!own && !dc && (dl ? cs : part != 2)) {
    const CleanShell sc(n3, n2, n1, ng, is, ie, js, je, ks, ke, md, td);
    int o = 0;
    for (auto &s : sg) {
      auto a = s.a;
      const int c0 = s.c0;
      if (sc.Empty()) {
        par_for("cln_hcopy", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, n3-1, 0, n2-1, 0, n1-1,
        KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
          h(m,o+n,k,j,i) = a(m,c0+n,k,j,i);
        });
      } else {
        par_for("cln_hcopy", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, sc.ncell-1,
        KOKKOS_LAMBDA(const int m, const int n, const int q) {
          int k, j, i;
          CleanShellCell(sc, q, k, j, i);
          h(m,o+n,k,j,i) = a(m,c0+n,k,j,i);
        });
      }
      o += s.nc;
    }
  }
  if (dc) {
    if constexpr (std::is_same<T, Real>::value) {
      for (auto &s : sg) {
        DvceArray5D<Real> a = s.a;
        pb->FillPanelCornersCC(a, false, 0, s.c0, s.nc);
      }
    }
  } else if (dl && cs) {
    pb->FillPanelCornersCC(ha);
  } else if (!dl) {
    if (part != 2) {
      while (pb->InitRecv(nc) == TaskStatus::incomplete) {}
      while (pb->PackAndSendCC(ha, hl.ca) == TaskStatus::incomplete) {}
      if (part == 1) return;
    }
    while (pb->RecvAndUnpackCC(ha, hl.ca) == TaskStatus::incomplete) {}
    while (pb->ClearSend() == TaskStatus::incomplete) {}
    while (pb->ClearRecv() == TaskStatus::incomplete) {}
  }
  const CleanShell sb(n3, n2, n1, 0, is, ie, js, je, ks, ke, md, td);
  int o = 0;
  for (auto &s : sg) {
    if (own || (dl && !cs) || dc) break;
    auto a = s.a;
    const int c0 = s.c0;
    par_for("cln_hback", DevExeSpace(), 0, nmb1, 0, s.nc-1, 0, sb.ncell-1,
    KOKKOS_LAMBDA(const int m, const int n, const int q) {
      int k, j, i;
      CleanShellCell(sb, q, k, j, i);
      a(m,c0+n,k,j,i) = h(m,o+n,k,j,i);
    });
    o += s.nc;
  }
  for (auto &s : sg) {
    if (s.phys) CleanGhostPhys(s.a, s.c0, s.nc);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloInit
//! \brief clean_halo_direct: build the ghost exchange's cell lists by sending cell ids
//! through cln_hpb: on-rank pairs make the copy kernel, off-rank ghosts are grouped by
//! rank (pack lists in the receiver's order, no headers).  Self-check: a test field must
//! come out bitwise equal through both paths, else the run stops

void RadiationVet::CleanHaloInit(const int n3, const int n2, const int n1) {
  CleanHaloList &q = cln_hdl;
  auto *pm = pmy_pack->pmesh;
  const int me = global_variable::my_rank;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const bool md = pm->multi_d, td = pm->three_d;
  const int nmb = pmy_pack->nmb_thispack, nrk = global_variable::nranks;
  int mx = 0;
  for (int r = 0; r < nrk; ++r) mx = std::max(mx, pm->nmb_eachrank[r]);
  const int n21 = n2*n1, n321 = n3*n21;
  const int64_t nid = static_cast<int64_t>(mx)*n321;
  q.n1 = n1; q.n2 = n2; q.n3 = n3;
  const int ng = indcs.ng;
  q.ng = ng;
  DvceArray5D<Real> x("cln_hdl_x", nmb, 1, n3, n2, n1),
                    y("cln_hdl_y", nmb, 1, n3, n2, n1);
  const std::vector<CleanGSeg<Real>> sx = {{x, 0, 1, false}}, sy = {{y, 0, 1, false}};
  // (pass 0: the ids; pass 1: other values, the check)
  for (int pass = 0; pass < 2; ++pass) {
    par_for("cln_hd_id", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real v = static_cast<Real>((me*nid + m*n321) + (k*n2 + j)*n1 + i);
      x(m,0,k,j,i) = (pass == 0) ? v : 2.0 + Kokkos::sin(1.0e-3*v);
      y(m,0,k,j,i) = x(m,0,k,j,i);
    });
    q.state = -1;
    CleanGhostExchange(sx, 0);
    if (pass == 1) {
      q.state = 1;
      CleanGhostExchange(sy, 0);
      break;
    }
    auto hx = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x);
    std::vector<int> ld, ls, lfm;
    std::vector<std::vector<int>> rg(nrk), rs(nrk), rf(nrk);   // off-rank ghosts, their
                                                                // sources and families
    // (family of a ghost copying id w: the axis along which its neighbour one layer in
    // differs by one cell; 3 at edges and corners)
    auto fam = [&](const int m, const int k, const int j, const int i, const int64_t w) {
      const int ox = (i < is || i > ie), oy = md && (j < js || j > je);
      const int oz = td && (k < ks || k > ke);
      if (ox + oy + oz != 1 || ng < 2) return 3;
      const int a = oy ? 1 : (oz ? 2 : 0), c = (a == 0) ? i : ((a == 1) ? j : k);
      const int lo = (a == 0) ? is : ((a == 1) ? js : ks);
      const int d = (c < lo) ? lo - c : c - ((a == 0) ? ie : ((a == 1) ? je : ke));
      const int t = ((c < lo) == (d > 1)) ? 1 : -1;
      const int64_t v = static_cast<int64_t>(hx(m,0,k+t*(a == 2),j+t*(a == 1),
                                                i+t*(a == 0)));
      const int64_t dw = (v > w) ? v - w : w - v;
      return (dw == 1) ? 0 : ((dw == n1) ? 1 : ((dw == n21) ? 2 : 3));
    };
    for (int m = 0; m < nmb; ++m) {
      for (int k = 0; k < n3; ++k) {
        for (int j = 0; j < n2; ++j) {
          for (int i = 0; i < n1; ++i) {
            const int g = (k*n2 + j)*n1 + i;
            const Real v = hx(m,0,k,j,i);
            const int64_t w = static_cast<int64_t>(v);
            if (v == static_cast<Real>(me*nid + m*n321 + g) || v < 0.0 ||
                static_cast<Real>(w) != v) continue;   // (not filled, or not a copy)
            const int r = static_cast<int>(w/nid);
            const int ms = static_cast<int>((w - r*nid)/n321);
            const int f = static_cast<int>(w - r*nid - ms*n321);
            const int kk = f/n21, jj = (f - kk*n21)/n1, ii = f - kk*n21 - jj*n1;
            if (r >= nrk || ms >= pm->nmb_eachrank[r] || ii < is || ii > ie ||
                (md && (jj < js || jj > je)) || (td && (kk < ks || kk > ke))) continue;
            if (r == me) {
              ld.push_back(m*n321 + g);
              ls.push_back(ms*n321 + f);
              lfm.push_back(fam(m, k, j, i, w));
            } else {
              rg[r].push_back(m*n321 + g);
              rs[r].push_back(ms*n321 + f);
              rf[r].push_back(fam(m, k, j, i, w));
            }
          }
        }
      }
    }
    // (each list in order of depth, the ghost layer of the receiving cell)
    auto dep = [&](const int f) {
      const int r = f%n321, k = r/n21, j = (r - k*n21)/n1, i = r - k*n21 - j*n1;
      const int d = std::max({is - i, i - ie, js - j, j - je, ks - k, k - ke});
      return std::min(std::max(d, 1), ng);
    };
    std::vector<int> nrd(static_cast<size_t>(nrk)*ng, 0), lrd, ldd;
    auto srt = [&](std::vector<int> &g, std::vector<int> &c, std::vector<int> &fv, int *n,
                   std::vector<int> &dv) {
      std::vector<int> o(g.size()), g2, c2, f2;
      for (size_t a = 0; a < o.size(); ++a) o[a] = static_cast<int>(a);
      std::stable_sort(o.begin(), o.end(),
                       [&](int a, int b) {return dep(g[a]) < dep(g[b]);});
      for (int a : o) {
        g2.push_back(g[a]);
        c2.push_back(c[a]);
        f2.push_back(fv[a]);
        dv.push_back(dep(g[a]));
        for (int d = dv.back(); d <= ng; ++d) ++n[d-1];
      }
      g.swap(g2);
      c.swap(c2);
      fv.swap(f2);
    };
    q.nld.assign(ng, 0);
    srt(ld, ls, lfm, q.nld.data(), ldd);
    for (int r = 0; r < nrk; ++r) srt(rg[r], rs[r], rf[r], nrd.data() + r*ng, lrd);
    std::vector<int> sc, rc;
#if MPI_PARALLEL_ENABLED
    MPI_Comm_dup(MPI_COMM_WORLD, &q.comm);
    std::vector<int> nr(nrk), ns(nrk);
    for (int r = 0; r < nrk; ++r) nr[r] = static_cast<int>(rs[r].size());
    MPI_Alltoall(nr.data(), 1, MPI_INT, ns.data(), 1, MPI_INT, q.comm);
    std::vector<std::vector<int>> ss(nrk), sf(nrk);
    std::vector<MPI_Request> rq;
    for (int r = 0; r < nrk; ++r) {
      if (nr[r] == 0 && ns[r] == 0) continue;
      q.rank.push_back(r);
      ss[r].resize(ns[r]);
      sf[r].resize(ns[r]);
      rq.resize(rq.size() + 4);
      MPI_Irecv(ss[r].data(), ns[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 4]);
      MPI_Isend(rs[r].data(), nr[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 3]);
      MPI_Irecv(sf[r].data(), ns[r], MPI_INT, r, 4243, q.comm, &rq[rq.size() - 2]);
      MPI_Isend(rf[r].data(), nr[r], MPI_INT, r, 4243, q.comm, &rq[rq.size() - 1]);
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    std::vector<int> nsd(nrd.size());
    MPI_Alltoall(nrd.data(), ng, MPI_INT, nsd.data(), ng, MPI_INT, q.comm);
    for (int r : q.rank) {
      q.rlend.insert(q.rlend.end(), nrd.begin() + r*ng, nrd.begin() + (r + 1)*ng);
      q.slend.insert(q.slend.end(), nsd.begin() + r*ng, nsd.begin() + (r + 1)*ng);
      q.soff.push_back(static_cast<int>(sc.size()));
      q.slen.push_back(ns[r]);
      sc.insert(sc.end(), ss[r].begin(), ss[r].end());
      q.hsf.insert(q.hsf.end(), sf[r].begin(), sf[r].end());
      q.roff.push_back(static_cast<int>(rc.size()));
      q.rlen.push_back(nr[r]);
      rc.insert(rc.end(), rg[r].begin(), rg[r].end());
      q.hrf.insert(q.hrf.end(), rf[r].begin(), rf[r].end());
    }
#endif
    auto mk = [](const char *nm, const std::vector<int> &v) {
      DvceArray1D<int> a(std::string(nm), std::max(v.size(), static_cast<size_t>(1)));
      Kokkos::deep_copy(Kokkos::subview(a, std::make_pair(static_cast<size_t>(0),
                                                          v.size())),
                        HostArray1D<int>(const_cast<int *>(v.data()), v.size()));
      return a;
    };
    q.nl = static_cast<int>(ld.size());
    q.ns = static_cast<int>(sc.size());
    q.nr = static_cast<int>(rc.size());
    q.ld = mk("cln_hdl_ld", ld);
    q.ls = mk("cln_hdl_ls", ls);
    q.sc = mk("cln_hdl_sc", sc);
    q.rc = mk("cln_hdl_rc", rc);
    q.rd = mk("cln_hdl_rd", lrd);   // (rg in rank order: rc's order)
    q.ldp = mk("cln_hdl_ldp", ldd);
    q.lf = mk("cln_hdl_lf", lfm);
    q.sf = mk("cln_hdl_sf", q.hsf);
    std::vector<int> sdv(sc.size(), ng);   // (each message's entries in order of depth)
    for (size_t s = 0; s < q.soff.size(); ++s) {
      for (int e = 0, d = 1; e < q.slen[s]; ++e) {
        while (d < ng && e >= q.slend[s*ng + d - 1]) ++d;
        sdv[q.soff[s] + e] = d;
      }
    }
    q.sdp = mk("cln_hdl_sdp", sdv);
    q.rf = mk("cln_hdl_rf", q.hrf);
    q.sbuf = decltype(q.sbuf)("cln_hdl_sbuf", std::max(q.ns, 1)*cln_hmax);
    q.rbuf = decltype(q.rbuf)("cln_hdl_rbuf", 2*std::max(q.nr, 1)*cln_hmax);
  }
  int bad = 0;
  Kokkos::parallel_reduce("cln_hd_chk", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          static_cast<int>(x.size())),
  KOKKOS_LAMBDA(const int p, int &b) {
    b += (x.data()[p] != y.data()[p]) ? 1 : 0;
  }, bad);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (bad != 0) {
    if (me == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "CleanHaloInit: the self-check of clean_halo_direct "
                << "failed: its cell lists give " << bad << " cells other than the "
                << "MeshBoundaryValuesCC exchange (array " << n3 << " x " << n2 << " x "
                << n1 << ", " << global_variable::nranks << " ranks)" << std::endl;
    }
    std::exit(EXIT_FAILURE);
  }
  q.state = 1;
  // (the lists replace cln_hpb's buffers from here on)
  Kokkos::fence();
  for (int b = 0; b < 56; ++b) {
    cln_hpb->sendbuf[b].vars = cln_hpb->recvbuf[b].vars = DvceArray2D<Real>();
    cln_hpb->sendbuf[b].flux = cln_hpb->recvbuf[b].flux = DvceArray2D<Real>();
  }
  if (me == 0) {
    std::cout << "<rad_vet> clean_halo_direct: ON, rank 0 copies " << q.nl
              << " ghosts, sends " << q.ns << " cells, receives " << q.nr << " in "
              << q.rank.size() << " messages; on-rank ghosts by depth";
    for (int d = 0; d < ng; ++d) std::cout << " " << q.nld[d];
    std::cout << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::CleanHaloKey, CleanHaloMask
//! \brief clean_halo_direct: the masks of an exchange (CleanGSeg::fm per family, then
//! nc) and their packed buffer layout

template <typename T>
std::array<uint32_t, 5> RadiationVet::CleanHaloKey(const std::vector<CleanGSeg<T>> &sg) {
  std::array<uint32_t, 5> k = {0, 0, 0, 0, 0};
  for (auto &s : sg) {
    const uint32_t a = (s.nc < 32) ? (1u << s.nc) - 1 : ~0u;
    for (int f = 0; f < 4; ++f) k[f] |= (s.fm[f] & a) << k[4];
    if (((s.fm[0] | s.fm[1] | s.fm[2] | s.fm[3]) & a) != a) {   // (no silent skip)
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "CleanHaloKey: a ghost fill whose masks leave a "
                << "component out of every face family" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    k[4] += s.nc;
  }
  return k;
}

RadiationVet::CleanHaloList::Mask &RadiationVet::CleanHaloMask(
    const std::array<uint32_t, 5> &key) {
  CleanHaloList &q = cln_hdl;
  auto &h = q.mask[key];
  if (h.so.extent(0) > 0) return h;
  const int nc = static_cast<int>(key[4]), nseg = static_cast<int>(q.rank.size());
  auto lay = [&](const std::vector<int> &fv, const std::vector<int> &off,
                 const std::vector<int> &len, const std::vector<int> &lend,
                 std::vector<int> &lm, const char *nm) {
    std::vector<int> o(std::max(fv.size(), static_cast<size_t>(1)), 0), c;
    lm.assign(static_cast<size_t>(nseg)*q.ng, 0);
    for (int s = 0; s < nseg; ++s) {
      c.assign(len[s] + 1, 0);
      for (int e = 0; e < len[s]; ++e) {
        o[off[s] + e] = off[s]*nc + c[e];
        c[e+1] = c[e] + static_cast<int>(std::bitset<32>(key[fv[off[s] + e]]).count());
      }
      for (int d = 0; d < q.ng; ++d) lm[s*q.ng + d] = c[lend[s*q.ng + d]];
    }
    DvceArray1D<int> a(std::string(nm), o.size());
    Kokkos::deep_copy(a, HostArray1D<int>(o.data(), o.size()));
    return a;
  };
  h.so = lay(q.hsf, q.soff, q.slen, q.slend, h.slm, "cln_hdl_so");
  h.ro = lay(q.hrf, q.roff, q.rlen, q.rlend, h.rlm, "cln_hdl_ro");
  return h;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloPost
//! \brief clean_halo_direct, first half: post receives (alternate buffer halves), pack,
//! fence, on-rank copy, one send per neighbour rank.  Entry p, component n at p nc + n

template <typename T>
void RadiationVet::CleanHaloPost(const std::vector<CleanGSeg<T>> &sg) {
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  q.par = 1 - q.par;
  // (fm: a ghost takes component b if bit b of its family's mask is set, msk)
  const std::array<uint32_t, 5> ky = CleanHaloKey(sg);
  const uint32_t k0 = ky[0], k1 = ky[1], k2 = ky[2], k3 = ky[3];
#if MPI_PARALLEL_ENABLED
  const uint32_t all = (nc < 32) ? (1u << nc) - 1 : ~0u;
  const bool msk = (k0 & k1 & k2 & k3) != all;
  CleanHaloList::Mask *mh = msk ? &CleanHaloMask(ky) : nullptr;
  const int nseg = static_cast<int>(q.rank.size());
  q.rq.assign(nseg, MPI_REQUEST_NULL);
  q.sq.assign(nseg, MPI_REQUEST_NULL);
  Real *rb = q.rbuf.data() + static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  for (int s = 0; s < nseg; ++s) {
    // (the entries of depth <= dcur)
    const int nrs = msk ? mh->rlm[s*q.ng + q.dcur - 1] : q.rlend[s*q.ng + q.dcur - 1]*nc;
    if (nrs > 0) {
      MPI_Irecv(rb + static_cast<size_t>(q.roff[s])*nc, nrs, MPI_ATHENA_REAL,
                q.rank[s], 4245, q.comm, &q.rq[s]);
    }
  }
  auto sb = q.sbuf;
  auto so = msk ? mh->so : q.sc, sf = q.sf, sdp = q.sdp;
  const int dsd = q.dcur;   // (entries deeper than the exchange are not sent: not packed)
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.sc;
    const int c0 = s.c0;
    if (q.ns > 0) {
      par_for("cln_hd_pack", DevExeSpace(), 0, q.ns-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        if (sdp(p) > dsd) return;
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        int b = p*nc + o + n;
        if (msk) {
          const int fp = sf(p);
          const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
          if (!((w >> (o + n)) & 1u)) return;
          b = so(p);
          for (int t = 0; t < o + n; ++t) b += (w >> t) & 1u;
        }
        sb(b) = a(m,c0+n,k,j,r - k*n21 - j*n1);
      });
    }
    o += s.nc;
  }
  if (q.ns > 0 || q.nr > 0) DevExeSpace().fence();
#endif
  // (on-rank copy runs while the sends are posted; ghosts deeper than dcur kept)
  const int dc = q.dcur;
  int oc = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto d = q.ld, c = q.ls, lf = q.lf;
    const int c0 = s.c0, o = oc;
    oc += s.nc;
    if (q.nl == 0) break;
    par_for("cln_hd_copy", DevExeSpace(), 0, q.nld[dc-1]-1, 0, s.nc-1,
    KOKKOS_LAMBDA(const int p, const int n) {
      const int f = d(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
      const int fp = lf(p);
      const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
      if (!((w >> (o + n)) & 1u)) return;
      const int g = c(p), mg = g/n321, t = g - mg*n321, kg = t/n21, jg = (t - kg*n21)/n1;
      a(m,c0+n,k,j,r - k*n21 - j*n1) = a(mg,c0+n,kg,jg,t - kg*n21 - jg*n1);
    });
  }
#if MPI_PARALLEL_ENABLED
  for (int s = 0; s < nseg; ++s) {
    const int nss = msk ? mh->slm[s*q.ng + q.dcur - 1] : q.slend[s*q.ng + q.dcur - 1]*nc;
    if (nss > 0) {
      MPI_Isend(sb.data() + static_cast<size_t>(q.soff[s])*nc, nss,
                MPI_ATHENA_REAL, q.rank[s], 4245, q.comm, &q.sq[s]);
    }
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanHaloFinish
//! \brief clean_halo_direct, second half: wait, unpack, wait for the sends

template <typename T>
void RadiationVet::CleanHaloFinish(const std::vector<CleanGSeg<T>> &sg) {
#if MPI_PARALLEL_ENABLED
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  const int nseg = static_cast<int>(q.rank.size());
  MPI_Waitall(nseg, q.rq.data(), MPI_STATUSES_IGNORE);
  const std::array<uint32_t, 5> ky = CleanHaloKey(sg);
  const uint32_t all = (nc < 32) ? (1u << nc) - 1 : ~0u;
  const bool msk = (ky[0] & ky[1] & ky[2] & ky[3]) != all;
  const uint32_t k0 = ky[0], k1 = ky[1], k2 = ky[2], k3 = ky[3];
  auto ro = msk ? CleanHaloMask(ky).ro : q.rc, rf = q.rf;
  auto rb = q.rbuf;
  const size_t b = static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  const int dc = q.dcur;
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.rc, dp = q.rd;
    const int c0 = s.c0;
    if (q.nr > 0) {
      par_for("cln_hd_unpack", DevExeSpace(), 0, q.nr-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        int u = p*nc + o + n;
        const int fp = rf(p);
        const uint32_t w = (fp == 0) ? k0 : ((fp == 1) ? k1 : ((fp == 2) ? k2 : k3));
        if (dp(p) > dc || !((w >> (o + n)) & 1u)) return;   // (not sent)
        if (msk) {
          u = ro(p);
          for (int t = 0; t < o + n; ++t) u += (w >> t) & 1u;
        }
        a(m,c0+n,k,j,r - k*n21 - j*n1) = static_cast<T>(rb(b + u));
      });
    }
    o += s.nc;
  }
  MPI_Waitall(nseg, q.sq.data(), MPI_STATUSES_IGNORE);
#else
  (void)sg;
#endif
}

template <typename T>
void RadiationVet::CleanGhostPhys(DvceArray5D<T> a, const int c0, const int nc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ng = indcs.ng;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = a.extent_int(4), n2 = a.extent_int(3), n3 = a.extent_int(2);
  auto bc = pmy_pack->pmb->mb_bcs.d_view;
  auto pf = [] KOKKOS_FUNCTION (const BoundaryFlag f) {
    return f != BoundaryFlag::block && f != BoundaryFlag::periodic &&
           f != BoundaryFlag::panel;
  };
  // (the cubed sphere's D is Cartesian: a radial wall mirrors it about n, not done here)
  const bool ts = (static_cast<const void *>(a.data()) ==
                   static_cast<const void *>(cln_s.data())) &&
                  !pmy_pack->pmesh->use_cubed_sphere;
  const int ib[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  int rf = 0;   // (bit s: side s reflects)
  for (int s = 0; s < 6; ++s) rf |= (ib[s] == VET_IBC_REFLECT) ? (1 << s) : 0;
  // (x1min without vet_bc_reflect = mirror: the rays do not mirror there, so the D_12
  // ghost stays even (no spurious shear); D_13 stays odd (no wall torque: L_z exact))
  const bool rm = vet_bc_mirror;
  auto sg = [=] KOKKOS_FUNCTION (const int c, const int s, const int m) {
    const int p = c - CLN_D11 - 3;   // (D12, D13, D23: axes 01, 02, 12; p != 2 - d)
    if (!ts || p < 0 || p > 2) return 1.0;
    if (bc(m,s) == BoundaryFlag::polar) return (p < 2) ? -1.0 : 1.0;
    if (s == 0 && !rm && p == 0) return 1.0;
    return (((rf >> s) & 1) && p != 2 - s/2) ? -1.0 : 1.0;
  };
  // (no launch for a direction without a physical side on this rank)
  auto bh = pmy_pack->pmb->mb_bcs.h_view;
  auto phys = [&](const int d) {
    for (int m = 0; m <= nmb1; ++m)
      if (pf(bh(m,2*d)) || pf(bh(m,2*d+1))) return true;
    return false;
  };
  // (x1, x2, x3 in turn: a corner takes the product of its sides' parities)
  if (phys(0)) {
    par_for("cln_ghost1", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, ng-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
      if (pf(bc(m,BoundaryFace::inner_x1))) {
        a(m,c0+n,k,j,is-1-g) = sg(c0+n, 0, m)*a(m,c0+n,k,j,is);
      }
      if (pf(bc(m,BoundaryFace::outer_x1))) {
        a(m,c0+n,k,j,ie+1+g) = sg(c0+n, 1, m)*a(m,c0+n,k,j,ie);
      }
    });
  }
  if (pmy_pack->pmesh->multi_d && phys(1)) {
    par_for("cln_ghost2", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, ng-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int g, const int i) {
      if (pf(bc(m,BoundaryFace::inner_x2))) {
        a(m,c0+n,k,js-1-g,i) = sg(c0+n, 2, m)*a(m,c0+n,k,js,i);
      }
      if (pf(bc(m,BoundaryFace::outer_x2))) {
        a(m,c0+n,k,je+1+g,i) = sg(c0+n, 3, m)*a(m,c0+n,k,je,i);
      }
    });
  }
  if (pmy_pack->pmesh->three_d && phys(2)) {
    par_for("cln_ghost3", DevExeSpace(), 0, nmb1, 0, nc-1, 0, ng-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int n, const int g, const int j, const int i) {
      if (pf(bc(m,BoundaryFace::inner_x3))) {
        a(m,c0+n,ks-1-g,j,i) = sg(c0+n, 4, m)*a(m,c0+n,ks,j,i);
      }
      if (pf(bc(m,BoundaryFace::outer_x3))) {
        a(m,c0+n,ke+1+g,j,i) = sg(c0+n, 5, m)*a(m,c0+n,ke,j,i);
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanReduce
//! \brief maximum (max) or sum of v over all ranks

Real RadiationVet::CleanReduce(Real v, const bool max) {
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_ATHENA_REAL, max ? MPI_MAX : MPI_SUM,
                MPI_COMM_WORLD);
#endif
  return v;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::CleanSumCells
//! \brief the sum of a over the active cells (CleanSums)

Real RadiationVet::CleanSumCells(const DvceArray4D<Real> &a) {
  auto a_ = a;
  return CleanSums(pmy_pack, "cln_sumc", 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &s1) {
    s0 += a_(m,k,j,i);
  })[0];
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> RadiationVet::CleanDots
//! \brief the n <= 2 dot products of w(ab[2q]) and w(ab[2q+1]) over nc components, one
//! pass and one reduction

std::array<Real, 2> RadiationVet::CleanDots(CleanW w, const int nc, const int n,
                                            const std::array<int, 4> ab) {
  const int a0 = ab[0], a1 = ab[1], a2 = ab[2], a3 = ab[3];   // (std::array is host-only)
  return CleanSums(pmy_pack, "cln_dot", n,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &s1) {
    for (int c = 0; c < nc; ++c) {
      s0 += w(m,a0+c,k,j,i)*w(m,a1+c,k,j,i);
      if (n > 1) s1 += w(m,a2+c,k,j,i)*w(m,a3+c,k,j,i);
    }
  });
}

}  // namespace radvet
