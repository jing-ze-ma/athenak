//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_solver.cpp
//! \brief Linear algebra of the clean implicit radiation solve: face coefficients, the
//! operator A E = E (1 + h c kE (1-r)) + h div F_lin(E), an x1-line tridiagonal
//! preconditioner (or the old production multigrid on it) and right-preconditioned
//! BiCGStab.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn CleanArgs RadiationM1::CleanMakeArgs
//! \brief collects the arrays and constants the device functions need

CleanArgs RadiationM1::CleanMakeArgs() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &size = pmy_pack->pmb->mb_size;
  CleanArgs x;
  x.cs = cln_s;
  x.w = cln_w;
  x.fc[0] = cln_fc1;
  x.fc[1] = cln_fc2;
  x.fc[2] = cln_fc3;
  x.beam = ibm_f;
  x.bstress = ibm_s;
  x.u = cln_x;
  x.h = dt_sub;
  x.c = c_light;
  x.idx[0] = 1.0/size.h_view(0).dx1;
  x.idx[1] = 1.0/size.h_view(0).dx2;
  x.idx[2] = 1.0/size.h_view(0).dx3;
  const int bc[6] = {ibc_x1min, ibc_x1max, ibc_x2min, ibc_x2max, ibc_x3min, ibc_x3max};
  // an outflow face is the open one without an entering bath
  for (int s = 0; s < 6; ++s) {
    x.bc[s] = (bc[s] == M1_IBC_X2FLUX) ? M1_IBC_FLUX :
              ((bc[s] == M1_IBC_OUTFLOW) ? M1_IBC_MARSHAK : bc[s]);
  }
  x.iflux[0] = iflux_x1min;
  x.iflux[1] = iflux_x1max;
  x.ifx2[0] = ifx2_lo;
  x.ifx2[1] = ifx2_hi;
  x.mbc = pmy_pack->pmb->mb_bcs.d_view;
  for (int n = 0; n < 6; ++n) x.ebath[n] = 0.0;
  x.ebath[0] = iebath_x1min;
  x.ebath[1] = iebath_x1max;
  for (int s = 0; s < 2; ++s) x.ebath[s] *= (bc[s] != M1_IBC_OUTFLOW);
  x.marshak_q = marshak_q;
  x.marshak = (cln_bdry_form == 1);
  x.use_beam = (ibm_n > 0);
  x.feedback = cln_feedback;
  x.vzero = cln_vzero;
  x.recon = cln_recon;
  x.all_low = cln_all_low;
  x.consist = cln_consist;
  x.milne = cln_milne;
  x.fvrays = (cln_rfv > 0);
  x.hi = false;
  x.mcl = false;
  x.fa = false;
  x.fwc = false;
  x.fs = false;
  x.tf = false;
  x.fd = cln_fsolve && cln_fkin;
  x.gf = cln_gflux;
  x.fre = cln_freal;
  x.drag = cln_vimp && cln_feedback;
  x.upw = cln_upw;
  x.sbi = cln_sbimp && cln_feedback;
  x.sar = !cln_gflux && !(vgd_on && cln_rfv > 0);
  x.op = opac;
  x.nf = pmy_pack->pmesh->three_d ? 3 : (pmy_pack->pmesh->multi_d ? 2 : 1);
  x.is = indcs.is;
  x.ie = indcs.ie;
  x.js = indcs.js;
  x.je = indcs.je;
  x.ks = indcs.ks;
  x.ke = indcs.ke;
  x.multi_d = pmy_pack->pmesh->multi_d;
  x.three_d = pmy_pack->pmesh->three_d;
  x.sph = pmy_pack->pmesh->use_spherical_polar;
  auto &ms = pmy_pack->pmesh->mesh_size;
  x.cseam = 0.5/tan(ms.x2min) + 0.5/tan(ms.x2max);
  x.cub = pmy_pack->pmesh->use_cubed_sphere;
  if (x.sph || (!x.cub && cln_rfv > 0)) {
    auto *pc = pmy_pack->pcoord;
    x.vol = pc->volume;
    x.area[0] = pc->area.x1f;
    x.area[1] = pc->area.x2f;
    x.area[2] = pc->area.x3f;
    x.dxf[0] = pc->dxface.x1f;
    x.dxf[1] = pc->dxface.x2f;
    x.dxf[2] = pc->dxface.x3f;
    x.x1v = pc->x1v;
    x.x1f = pc->xx1f;
    x.x2v = pc->x2v;
    x.x2f = pc->xx2f;
  }
  if (x.cub) {
    x.vol = pmy_pack->pcoord->volume;
    for (int d = 0; d < 3; ++d) {
      x.area[d] = cln_ga[d];
      x.dxf[d] = cln_gd[d];
    }
    x.g = cln_gn;
    x.x1v = pmy_pack->pcoord->x1v;
    x.x1f = pmy_pack->pcoord->xx1f;
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  x.etg = fl.on && fl.use_etotgrav;
  x.mhd = fl.on && fl_mhd;
  x.phi = fl.phicc0;
  x.emag = emag0;
  return x;
}

//----------------------------------------------------------------------------------------
//! \fn const void *RadiationM1::CleanDevCopy
//! \brief a device copy (slot s < 8) of an argument struct of n <= 8 KB bytes (CleanArgs,
//! the EOS) for the device calls that are not inlined: a reference to a kernel's capture
//! is copied to the stack of every thread (~3 KB for CleanArgs), one to device memory is
//! read in place (CleanDev, CleanDevArgs).

const void *RadiationM1::CleanDevCopy(const void *p, const size_t n, const int s) {
  // (the host too: a copy that outlives the struct, as cln_gdx is kept across calls)
  constexpr size_t slot = 8192;
  if (n > slot || s < 0 || s >= 8) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "CleanDevCopy: " << n << " bytes in slot " << s
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (cln_xdev.extent(0) == 0) Kokkos::realloc(cln_xdev, 8*slot);
  Kokkos::View<const char *, HostMemSpace, Kokkos::MemoryUnmanaged>
      h(static_cast<const char *>(p), n);
  Kokkos::deep_copy(Kokkos::subview(cln_xdev, std::make_pair(s*slot, s*slot + n)), h);
  return cln_xdev.data() + s*slot;
}

const CleanArgs *RadiationM1::CleanDevArgs(const CleanArgs &x, const int s) {
  return CleanDev(x, s);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanGhostFill
//! \brief ghost cells of components [c0, c0+nc) of a work array: the neighbours' cells
//! (other MeshBlocks, periodic wrap; across a cubed-sphere seam the neighbour's own
//! cells, without the along-seam resample) by the cell-centred boundary exchange on a
//! scratch array of nc components (only the ghosts and the ng layers the exchange reads
//! are copied in, only the ghosts back), then (phys) a copy of the edge cell on physical
//! boundaries.  part = 1 only starts the exchange (nc <= clean_halo_vars), part = 2
//! finishes it: work on the cells CleanInner by ng can run in between.  The fills
//! deferred by CleanGhostDefer (their ghosts are read only after the next exchange) go
//! in the same exchange where they fit (one exchange of all their components).  The
//! closure's off-diagonal D12, D13, D23 take the edge cell's with their parity: odd
//! across a reflecting radiation wall normal to one of their axes, and D_r,theta,
//! D_r,phi across the pole (e_theta, e_phi flip), so that d_b(D E) puts no shear
//! traction on the wall (as ce91f387d)

void RadiationM1::CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                                 const bool phys, const int part) {
  CleanGhostFillT(a, c0, nc, phys, part);
}

void RadiationM1::CleanGhostFill(DvceArray5D<float> a, const int c0, const int nc) {
  CleanGhostFillT(a, c0, nc, true, 0);
}

void RadiationM1::CleanGhostDefer(DvceArray5D<Real> a, const int c0, const int nc) {
  int np = nc;
  for (auto &s : cln_hpend) np += s.nc;
  if (np > cln_hmax) CleanGhostFlush();
  cln_hpend.push_back({a, c0, nc, true});
}

// (a range lies in one of the two arrays: CLW_NS)
void RadiationM1::CleanGhostFill(const CleanW &a, const int c0, const int nc) {
  CleanGhostFill((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc);
}

void RadiationM1::CleanGhostDefer(const CleanW &a, const int c0, const int nc) {
  CleanGhostDefer((c0 < a.ns) ? a.s : a.o, (c0 < a.ns) ? c0 : c0 - a.ns, nc);
}

void RadiationM1::CleanGhostFlush() {
  if (cln_hpend.empty()) return;
  std::vector<CleanGSeg<Real>> sg;
  sg.swap(cln_hpend);
  CleanGhostExchange(sg, 0);
}

template <typename T>
void RadiationM1::CleanGhostFillT(DvceArray5D<T> a, const int c0, const int nc,
                                  const bool phys, const int part) {
  if (nc > cln_hmax) {   // (the exchange moves at most clean_halo_vars variables at once)
    CleanGhostFillT(a, c0, cln_hmax, phys, 0);
    CleanGhostFillT(a, c0 + cln_hmax, nc - cln_hmax, phys, 0);
    return;
  }
  std::vector<CleanGSeg<T>> sg;
  if constexpr (std::is_same<T, Real>::value) {
    int np = nc;
    for (auto &s : cln_hpend) np += s.nc;
    if (part == 0 && np <= cln_hmax) sg.swap(cln_hpend);
  }
  if (part != 2) CleanGhostFlush();   // (those that do not fit: on their own first)
  sg.push_back({a, c0, nc, phys});
  CleanGhostExchange(sg, part);
}

template <typename T>
void RadiationM1::CleanGhostExchange(const std::vector<CleanGSeg<T>> &sg,
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
    // (the exchanges run one at a time: their scratch arrays share one allocation)
    Kokkos::realloc(cln_hbuf, nmb1 + 1, cln_hmax, n3, n2, n1);
    cln_hpb.reset(new MeshBoundaryValuesCC(pmy_pack, nullptr, false));
    cln_hpb->InitializeBuffers(cln_hmax);
    cln_hpb->SetVectorPairs(cln_hmax, {});
    cln_hpb->cs_noresample = true;
    // (16 teams per buffer: a MeshBlock face is one buffer, which one team packed alone,
    // ~0.5 ms per exchange on an A100 at 831 x 64 x 32; bitwise, copies)
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
  if (dl && part != 2) CleanHaloPost(sg);
  if (dl && part == 1) return;
  if (dl) CleanHaloFinish(sg);
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  CleanHalo &hl = cln_halo[nc];
  if (hl.ca.extent_int(0) == 0) {
    hl.a = DvceArray5D<Real>(cln_hbuf.data(), nmb1 + 1, nc, n3, n2, n1);
    Kokkos::realloc(hl.ca, nmb1 + 1, nc, 1, 1, 1);
  }
  auto h = hl.a;
  const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
  MeshBoundaryValuesCC *pb = cln_hpb.get();
  // one MeshBlock and one array of doubles: the exchange runs on the array's own
  // components [c0, c0 + nc), no copies (the copies took every cell outside CleanInner,
  // ghosts included: the same values)
  DvceArray5D<Real> ha = hl.a;
  bool own = false;
  if constexpr (std::is_same<T, Real>::value) {
    own = (nmb1 == 0 && sg.size() == 1);
    const size_t off = static_cast<size_t>(sg[0].c0)*n3*n2*n1;
    if (own) ha = DvceArray5D<Real>(sg[0].a.data() + off, 1, nc, n3, n2, n1);
  }
  // (only the cells outside CleanInner are copied: CleanShell.  clean_halo_direct: only
  // on the cubed sphere, after the exchange, for the cube-vertex corners of the old one)
  if (!own && (dl ? cs : part != 2)) {
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
  if (dl && cs) {
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
    if (own || (dl && !cs)) break;
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
//! \fn void RadiationM1::CleanHaloInit
//! \brief clean_halo_direct: the cell lists of the ghost exchange, read off the
//! exchange itself: an array of cell ids (rank, MeshBlock, cell) through cln_hpb gives
//! every ghost the cell it copies; the on-rank pairs make the copy kernel, the off-rank
//! ghosts are grouped by rank and their source cells sent to that rank once (its pack
//! list, in the receiver's order: the messages carry no header).  Checked once: a field
//! of other values through both paths must give the same array bitwise, else (and on a
//! refined mesh) the lists stay off and the exchange is cln_hpb's, with a log line.

void RadiationM1::CleanHaloInit(const int n3, const int n2, const int n1) {
  CleanHaloList &q = cln_hdl;
  auto *pm = pmy_pack->pmesh;
  const int me = global_variable::my_rank;
  if (pm->multilevel) {
    q.state = -1;
    if (me == 0) {
      std::cout << "<rad_m1> clean_halo_direct: refined mesh, the clean ghost exchange "
                << "is MeshBoundaryValuesCC's" << std::endl;
    }
    return;
  }
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
  DvceArray5D<Real> x("cln_hdl_x", nmb, 1, n3, n2, n1), y("cln_hdl_y", nmb, 1, n3, n2, n1);
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
    std::vector<int> ld, ls;
    std::vector<std::vector<int>> rg(nrk), rs(nrk);   // off-rank ghosts, their sources
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
            } else {
              rg[r].push_back(m*n321 + g);
              rs[r].push_back(ms*n321 + f);
            }
          }
        }
      }
    }
    std::vector<int> sc, rc;
#if MPI_PARALLEL_ENABLED
    MPI_Comm_dup(MPI_COMM_WORLD, &q.comm);
    std::vector<int> nr(nrk), ns(nrk);
    for (int r = 0; r < nrk; ++r) nr[r] = static_cast<int>(rs[r].size());
    MPI_Alltoall(nr.data(), 1, MPI_INT, ns.data(), 1, MPI_INT, q.comm);
    std::vector<std::vector<int>> ss(nrk);
    std::vector<MPI_Request> rq;
    for (int r = 0; r < nrk; ++r) {
      if (nr[r] == 0 && ns[r] == 0) continue;
      q.rank.push_back(r);
      ss[r].resize(ns[r]);
      rq.resize(rq.size() + 2);
      MPI_Irecv(ss[r].data(), ns[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 2]);
      MPI_Isend(rs[r].data(), nr[r], MPI_INT, r, 4244, q.comm, &rq[rq.size() - 1]);
    }
    MPI_Waitall(static_cast<int>(rq.size()), rq.data(), MPI_STATUSES_IGNORE);
    for (int r : q.rank) {
      q.soff.push_back(static_cast<int>(sc.size()));
      q.slen.push_back(ns[r]);
      sc.insert(sc.end(), ss[r].begin(), ss[r].end());
      q.roff.push_back(static_cast<int>(rc.size()));
      q.rlen.push_back(nr[r]);
      rc.insert(rc.end(), rg[r].begin(), rg[r].end());
    }
#endif
    auto mk = [](const char *nm, const std::vector<int> &v) {
      DvceArray1D<int> a(std::string(nm), std::max(v.size(), static_cast<size_t>(1)));
      Kokkos::deep_copy(Kokkos::subview(a, std::make_pair(static_cast<size_t>(0), v.size())),
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
  q.state = (bad == 0) ? 1 : -1;
  if (me == 0) {
    if (bad == 0) {
      std::cout << "<rad_m1> clean_halo_direct: ON, rank 0 copies " << q.nl
                << " ghosts, sends " << q.ns << " cells, receives " << q.nr << " in "
                << q.rank.size() << " messages" << std::endl;
    } else {
      std::cout << "<rad_m1> clean_halo_direct: the lists miss " << bad << " cells of "
                << "the exchange, the clean ghost exchange is MeshBoundaryValuesCC's"
                << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanHaloPost
//! \brief clean_halo_direct, first half: receives posted (into the half of the receive
//! buffer the previous exchange's unpack does not read), pack, wait for the pack (the
//! stream holds nothing behind it), one send per neighbour rank, then the on-rank copy.
//! Buffer entry p of component n (over the segments in order) at p nc + n

template <typename T>
void RadiationM1::CleanHaloPost(const std::vector<CleanGSeg<T>> &sg) {
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  q.par = 1 - q.par;
#if MPI_PARALLEL_ENABLED
  const int nseg = static_cast<int>(q.rank.size());
  q.rq.assign(nseg, MPI_REQUEST_NULL);
  q.sq.assign(nseg, MPI_REQUEST_NULL);
  Real *rb = q.rbuf.data() + static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  for (int s = 0; s < nseg; ++s) {
    if (q.rlen[s] > 0) {
      MPI_Irecv(rb + static_cast<size_t>(q.roff[s])*nc, q.rlen[s]*nc, MPI_ATHENA_REAL,
                q.rank[s], 4245, q.comm, &q.rq[s]);
    }
  }
  auto sb = q.sbuf;
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.sc;
    const int c0 = s.c0;
    if (q.ns > 0) {
      par_for("cln_hd_pack", DevExeSpace(), 0, q.ns-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        sb(p*nc + o + n) = a(m,c0+n,k,j,r - k*n21 - j*n1);
      });
    }
    o += s.nc;
  }
  if (q.ns > 0 || q.nr > 0) DevExeSpace().fence();
  for (int s = 0; s < nseg; ++s) {
    if (q.slen[s] > 0) {
      MPI_Isend(sb.data() + static_cast<size_t>(q.soff[s])*nc, q.slen[s]*nc,
                MPI_ATHENA_REAL, q.rank[s], 4245, q.comm, &q.sq[s]);
    }
  }
#endif
  for (auto &s : sg) {
    auto a = s.a;
    auto d = q.ld, c = q.ls;
    const int c0 = s.c0;
    if (q.nl == 0) break;
    par_for("cln_hd_copy", DevExeSpace(), 0, q.nl-1, 0, s.nc-1,
    KOKKOS_LAMBDA(const int p, const int n) {
      const int f = d(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
      const int g = c(p), mg = g/n321, t = g - mg*n321, kg = t/n21, jg = (t - kg*n21)/n1;
      a(m,c0+n,k,j,r - k*n21 - j*n1) = a(mg,c0+n,kg,jg,t - kg*n21 - jg*n1);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanHaloFinish
//! \brief clean_halo_direct, second half: wait for the messages, unpack, wait for the
//! sends

template <typename T>
void RadiationM1::CleanHaloFinish(const std::vector<CleanGSeg<T>> &sg) {
#if MPI_PARALLEL_ENABLED
  CleanHaloList &q = cln_hdl;
  int nc = 0;
  for (auto &s : sg) nc += s.nc;
  const int n1 = q.n1, n21 = q.n2*q.n1, n321 = q.n3*n21;
  const int nseg = static_cast<int>(q.rank.size());
  MPI_Waitall(nseg, q.rq.data(), MPI_STATUSES_IGNORE);
  auto rb = q.rbuf;
  const size_t b = static_cast<size_t>(q.par)*std::max(q.nr, 1)*cln_hmax;
  int o = 0;
  for (auto &s : sg) {
    auto a = s.a;
    auto e = q.rc;
    const int c0 = s.c0;
    if (q.nr > 0) {
      par_for("cln_hd_unpack", DevExeSpace(), 0, q.nr-1, 0, s.nc-1,
      KOKKOS_LAMBDA(const int p, const int n) {
        const int f = e(p), m = f/n321, r = f - m*n321, k = r/n21, j = (r - k*n21)/n1;
        a(m,c0+n,k,j,r - k*n21 - j*n1) = static_cast<T>(rb(b + p*nc + o + n));
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
void RadiationM1::CleanGhostPhys(DvceArray5D<T> a, const int c0, const int nc) {
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
  for (int s = 0; s < 6; ++s) rf |= (ib[s] == M1_IBC_REFLECT) ? (1 << s) : 0;
  // (D is the rays'. The x1min end of their sweep keeps the diffusion intensity unless
  // vet_bc_reflect = mirror: there the wall cell's D_12 is that of a field the rays do
  // not mirror, and an odd ghost turned it into a shear gradient D_12 E/dx in the theta
  // faces of the wall row that moved the gas at the flz q wall outward. Its D_12 ghost
  // stays even; D_13 stays odd (no wall torque: L_z exact))
  const bool rm = vet_bc_mirror;
  auto sg = [=] KOKKOS_FUNCTION (const int c, const int s, const int m) {
    const int p = c - CLN_D11 - 3;   // (D12, D13, D23: axes 01, 02, 12; p != 2 - d)
    if (!ts || p < 0 || p > 2) return 1.0;
    if (bc(m,s) == BoundaryFlag::polar) return (p < 2) ? -1.0 : 1.0;
    if (s == 0 && !rm && p == 0) return 1.0;
    return (((rf >> s) & 1) && p != 2 - s/2) ? -1.0 : 1.0;
  };
  // (x1, x2, x3 in turn: a corner takes the product of its sides' parities)
  par_for("cln_ghost1", DevExeSpace(), 0, nmb1, 0, nc-1, 0, n3-1, 0, n2-1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int g) {
    if (pf(bc(m,BoundaryFace::inner_x1))) {
      a(m,c0+n,k,j,is-1-g) = sg(c0+n, 0, m)*a(m,c0+n,k,j,is);
    }
    if (pf(bc(m,BoundaryFace::outer_x1))) {
      a(m,c0+n,k,j,ie+1+g) = sg(c0+n, 1, m)*a(m,c0+n,k,j,ie);
    }
  });
  if (pmy_pack->pmesh->multi_d) {
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
  if (pmy_pack->pmesh->three_d) {
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
//! \fn Real RadiationM1::CleanReduce
//! \brief maximum (max) or sum of v over all ranks

Real RadiationM1::CleanReduce(Real v, const bool max) {
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_ATHENA_REAL, max ? MPI_MAX : MPI_SUM,
                MPI_COMM_WORLD);
#endif
  return v;
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> CleanSums
//! \brief the n <= 2 sums s0, s1 that f(m, k, j, i, s0, s1) adds over the active cells
//! (it may also update the cell's vectors: a BiCGStab step and its dots in one pass), one
//! reduction over the ranks

template <typename F>
std::array<Real, 2> CleanSums(MeshBlockPack *pp, const char *name, const int n,
                              const F &f) {
  auto &indcs = pp->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pp->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  std::array<Real, 2> d = {0.0, 0.0};
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &s0, Real &s1) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    f(m, k+ks, j+js, i+is, s0, s1);
  }, Kokkos::Sum<Real>(d[0]), Kokkos::Sum<Real>(d[1]));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, d.data(), n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  return d;
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> RadiationM1::CleanDots
//! \brief the n <= 2 sums over active cells and the nc components of w(ab[2q])
//! w(ab[2q+1]) (q < n): one pass, one reduction over the ranks

std::array<Real, 2> RadiationM1::CleanDots(CleanW w, const int nc, const int n,
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

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanFaceCoefficients
//! \brief per face: the blend weight w(C1), theta, the old-vector flux G, and the
//! half-range factors h+-_f.  h+_f = H+_f/J_L and h-_f = H-_f/J_R, with H+-_f the rays'
//! half-range fluxes reconstructed to the face from the upwind side (minmod-limited);
//! on the cubed sphere with FV rays their half ranges about the face normal.
//! |h+-_f| <= 1 keeps the operator an M-matrix.

void RadiationM1::CleanFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  auto cs = cln_s;
  auto op = opac;
  auto inc = t2inc;
  const bool fsk = cln_fsolve;
  // (cubed sphere: the opacities of seam ghosts are the neighbour's cells, not of the
  // resampled gas ghosts, as CleanCellFaceCoefficients: a seam face is single-valued)
  if (pmy_pack->pmesh->use_cubed_sphere) {
    CleanGhostFill(opac, 0, opac.extent_int(1), false);
  }
  auto gn = cln_gn;
  auto fi = cln_fvi;
  auto dv = cln_fvd;
  const int nd = pmy_pack->pmesh->use_cubed_sphere ? cln_fvn : 0;

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = M1_T2_F1 + dir;

    par_for("cln_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real c1 = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      fc(m,CLF_W,k,j,i) = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_G,k,j,i) = f0(m,k,j,i) + (stage ? inc(m,tf,k,j,i) : 0.0);
      if (fsk) fc(m,CLF_N+2,k,j,i) = fc(m,CLF_G,k,j,i);   // (fsolve: start of its F)

      // rays' half-range fluxes H+ = h+ J (from the left) and H- = h- J (from the right)
      auto hp = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HP+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto hm = [&](const int kk, const int jj, const int ii) {
        return cs(m,CLN_HM+dir,kk,jj,ii)*cs(m,CLN_J,kk,jj,ii);
      };
      auto minmod = [](const Real a, const Real b) {
        return (a*b <= 0.0) ? 0.0 : ((fabs(a) < fabs(b)) ? a : b);
      };
      const Real hpl = hp(kl,jl,il);
      const Real hpf = hpl + 0.5*minmod(hpl - hp(kl-dk,jl-dj,il-di), hp(k,j,i) - hpl);
      const Real hmr = hm(k,j,i);
      const Real hmf = hmr - 0.5*minmod(hmr - hm(kl,jl,il), hm(k+dk,j+dj,i+di) - hmr);
      const Real jl_ = cs(m,CLN_J,kl,jl,il), jr_ = cs(m,CLN_J,k,j,i);
      fc(m,CLF_HP,k,j,i) = (jl_ > 0.0) ? fmin(fmax(hpf/jl_, 0.0), 1.0)
                                       : cs(m,CLN_HP+dir,kl,jl,il);
      fc(m,CLF_HM,k,j,i) = (jr_ > 0.0) ? fmax(fmin(hmf/jr_, 0.0), -1.0)
                                       : cs(m,CLN_HM+dir,k,j,i);
      if (nd == 0) return;
      // (cubed sphere: about the face's own normal from the FV intensities of L and R,
      // as CleanCellConsistency, so that a seam face takes the same h+- on both sides)
      // (J_L, J_R, H+_L, H-_R as scalars: nvcc -O3 miscompiles a runtime-indexed local)
      Real jl0 = 0.0, jr0 = 0.0, hp0 = 0.0, hm0 = 0.0;
      for (int d = 0; d < nd; ++d) {
        Real mu = 0.0;
        for (int b = 0; b < 3; ++b) mu += dv(d,b)*gn(m,CLG_N+3*dir+b,k,j,i);
        const Real wl = dv(d,3)*fi(m,d,kl,jl,il), wr = dv(d,3)*fi(m,d,k,j,i);
        jl0 += wl;
        jr0 += wr;
        if (mu > 0.0) hp0 += wl*mu;
        if (!(mu > 0.0)) hm0 += wr*mu;
      }
      fc(m,CLF_HP,k,j,i) = (jl0 > 0.0) ? hp0/jl0 : 0.25;
      fc(m,CLF_HM,k,j,i) = (jr0 > 0.0) ? hm0/jr0 : -0.25;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanApply
//! \brief w(dst) = A w(src) over the active cells (the linear part of the operator)

void RadiationM1::CleanApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanGhostFill(cln_w, src, 1);
  CleanStagCsLateral(src, cln_stst);
  CleanArgs x = CleanMakeArgs(), xh = x;
  x.fwc = true;
  xh.hi = cln_split;   // (clean_stag_sbimp: the face momentum, not from CLF_QL)
  auto op = opac;
  if (cln_stst) {   // (clean_stag_stencil: from the rows of the pass)
    auto w = cln_w;
    auto st = cln_st;
    const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
    // (cubed sphere: the lateral stress of every face from the sums Y, CleanFaceCross)
    const bool cub = x.cub;
    const CleanArgs *px = cub ? CleanDevArgs(x, 0) : nullptr;
    par_for("cln_apply_st", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real a = 0.0;
#pragma unroll
      for (int s = 0; s < CLN_ST; ++s) {
        int dk, dj, di;
        CleanStOff(s, dk, dj, di);
        if ((dj != 0 && !md) || (dk != 0 && !td) || (cub && dk*dk + dj*dj + di*di > 1)) {
          continue;   // (no such neighbour)
        }
        a += st(m,s,k-ks,j-js,i-is)*w(m,src,k+dk,j+dj,i+di);
      }
      if (cub) {
        const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
        for (int d = 0; d < 3; ++d) {
          for (int sd = 0; sd < 2 && d < x.nf; ++sd) {
            const int kf = k + sd*(d == 2), jf = j + sd*(d == 1), iff = i + sd*(d == 0);
            a += (sd ? x.h : -x.h)*CleanDivW(x, d, m, k, j, i, sd)*
                 CleanFaceCross(x, src, d, m, kf, jf, iff);
          }
        }
      }
      w(m,dst,k,j,i) = a;
    });
    return;
  }
  const CleanArgs *px = CleanDevArgs(x, 0), *pxh = CleanDevArgs(xh, 1);

  par_for("cln_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const CleanArgs &x = *px, &xh = *pxh;   // (in device memory: no per-thread copy)
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    x.w(m,dst,k,j,i) = x.w(m,src,k,j,i)*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr)
                       + CleanFluxDivergence(x, src, m, k, j, i, false);
    // (clean_stag_sbimp: - (1 - r) sb(E), CleanRightHandSide)
    if (x.sbi) x.w(m,dst,k,j,i) -= qr*CleanStagBetaE(xh, src, m, k, j, i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPreconditionerSetup
//! \brief coefficients of the preconditioner: the x1-line tridiagonal part of A plus the
//! diagonal of its x2 and x3 parts (the off-diagonal Eddington terms are left to the
//! Krylov iteration); clean_stag_pc = mg also the x2 and x3 couplings and the coarse
//! levels

void RadiationM1::CleanPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto op = opac;
  const bool mg = cln_stagmg;
  // the face weights of the pass, once per face for CleanApply (x.fwc)
  CLEAN_HIPX(const CleanArgs *pxw = CleanDevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    auto fc = x.fc[dir];
    par_for("cln_fwc", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const CleanArgs &x = *pxw;)   // (in device memory: no per-thread copy)
      CleanFaceWeights(x, dir, m, k, j, i, fc(m,CLF_QL,k,j,i), fc(m,CLF_QL+1,k,j,i));
      const Real t0 = fc(m,CLF_TH,k,j,i);   // (clean_stag_vimp: theta_d, CleanFaceCross)
      const Real th = t0 + (x.drag ? fc(m,CLF_DB,k,j,i) : 0.0)*(1.0 - t0);
      fc(m,CLF_QL+2,k,j,i) = -(1.0 - fc(m,CLF_W,k,j,i))*th*x.h*x.c*x.c;
    });
  }
  x.fwc = true;
  const CleanArgs *px = CleanDevArgs(x, 0);
  if (cln_stst) {   // (clean_stag_stencil: the operator's rows, CleanApply)
    auto st = cln_st;
    par_for("cln_stencil", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
      Real c[CLN_ST];
      CleanStencilRow(x, m, k, j, i, c);
      for (int s = 0; s < CLN_ST; ++s) st(m,s,k-ks,j-js,i-is) = c[s];
    });
    if (x.cub) {   // (the lateral stress sums' coefficients, CleanStagCsLateral)
      auto cy = cln_csy;
      par_for("cln_csyb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const CleanArgs &x = *px;   // (in device memory: no per-thread copy)
        for (int q = 0; q < 3; ++q) {
          const int dk = (q == 2), dj = (q == 1), di = (q == 0);
          Real y0[3] = {0.0, 0.0, 0.0};
          for (int s = -1; s <= 1; s += 2) {
            const int kf = k + (s > 0)*dk, jf = j + (s > 0)*dj, iff = i + (s > 0)*di;
            const int o = 6*q + 3*(s > 0);
            Real yn[3] = {0.0, 0.0, 0.0};
            if (!CleanBoundaryFace(x, q, m, kf, jf, iff)) {
              const Real a = 0.25*s*x.area[q](m,kf,jf,iff);
              for (int b = 0; b < 3; ++b) {
                for (int c = 0; c < 3; ++c) {
                  const int n = CLN_D11 + ((b == c) ? b : 2 + b + c);
                  const Real ag = a*x.g(m,CLG_N+3*q+c,kf,jf,iff);
                  yn[b] += ag*x.cs(m,n,k+s*dk,j+s*dj,i+s*di);
                  y0[b] -= ag*x.cs(m,n,k,j,i);
                }
              }
            }
            for (int b = 0; b < 3; ++b) cy(m,o+b,k-ks,j-js,i-is) = yn[b];
          }
          for (int b = 0; b < 3; ++b) cy(m,18+3*q+b,k-ks,j-js,i-is) = y0[b];
        }
      });
    }
  }
  par_for("cln_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    CLEAN_HIPX(const CleanArgs &x = *px;)   // (in device memory: no per-thread copy)
    Real al_p, ar_p, al_m, ar_m;
    CleanFaceWeights(x, 0, m, k, j, i+1, al_p, ar_p);
    CleanFaceWeights(x, 0, m, k, j, i, al_m, ar_m);
    const Real qr = x.feedback ? x.w(m,CLW_RY,k,j,i) : 1.0;   // (1 - r)
    Real a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr + x.h*x.idx[0]*(al_p - ar_m);
    x.w(m,CLW_AM,k,j,i) = -x.h*x.idx[0]*al_m;
    x.w(m,CLW_AP,k,j,i) = x.h*x.idx[0]*ar_p;
    if (x.sph || x.cub) {   // (A_f/V weights, CleanDivW)
      const Real wp = CleanDivW(x, 0, m, k, j, i, 1), wm = CleanDivW(x, 0, m, k, j, i, 0);
      a0 = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*qr + x.h*(wp*al_p - wm*ar_m);
      x.w(m,CLW_AM,k,j,i) = -x.h*wm*al_m;
      x.w(m,CLW_AP,k,j,i) = x.h*wp*ar_p;
    }
    if (x.multi_d) {
      CleanFaceWeights(x, 1, m, k, j+1, i, al_p, ar_p);
      CleanFaceWeights(x, 1, m, k, j, i, al_m, ar_m);
      a0 += (x.sph || x.cub) ? x.h*(CleanDivW(x, 1, m, k, j, i, 1)*al_p -
                                    CleanDivW(x, 1, m, k, j, i, 0)*ar_m)
                             : x.h*x.idx[1]*(al_p - ar_m);
      if (mg) {
        const bool cv = x.sph || x.cub;
        x.w(m,CLW_TM,k,j,i) = -x.h*al_m*(cv ? CleanDivW(x, 1, m, k, j, i, 0) : x.idx[1]);
        x.w(m,CLW_TM+1,k,j,i) = x.h*ar_p*(cv ? CleanDivW(x, 1, m, k, j, i, 1) : x.idx[1]);
      }
    }
    if (x.three_d) {
      CleanFaceWeights(x, 2, m, k+1, j, i, al_p, ar_p);
      CleanFaceWeights(x, 2, m, k, j, i, al_m, ar_m);
      const bool cv = x.sph || x.cub;
      a0 += cv ? x.h*(CleanDivW(x, 2, m, k, j, i, 1)*al_p -
                      CleanDivW(x, 2, m, k, j, i, 0)*ar_m)
               : x.h*x.idx[2]*(al_p - ar_m);
      if (mg) {
        x.w(m,CLW_TM+2,k,j,i) = -x.h*al_m*(cv ? CleanDivW(x, 2, m, k, j, i, 0)
                                               : x.idx[2]);
        x.w(m,CLW_TM+3,k,j,i) = x.h*ar_p*(cv ? CleanDivW(x, 2, m, k, j, i, 1) : x.idx[2]);
      }
    }
    x.w(m,CLW_A0,k,j,i) = a0;
  });
  if (mg) CleanStagCoarseSetup();
}

//----------------------------------------------------------------------------------------
//! \fn CleanLevel RadiationM1::CleanStagLevel
//! \brief level l of the staggered mg: l = 0 is the mesh (rows in cln_w), l > 0 the
//! coarse levels (cln_mg[l-1]).  The rows are scalar (CleanLevel with 1x1 blocks) and
//! keep the slot order AM, A0, AP, Thomas scratch of cln_w on every level.

CleanLevel RadiationM1::CleanStagLevel(const int l) {
  CleanLevel lv;
  if (l == 0) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    lv.u = cln_w;
    lv.o[0] = indcs.is;
    lv.o[1] = indcs.js;
    lv.o[2] = indcs.ks;
    lv.n[0] = indcs.nx1;
    lv.n[1] = indcs.nx2;
    lv.n[2] = indcs.nx3;
    lv.a0 = CLW_A0;
  } else {
    lv.u = cln_mg[l-1];
    for (int d = 0; d < 3; ++d) {
      lv.o[d] = 0;
      lv.n[d] = cln_mgn[l-1][d];
    }
    lv.a0 = CLS_A0;
  }
  const int tm = (l == 0) ? CLW_TM : CLS_TM;
  lv.sm[0] = lv.a0 - 1;
  lv.sp[0] = lv.a0 + 1;
  for (int d = 1; d < 3; ++d) {
    lv.sm[d] = tm + 2*(d-1);
    lv.sp[d] = tm + 2*(d-1) + 1;
  }
  return lv;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagCoarseSetup
//! \brief rows of the coarse levels of the staggered mg (Galerkin R A7 P, the old
//! production mg): P piecewise constant over 2 x 2 (x2,x3) aggregates of the same x1
//! index inside one MeshBlock, R = P^T.  A coupling inside the aggregate goes onto the
//! diagonal, one to another aggregate stays, one across the block face is never used.
//! On curvilinear meshes the mesh rows are summed with their cell volumes (as
//! CleanCellCoarseSetup).

void RadiationM1::CleanStagCoarseSetup() {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nlev = static_cast<int>(cln_mg.size()) + 1;
  const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                    pmy_pack->pmesh->use_cubed_sphere;
  auto vol = pmy_pack->pcoord->volume;
  for (int lc = 1; lc < nlev; ++lc) {
    const CleanLevel f = CleanStagLevel(lc-1), c = CleanStagLevel(lc);
    const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
    const bool vw = (lc == 1) && curv;

    par_for("cls_galerkin", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1,
            0, c.n[0]-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
      Real a0 = 0.0, am[3] = {0.0, 0.0, 0.0}, ap[3] = {0.0, 0.0, 0.0};
      for (int a2 = 0; a2 < r[2]; ++a2) {
        for (int a1 = 0; a1 < r[1]; ++a1) {
          const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
          const int i = f.o[0] + ic;
          const int a[3] = {0, a1, a2};
          const Real v = vw ? vol(m,k,j,i) : 1.0;
          a0 += v*f.u(m,f.a0,k,j,i);
          for (int d = 0; d < 3; ++d) {
            if (f.n[d] == 1) continue;
            const Real lo = v*f.u(m,f.sm[d],k,j,i), hi = v*f.u(m,f.sp[d],k,j,i);
            if (a[d] > 0) {
              a0 += lo;
            } else {
              am[d] += lo;
            }
            if (a[d] < r[d] - 1) {
              a0 += hi;
            } else {
              ap[d] += hi;
            }
          }
        }
      }
      c.u(m,c.a0,kc,jc,ic) = a0;
      for (int d = 0; d < 3; ++d) {
        c.u(m,c.sm[d],kc,jc,ic) = am[d];
        c.u(m,c.sp[d],kc,jc,ic) = ap[d];
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagLineSolve
//! \brief u(dst) = L^{-1} u(src) on a level of the staggered preconditioner, L its
//! x1-line tridiagonal part: Thomas solve along the x1 lines of colour col (-1: every
//! line; 0: the lines of even j + k; 1: the odd ones, with the couplings to u(dst) of the
//! even neighbours inside the level on the right-hand side).  On a GPU one team per line
//! solves it by parallel cyclic reduction in team scratch, coalesced along x1, as the old
//! production line solves (M1PCRMG); the host keeps Thomas.

void RadiationM1::CleanStagLineSolve(const CleanLevel &lv, const int src, const int dst,
                                     const int col) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel l = lv;
  const int is = l.o[0], ie = l.o[0] + l.n[0] - 1, cp = l.a0 + 2;
  auto w = l.u;

  if (!std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace> && cln_pcr) {
    const int nx = l.n[0], n2 = l.n[1], n32 = l.n[2]*l.n[1];
    int nround = 0;
    while ((1 << nround) < nx) ++nround;
    const size_t scr = ScrArray1D<Real>::shmem_size(8*nx);
    Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), (nmb1 + 1)*n32, 256);
    const int lv0 = TeamScratchLevel(scr, 0, pol.team_size());
    Kokkos::parallel_for("cln_pcr", pol.set_scratch_size(lv0, Kokkos::PerTeam(scr)),
    KOKKOS_LAMBDA(TeamMember_t tm) {
      const int m = tm.league_rank()/n32, pk = (tm.league_rank() - m*n32)/n2;
      const int pj = tm.league_rank() - m*n32 - pk*n2, k = l.o[2] + pk, j = l.o[1] + pj;
      if (col >= 0 && (pj + pk) % 2 != col) return;   // (team-uniform)
      ScrArray1D<Real> sw(tm.team_scratch(lv0), 8*nx);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
        const int i = is + q;
        Real r = w(m,src,k,j,i);
        if (col == 1) {
          if (pj > 0) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
          if (pj < l.n[1] - 1) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
          if (pk > 0) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
          if (pk < l.n[2] - 1) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
        }
        sw(q) = (q > 0) ? w(m,l.sm[0],k,j,i) : 0.0;
        sw(nx + q) = w(m,l.a0,k,j,i);
        sw(2*nx + q) = (q < nx - 1) ? w(m,l.sp[0],k,j,i) : 0.0;
        sw(3*nx + q) = r;
      });
      tm.team_barrier();
      int o = 0;
      for (int rd = 0, s = 1; rd < nround; ++rd, s *= 2) {
        const int d = 4*nx - o;
        Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
          Real a = 0.0, b = sw(o + nx + q), c = 0.0, r = sw(o + 3*nx + q);
          if (q >= s) {
            const Real f = -sw(o + q)/sw(o + nx + q - s);
            a = f*sw(o + q - s);
            b += f*sw(o + 2*nx + q - s);
            r += f*sw(o + 3*nx + q - s);
          }
          if (q + s < nx) {
            const Real g = -sw(o + 2*nx + q)/sw(o + nx + q + s);
            c = g*sw(o + 2*nx + q + s);
            b += g*sw(o + q + s);
            r += g*sw(o + 3*nx + q + s);
          }
          sw(d + q) = a;
          sw(d + nx + q) = b;
          sw(d + 2*nx + q) = c;
          sw(d + 3*nx + q) = r;
        });
        tm.team_barrier();
        o = d;
      }
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int q) {
        w(m,dst,k,j,is+q) = sw(o + 3*nx + q)/sw(o + nx + q);
      });
    });
    return;
  }
  par_for("cln_thomas", DevExeSpace(), 0, nmb1, l.o[2], l.o[2]+l.n[2]-1,
          l.o[1], l.o[1]+l.n[1]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int pj = j - l.o[1], pk = k - l.o[2];
    if (col >= 0 && (pj + pk) % 2 != col) return;
    Real cprev = 0.0, dprev = 0.0;
    for (int i = is; i <= ie; ++i) {
      const Real am = (i > is) ? w(m,l.sm[0],k,j,i) : 0.0;
      const Real ap = (i < ie) ? w(m,l.sp[0],k,j,i) : 0.0;
      Real r = w(m,src,k,j,i);
      if (col == 1) {
        if (pj > 0) r -= w(m,l.sm[1],k,j,i)*w(m,dst,k,j-1,i);
        if (pj < l.n[1] - 1) r -= w(m,l.sp[1],k,j,i)*w(m,dst,k,j+1,i);
        if (pk > 0) r -= w(m,l.sm[2],k,j,i)*w(m,dst,k-1,j,i);
        if (pk < l.n[2] - 1) r -= w(m,l.sp[2],k,j,i)*w(m,dst,k+1,j,i);
      }
      const Real den = w(m,l.a0,k,j,i) - am*cprev;
      cprev = ap/den;
      dprev = (r - am*dprev)/den;
      w(m,cp,k,j,i) = cprev;
      w(m,dst,k,j,i) = dprev;
    }
    for (int i = ie-1; i >= is; --i) {
      w(m,dst,k,j,i) -= w(m,cp,k,j,i)*w(m,dst,k,j,i+1);
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanStagVCycle
//! \brief u(dst) = M^{-1} u(src) on level l of the staggered mg (the old production mg,
//! rad_m1_precond.cpp at 391cd4c15): z0 = S r, one forward red-black sweep of x1-line
//! solves (the even lines, then the odd ones with the even z0 on the right-hand side);
//! r1 = R (r - A7 z0), the residual of the level's rows inside the MeshBlock summed over
//! each aggregate (with the cell volumes from a curvilinear mesh); z = z0 + P z1, z1 the
//! same on level l+1.  The coarsest level takes the sweep alone; no post-smoothing.

void RadiationM1::CleanStagVCycle(const int l, const int src, const int dst) {
  const CleanLevel f = CleanStagLevel(l);
  CleanStagLineSolve(f, src, dst, 0);
  CleanStagLineSolve(f, src, dst, 1);
  if (l == static_cast<int>(cln_mg.size())) return;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel c = CleanStagLevel(l+1);
  const int r[3] = {1, f.n[1]/c.n[1], f.n[2]/c.n[2]};
  const bool vw = (l == 0) && (pmy_pack->pmesh->use_spherical_polar ||
                                pmy_pack->pmesh->use_cubed_sphere);
  auto vol = pmy_pack->pcoord->volume;

  par_for("cls_restrict", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1, 0, c.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
    Real s = 0.0;
    for (int a2 = 0; a2 < r[2]; ++a2) {
      for (int a1 = 0; a1 < r[1]; ++a1) {
        const int k = f.o[2] + r[2]*kc + a2, j = f.o[1] + r[1]*jc + a1;
        const int i = f.o[0] + ic;
        const int pos[3] = {ic, j - f.o[1], k - f.o[2]};
        Real y = f.u(m,f.a0,k,j,i)*f.u(m,dst,k,j,i);
        for (int d = 0; d < 3; ++d) {
          const int dk = (d == 2), dj = (d == 1), di = (d == 0);
          if (pos[d] > 0) y += f.u(m,f.sm[d],k,j,i)*f.u(m,dst,k-dk,j-dj,i-di);
          if (pos[d] < f.n[d] - 1) y += f.u(m,f.sp[d],k,j,i)*f.u(m,dst,k+dk,j+dj,i+di);
        }
        s += (vw ? vol(m,k,j,i) : 1.0)*(f.u(m,src,k,j,i) - y);
      }
    }
    c.u(m,CLS_B,kc,jc,ic) = s;
  });
  CleanStagVCycle(l+1, CLS_B, CLS_X);
  par_for("cls_prolong", DevExeSpace(), 0, nmb1, f.o[2], f.o[2]+f.n[2]-1,
          f.o[1], f.o[1]+f.n[1]-1, f.o[0], f.o[0]+f.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f.u(m,dst,k,j,i) += c.u(m,CLS_X,(k - f.o[2])/r[2],(j - f.o[1])/r[1],i - f.o[0]);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanPrecondition
//! \brief w(dst) = M^{-1} w(src): tridiagonal (Thomas) solve along every x1 line
//! (clean_stag_pc = line), or the old production mg (mg, CleanStagVCycle)

void RadiationM1::CleanPrecondition(const int src, const int dst) {
  if (cln_stagmg) {
    CleanStagVCycle(0, src, dst);
  } else {
    CleanStagLineSolve(CleanStagLevel(0), src, dst, -1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanBiCGStab
//! \brief right-preconditioned BiCGStab for A X = B, starting from the X it is given:
//! scalar E vectors of cln_w (staggered layout) or 4-component vectors of cln_x (cell
//! layout).  Returns the number of iterations, negative if not converged.

int RadiationM1::CleanBiCGStab() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool cell = cln_cell;
  CleanW w = cell ? CleanW(cln_x) : cln_w;
  const int nc = cell ? 4 : 1;
  // clean_conv = budget, in the Newton passes (a single solve keeps clean_tol in its own
  // norm): the rows over their gross terms S (CleanRightHandSide), W A X = W B with
  // W = 1/S, right-preconditioned by M^-1 W^-1, so that every norm below is the budget's
  const bool bud = cln_bon && !cell && cln_ewr > -1.5;
  const int sx = cell ? CLC_X : CLW_X, sb = cell ? CLC_B : (bud ? CLW_BW : CLW_B);
  const int sr = cell ? CLC_R : CLW_R, sr0 = cell ? CLC_R0 : CLW_R0;
  const int sp = cell ? CLC_P : CLW_P, sv = cell ? CLC_V : CLW_V;
  const int ss = cell ? CLC_R : CLW_S, sq = cell ? CLC_Q : CLW_Q;
  const int sy = cell ? CLC_Y : CLW_Y, sz = cell ? CLC_Z : CLW_Z;

  // w(dst) = sa w(a) + sb w(b)
  auto axpby = [&](const int dst, const int a, const Real sa, const int b,
                   const Real sb_) {
    par_for("cln_axpby", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,dst+c,k,j,i) = sa*w(m,a+c,k,j,i) + sb_*w(m,b+c,k,j,i);
      }
    });
  };
  // w(dst) = w(a) * w(s)^p
  auto scale = [&](const int dst, const int a, const int s, const Real p) {
    par_for("cln_scale", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,dst,k,j,i) = (p > 0.0) ? w(m,a,k,j,i)*w(m,s,k,j,i) : w(m,a,k,j,i)/w(m,s,k,j,i);
    });
  };
  auto apply = [&](const int src, const int dst) {
    if (cell) {
      CleanCellApply(src, dst);
    } else {
      CleanApply(src, dst);
      if (bud) scale(dst, dst, CLW_SE, -1.0);
    }
  };
  auto precondition = [&](const int src, const int dst) {
    if (cell) {
      CleanCellPrecondition(src, dst);
    } else if (bud) {
      scale(CLW_PW, src, CLW_SE, 1.0);
      CleanPrecondition(CLW_PW, dst);
    } else {
      CleanPrecondition(src, dst);
    }
  };
  if (bud) scale(CLW_BW, CLW_B, CLW_SE, -1.0);
  auto dots = [&](const int n, const std::array<int, 4> ab) {
    return CleanDots(w, nc, n, ab);
  };

  // (re)start from the true residual of X with the shadow residual r0 = r.  A breakdown
  // (r0.r at the rounding level of |r0| |r|, or r0.v or omega zero) restarts; a second
  // breakdown right after a restart fails the solve.  (A step's vector updates and the
  // dots after them take one pass, CleanSums, in the arithmetic and summation order of
  // separate update kernels and CleanDots.)  The dots of an iteration take four
  // reductions: r0.v, s.s, (q.q, q.s), and (r.r, r0.r) for the next one.  Convergence of
  // the recursive residual is checked on the true one, which restarts if it is above
  // 10 tol (rows of very different scale, e.g. the empty side cells of a diagonal beam
  // whose flux factors are uncapped, can leave the recursion far from it)
  Real rho = 1.0, alpha = 1.0, omega = 1.0, rr = 0.0, r0r0 = 0.0, rho_new = 0.0, bb = 0.0;
  int it0 = 0;
  auto restart = [&](const int it, const bool force = false) {
    if (it == it0 + 1 && !force) return false;
    apply(sx, sr);
    const int n = (it < 0) ? 2 : 1;
    const auto d = CleanSums(pmy_pack, "cln_dot_rst", n,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {   // (r = b - A x, r0 = r, p = v = 0)
        w(m,sr+c,k,j,i) = 1.0*w(m,sb+c,k,j,i) + -1.0*w(m,sr+c,k,j,i);
        w(m,sr0+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sp+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        w(m,sv+c,k,j,i) = 0.0*w(m,sr+c,k,j,i) + 0.0*w(m,sr+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        if (n > 1) s1 += w(m,sb+c,k,j,i)*w(m,sb+c,k,j,i);
      }
    });
    rho = alpha = omega = 1.0;
    rr = r0r0 = rho_new = d[0];
    if (it < 0) bb = d[1];
    it0 = it;
    return true;
  };
  restart(-1);
  Real tol = cln_tol*fmax(std::sqrt(bb), 1.0e-300);
  // Eisenstat-Walker forcing (choice 2, safeguarded) on |r0|, the Newton residual of the
  // pass in the budget's norm
  if (bud) {
    const Real r0 = std::sqrt(rr);
    Real eta = 0.1;
    if (cln_ewr > 0.0) {
      eta = 0.9*SQR(r0/cln_ewr);
      if (0.9*SQR(cln_eta) > 0.1) eta = fmax(eta, 0.9*SQR(cln_eta));
      eta = fmin(eta, 0.1);
    }
    cln_ewr = r0;
    cln_eta = eta;
    tol = fmax(eta*r0, 0.1*cln_tolb);
  }
  // clean_init_cnorm > 0, the steady start's solves, and clean_lin_cnorm > 0, the steps'
  // strict ones (as the old production's implicit_lin_cnorm; the budget's passes keep
  // their own per-row test): besides the test relative to |b|, which leaves an envelope
  // far below the core at its first preconditioned iterate (MeshBlock-dependent, rule
  // 5), max_i |r_i|/(a_i E_i) <= cnorm, a_i the row's diagonal (CLW_A0), E_i the
  // largest of the Newton iterate E^k and the solve's |X| over the row's stencil: an
  // M-matrix row's off-diagonals sum to at most a_i, so a_i E_i bounds its round-off
  // (Oettli-Prager), also in a cell next to a beam whose own E is ~0 (with E^k alone the
  // test was unattainable there, and every solve ran to the stall window).
  // (Production's s_i = 1 + h c chi_E is below the round-off of a transport-dominated
  // row.)  Where that measure stops falling (stall window below) with the global test
  // met, it is at the round-off of an empty or free-streaming row, and the solve is
  // converged
  const Real cnt = cln_insteady ? cln_cnorm : cln_lcnorm;
  const bool cn = (cnt > 0.0) && !bud && !cell;
  const Real efl = e_floor, inf = std::numeric_limits<Real>::infinity();
  Real cmin = inf, cwin = inf;
  bool cmet = false;
  auto conv = [&](const int s, const Real r2, const Real t) {
    if (!(std::sqrt(r2) <= t)) return false;
    if (!cn) return true;
    Real a = 0.0;
    const int dj = pmy_pack->pmesh->multi_d, dk = pmy_pack->pmesh->three_d;
    Kokkos::parallel_reduce("cln_cnorm", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
        DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &mx) {
      Real e = fmax(w(m,CLW_EK,k,j,i), efl);
      for (int o = -1; o <= 1; o += 2) {
        e = fmax(e, fmax(fabs(w(m,sx,k,j,i+o)), fabs(w(m,sx,k,j+o*dj,i))));
        e = fmax(e, fmax(fabs(w(m,sx,k+o*dk,j,i)), fabs(w(m,sx,k,j,i))));
      }
      mx = fmax(mx, fabs(w(m,s,k,j,i))/fmax(cnt*fabs(w(m,CLW_A0,k,j,i))*e, 1.0e-300));
    }, Kokkos::Max<Real>(a));
    a = CleanReduce(a, true);
    cmet = true;
    cmin = fmin(cmin, a);
    return a <= 1.0;
  };
  if (conv(sr, rr, tol)) return 0;
  auto done = [&](const int it) {return restart(it, true) && conv(sr, rr, 10.0*tol);};

  const Real eps = std::numeric_limits<Real>::epsilon();
  // clean_bicg_stall_its N: a solve whose least residual (the recursive one, or the true
  // one of a restart) fell by less than clean_bicg_stall_fac over the last N iterations
  // has stalled, and fails as at the iteration cap
  Real rmin = std::sqrt(rr), rwin = rmin;
  for (int it = 1; it <= cln_maxit; ++it) {
    rmin = fmin(rmin, std::sqrt(rr));
    if (cln_stall_n > 0 && it > 1 && (it - 1) % cln_stall_n == 0) {
      if (cmet) {   // (cnorm with the global test met: the per-cell measure's window)
        if (!(cmin < cln_stall_f*cwin)) {
          const bool ok = restart(it, true) && std::sqrt(rr) <= 10.0*tol;
          cln_ncwin += ok ? 1 : 0;   // (counted: the CLEAN line's cwin)
          return ok ? it : -it;
        }
        cwin = cmin;
      } else if (!(rmin < cln_stall_f*rwin)) {
        if (global_variable::my_rank == 0) {
          std::printf("CLEAN BICGSTAB stall cycle %d it %d |r|/tol %.3e\n",
                      pmy_pack->pmesh->ncycle, it - 1, rmin/tol);
        }
        return -it;
      }
      rwin = rmin;
    }
    if (!(fabs(rho_new) > eps*std::sqrt(r0r0*rr))) {
      if (restart(it)) continue;
      return -it;
    }
    const Real beta = (rho_new/rho)*(alpha/omega);
    rho = rho_new;
    par_for("cln_pupdate", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < nc; ++c) {
        w(m,sp+c,k,j,i) = w(m,sr+c,k,j,i)
                          + beta*(w(m,sp+c,k,j,i) - omega*w(m,sv+c,k,j,i));
      }
    });
    precondition(sp, sy);
    apply(sy, sv);
    const Real r0v = dots(1, {sr0, sv, 0, 0})[0];
    if (!std::isfinite(r0v) || r0v == 0.0) {
      if (restart(it)) continue;
      return -it;
    }
    alpha = rho/r0v;
    const Real ma = -alpha;
    const Real s2 = CleanSums(pmy_pack, "cln_dot_s", 1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0, Real &) {
      for (int c = 0; c < nc; ++c) {   // (s = r - alpha v)
        w(m,ss+c,k,j,i) = 1.0*w(m,sr+c,k,j,i) + ma*w(m,sv+c,k,j,i);
        s0 += w(m,ss+c,k,j,i)*w(m,ss+c,k,j,i);
      }
    })[0];
    if (conv(ss, s2, tol)) {
      axpby(sx, sx, 1.0, sy, alpha);
      if (done(it)) return it;
      continue;
    }
    precondition(ss, sz);
    apply(sz, sq);
    const auto d = dots(2, {sq, sq, sq, ss});
    omega = (d[0] > 0.0) ? d[1]/d[0] : 0.0;
    const Real al_ = alpha, om = omega;
    const auto e = CleanSums(pmy_pack, "cln_dot_x", 2,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &s0,
                  Real &s1) {
      for (int c = 0; c < nc; ++c) {
        w(m,sx+c,k,j,i) += al_*w(m,sy+c,k,j,i) + om*w(m,sz+c,k,j,i);
        w(m,sr+c,k,j,i) = w(m,ss+c,k,j,i) - om*w(m,sq+c,k,j,i);
        s0 += w(m,sr+c,k,j,i)*w(m,sr+c,k,j,i);
        s1 += w(m,sr0+c,k,j,i)*w(m,sr+c,k,j,i);
      }
    });
    rr = e[0];
    rho_new = e[1];
    if (conv(sr, rr, tol) && done(it)) return it;
    if (omega == 0.0 && !restart(it)) return -it;
  }
  return -cln_maxit;
}

}  // namespace radm1
