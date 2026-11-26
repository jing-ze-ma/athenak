//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cs_halo.cpp
//! \brief cubed-sphere rays: band halos, full and compact.
//!
//!  called from: rays/sc_cs_setup.cpp, rays/sc_cs_sweep.cpp, rays/sc_sph_halo.cpp

#include <math.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <thread>  // NOLINT(build/c++11)
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_cs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

extern int64_t g_vgd_csgen[16];       // rays/sc_cs_setup.cpp; VET_VGD_CSCL_MAX
extern int64_t g_vgd_csgen_save[16];

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsCompact
//! \brief compact lists of the current direction set, LRU cache of
//! max(2, vet_rotate_avg) entries. (1) read mask: replay the sweep's upwind geometry and
//! mark every remote band value read (1 inward, 2 outward; both cut variants, so the
//! lists depend on the set alone); (2) receivers send marks to owners; (3) both sides
//! build lists in the same (pass, shell, partner, entry, direction) order.

void RadiationVet::VetGdCsCompact() {
  int c = -1;
  for (int q = 0; q < vgd_csncl; ++q) {
    if (vgd_cscl[q].alpha == vgd_alpha) {c = q;}
  }
  if (c >= 0 && vgd_cscl[c].msk.data() != nullptr) {
    vgd_cscur = c;
    vgd_cscl[c].use = ++vgd_csuse;
    return;
  }
  // a cached set whose mask was reused gets its mask again (step 1 only)
  const bool remask = (c >= 0);
  Kokkos::Timer tb;
  if (!remask) {
    c = 0;   // least recently used entry (first of equal ones)
    for (int q = 1; q < vgd_csncl; ++q) {
      if (vgd_cscl[q].use < vgd_cscl[c].use) {c = q;}
    }
  }
  VgdCsCl &E = vgd_cscl[c];
  // two read masks only (vgd_csmb); the mask of the entry in use is kept
  if (E.msk.data() == nullptr) {
    const auto *k = vgd_cscl[vgd_cscur].msk.data();
    const int b = (k != nullptr && k == vgd_csmb[0].data()) ? 1 : 0;
    for (auto &e : vgd_cscl) {
      if (e.msk.data() == vgd_csmb[b].data()) {e.msk = decltype(e.msk)();}
    }
    const int64_t nmv = static_cast<int64_t>(std::max(vgd_csnrcv, 1))*vgd_csper;
    if (vgd_csmb[b].extent(0) != static_cast<size_t>((nmv + 15)/16)) {
      vgd_csmb[b] = decltype(E.msk)();
      Kokkos::realloc(vgd_csmb[b], (nmv + 15)/16);
    }
    E.msk = vgd_csmb[b];
  }
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int c2b = indcs.nx2 + 2*wb, c3b = indcs.nx3 + 2*wb;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  auto msk_ = E.msk;
  Kokkos::deep_copy(msk_, static_cast<uint32_t>(0));
  auto dir_ = vgd_dir;
  auto mr_ = vgd_mr;
  auto cr_ = vgd_csr;
  auto bl_ = vgd_csblk;
  auto cm_ = vgd_cscm;
  auto nb_ = vgd_csnb;
  auto rx_ = vgd_csrix;
  const int nc = vgd_csn;
  const Real h = 0.5*M_PI/nc;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const int64_t per_ = vgd_csper;   // n*n1 padded to 16
  Kokkos::Timer tsub;
  // (1) read mask (VetGdCsSweep's geometry; values never matter)
  par_for("vet_vgdcs_rmask", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, n - 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
    const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
    const Real mr = mr_(m,k,j,d);
    const bool inw = !(mr >= 0.0);
    const Real rx = cr_(m,k,j,0), ry = cr_(m,k,j,1), rz = cr_(m,k,j,2);
    // mark the reads of one upwind segment (ds, upwind shell iu, branch want_in)
    auto mark = [&](const Real r, const Real ds, const int iu, const bool want_in) {
      Real xu[3] = {r*rx - ds*nx, r*ry - ds*ny, r*rz - ds*nz};
      const Real ru = sqrt(xu[0]*xu[0] + xu[1]*xu[1] + xu[2]*xu[2]);
      xu[0] /= ru;
      xu[1] /= ru;
      xu[2] /= ru;
      VgdCsSt st;
      VgdCsStencil(xu, bl_(m,0), bl_(m,1), bl_(m,2), wb, c2b, c3b, nc, h, cm_, nb_, st);
      if (st.n == 0) {return;}   // the own column (interior)
      bool ok[4];
      bool all = true;
      for (int e = 0; e < st.n; ++e) {
        const Real mq = mr_(m,st.k[e],st.j[e],d);
        ok[e] = want_in ? (mq < 0.0) : (mq >= 0.0);
        all = all && ok[e];
      }
      for (int e = 0; e < st.n; ++e) {
        if (!all && !(ok[e] && st.w[e] > 0.0)) {continue;}
        const int x = rx_(m,st.k[e],st.j[e]);
        if (x >= 0) {
          const int64_t t = static_cast<int64_t>(x)*per_ + d*n1 + (iu - is);
          Kokkos::atomic_or(&msk_(t >> 4), (want_in ? 1u : 2u) << (2*(t & 15)));
        }
      }
    };
    for (int l = 0; l < n1; ++l) {
      const int i = is + l;
      const Real r = cx1v(m,i);
      const Real zr = r*fabs(mr);
      const Real p2 = fmax(r*r - zr*zr, 0.0);
      if (inw) {
        if (l == n1 - 1) {continue;}
        const Real ru = cx1v(m,i+1);
        mark(r, sqrt(ru*ru - p2) - zr, i + 1, true);
        continue;
      }
      if (zr < 1.0e-12*r) {continue;}
      // as the shell above the cut (l > lcut)
      if (l > 0) {
        const Real rd = cx1v(m,i-1);
        if (p2 <= rd*rd*(1.0 + 1.0e-13)) {
          mark(r, zr - sqrt(fmax(rd*rd - p2, 0.0)), i - 1, false);
        } else {
          mark(r, 2.0*zr, i, true);
        }
      }
      // as the cut shell (l == lcut): the inner face (no intensity read) or turned
      const Real rf = cx1f(m,i);
      if (!(p2 <= rf*rf)) {mark(r, 2.0*zr, i, true);}
    }
  });
  Kokkos::fence();
  vgd_cstb_k += tsub.seconds();
  tsub.reset();
  if (remask) {
    vgd_cscur = c;
    E.use = ++vgd_csuse;
    return;
  }
  const int np = static_cast<int>(vgd_csprk.size());
  const size_t per = static_cast<size_t>(per_);
  // (3) the lists: counts per (pass, shell, partner), offsets, fill in (entry, d) order
  const size_t nsl = static_cast<size_t>(2)*n1;
  std::vector<int> sl_h, rl_h;
  int64_t ts = 0, tr = 0;
  if (per % 16 != 0) {VgdCsFatal("compact lists: n*n1 not padded to 16");}
  // lists from packed 2-bit words (partner blocks word aligned)
  const size_t wpe = per/16;   // words per entry
  // only the nonzero (word index, word) pairs, in increasing word order, are downloaded,
  // exchanged and scanned. Host copies use persistent buffers (grown, never freed): a
  // fresh buffer per build let UCX's registration cache send stale pages under
  // LD_PRELOAD profilers (MPI_ERR_TRUNCATE or hang).
  const int64_t nwr = static_cast<int64_t>(vgd_csnrcv)*static_cast<int64_t>(wpe);
  using RpL = Kokkos::RangePolicy<DevExeSpace, Kokkos::IndexType<int64_t>>;
  int64_t nz = 0;
  Kokkos::parallel_scan("vet_vgdcs_nzn", RpL(DevExeSpace(), 0, nwr),
  KOKKOS_LAMBDA(const int64_t w, int64_t &u, const bool) {
    if (msk_(w) != 0u) {u += 1;}
  }, nz);
  DvceArray1D<int64_t> zi("vet_vgdcs_zi", std::max<int64_t>(nz, 1));
  DvceArray1D<uint32_t> zv("vet_vgdcs_zv", std::max<int64_t>(nz, 1));
  Kokkos::parallel_scan("vet_vgdcs_nzc", RpL(DevExeSpace(), 0, nwr),
  KOKKOS_LAMBDA(const int64_t w, int64_t &u, const bool fin) {
    const uint32_t x = msk_(w);
    if (x == 0u) {return;}
    if (fin) {zi(u) = w; zv(u) = x;}
    u += 1;
  });
  std::vector<int64_t> hzi(std::max<int64_t>(nz, 1));
  std::vector<uint32_t> hzv(std::max<int64_t>(nz, 1));
  Kokkos::deep_copy(Kokkos::View<int64_t *, Kokkos::HostSpace,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(hzi.data(), hzi.size()), zi);
  Kokkos::deep_copy(Kokkos::View<uint32_t *, Kokkos::HostSpace,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(hzv.data(), hzv.size()), zv);
  bool both = false;
  for (int64_t q = 0; q < nz; ++q) {
    const uint32_t x = hzv[q];
    if ((x & (x >> 1) & 0x55555555u) != 0u) {both = true;}
  }
  if (both) {VgdCsFatal("read mask: a value marked both inward and outward");}
  // own pairs per partner, into the persistent buffer
  std::vector<size_t> ho(np + 1, 0), hn(np, 0);
  std::vector<int64_t> hlo(np, 0);
  for (int p = 0; p < np; ++p) {
    const int64_t w0 = static_cast<int64_t>(vgd_cspr[p])*wpe;
    const int64_t w1 = w0 + static_cast<int64_t>(vgd_cspm[p])*wpe;
    if (w1 - w0 > static_cast<int64_t>(UINT32_MAX)) {
      VgdCsFatal("compact lists: a partner block of more than 2^32 mask words");
    }
    hlo[p] = std::lower_bound(hzi.begin(), hzi.begin() + nz, w0) - hzi.begin();
    hn[p] = (vgd_cspm[p] > 0)
            ? (std::lower_bound(hzi.begin(), hzi.begin() + nz, w1) - hzi.begin()) - hlo[p]
            : 0;
    ho[p+1] = ho[p] + 2*hn[p];
  }
  if (vgd_csmw_h.size() < std::max<size_t>(ho[np], 1)) {vgd_csmw_h.resize(ho[np]);}
  uint32_t *mp = vgd_csmw_h.data();
  for (int p = 0; p < np; ++p) {
    const int64_t w0 = static_cast<int64_t>(vgd_cspr[p])*wpe;
    for (size_t q = 0; q < hn[p]; ++q) {
      mp[ho[p] + 2*q] = static_cast<uint32_t>(hzi[hlo[p] + q] - w0);
      mp[ho[p] + 2*q + 1] = hzv[hlo[p] + q];
    }
  }
  vgd_cstb_d += tsub.seconds();
  tsub.reset();
  // owners' pairs (marks of the entries they send me), per partner
  std::vector<uint64_t> rn(np, 0), sn(np, 0);
  for (int p = 0; p < np; ++p) {sn[p] = hn[p];}
  std::vector<size_t> ro(np + 1, 0);
#if MPI_PARALLEL_ENABLED
  {
    // tags 7109 pair counts, 7107 pairs, 7108 checksums (7105: fused twin exchange);
    // a checksum mismatch is fatal
    auto hsum = [](const uint32_t *a, const size_t nw) {
      uint64_t h = 1469598103934665603ull;
      for (size_t w = 0; w < nw; ++w) {h = (h ^ a[w])*1099511628211ull;}
      return h;
    };
    std::vector<MPI_Request> req;
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_csprk[p];
      if (vgd_cspn[p] > 0) {
        req.emplace_back();
        MPI_Irecv(&rn[p], 1, MPI_UINT64_T, rk, 7109, MPI_COMM_WORLD, &req.back());
      }
      if (vgd_cspm[p] > 0) {
        req.emplace_back();
        MPI_Isend(&sn[p], 1, MPI_UINT64_T, rk, 7109, MPI_COMM_WORLD, &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    req.clear();
    for (int p = 0; p < np; ++p) {ro[p+1] = ro[p] + 2*rn[p];}
    if (vgd_csmw_r.size() < std::max<size_t>(ro[np], 1)) {vgd_csmw_r.resize(ro[np]);}
    std::vector<uint64_t> hs(np, 0), hr(np, 0);
    for (int p = 0; p < np; ++p) {
      const int rk = vgd_csprk[p];
      if (vgd_cspn[p] > 0) {
        req.emplace_back();
        MPI_Irecv(vgd_csmw_r.data() + ro[p], static_cast<int>(2*rn[p]), MPI_UINT32_T,
                  rk, 7107, MPI_COMM_WORLD, &req.back());
        req.emplace_back();
        MPI_Irecv(&hr[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
      }
      if (vgd_cspm[p] > 0) {
        hs[p] = hsum(mp + ho[p], 2*hn[p]);
        req.emplace_back();
        MPI_Isend(mp + ho[p], static_cast<int>(2*hn[p]), MPI_UINT32_T, rk, 7107,
                  MPI_COMM_WORLD, &req.back());
        req.emplace_back();
        MPI_Isend(&hs[p], 1, MPI_UINT64_T, rk, 7108, MPI_COMM_WORLD, &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    for (int p = 0; p < np; ++p) {
      if (vgd_cspn[p] <= 0) {continue;}
      if (hsum(vgd_csmw_r.data() + ro[p], 2*rn[p]) != hr[p]) {
        VgdCsFatal("compact lists: the read marks from rank " +
                   std::to_string(vgd_csprk[p]) + " arrived corrupted (rank " +
                   std::to_string(global_variable::my_rank) + ", checksum mismatch): "
                   "the MPI transport delivered other bytes than were sent (seen with "
                   "UCX's registration cache under an LD_PRELOAD profiler such as nsys "
                   "-t osrt / -s cpu: export UCX_RCACHE_ENABLE=n for such runs)");
      }
    }
  }
#endif
  vgd_cstb_m += tsub.seconds();
  tsub.reset();
  // each partner's counts and fills on its own host thread (disjoint ranges)
  auto parts = [np](const auto &f) {
    std::vector<std::thread> th;
    for (int p = 1; p < np; ++p) th.emplace_back(f, p);
    if (np > 0) f(0);
    for (auto &t : th) t.join();
  };
  // lists from partner p's pairs (pr + po[p], pn[p] of them), entries from st0[p]
  auto buildw = [&](const uint32_t *pr, const std::vector<size_t> &po,
                    const std::vector<int> &st0, std::vector<int64_t> &off,
                    std::vector<int> &lst) {
    std::vector<int64_t> nn(nsl*np, 0);
    parts([&](const int p) {
      const size_t q1 = (po[p+1] - po[p])/2;
      for (size_t q = 0; q < q1; ++q) {
        uint32_t x = pr[po[p] + 2*q + 1];
        const size_t r0 = static_cast<size_t>(pr[po[p] + 2*q])*16;   // in-block t
        while (x != 0u) {
          const int b = __builtin_ctz(x) >> 1;
          const uint32_t v = (x >> (2*b)) & 3u;
          x &= ~(3u << (2*b));
          const size_t r = (r0 + b) % per;
          const int sh = static_cast<int>(r % n1);
          nn[(static_cast<size_t>(v - 1)*n1 + sh)*np + p] += 1;
        }
      }
    });
    off.assign(nsl*(np + 1), 0);
    int64_t t = 0;
    for (size_t sl = 0; sl < nsl; ++sl) {
      for (int p = 0; p < np; ++p) {
        off[sl*(np + 1) + p] = t;
        t += nn[sl*np + p];
      }
      off[sl*(np + 1) + np] = t;
    }
    lst.assign(std::max<int64_t>(t, 1), 0);
    std::vector<int64_t> cur(nsl*np);
    for (size_t sl = 0; sl < nsl; ++sl) {
      for (int p = 0; p < np; ++p) {cur[sl*np + p] = off[sl*(np + 1) + p];}
    }
    parts([&](const int p) {
      const size_t q1 = (po[p+1] - po[p])/2;
      for (size_t q = 0; q < q1; ++q) {
        uint32_t x = pr[po[p] + 2*q + 1];
        const size_t r0 = static_cast<size_t>(pr[po[p] + 2*q])*16;
        while (x != 0u) {   // bits low to high = t increasing
          const int b = __builtin_ctz(x) >> 1;
          const uint32_t v = (x >> (2*b)) & 3u;
          x &= ~(3u << (2*b));
          const size_t rr = r0 + b;
          const int e = st0[p] + static_cast<int>(rr/per);
          const size_t r = rr % per;
          const int d = static_cast<int>(r/n1), sh = static_cast<int>(r % n1);
          lst[cur[(static_cast<size_t>(v - 1)*n1 + sh)*np + p]++] = e*n + d;
        }
      }
    });
    return t;
  };
  // send and receive lists concurrently
  std::thread tl([&]() {tr = buildw(mp, ho, vgd_cspr, E.ro, rl_h);});
  ts = buildw(vgd_csmw_r.data(), ro, vgd_csps, E.so, sl_h);
  tl.join();
  vgd_cstb_h += tsub.seconds();
  tsub.reset();
  // (entry, direction) order is the same on both sides; counts per (pass, shell,
  // partner) are checked by the message sizes below
  Kokkos::realloc(E.si, sl_h.size());
  Kokkos::realloc(E.ri, rl_h.size());
  {
    auto a = Kokkos::create_mirror_view(E.si);
    for (size_t q = 0; q < sl_h.size(); ++q) {a(q) = sl_h[q];}
    Kokkos::deep_copy(E.si, a);
    auto b = Kokkos::create_mirror_view(E.ri);
    for (size_t q = 0; q < rl_h.size(); ++q) {b(q) = rl_h[q];}
    Kokkos::deep_copy(E.ri, b);
  }
  vgd_cstb_u += tsub.seconds();
  E.alpha = vgd_alpha;
  E.use = ++vgd_csuse;
  vgd_cscur = c;
  g_vgd_csgen[c] += 1;
  vgd_csnbuild += 1.0;
  vgd_cstbuild += tb.seconds();
  if (global_variable::my_rank == 0 && vgd_csnbuild <= 1.0) {
    const double dn = static_cast<double>(vgd_csnsnd)*n*n1*2.0;
    std::cout << "<rad_vet> vet_gd (cubed sphere) compact lists: rank 0 sends " << ts
              << " values per sweep (dense " << dn << ": "
              << ((dn > 0.0) ? ts/dn : 0.0) << "), receives " << tr << "; built in "
              << tb.seconds() << " s (once per direction set)" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsHaloC
//! \brief compact exchange of shell i after pass `pass`: the listed (entry, direction)
//! values of every partner (VetGdCsCompact)

void RadiationVet::VetGdCsHaloC(const int pass, const int i) {
#if MPI_PARALLEL_ENABLED
  const VgdCsCl &E = vgd_cscl[vgd_cscur];
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1;
  const int n = vgd_n;
  const int np = static_cast<int>(vgd_csprk.size());
  const size_t sl = static_cast<size_t>(pass)*n1 + (i - indcs.is);
  const int64_t *so = &E.so[sl*(np + 1)];
  const int64_t *ro = &E.ro[sl*(np + 1)];
  const int64_t s0 = so[0], ns = so[np] - so[0];
  const int64_t r0 = ro[0], nr = ro[np] - ro[0];
  auto a_ = vgd_i;
  auto sb_ = vgd_cssb;
  auto rb_ = vgd_csrb;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto si_ = E.si;
  auto ri_ = E.ri;
  if (static_cast<int64_t>(sb_.extent(0)) < ns ||
      static_cast<int64_t>(rb_.extent(0)) < nr) {
    VgdCsFatal("VetGdCsHaloC: message buffer too small");
  }
  Real *sbp = sb_.data();
  Real *rbp = rb_.data();
  if (ns > 0) {
    Kokkos::parallel_for("vet_vgdcs_cpack",
                         VgdLight(DevExeSpace(), 0, static_cast<int>(ns)),
    KOKKOS_LAMBDA(const int o) {
      const int x = si_(s0 + o);
      const int e = x/n, d = x - e*n;
      sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
    });
  }
  Kokkos::fence();
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    const int64_t cr = ro[p+1] - ro[p], cs = so[p+1] - so[p];
    if (cr > 0) {
      req.emplace_back();
      MPI_Irecv(rb_.data() + (ro[p] - r0), static_cast<int>(cr), MPI_ATHENA_REAL, rk,
                7103, MPI_COMM_WORLD, &req.back());
    }
    if (cs > 0) {
      req.emplace_back();
      MPI_Isend(sb_.data() + (so[p] - s0), static_cast<int>(cs), MPI_ATHENA_REAL, rk,
                7103, MPI_COMM_WORLD, &req.back());
      vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  vgd_tmpi += tq.seconds();
  vgd_nexch += 1.0;
  vgd_csvsent += static_cast<double>(ns);
  vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;
  if (nr > 0) {
    Kokkos::parallel_for("vet_vgdcs_cunpack",
                         VgdLight(DevExeSpace(), 0, static_cast<int>(nr)),
    KOKKOS_LAMBDA(const int o) {
      const int x = ri_(r0 + o);
      const int e = x/n, d = x - e*n;
      a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rbp[o];
    });
  }
#else
  (void) pass;
  (void) i;
#endif
}

namespace {

// per-partner offsets of the merged two-set exchange (VetGdCsHaloCX)
constexpr int VGD_CSX_NS = 32;

struct VgdCsXIdx {int v[VGD_CSX_NS + 1];};

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::VetGdCsHaloCX
//! \brief both compact exchanges of shell i of a fused two-set sweep in one round (one
//! message per partner: current set's piece, then the next set's; tag 7107); bitwise two
//! VetGdCsHaloC calls. false (nothing done): the caller runs the two exchanges.

bool RadiationVet::VetGdCsHaloCX(const int pass, const int i, const int eb) {
#if MPI_PARALLEL_ENABLED
  if (eb < 0) {return false;}
  const VgdCsCl &EA = vgd_cscl[vgd_cscur];
  const VgdCsCl &EB = vgd_cscl[eb];
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int n1 = indcs.nx1;
  const int n = vgd_n;
  const int np = static_cast<int>(vgd_csprk.size());
  if (np > VGD_CSX_NS) {return false;}
  const size_t sl = static_cast<size_t>(pass)*n1 + (i - indcs.is);
  const int64_t *soA = &EA.so[sl*(np + 1)];
  const int64_t *roA = &EA.ro[sl*(np + 1)];
  const int64_t *soB = &EB.so[sl*(np + 1)];
  const int64_t *roB = &EB.ro[sl*(np + 1)];
  const int64_t s0A = soA[0], nsA = soA[np] - soA[0];
  const int64_t r0A = roA[0], nrA = roA[np] - roA[0];
  const int64_t s0B = soB[0], nsB = soB[np] - soB[0];
  const int64_t r0B = roB[0], nrB = roB[np] - roB[0];
  if (nsA + nsB >= static_cast<int64_t>(std::numeric_limits<int>::max()) ||
      nrA + nrB >= static_cast<int64_t>(std::numeric_limits<int>::max())) {
    return false;
  }
  // buffers hold both sets (grown once if needed)
  if (static_cast<int64_t>(vgd_cssb.extent(0)) < nsA + nsB) {
    Kokkos::realloc(vgd_cssb, static_cast<size_t>(nsA + nsB));
  }
  if (static_cast<int64_t>(vgd_csrb.extent(0)) < nrA + nrB) {
    Kokkos::realloc(vgd_csrb, static_cast<size_t>(nrA + nrB));
  }
  // per-set list offsets relative to the shell's first value, and the merged offsets
  VgdCsXIdx sa{}, sb{}, ra{}, rb{}, sm{}, rm{};
  for (int p = 0; p <= np; ++p) {
    sa.v[p] = static_cast<int>(soA[p] - s0A);
    sb.v[p] = static_cast<int>(soB[p] - s0B);
    ra.v[p] = static_cast<int>(roA[p] - r0A);
    rb.v[p] = static_cast<int>(roB[p] - r0B);
    sm.v[p] = sa.v[p] + sb.v[p];
    rm.v[p] = ra.v[p] + rb.v[p];
  }
  auto a_ = vgd_i;
  auto b_ = vgd_i_xb;
  auto S_ = vgd_cssl;
  auto R_ = vgd_csrl;
  auto siA_ = EA.si;
  auto riA_ = EA.ri;
  auto siB_ = EB.si;
  auto riB_ = EB.ri;
  Real *sbp = vgd_cssb.data();
  Real *rbp = vgd_csrb.data();
  const int nsa = static_cast<int>(nsA), nsb = static_cast<int>(nsB);
  const int nra = static_cast<int>(nrA), nrb = static_cast<int>(nrB);
  const int npk = np;
  if (nsa + nsb > 0) {
    Kokkos::parallel_for("vet_vgdcs_cpackx", VgdLight(DevExeSpace(), 0, nsa + nsb),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nsa);
      const int o = b ? (c - nsa) : c;
      const VgdCsXIdx &sx = b ? sb : sa;
      int p = 0;
      while (p + 1 < npk && sx.v[p+1] <= o) {++p;}
      const int x = b ? siB_(s0B + o) : siA_(s0A + o);
      const int e = x/n, d = x - e*n;
      const int pos = sm.v[p] + (b ? (sa.v[p+1] - sa.v[p]) : 0) + (o - sx.v[p]);
      sbp[pos] = b ? b_(S_(e,0),d,S_(e,1),S_(e,2),i) : a_(S_(e,0),d,S_(e,1),S_(e,2),i);
    });
  }
  Kokkos::fence();
  Kokkos::Timer tq;
  std::vector<MPI_Request> req;
  req.reserve(2*np);
  for (int p = 0; p < np; ++p) {
    const int rk = vgd_csprk[p];
    const int cr = rm.v[p+1] - rm.v[p], cs = sm.v[p+1] - sm.v[p];
    if (cr > 0) {
      req.emplace_back();
      MPI_Irecv(rbp + rm.v[p], cr, MPI_ATHENA_REAL, rk, 7107, MPI_COMM_WORLD,
                &req.back());
    }
    if (cs > 0) {
      req.emplace_back();
      MPI_Isend(sbp + sm.v[p], cs, MPI_ATHENA_REAL, rk, 7107, MPI_COMM_WORLD,
                &req.back());
      vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  vgd_tmpi += tq.seconds();
  vgd_nexch += 2.0;
  vgd_csvsent += static_cast<double>(nsA + nsB);
  vgd_csvdense += static_cast<double>(vgd_csnsnd)*n*2;
  if (nra + nrb > 0) {
    Kokkos::parallel_for("vet_vgdcs_cunpackx", VgdLight(DevExeSpace(), 0, nra + nrb),
    KOKKOS_LAMBDA(const int c) {
      const bool b = (c >= nra);
      const int o = b ? (c - nra) : c;
      const VgdCsXIdx &rx = b ? rb : ra;
      int p = 0;
      while (p + 1 < npk && rx.v[p+1] <= o) {++p;}
      const int x = b ? riB_(r0B + o) : riA_(r0A + o);
      const int e = x/n, d = x - e*n;
      const Real v = rbp[rm.v[p] + (b ? (ra.v[p+1] - ra.v[p]) : 0) + (o - rx.v[p])];
      if (b) {
        b_(R_(e,0),d,R_(e,1),R_(e,2),i) = v;
      } else {
        a_(R_(e,0),d,R_(e,1),R_(e,2),i) = v;
      }
    });
  }
  return true;
#else
  (void) pass;
  (void) i;
  (void) eb;
  return false;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsHalo
//! \brief fill every band position of a(m, v, k, j, i), v < nv, shells [i0, i1] with its
//! owner's interior value (exact copies, chunked). ws >= 0 (vet_gd_band_global,
//! depth-sorted lists): only positions within depth ws; ws < 0 (cubed sphere): all

template <class V>
void RadiationVet::VetGdCsHalo(V &a, const int nv, const int i0, const int i1,
                               const bool loc, const int ws) {
  auto a_ = a;
  const bool dp = (ws >= 0) && (vgd_gband || vgd_raw) && !vgd_gdnl.empty();
  const int wq = dp ? std::min(ws, static_cast<int>(vgd_gdnl.size()) - 1) : 0;
  const int nl = loc ? (dp ? vgd_gdnl[wq] : vgd_csnloc) : 0;
  {
    const int ni = i1 - i0 + 1;
    if (nl > 0 && ni > 0) {
      auto L_ = vgd_csloc;
      par_for("vet_vgdcs_loc", DevExeSpace(), 0, nl - 1, 0, nv - 1, 0, ni - 1,
      KOKKOS_LAMBDA(const int e, const int v, const int ii) {
        a_(L_(e,0),v,L_(e,1),L_(e,2),i0+ii) = a_(L_(e,3),v,L_(e,4),L_(e,5),i0+ii);
      });
    }
  }
#if MPI_PARALLEL_ENABLED
  if (vgd_hmpi) {
    const int nent = std::max(std::max(vgd_csnsnd, vgd_csnrcv), 1);
    // full-sphere global band: device buffers (CUDA-aware MPI); cubed sphere: pinned host
    const int cap = static_cast<int>((dp ? vgd_sbuf.extent(0) : vgd_cssb.extent(0))/nent);
    const int nich = std::max(1, cap/nv);
    for (int c0 = i0; c0 <= i1; c0 += nich) {
      const int c1 = std::min(c0 + nich - 1, i1);
      const int ni = c1 - c0 + 1;
      const int nvi = nv*ni;
      auto sb_ = vgd_cssb;
      auto rb_ = vgd_csrb;
      auto S_ = vgd_cssl;
      auto R_ = vgd_csrl;
      auto DS_ = vgd_gdps;
      auto DR_ = vgd_gdpr;
      auto gsb_ = vgd_sbuf;
      auto grb_ = vgd_rbuf;
      Real *sbp = dp ? gsb_.data() : sb_.data();
      Real *rbp = dp ? grb_.data() : rb_.data();
      if (vgd_csnsnd > 0 && dp) {
        // depth-ws prefix of every partner segment (same offsets)
        par_for("vet_vgdgb_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DS_(e) > wq) {return;}
          gsb_(static_cast<size_t>(e)*nvi + v*ni + ii) =
              a_(S_(e,0),v,S_(e,1),S_(e,2),c0+ii);
        });
      } else if (vgd_csnsnd > 0) {
        par_for("vet_vgdcs_pack", DevExeSpace(), 0, vgd_csnsnd - 1, 0, nv - 1, 0, ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          sb_(static_cast<size_t>(e)*nvi + v*ni + ii) =
              a_(S_(e,0),v,S_(e,1),S_(e,2),c0+ii);
        });
      }
      Kokkos::fence();
      Kokkos::Timer tq;
      std::vector<MPI_Request> req;
      for (size_t p = 0; p < vgd_csprk.size(); ++p) {
        const int rk = vgd_csprk[p];
        const int cm = dp ? vgd_gdrn[wq][p] : vgd_cspm[p];
        const int cn = dp ? vgd_gdsn[wq][p] : vgd_cspn[p];
        if (cm > 0) {
          req.emplace_back();
          MPI_Irecv(rbp + static_cast<size_t>(vgd_cspr[p])*nvi, cm*nvi,
                    MPI_ATHENA_REAL, rk, 7101, MPI_COMM_WORLD, &req.back());
        }
        if (cn > 0) {
          req.emplace_back();
          MPI_Isend(sbp + static_cast<size_t>(vgd_csps[p])*nvi, cn*nvi,
                    MPI_ATHENA_REAL, rk, 7101, MPI_COMM_WORLD, &req.back());
          vgd_nbyte += static_cast<Real>(cn)*nvi*sizeof(Real);
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
      vgd_tmpi += tq.seconds();
      vgd_nexch += 1.0;
      if (!loc) {
        vgd_csvsent += static_cast<double>(vgd_csnsnd)*nvi;
        vgd_csvdense += static_cast<double>(vgd_csnsnd)*nvi;
      }
      if (vgd_csnrcv > 0 && dp) {
        par_for("vet_vgdgb_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0,
                ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          if (DR_(e) > wq) {return;}
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              grb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      } else if (vgd_csnrcv > 0) {
        par_for("vet_vgdcs_unpack", DevExeSpace(), 0, vgd_csnrcv - 1, 0, nv - 1, 0,
                ni - 1,
        KOKKOS_LAMBDA(const int e, const int v, const int ii) {
          a_(R_(e,0),v,R_(e,1),R_(e,2),c0+ii) =
              rb_(static_cast<size_t>(e)*nvi + v*ni + ii);
        });
      }
    }
  }
#endif
}

template void RadiationVet::VetGdCsHalo<DvceArray5D<Real>>(DvceArray5D<Real> &a,
                                                         const int nv, const int i0,
                                                         const int i1, const bool loc,
                                                         const int ws);

template void RadiationVet::VetGdCsHalo<RadiationVet::VgdRag>(RadiationVet::VgdRag &a,
                                                              const int nv, const int i0,
                                                              const int i1,
                                                              const bool loc,
                                                              const int ws);

}  // namespace radvet
