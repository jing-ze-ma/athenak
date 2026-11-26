//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cs_sweep.cpp
//! \brief <rad_vet>/vet_gd on the cubed sphere: the sp sweep of sc_sph_sweep.cpp
//! (shells, passes, segment and boundary rules) with the cs lateral stencil and band.
//! Band: the own panel's index space extended by vgd_w cells; beyond a seam the
//! neighbour's cells are stored unfolded (mirror image across the seam). No cell beyond
//! two edges or at depth >= N. Filled per shell by exact copies keyed on the global cell.
//! Read: the upwind point is located on its panel at fractional native indices:
//! bilinear inside the centre lattice; linear across a seam strip, each edge column
//! along the seam in its own chart; barycentric on the 3 cells at a cube vertex.
//! Weights >= 0, sum 1. Hand-over: D_rr - f_K into LAT0, the Cartesian tangential
//! vector D.rhat - D_rr rhat into LAT1..LAT3 (no tangential anisotropy on cs).
//!
//!  called from: rays/sc_sph_sweep.cpp
//!  calls into: rad_vet_params_closure.cpp, rays/sc_bandring.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cs_halo.cpp, rays/sc_setup.cpp

#include <math.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_cs.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetGdCsSweep
//! \brief one sweep (inward then outward) of every direction on the cubed sphere:
//! the VetGdSweep body with the cs upwind stencil

void RadiationVet::VetGdCsSweep() {
  if (vgd_xff_prep) {return;}   // nothing to prepare (lists made below)
  // vet_gd_rotate_xfade_fuse: the next set (own intensities vgd_i_xb, tables, lists)
  // is swept in the same shell kernels
  const bool xffs = vgd_xff_now;
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int wb = vgd_w;
  const int js = wb, je = wb + indcs.nx2 - 1;
  const int ks = wb, ke = wb + indcs.nx3 - 1;
  const int c2b = indcs.nx2 + 2*wb, c3b = indcs.nx3 + 2*wb;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vgd_scut;
  auto cs_ = vgd_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vray_cnt;
  auto mr_ = vgd_mr;
  auto cr_ = vgd_csr;
  auto bl_ = vgd_csblk;
  auto cm_ = vgd_cscm;
  auto nb_ = vgd_csnb;
  auto ow_ = vgd_csown;
  const int nc = vgd_csn;
  const Real h = 0.5*M_PI/nc;
  const bool bfc = vet_bc_face;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  const bool emi = VetBcEmit() && (lcut == 0) &&
                   (ibc_x1min != VET_IBC_REFLECT);
  const Real iemit = VetBathLo();
  const bool tbh = vgd_tbath && (ibc_x1max == VET_IBC_MARSHAK) && (iebath_x1max != 0.0);
  const Real itop = tbh ? iebath_x1max : 0.0;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  const Real lcr = log(kVetCornerRatio);
  vgd_hsweep += xffs ? 2 : 1;
  // compact lists of this set (built on first sweep); a remote read missing from the
  // read mask is fatal below
  const bool ccl = vgd_cscomp;
  int xcurb = -1;   // next set's compact-list entry
  if (ccl) {
    VetGdCsCompact();
    if (xffs) {
      // current set first, then the next (the order of two sweeps: same LRU entry)
      const int ca = vgd_cscur;
      std::swap(vgd_dir, vgd_tabd.dir);
      std::swap(vgd_mr, vgd_tabd.mr);
      std::swap(vgd_alpha, vgd_xf_alt);
      VetGdCsCompact();
      xcurb = vgd_cscur;
      std::swap(vgd_alpha, vgd_xf_alt);
      std::swap(vgd_mr, vgd_tabd.mr);
      std::swap(vgd_dir, vgd_tabd.dir);
      vgd_cscur = ca;
    }
    Kokkos::deep_copy(vgd_csmiss, 0);
  }
  auto msk_ = vgd_cscl[vgd_cscur].msk;
  auto rx_ = vgd_csrix;
  auto miss_ = vgd_csmiss;
  const int64_t per_ = vgd_csper;   // mask stride per entry
  // vet_gd_band_ring (needs compact lists): the sp ring; the saved inward band values of
  // shell i are the receive list of its inward compact exchange (record vgd_hxe[i];
  // all addresses in the unused entry vgd_hxe[n1c]). Pure copies (bitwise).
  const bool rga = vgd_ring && ccl;
  if (vgd_ring && !ccl) {VgdCsFatal("vet_gd_band_ring without the compact lists");}
  // next set's records in vgd_hxe_alt; the ring's second array is its intensities
  auto ringrec = [&](std::vector<VgdHlEntry> &hxe, const VgdCsCl &E, const Real alr) {
    const int n1c = indcs.nx1 + 2*indcs.ng;
    if (static_cast<int>(hxe.size()) != 2*n1c) {hxe.resize(2*n1c);}
    const int np = static_cast<int>(vgd_csprk.size());
    const int64_t a0 = E.ro[0];
    const int64_t a1 = E.ro[static_cast<size_t>(n1 - 1)*(np + 1) + np];
    const int64_t na = a1 - a0;
    if (static_cast<int64_t>(pmy_pack->nmb_thispack)*n*c3b*c2b >= INT32_MAX ||
        na >= INT32_MAX) {
      VgdCsFatal("vet_gd_band_ring: band addresses >= 2^31 (32-bit records)");
    }
    VgdHlEntry &B = hxe[n1c];
    B.gen = -1;
    B.exact = false;
    if (B.dl.extent_int(0) < na + 1) {
      B.dl = DvceArray1D<int>();
      Kokkos::realloc(B.dl, na + na/8 + 1024);
    }
    auto bd_ = B.dl;
    auto R_ = vgd_csrl;
    auto ri_ = E.ri;
    if (na > 0) {
      par_for("vet_vgdcs_ringrec", DevExeSpace(), 0, static_cast<int>(na) - 1,
      KOKKOS_LAMBDA(const int o) {
        const int x = ri_(a0 + o);
        const int e = x/n, d = x - e*n;
        bd_(o) = ((R_(e,0)*n + d)*c3b + R_(e,1))*c2b + R_(e,2);
      });
    }
    for (int l = 0; l < n1; ++l) {
      VgdHlEntry &X = hxe[is + l];
      const int64_t r0 = E.ro[static_cast<size_t>(l)*(np + 1)];
      const int64_t r1 = E.ro[static_cast<size_t>(l)*(np + 1) + np];
      X.dl = DvceArray1D<int>(bd_.data() + (r0 - a0), r1 - r0);
      X.dw = DvceArray1D<int>();
      X.nr = static_cast<int>(r1 - r0);
      X.nw = 0;
      X.exact = true;
      X.gen = 1;
      X.alpha = alr;
      X.scut = vgd_scut;
    }
  };
  if (rga) {
    ringrec(vgd_hxe, vgd_cscl[vgd_cscur], vgd_alpha);
    if (xffs) {ringrec(vgd_hxe_alt, vgd_cscl[xcurb], vgd_xf_alt);}
  }
  if (rga) {vgd_cur = DevExeSpace();}   // ring copy kernel on this sweep's stream
  // per pass: codes ((m k) j) d of the (cell, direction) pairs of that pass's branch
  // (inward: n.r_hat < 0), in flat order, made each sweep from vgd_mr; the shell kernel
  // runs over this list only
  const int nkl = ke - ks + 1, njl = je - js + 1;
  const int nlt = (nmb1 + 1)*nkl*njl*n;
  int nlist[2] = {0, 0};
  if (vgd_cslst.extent_int(0) < 2*nlt) {
    vgd_cslst = DvceArray1D<int>();
    Kokkos::realloc(vgd_cslst, 2*nlt);
  }
  auto L_ = vgd_cslst;
  for (int p = 0; p < 2; ++p) {
    const bool inp = (p == 0);
    const int off = p*nlt;
    int cnt = 0;
    Kokkos::parallel_scan("vet_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                                  nlt),
    KOKKOS_LAMBDA(const int c, int &u, const bool fin) {
      int t = c;
      const int d = t % n;
      t /= n;
      const int j = js + t % njl;
      t /= njl;
      const int k = ks + t % nkl;
      const int m = t/nkl;
      if (inp != (mr_(m,k,j,d) >= 0.0)) {
        if (fin) {L_(off + u) = c;}
        u += 1;
      }
    }, cnt);
    nlist[p] = cnt;
  }
  auto lst_ = vgd_cslst;
  // next set's pass lists (same scan over its own n . r_hat)
  int nlistb[2] = {0, 0};
  if (xffs) {
    if (vgd_cslst_xb.extent_int(0) < 2*nlt) {
      vgd_cslst_xb = DvceArray1D<int>();
      Kokkos::realloc(vgd_cslst_xb, 2*nlt);
    }
    auto L_ = vgd_cslst_xb;
    auto mrb_ = vgd_tabd.mr;
    for (int p = 0; p < 2; ++p) {
      const bool inp = (p == 0);
      const int off = p*nlt;
      int cnt = 0;
      Kokkos::parallel_scan("vet_vgdcs_list", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                                                                    nlt),
      KOKKOS_LAMBDA(const int c, int &u, const bool fin) {
        int t = c;
        const int d = t % n;
        t /= n;
        const int j = js + t % njl;
        t /= njl;
        const int k = ks + t % nkl;
        const int m = t/nkl;
        if (inp != (mrb_(m,k,j,d) >= 0.0)) {
          if (fin) {L_(off + u) = c;}
          u += 1;
        }
      }, cnt);
      nlistb[p] = cnt;
    }
  }
  auto lstb_ = vgd_cslst_xb;
  // cs early send (see the shell loop)
  bool early = ccl;
#if !MPI_PARALLEL_ENABLED
  early = false;
#endif
  if (xffs) {early = false;}   // early send: one set per shell kernel
#if defined(KOKKOS_ENABLE_HIP) && MPI_PARALLEL_ENABLED
  // HIP: > 2 ranks per GPU hang in the early send's event wait unless
  // GPU_MAX_HW_QUEUES <= 2 (same bits)
  const char *hq = std::getenv("GPU_MAX_HW_QUEUES");
  static bool s_hqchk = false;
  if (early && !s_hqchk && (hq == nullptr || std::atoi(hq) > 2)) {
    s_hqchk = true;
    char bus[32] = {0};
    (void) hipDeviceGetPCIBusId(bus, sizeof(bus), Kokkos::device_id());
    MPI_Comm nc;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nc);
    int nn = 1, sh = 0;
    MPI_Comm_size(nc, &nn);
    std::vector<char> all(32*static_cast<size_t>(nn));
    MPI_Allgather(bus, 32, MPI_CHAR, all.data(), 32, MPI_CHAR, nc);
    MPI_Comm_free(&nc);
    for (int q = 0; q < nn; ++q) {sh += (std::strncmp(bus, &all[32*q], 32) == 0);}
    MPI_Allreduce(MPI_IN_PLACE, &sh, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (sh > 2) {
      VgdCsFatal(std::to_string(sh) + " ranks share one GPU: the early send hangs "
                 "there (HIP hardware queues); set GPU_MAX_HW_QUEUES=2 (same bits) or "
                 "use at most 2 ranks per GPU");
    }
  }
#endif
  // early-send cells are computed once: vet_vgdcs_estamp stamps the (m, k, j, d) that
  // vet_vgdcs_eshell computes (stamp vgd_esid); the shell kernel skips them (bitwise)
  const int nkb = ke - ks + 1, njb = je - js + 1;
  if (early && (vgd_est[0].extent_int(0) < nmb1 + 1 || vgd_est[0].extent_int(1) < nkb ||
                vgd_est[0].extent_int(2) < njb || vgd_est[0].extent_int(3) < n)) {
    vgd_est[0] = DvceArray4D<int>("vgd_est", nmb1 + 1, nkb, njb, n);   // zeroed
    vgd_esid = 0;
  }
  auto est_ = vgd_est[0];
  int esid = 0;   // 0: no early cells (never a stamp)
  VgdCsEv eev, efork;
  // shell kernel LaunchBounds<128, kLb> (register cap): HIP 2 (min waves per EU; more
  // spills), CUDA 8
#if defined(KOKKOS_ENABLE_HIP)
  constexpr unsigned kLb = 2;
#else
  constexpr unsigned kLb = 8;
#endif
  VetGdRingSweepBegin(xffs);   // the other set's array second
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    const int nsh = n1 - lcut;
    for (int q = 0; q < nsh; ++q) {
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      if (rga) {
        // this shell's wide slab: inward from scratch, outward with its interior and
        // its inward-exchange band values
        VetGdRingFlush();
        if (vgd_rg_state[i] < 0) {
          VgdRingOp o;
          o.kind = 3;
          o.i = i;
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
        VetGdRingFlush();
      }
      // flat (m, k, j, d) range, d fastest; with the fused twin one thread computes the
      // geometry once and the values of both fields
      const int nkk = ke - ks + 1, njj = je - js + 1;
      const int nper = nlist[pass];
      const int lo = pass*nlt;
      // shell body for one (cell, direction): run over the pass's list and (early send)
      // first over the compact send list with cnt = false (counters counted once)
#define VET_VGDCS_SHELL_NAME shell_
#include "rad_vet/rays/sc_cs_shell.inc"
#undef VET_VGDCS_SHELL_NAME
      // early send: compute and pack this shell's send values first, record an event;
      // the full shell kernel overlaps the host exchange (finished at the halo point)
      if (early) {
        const VgdCsCl &E = vgd_cscl[vgd_cscur];
        const int np = static_cast<int>(vgd_csprk.size());
        const size_t sl = static_cast<size_t>(pass)*n1 + l;
        const int64_t *so = &E.so[sl*(np + 1)];
        const int64_t es0 = so[0], ens = so[np] - so[0];
        if (ens > 0) {
          auto si_ = E.si;
          auto S_ = vgd_cssl;
          vgd_esid = (vgd_esid >= std::numeric_limits<int>::max() - 1) ? 1
                                                                       : (vgd_esid + 1);
          if (vgd_esid == 1) {Kokkos::deep_copy(est_, 0);}   // first use or wrapped
          esid = vgd_esid;
          const int sid_ = esid;
          Kokkos::parallel_for("vet_vgdcs_estamp", VgdLight(DevExeSpace(), 0,
                                                            static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            est_(S_(e,0), S_(e,1) - ks, S_(e,2) - js, d) = sid_;
          });
          // send cells + pack on the second stream, concurrent with the shell kernel;
          // the sweep stream waits for them after it (efork / eev)
          VgdEs2Make();
          efork.Record();
          efork.StreamWait(vgd_es2);
          Kokkos::parallel_for("vet_vgdcs_eshell",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            shell_(S_(e,0), S_(e,1), S_(e,2), d, true);
          });
          auto a_ = vgd_i;
          Real *sbp = vgd_cssb.data();
          Kokkos::parallel_for("vet_vgdcs_cpack",
                               VgdLight(vgd_es2, 0, static_cast<int>(ens)),
          KOKKOS_LAMBDA(const int o) {
            const int x = si_(es0 + o);
            const int e = x/n, d = x - e*n;
            sbp[o] = a_(S_(e,0),d,S_(e,1),S_(e,2),i);
          });
        }
        eev.Record((esid != 0) ? vgd_es2 : DevExeSpace());
      }
      // fused two-set shell kernel: this pass's list of the current set, then the next
      // set's, one launch (next set's body: the include with its names)
#define VET_VGDCS_LAUNCHN(F, N) \
      Kokkos::parallel_for("vet_vgdcs_shell_xff", Kokkos::Experimental::require(        \
          Kokkos::RangePolicy<Kokkos::LaunchBounds<128, kLb>>(DevExeSpace(), 0, N),    \
          Kokkos::Experimental::WorkItemProperty::HintLightWeight), F);
      const int esid_ = esid;
      esid = 0;
      if (xffs) {
        auto xsb_ = [&]() {
          auto dir_ = vgd_tabd.dir;
          auto mr_ = vgd_tabd.mr;
          auto vi_ = vgd_i_xb;
          auto msk_ = vgd_cscl[(xcurb >= 0) ? xcurb : vgd_cscur].msk;
#define VET_VGDCS_SHELL_NAME shellx_
#include "rad_vet/rays/sc_cs_shell.inc"
#undef VET_VGDCS_SHELL_NAME
          return shellx_;
        }();
        const int nA = nlist[pass], nB = nlistb[pass];
        auto mainx_ = KOKKOS_LAMBDA(const int c) {
          const bool sb = (c >= nA);
          int t = sb ? lstb_(lo + c - nA) : lst_(lo + c);
          const int d = t % n;
          t /= n;
          const int j = js + t % njj;
          t /= njj;
          const int k = ks + t % nkk;
          const int m = t/nkk;
          if (sb) {
            xsb_(m, k, j, d, true);
          } else {
            shell_(m, k, j, d, true);
          }
        };
        if (nA + nB > 0) {VET_VGDCS_LAUNCHN(mainx_, nA + nB)}
      } else {
      auto mainf_ = KOKKOS_LAMBDA(const int c) {
        int t = lst_(lo + c);
        const int d = t % n;
        t /= n;
        const int j = js + t % njj;
        t /= njj;
        const int k = ks + t % nkk;
        const int m = t/nkk;
        if (esid_ != 0 && est_(m, k - ks, j - js, d) == esid_) {return;}
        shell_(m, k, j, d, true);
      };
      Kokkos::parallel_for("vet_vgdcs_shell", Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::LaunchBounds<128, kLb>>(DevExeSpace(), 0, nper),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight), mainf_);
      }   // !xffs
      if (esid_ != 0) {eev.StreamWait(DevExeSpace());}   // the send cells
#undef VET_VGDCS_LAUNCHN
      if (rga) {
        // previous shell read: save its inward band values (inward pass), interior back
        // to the compact slab; the next shell gets a slab after this exchange
        if (q > 0) {
          const int ip = inw ? (i + 1) : (i - 1);
          VgdRingOp o;
          if (inw) {
            o.kind = 1;
            o.i = ip;
            vgd_rg_pend.push_back(o);
          }
          o.kind = 2;
          o.i = ip;
          vgd_rg_pend.push_back(o);
        }
        if (q + 1 < nsh) {
          VgdRingOp o;
          o.kind = 3;
          o.i = inw ? (i - 1) : (i + 1);
          o.copyin = !inw;
          o.restore = !inw;
          vgd_rg_pend.push_back(o);
        }
      }
      Kokkos::Timer th;
      if (early) {
#if MPI_PARALLEL_ENABLED
        // early send: pack done (event); messages of VetGdCsHaloC, then the unpack
        // (stream-ordered after the shell kernel)
        const VgdCsCl &E = vgd_cscl[vgd_cscur];
        const int np = static_cast<int>(vgd_csprk.size());
        const size_t sl = static_cast<size_t>(pass)*n1 + l;
        const int64_t *so = &E.so[sl*(np + 1)];
        const int64_t *ro = &E.ro[sl*(np + 1)];
        const int64_t s0 = so[0], ns = so[np] - so[0];
        const int64_t r0 = ro[0], nr = ro[np] - ro[0];
        if (static_cast<int64_t>(vgd_cssb.extent(0)) < ns ||
            static_cast<int64_t>(vgd_csrb.extent(0)) < nr) {
          VgdCsFatal("cs early send: message buffer too small");
        }
        eev.Wait();
        Kokkos::Timer tq;
        std::vector<MPI_Request> req;
        for (int p = 0; p < np; ++p) {
          const int rk = vgd_csprk[p];
          const int64_t cr = ro[p+1] - ro[p], cs = so[p+1] - so[p];
          if (cr > 0) {
            req.emplace_back();
            MPI_Irecv(vgd_csrb.data() + (ro[p] - r0), static_cast<int>(cr),
                      MPI_ATHENA_REAL, rk, 7103, MPI_COMM_WORLD, &req.back());
          }
          if (cs > 0) {
            req.emplace_back();
            MPI_Isend(vgd_cssb.data() + (so[p] - s0), static_cast<int>(cs),
                      MPI_ATHENA_REAL, rk, 7103, MPI_COMM_WORLD, &req.back());
            vgd_nbyte += static_cast<Real>(cs)*sizeof(Real);
          }
        }
        MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
        vgd_tmpi += tq.seconds();
        vgd_nexch += 1.0;
        vgd_csvsent += static_cast<double>(ns);
        vgd_csvdense += static_cast<double>(vgd_csnsnd)*n;
        if (nr > 0) {
          auto a_ = vgd_i;
          auto R_ = vgd_csrl;
          auto ri_ = E.ri;
          const Real *rbp = vgd_csrb.data();
          Kokkos::parallel_for("vet_vgdcs_cunpack",
                               VgdLight(DevExeSpace(), 0, static_cast<int>(nr)),
          KOKKOS_LAMBDA(const int o) {
            const int x = ri_(r0 + o);
            const int e = x/n, d = x - e*n;
            a_(R_(e,0),d,R_(e,1),R_(e,2),i) = rbp[o];
          });
        }
#endif
      } else if (ccl) {
        // fused sweep: both sets in one round (VetGdCsHaloCX), else two exchanges
        if (!(xffs && VetGdCsHaloCX(pass, i, xcurb))) {
          VetGdCsHaloC(pass, i);          // read values only
          if (xffs) {
            // next set's exchange (own compact lists and array)
            const int ca = vgd_cscur;
            vgd_cscur = xcurb;
            std::swap(vgd_i, vgd_i_xb);
            VetGdCsHaloC(pass, i);
            std::swap(vgd_i, vgd_i_xb);
            vgd_cscur = ca;
          }
        }
      } else {
        VetGdCsHalo(vgd_i, n, i, i, false);  // remote band only (on-rank: read direct)
        if (xffs) {VetGdCsHalo(vgd_i_xb, n, i, i, false);}
      }
      if (rga && inw) {
        // values to save from this inward band: its compact receive list (vgd_hxe[i])
        vgd_rg_sk[i] = 1;
        vgd_rg_se[i] = i;
      }
      if (rga) {VetGdRingFlush();}
      vgd_thalo += th.seconds();
    }
    if (rga) {
      // last shell of the pass: saved (inward), back to its compact slab
      const int il = is + (inw ? lcut : (n1 - 1));
      if (vgd_rg_state[il] >= 0) {
        VgdRingOp o;
        if (inw) {
          o.kind = 1;
          o.i = il;
          vgd_rg_pend.push_back(o);
        }
        o.kind = 2;
        o.i = il;
        vgd_rg_pend.push_back(o);
      }
      VetGdRingFlush();
    }
  }
  VetGdRingSweepEnd();
  Kokkos::fence();
  if (ccl) {
    auto mh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_csmiss);
    if (mh(0) != 0) {
      VgdCsFatal("compact lists: " + std::to_string(mh(0)) + " remote band reads not "
                 "in the read mask (rank " + std::to_string(global_variable::my_rank) +
                 "); run with vet_gd_halo_compact = 0 (the dense band) and report");
    }
  }
}

}  // namespace radvet
