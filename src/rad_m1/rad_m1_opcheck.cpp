//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_opcheck.cpp
//! \brief implicit M1: the operator-equivalence check (tests_m1/gates).
//!
//! <rad_m1>/implicit_op_check = K (default 0 = off, read only when named).  At the first
//! |K| implicit solves (the top of ImplicitBiCGStab), before the solve touches its
//! vectors, every way the code has of forming y = A x for the CURRENT configuration is
//! applied to the same pseudo-random x (a fixed hash of gid, k, j, i: no state, the same
//! on every rank count), and compared with the first one:
//!   * halo paths: the ordinary exchange (pbval_kr), implicit_halo_direct (every
//!     neighbour on the rank), implicit_halo_mpi (on-rank copy + one message per rank);
//!   * operators: the stored stencil (ImplicitStencilOp, with and without
//!     implicit_op_split_red), the overlapped interior + shell (ImplicitStencilOpPart,
//!     with implicit_halo_ovl_faces off and on), the legacy 7-point row + M1VimpRow +
//!     off-diagonal kernels (with and without implicit_od_cache), and the device-scalar
//!     operator of implicit_krylov_dev (halo read in place, and halo kernel);
//!   * under implicit_vimp the stencil is also rebuilt with implicit_vimp_fold flipped,
//!     and the stencil, overlap and device variants are run again on it.
//! Every variant that the configuration allows is run whatever the input selects.  The
//! reductions of red = 3 ((y,s), (y,y), (rhat,y)) are compared too, where the variant
//! forms them.  Before the operator, the ghost zones of x are POISONED (1e20): the ghosts
//! each halo path fills are checked against the neighbour's active cells (exact), and a
//! reference y beyond 1e12 means an unfilled ghost is read with a non-zero coefficient.
//! A variant FAILS when max|y - y_ref|/max|y_ref| or a reduction differs by more than
//! implicit_op_check_tol (default 1e-12).  K > 0: fatal after the report of a failing
//! solve; K < 0: report only.  The whole work array (and the stored stencil, kd) is
//! restored afterwards, so the solve itself is unchanged (bitwise).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/nghbr_index.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

namespace {
//! a fixed pseudo-random number in [-1, 1) of (gid, k, j, i, salt): splitmix64
KOKKOS_INLINE_FUNCTION
Real M1ChkHash(const int gid, const int k, const int j, const int i, const int salt) {
  uint64_t z = static_cast<uint64_t>(gid) * 0x100000001B3ULL;
  z = (z + static_cast<uint64_t>(k + 64)) * 0x9E3779B97F4A7C15ULL;
  z = (z + static_cast<uint64_t>(j + 64)) * 0xBF58476D1CE4E5B9ULL;
  z = (z + static_cast<uint64_t>(i + 64)) * 0x94D049BB133111EBULL;
  z += static_cast<uint64_t>(salt) * 0xD6E8FEB86659FD93ULL + 0x9E3779B97F4A7C15ULL;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z = z ^ (z >> 31);
  return static_cast<Real>(z >> 11) * (2.0/9007199254740992.0) - 1.0;
}
constexpr Real M1CHK_POISON = 1.0e20;

//! global max over ranks (in place)
void M1ChkMax(Real *v, int n) {
#if MPI_PARALLEL_ENABLED
  std::vector<Real> g(n);
  MPI_Allreduce(v, g.data(), n, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  for (int q = 0; q < n; ++q) {v[q] = g[q];}
#else
  (void)v;
  (void)n;
#endif
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitOpCheck
//! \brief implicit_op_check (see the file header)

void RadiationM1::ImplicitOpCheck() {
  const int nchk = (impl_opchk > 0) ? impl_opchk : -impl_opchk;
  if (opchk_n >= nchk) {return;}
  ++opchk_n;
  auto *pm = pmy_pack->pmesh;
  if (pm->multilevel || pm->use_cubed_sphere || pm->use_polar_boundary) {
    if (global_variable::my_rank == 0) {
      std::cout << "M1OPCHK solve " << opchk_n << ": SKIPPED (multilevel, cubed sphere "
                << "or polar mesh)" << std::endl;
    }
    return;
  }
  DevExeSpace().fence();
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int ng = indcs.ng;
  const int n1 = nx1 + 2*ng;
  const int n2 = (nx2 > 1) ? (nx2 + 2*ng) : 1;
  const int n3 = (nx3 > 1) ? (nx3 + 2*ng) : 1;
  const bool md = (nx2 > 1), td = (nx3 > 1);
  const int nmb = pmy_pack->nmb_thispack;
  const int nmb1 = nmb - 1;
  const int gid0 = pmy_pack->gids;
  const int reach = vimp_now ? 2 : 1;
  const Real tol = impl_opchk_tol;
  const bool me0 = (global_variable::my_rank == 0);

  // ---- save everything the variants may touch
  DvceArray5D<Real> iw_bak("m1_chk_iw", iw.extent(0), iw.extent(1), iw.extent(2),
                           iw.extent(3), iw.extent(4));
  Kokkos::deep_copy(iw_bak, iw);
  const bool s_stencil = impl_stencil, s_odc = impl_odc, s_vfold = impl_vfold;
  const bool s_ovl = impl_halo_ovl;
  const bool s_faces = impl_ovl_faces;   // OVLF
  const bool s_hdon = halo_direct_on, s_hmpi = impl_halo_mpi, s_split = impl_opsplit;
  const bool s_edges = st_edges;
  const int s_kdev = impl_kdev, s_kdevh = impl_kdev_halo;
  auto s_ost = ost;
  auto s_kdv = kdv;
  if (impl_stencil) {ImplicitStencilBuild();}   // the operator of this pass (idempotent)

  // ---- the neighbour gid of each block in each of the 27 directions (same level)
  DualArray2D<int> nbg("m1_chk_nbg", nmb, 27);
  {
    auto &nb = pmy_pack->pmb->nghbr;
    for (int m = 0; m < nmb; ++m) {
      for (int d = 0; d < 27; ++d) {nbg.h_view(m,d) = -1;}
      for (int o3 = (td ? -1 : 0); o3 <= (td ? 1 : 0); ++o3) {
        for (int o2 = (md ? -1 : 0); o2 <= (md ? 1 : 0); ++o2) {
          for (int o1 = -1; o1 <= 1; ++o1) {
            if (o1 == 0 && o2 == 0 && o3 == 0) continue;
            int n = NeighborIndex(o1, o2, o3, 0, 0);
            if (n < 0 || n >= pmy_pack->pmb->nnghbr) continue;
            const NeighborBlock &q = nb.h_view(m,n);
            if (q.gid >= 0) {nbg.h_view(m, (o1+1) + 3*(o2+1) + 9*(o3+1)) = q.gid;}
          }
        }
      }
    }
    nbg.modify_host();
    nbg.sync_device();
  }
  auto nbg_ = nbg.d_view;

  constexpr int XC = M1_IW_KY, YC = M1_IW_KV;
  auto iw_ = iw;
  // x: the hash in the active cells, the poison in every ghost; s and rhat: hashes; y:
  // the poison (a cell a variant does not write shows up as a huge difference)
  auto fill = [&]() {
    auto a = iw;
    par_for("m1_chk_fill", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie) && (!md || (j >= js && j <= je)) &&
                       (!td || (k >= ks && k <= ke));
      const int g = gid0 + m;
      a(m,XC,k,j,i) = act ? M1ChkHash(g, k, j, i, 1) : M1CHK_POISON;
      a(m,M1_IW_KS,k,j,i) = M1ChkHash(g, k, j, i, 2);
      a(m,M1_IW_KRH,k,j,i) = M1ChkHash(g, k, j, i, 3);
      a(m,YC,k,j,i) = M1CHK_POISON;
    });
  };
  // the ghosts of x within `reach` of the active box whose owner is a neighbour block:
  // they must equal that block's active cell exactly.  Returns (checked, wrong).
  auto halo_check = [&](Real *res) {
    auto a = iw;
    const int rch = reach;
    Real nchecked = 0.0, nbad = 0.0;
    const int nkji = n3*n2*n1;
    Kokkos::parallel_reduce("m1_chk_halo",
    Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &nc, Real &nw) {
      const int m = idx/nkji;
      int r = idx - m*nkji;
      const int k = r/(n2*n1);
      r -= k*n2*n1;
      const int j = r/n1;
      const int i = r - j*n1;
      const int o1 = (i < is) ? -1 : ((i > ie) ? 1 : 0);
      const int o2 = md ? ((j < js) ? -1 : ((j > je) ? 1 : 0)) : 0;
      const int o3 = td ? ((k < ks) ? -1 : ((k > ke) ? 1 : 0)) : 0;
      if (o1 == 0 && o2 == 0 && o3 == 0) return;
      if (i < is - rch || i > ie + rch) return;
      if (md && (j < js - rch || j > je + rch)) return;
      if (td && (k < ks - rch || k > ke + rch)) return;
      const int g = nbg_(m, (o1+1) + 3*(o2+1) + 9*(o3+1));
      if (g < 0) return;
      nc += 1.0;
      const Real want = M1ChkHash(g, k - o3*nx3, j - o2*nx2, i - o1*nx1, 1);
      if (a(m,XC,k,j,i) != want) {nw += 1.0;}
    }, nchecked, nbad);
    res[0] = nchecked;
    res[1] = nbad;
  };

  // the result of the reference and of each variant
  DvceArray4D<Real> yref("m1_chk_yref", nmb, n3, n2, n1);
  // the scale of each row: sum_o |a_o x_o| over the stored stencil of the production
  // configuration (1 where there is none).  A difference is measured against the row it
  // is in, so a wrong term that is small against max|y| (a vimp coefficient against the
  // transport diagonal) is not hidden by the rows where |y| is large; round-off gives a
  // few 1e-16 in this measure.
  DvceArray4D<Real> rabs("m1_chk_rabs", nmb, n3, n2, n1);
  Kokkos::deep_copy(rabs, 1.0);
  auto rowabs = [&]() {
    auto a = iw;
    auto st_ = ost;
    auto ra = rabs;
    const bool edg = st_edges, thrd = trans_x3, vf = impl_vfold && vimp_now;
    par_for("m1_chk_rabs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real r = fabs(st_(m,0,k,j,i)*a(m,XC,k,j,i))
               + fabs(st_(m,1,k,j,i)*a(m,XC,k,j,i-1))
               + fabs(st_(m,2,k,j,i)*a(m,XC,k,j,i+1))
               + fabs(st_(m,3,k,j,i)*a(m,XC,k,j-1,i))
               + fabs(st_(m,4,k,j,i)*a(m,XC,k,j+1,i));
      if (edg) {
        r += fabs(st_(m,7,k,j,i)*a(m,XC,k,j-1,i-1))
             + fabs(st_(m,8,k,j,i)*a(m,XC,k,j-1,i+1))
             + fabs(st_(m,9,k,j,i)*a(m,XC,k,j+1,i-1))
             + fabs(st_(m,10,k,j,i)*a(m,XC,k,j+1,i+1));
      }
      if (thrd) {
        r += fabs(st_(m,5,k,j,i)*a(m,XC,k-1,j,i)) + fabs(st_(m,6,k,j,i)*a(m,XC,k+1,j,i));
        if (edg) {
          r += fabs(st_(m,11,k,j,i)*a(m,XC,k-1,j,i-1))
               + fabs(st_(m,12,k,j,i)*a(m,XC,k-1,j,i+1))
               + fabs(st_(m,13,k,j,i)*a(m,XC,k+1,j,i-1))
               + fabs(st_(m,14,k,j,i)*a(m,XC,k+1,j,i+1))
               + fabs(st_(m,15,k,j,i)*a(m,XC,k-1,j-1,i))
               + fabs(st_(m,16,k,j,i)*a(m,XC,k-1,j+1,i))
               + fabs(st_(m,17,k,j,i)*a(m,XC,k+1,j-1,i))
               + fabs(st_(m,18,k,j,i)*a(m,XC,k+1,j+1,i));
        }
      }
      if (vf) {
        r += fabs(st_(m,19,k,j,i)*a(m,XC,k,j,i-2)) + fabs(st_(m,20,k,j,i)*a(m,XC,k,j,i+2))
             + fabs(st_(m,21,k,j,i)*a(m,XC,k,j-2,i))
             + fabs(st_(m,22,k,j,i)*a(m,XC,k,j+2,i));
        if (thrd) {
          r += fabs(st_(m,23,k,j,i)*a(m,XC,k-2,j,i))
               + fabs(st_(m,24,k,j,i)*a(m,XC,k+2,j,i));
        }
      }
      // a poisoned ghost under a zero coefficient gives 0*1e20 = 0; under a non-zero
      // one the reference fails anyway
      ra(m,k,j,i) = (r > 0.0 && r < 1.0e12) ? r : 1.0;
    });
  };
  Real sref[3] = {0.0, 0.0, 0.0};
  bool sref_ok = false;
  Real yscale = 0.0, sscale = 1.0;
  int nfail = 0, nrun = 0;
  std::string refname;
  std::vector<std::string> lines;

  // y := yref (first call) or the differences to it
  auto compare = [&](const std::string &name, bool has_red, const Real *out) {
    DevExeSpace().fence();
    auto a = iw;
    auto yr = yref;
    const bool first = refname.empty();
    if (first && impl_stencil) {rowabs();}
    auto ra = rabs;
    Real dmax = 0.0, ymax = 0.0, nbad = 0.0;
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    Kokkos::parallel_reduce("m1_chk_cmp",
    Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &ld, Real &ly, Real &lb) {
      const int m = idx/nkji;
      int r = idx - m*nkji;
      const int k = ks + r/nji;
      r -= (k - ks)*nji;
      const int j = js + r/ni;
      const int i = is + (r - (j - js)*ni);
      const Real y = a(m,YC,k,j,i);
      if (first) {yr(m,k,j,i) = y;}
      if (!(fabs(y) < 1.0e12)) {lb += 1.0;}   // NaN, inf or a poisoned ghost read
      ly = fmax(ly, fabs(y));
      ld = fmax(ld, fabs(y - yr(m,k,j,i))/ra(m,k,j,i));
    }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(ymax), nbad);
    Real v[4] = {dmax, ymax, 0.0, 0.0};
    if (has_red && sref_ok) {
      v[2] = fmax(fabs(out[0] - sref[0]), fabs(out[2] - sref[2]))/sscale;
      v[3] = fabs(out[1] - sref[1])/fmax(sref[1], 1.0e-300);
    }
#if MPI_PARALLEL_ENABLED
    {Real g = nbad;
    MPI_Allreduce(&nbad, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    nbad = g;}
#endif
    M1ChkMax(v, 4);
    char buf[256];
    if (first) {
      refname = name;
      yscale = fmax(v[1], 1.0e-300);
      if (has_red) {
        for (int q = 0; q < 3; ++q) {sref[q] = out[q];}
        sref_ok = true;
        // |(y,s)| and |(rhat,y)| <= sqrt((y,y) n) (|s|, |rhat| < 1), per rank
        sscale = std::sqrt(fmax(sref[1], 0.0)*nmb*nx1*nx2*nx3) + 1.0e-300;
      }
      const bool bad = (nbad > 0.0);
      if (bad) {++nfail;}
      std::snprintf(buf, sizeof(buf), "  ref  %-34s max|y| %.3e  unfilled-ghost reads "
                    "%.0f  %s", name.c_str(), v[1], nbad, bad ? "FAIL" : "ok");
      lines.push_back(buf);
      return;
    }
    ++nrun;
    const Real dy = s_stencil ? v[0] : v[0]/yscale;   // per row, or against max|y|
    const bool bad = !(dy <= tol) || (nbad > 0.0) || (has_red && sref_ok &&
                                                        !(v[2] <= tol && v[3] <= tol));
    if (bad) {++nfail;}
    if (has_red && sref_ok) {
      std::snprintf(buf, sizeof(buf), "  var  %-34s dy/row %.2e  d(y,s)/(rhat,y) %.2e  "
                    "d(y,y) %.2e  %s", name.c_str(), dy, v[2], v[3],
                    bad ? "FAIL" : "PASS");
    } else {
      std::snprintf(buf, sizeof(buf), "  var  %-34s dy/row %.2e  %s", name.c_str(), dy,
                    bad ? "FAIL" : "PASS");
    }
    lines.push_back(buf);
  };

  // ---- the halo paths available here
  struct HPath {const char *name; bool hdon, hmpi;};
  std::vector<HPath> hp;
  const bool gen_ok = (pbval_kr != nullptr);
  const bool dir_ok = s_hdon;
  bool mpi_ok = false;
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1 && impl_halo_direct && hd_src.extent_int(0) == nmb) {
    if (hm_state == 0) {ImplicitHaloMPIInit();}
    mpi_ok = (hm_state == 1);
  }
#endif
  // the production path first: it gives the reference
  if (s_hdon) {
    hp.push_back({"direct", true, false});
  } else if (s_hmpi && mpi_ok) {
    hp.push_back({"mpi", false, true});
  } else if (gen_ok) {
    hp.push_back({"exch", false, false});
  }
  if (gen_ok && (s_hdon || (s_hmpi && mpi_ok))) {hp.push_back({"exch", false, false});}
  if (dir_ok && !s_hdon) {hp.push_back({"direct", true, false});}
  if (mpi_ok && !(s_hmpi && !s_hdon)) {hp.push_back({"mpi", false, true});}

  // the overlap needs an interior box in every non-degenerate direction
  const int w = reach;
  const bool ovl_geom = (nx1 > 2*w) && (!md || (nx2 > 2*w)) && (!td || (nx3 > 2*w));
  bool devok = false;
  if (s_hdon && global_variable::nranks == 1 && impl_stencil) {
    impl_kdev = 1;
    impl_vfold = true;   // (only asks whether the rest of the conditions hold)
    devok = ImplicitKrylovDevOK();
    impl_vfold = s_vfold;
    impl_kdev = s_kdev;
  }
  DvceArray1D<Real> kd_tmp("m1_chk_kd", M1_KD_SIZE);

  // one pass over the operator variants for the current fold state
  auto run_ops = [&](const std::string &tag) {
    Real out[4];
    for (const HPath &h : hp) {
      halo_direct_on = h.hdon;
      impl_halo_mpi = h.hmpi;
      const std::string hn = std::string(h.name) + tag;
      impl_halo_ovl = false;
      impl_opsplit = false;
      // the halo of this path, checked against the neighbours
      fill();
      ImplicitKrylovHalo(XC);
      if (tag.empty()) {
        Real hr[2];
        halo_check(hr);
        Real hv[2] = {hr[0], hr[1]};
#if MPI_PARALLEL_ENABLED
        MPI_Allreduce(hr, hv, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
        const bool bad = (hv[1] > 0.0);
        if (bad) {++nfail;}
        char buf[256];
        std::snprintf(buf, sizeof(buf), "  halo %-34s ghosts checked %.0f  wrong %.0f"
                      "  %s", hn.c_str(), hv[0], hv[1], bad ? "FAIL" : "PASS");
        lines.push_back(buf);
      }
      if (impl_stencil) {
        ImplicitStencilOp(XC, YC, 3, out);
        compare(hn + "/stencil", true, out);
        impl_opsplit = true;
        fill();
        ImplicitKrylovHalo(XC);
        ImplicitStencilOp(XC, YC, 3, out);
        compare(hn + "/stencil_split_red", true, out);
        impl_opsplit = false;
        if (h.hmpi && ovl_geom) {
          impl_halo_ovl = true;
          impl_ovl_faces = false;   // OVLF
          fill();
          ImplicitHaloOp(XC, YC, 3, out);
          compare(hn + "/overlap", true, out);
          impl_ovl_faces = true;                                      // OVLF
          fill();                                                     // OVLF
          ImplicitHaloOp(XC, YC, 3, out);                             // OVLF
          compare(hn + "/overlap_faces", true, out);                  // OVLF
          impl_ovl_faces = s_faces;                                   // OVLF
          impl_halo_ovl = false;
        }
        if (h.hdon && devok && (impl_vfold || !vimp_now)) {
          for (int allin = 0; allin <= 1; ++allin) {
            kdv = kd_tmp;
            Kokkos::deep_copy(kdv, 0.0);
            impl_kdev_halo = allin;
            fill();   // allin = 0: the ghosts stay poisoned, the operator reads in place
            ImplicitOpXD(XC, YC, 3);
            DevExeSpace().fence();
            auto kh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), kdv);
            // M1DVal (s0, s1, s2, mx) sits at KD_N = 16 (rad_m1_launch.cpp)
            for (int q = 0; q < 3; ++q) {out[q] = kh(16 + q);}
            compare(hn + (allin ? "/krylov_dev_halo" : "/krylov_dev"), true, out);
            impl_kdev_halo = s_kdevh;
            kdv = s_kdv;
          }
        }
      }
      if (tag.empty()) {
        // the legacy operators (no stored stencil): 7-point row + M1VimpRow + od terms
        impl_stencil = false;
        if (impl_odc && odc.extent_int(0) == nmb) {
          fill();
          ImplicitKrylovHalo(XC);
          ImplicitOffDiagOpC(XC, YC, 1.0, true, 3, out);
          compare(hn + "/legacy_od_cache", true, out);
        }
        impl_odc = false;
        fill();
        ImplicitApplyOp(XC, YC);
        compare(hn + "/legacy_7pt", false, out);
        impl_odc = s_odc;
        impl_stencil = s_stencil;
      }
    }
    halo_direct_on = s_hdon;
    impl_halo_mpi = s_hmpi;
  };

  run_ops("");
  if (impl_stencil && vimp_now) {
    // the other fold state, on a stencil array of its own
    impl_vfold = !s_vfold;
    ost = DvceArray5D<Real>("m1_chk_ost", nmb, 25, n3, n2, n1);
    ImplicitStencilBuild();
    run_ops(impl_vfold ? "+fold" : "+nofold");
    impl_vfold = s_vfold;
    ost = s_ost;
  }

  // ---- restore
  DevExeSpace().fence();
  Kokkos::deep_copy(iw, iw_bak);
  impl_stencil = s_stencil; impl_odc = s_odc; impl_vfold = s_vfold;
  impl_halo_ovl = s_ovl;
  impl_ovl_faces = s_faces;   // OVLF
  halo_direct_on = s_hdon; impl_halo_mpi = s_hmpi; impl_opsplit = s_split;
  st_edges = s_edges;
  impl_kdev = s_kdev; impl_kdev_halo = s_kdevh;
  ost = s_ost;
  kdv = s_kdv;
  if (impl_stencil) {ImplicitStencilBuild();}
  DevExeSpace().fence();

  if (me0) {
    std::cout << "M1OPCHK solve " << opchk_n << " (" << global_variable::nranks
              << " rank(s), vimp " << (vimp_now ? 1 : 0) << ", vimp_fold "
              << (s_vfold ? 1 : 0) << ", reach " << reach << ", tol " << tol
              << "): reference = " << refname << std::endl;
    for (const auto &l : lines) {std::cout << l << std::endl;}
    std::cout << "M1OPCHK solve " << opchk_n << ": " << (nfail ? "FAIL" : "PASS")
              << " (" << nrun << " variants, " << nfail << " failed)" << std::endl;
  }
  if (nfail > 0 && impl_opchk > 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << "<rad_m1>/"
              << "implicit_op_check: operator variants disagree (see M1OPCHK)"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

} // namespace radm1
