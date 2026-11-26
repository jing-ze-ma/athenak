#ifndef RAD_VET_RAD_VET_REDUCE_HPP_
#define RAD_VET_RAD_VET_REDUCE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_reduce.hpp
//! \brief global sums over the active cells, layout-invariant on request.

#include <array>
#include <utility>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/repro_sum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//! up to three sums and one max (VetLayReduce)
struct VetBcgVal {
  Real s0, s1, s2, mx;
};

constexpr int VET_DR_NT = 256;    // slots per team (a power of 2)

// block mode: level-1 slots per MeshBlock (fills the GPU)
constexpr int VET_LB_NS = 32768;

//----------------------------------------------------------------------------------------
//! \fn VetGlobalMax
inline void VetGlobalMax(Real &a) {
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&a, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  a = g;
#else
  (void)a;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn VetLayReduce
//! \brief <mesh>/layout_invariant: the VetBcgVal reduction of fn(idx, v) over idx in
//! [0, n) and all ranks with exact order-independent sums (repro::ReproSum).
//! Cell mode: every cell term deposited exactly (independent of layout and backend).
//! Block mode (default): each MeshBlock's (or tile's) partial in a fixed order
//! (VET_LB_NS strided slots, then a fixed pairwise tree), partials summed exactly on the
//! host and merged over ranks.  bpart / hpart: device / host scratch.

template <class F>
void VetLayReduce(const char *name, const int n, const int nkji, DvceArray1D<Real> &bpart,
                 Kokkos::View<Real*, Kokkos::HostSpace> &hpart, const F &fn,
                 VetBcgVal &out, const int *nxb = nullptr) {
  if (repro::block_mode) {
    // level 1: slot s of VET_LB_NS sums q = s, s + VET_LB_NS, ...; level 2: one team per
    // unit sums slots t, t + VET_DR_NT, ..., then the pairwise tree.  Tile mode (nxb):
    // the units are N x N (x2, x3) tiles (a MeshBlock-sized tile is bitwise the block)
    const bool tl = repro::tile_n > 0 && nxb != nullptr;
    const int t2 = (tl && nxb[1] > 1) ? repro::tile_n : 1;
    const int t3 = (tl && nxb[2] > 1) ? repro::tile_n : 1;
    const int n1 = tl ? nxb[0] : 1, n2 = tl ? nxb[1] : 1, nt2 = n2/t2;
    const int ntb = tl ? nt2*(nxb[2]/t3) : 1;
    const int nmb = (n/nkji)*ntb;   // (the units)
    const int ntc = n1*t2*t3;   // (tile mode: the tile's cells)
    const int ns = (tl && ntc < VET_LB_NS) ? ntc : VET_LB_NS;
    const int off2 = 4*nmb*ns;
    if (bpart.extent_int(0) < off2 + 4*nmb) {
      bpart = DvceArray1D<Real>(
          Kokkos::view_alloc("vet_lay_part", Kokkos::WithoutInitializing), off2 + 4*nmb);
    }
    if (hpart.extent_int(0) < 4*nmb) {
      hpart = Kokkos::View<Real*, Kokkos::HostSpace>("vet_lay_part_h", 4*nmb);
    }
    auto part = bpart;
    Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nmb*ns),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/ns, sl = idx - m*ns;
      VetBcgVal v;
      v.s0 = v.s1 = v.s2 = v.mx = 0.0;
      if (tl) {
        const int b = m/ntb, u = m - b*ntb, k0 = (u/nt2)*t3, j0 = (u%nt2)*t2;
#pragma unroll 4
        for (int q = sl; q < ntc; q += ns) {
          const int r = q/n1, k = k0 + r/t2, j = j0 + r%t2;
          fn(b*nkji + (k*n2 + j)*n1 + q - r*n1, v);
        }
      } else {
#pragma unroll 4
        for (int q = sl; q < nkji; q += ns) {fn(m*nkji + q, v);}
      }
      part(4*idx) = v.s0;
      part(4*idx + 1) = v.s1;
      part(4*idx + 2) = v.s2;
      part(4*idx + 3) = v.mx;
    });
    using TP = Kokkos::TeamPolicy<DevExeSpace>;
    using Mem = Kokkos::View<Real*, ScratchMemSpace,
                             Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    const size_t sb = Mem::shmem_size(4*VET_DR_NT);
    Kokkos::parallel_for("vet_lay_red2", TP(DevExeSpace(), nmb, Kokkos::AUTO)
                         .set_scratch_size(0, Kokkos::PerTeam(sb)),
    KOKKOS_LAMBDA(const TP::member_type &tm) {
      Mem sh(tm.team_scratch(0), 4*VET_DR_NT);
      const int m = tm.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, VET_DR_NT), [&](const int t) {
        Real a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        // (kU slots' values loaded at once, then added in the slot order: bitwise)
        constexpr int kU = 8;
        int q = t;
        for (; q + (kU - 1)*VET_DR_NT < ns; q += kU*VET_DR_NT) {
          Real p[kU][4];
          for (int u = 0; u < kU; ++u) {
            for (int c = 0; c < 4; ++c) p[u][c] = part(4*(m*ns + q + u*VET_DR_NT) + c);
          }
          for (int u = 0; u < kU; ++u) {
            a0 += p[u][0];
            a1 += p[u][1];
            a2 += p[u][2];
            a3 = (p[u][3] > a3) ? p[u][3] : a3;
          }
        }
        for (; q < ns; q += VET_DR_NT) {
          const int b = 4*(m*ns + q);
          a0 += part(b);
          a1 += part(b + 1);
          a2 += part(b + 2);
          a3 = (part(b + 3) > a3) ? part(b + 3) : a3;
        }
        sh(4*t) = a0;
        sh(4*t + 1) = a1;
        sh(4*t + 2) = a2;
        sh(4*t + 3) = a3;
      });
      tm.team_barrier();
      for (int w = VET_DR_NT/2; w > 0; w /= 2) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, w), [&](const int t) {
          sh(4*t) += sh(4*(t + w));
          sh(4*t + 1) += sh(4*(t + w) + 1);
          sh(4*t + 2) += sh(4*(t + w) + 2);
          sh(4*t + 3) = (sh(4*(t + w) + 3) > sh(4*t + 3)) ? sh(4*(t + w) + 3)
                                                          : sh(4*t + 3);
        });
        tm.team_barrier();
      }
      Kokkos::single(Kokkos::PerTeam(tm), [&]() {
        for (int c = 0; c < 4; ++c) {part(off2 + 4*m + c) = sh(c);}
      });
    });
    auto hp = Kokkos::subview(hpart, Kokkos::make_pair(0, 4*nmb));
    // (stream-ordered copy + one fence)
    Kokkos::deep_copy(DevExeSpace(), hp, Kokkos::subview(part, Kokkos::make_pair(off2,
                      off2 + 4*nmb)));
    DevExeSpace().fence();
    repro::ReproSum rs[3];
    for (int c = 0; c < 3; ++c) {rs[c].Init();}
    Real mx = 0.0;
    for (int m = 0; m < nmb; ++m) {
      for (int c = 0; c < 3; ++c) {rs[c].Add(hp(4*m + c));}
      mx = (hp(4*m + 3) > mx) ? hp(4*m + 3) : mx;
    }
    repro::AllMergeMax(rs, 3, mx);
    out.s0 = rs[0].Value();
    out.s1 = rs[1].Value();
    out.s2 = rs[2].Value();
    out.mx = mx;
    return;
  }
  repro::ReproVal<3> rv;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, n),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<3> &a) {
    VetBcgVal t;
    t.s0 = t.s1 = t.s2 = t.mx = 0.0;
    fn(idx, t);
    a.s[0].Add(t.s0);
    a.s[1].Add(t.s1);
    a.s[2].Add(t.s2);
    a.mx = (t.mx > a.mx) ? t.mx : a.mx;
  }, repro::ReproRed<3, Kokkos::HostSpace>(rv));
  repro::AllMerge(rv.s, 3);
  VetGlobalMax(rv.mx);
  out.s0 = rv.s[0].Value();
  out.s1 = rv.s[1].Value();
  out.s2 = rv.s[2].Value();
  out.mx = rv.mx;
}

//----------------------------------------------------------------------------------------
//! \fn std::array<Real, 2> GlobalSums
//! \brief the n <= 2 sums f(m, k, j, i, s0, s1) adds over the active cells and ranks
//! (f may also update vectors); layout_invariant: via VetLayReduce

template <typename F>
std::array<Real, 2> GlobalSums(MeshBlockPack *pp, const char *name, const int n,
                               const F &f) {
  auto &indcs = pp->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pp->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  std::array<Real, 2> d = {0.0, 0.0};
  if (pp->pmesh->layout_invariant) {
    VetBcgVal v;
    auto *r = pp->pradvet;
    const int nxb[3] = {nx1, nx2, nx3};   // (layout_invariant_mode = tile)
    VetLayReduce(name, nmb*nkji, nkji, r->lay_part, r->lay_part_h,
    KOKKOS_LAMBDA(const int idx, VetBcgVal &a) {
      int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
      f(m, k+ks, j+js, idx - m*nkji - k*nji - j*nx1 + is, a.s0, a.s1);
    }, v, nxb);
    return {v.s0, v.s1};
  }
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
//! \fn std::array<Real, 3> GlobalSumsMax
//! \brief GlobalSums' n <= 2 sums of f and the max over the active cells and ranks of
//! g(m, k, j, i) >= 0, read after f updated the cell.  layout_invariant: one pass and one
//! reduction (VetLayReduce's max slot rides in the same merge; the sums are those of
//! GlobalSums); else GlobalSums, then the max in a pass of its own

template <typename F, typename G>
std::array<Real, 3> GlobalSumsMax(MeshBlockPack *pp, const char *name, const int n,
                                  const F &f, const G &g) {
  auto &indcs = pp->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pp->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  if (pp->pmesh->layout_invariant) {
    VetBcgVal v;
    auto *r = pp->pradvet;
    const int nxb[3] = {nx1, nx2, nx3};   // (layout_invariant_mode = tile)
    VetLayReduce(name, nmb*nkji, nkji, r->lay_part, r->lay_part_h,
    KOKKOS_LAMBDA(const int idx, VetBcgVal &a) {
      int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
      const int i = idx - m*nkji - k*nji - j*nx1 + is;
      f(m, k+ks, j+js, i, a.s0, a.s1);
      a.mx = fmax(a.mx, g(m, k+ks, j+js, i));
    }, v, nxb);
    return {v.s0, v.s1, v.mx};
  }
  const auto d = GlobalSums(pp, name, n, f);
  Real mx = 0.0;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/nx1;
    a = fmax(a, g(m, k+ks, j+js, idx - m*nkji - k*nji - j*nx1 + is));
  }, Kokkos::Max<Real>(mx));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  return {d[0], d[1], mx};
}

//----------------------------------------------------------------------------------------
//! \fn Real SumExact
//! \brief GlobalSums' s0 by the exact accumulator (utils/repro_sum.hpp): no round-off of
//! the sum, the same for every layout and rank count

template <typename F>
Real SumExact(MeshBlockPack *pp, const char *name, const F &f) {
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx1 = indcs.nx1, nji = indcs.nx2*nx1, nkji = indcs.nx3*nji;
  repro::ReproVal<1> v;
  Kokkos::parallel_reduce(name, Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pp->nmb_thispack*nkji),
  KOKKOS_LAMBDA(const int idx, repro::ReproVal<1> &a) {
    const int m = idx/nkji, k = (idx%nkji)/nji, j = (idx%nji)/nx1;
    Real s0 = 0.0, s1 = 0.0;
    f(m, k+ks, j+js, idx%nx1 + is, s0, s1);
    a.s[0].Add(s0);
  }, repro::ReproRed<1, Kokkos::HostSpace>(v));
  repro::AllMerge(v.s, 1);
  return v.s[0].Value();
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_REDUCE_HPP_
