#ifndef RAD_VET_RAD_VET_LAYREDUCE_HPP_
#define RAD_VET_RAD_VET_LAYREDUCE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_layreduce.hpp
//! \brief layinv-1105 <mesh>/layout_invariant: the order-independent reduction of the
//! implicit BiCGStab (M1LayReduce), shared with the clean solver's sums (CleanSums)

#include "athena.hpp"
#include "utils/repro_sum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//! up to three sums and one max (M1BcgRed, M1LayReduce)
struct M1BcgVal {
  Real s0, s1, s2, mx;
};

constexpr int M1_DR_NT = 256;    // slots per team (a power of 2)
// block mode: level-1 slots per MeshBlock (fills the GPU)
constexpr int M1_LB_NS = 32768;

//----------------------------------------------------------------------------------------
//! \fn M1GlobalMax
inline void M1GlobalMax(Real &a) {
#if MPI_PARALLEL_ENABLED
  Real g;
  MPI_Allreduce(&a, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  a = g;
#else
  (void)a;
#endif
}

//----------------------------------------------------------------------------------------
//! \fn M1LayReduce
//! \brief layinv-1105 <mesh>/layout_invariant: the M1BcgVal reduction of fn(idx, v) over
//! idx in [0, n) AND over all ranks with ORDER-INDEPENDENT sums (repro::ReproSum): each
//! cell's contribution is taken from fn on a zeroed M1BcgVal (0 + term = term exactly;
//! fn runs once per idx, so its side effects are those of the plain reduction) and
//! deposited exactly.  The result depends only on the multiset of cell terms: the same
//! for every MeshBlock layout, rank count, thread schedule and backend.

//!
//! block mode (<mesh>/layout_invariant_mode = block, the default;
//! repro::block_mode): idx = m*nkji + q (M1BcgIdx).  Each MeshBlock's partial is taken in
//! a fixed code-set order (M1_LB_NS strided slots, then one team per MeshBlock with the
//! fixed pairwise tree of M1DetReduce), a function of the MeshBlock alone (any team size,
//! any number of MeshBlocks on the rank).  The nmb partials are copied to the host and
//! deposited exactly (repro::ReproSum), then merged over ranks with the max in one
//! Allgather: the result depends only on the set of MeshBlock partials.
//! bpart / hpart: device / host scratch (sized on first use).

template <class F>
void M1LayReduce(const char *name, const int n, const int nkji, DvceArray1D<Real> &bpart,
                 Kokkos::View<Real*, Kokkos::HostSpace> &hpart, const F &fn,
                 M1BcgVal &out, const int *nxb = nullptr) {
  if (repro::block_mode) {
    // level 1: M1_LB_NS slots per MeshBlock (independent of nmb), slot s accumulates
    // q = s, s + M1_LB_NS, ... (coalesced over s); level 2: one team per MeshBlock, slot
    // t of M1_DR_NT accumulates the level-1 slots t, t + M1_DR_NT, ... in order, then
    // the fixed pairwise tree.  Every order is set by the code.
    // (tile mode, nxb the MeshBlock's cells: the units are the N x N (x2, x3) tiles
    // instead, q the cell's index in its tile: a tile of the MeshBlock's size is the
    // MeshBlock's partial, bitwise)
    const bool tl = repro::tile_n > 0 && nxb != nullptr;
    const int t2 = (tl && nxb[1] > 1) ? repro::tile_n : 1;
    const int t3 = (tl && nxb[2] > 1) ? repro::tile_n : 1;
    const int n1 = tl ? nxb[0] : 1, n2 = tl ? nxb[1] : 1, nt2 = n2/t2;
    const int ntb = tl ? nt2*(nxb[2]/t3) : 1;
    const int nmb = (n/nkji)*ntb;   // (the units)
    const int ntc = n1*t2*t3;   // (tile mode: the tile's cells)
    const int ns = (tl && ntc < M1_LB_NS) ? ntc : M1_LB_NS;
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
      M1BcgVal v;
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
    const size_t sb = Mem::shmem_size(4*M1_DR_NT);
    Kokkos::parallel_for("vet_lay_red2", TP(DevExeSpace(), nmb, Kokkos::AUTO)
                         .set_scratch_size(0, Kokkos::PerTeam(sb)),
    KOKKOS_LAMBDA(const TP::member_type &tm) {
      Mem sh(tm.team_scratch(0), 4*M1_DR_NT);
      const int m = tm.league_rank();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, M1_DR_NT), [&](const int t) {
        Real a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        // (kU slots' values loaded at once, then added in the slot order: bitwise)
        constexpr int kU = 8;
        int q = t;
        for (; q + (kU - 1)*M1_DR_NT < ns; q += kU*M1_DR_NT) {
          Real p[kU][4];
          for (int u = 0; u < kU; ++u) {
            for (int c = 0; c < 4; ++c) p[u][c] = part(4*(m*ns + q + u*M1_DR_NT) + c);
          }
          for (int u = 0; u < kU; ++u) {
            a0 += p[u][0];
            a1 += p[u][1];
            a2 += p[u][2];
            a3 = (p[u][3] > a3) ? p[u][3] : a3;
          }
        }
        for (; q < ns; q += M1_DR_NT) {
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
      for (int w = M1_DR_NT/2; w > 0; w /= 2) {
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
    // (stream-ordered copy + one stream fence, not deep_copy's two device syncs; h200)
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
    M1BcgVal t;
    t.s0 = t.s1 = t.s2 = t.mx = 0.0;
    fn(idx, t);
    a.s[0].Add(t.s0);
    a.s[1].Add(t.s1);
    a.s[2].Add(t.s2);
    a.mx = (t.mx > a.mx) ? t.mx : a.mx;
  }, repro::ReproRed<3, Kokkos::HostSpace>(rv));
  repro::AllMerge(rv.s, 3);
  M1GlobalMax(rv.mx);
  out.s0 = rv.s[0].Value();
  out.s1 = rv.s[1].Value();
  out.s2 = rv.s[2].Value();
  out.mx = rv.mx;
}

} // namespace radvet

#endif // RAD_VET_RAD_VET_LAYREDUCE_HPP_
