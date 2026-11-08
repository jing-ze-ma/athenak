#ifndef UTILS_REPRO_SUM_HPP_
#define UTILS_REPRO_SUM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file repro_sum.hpp
//! \brief layinv-1105 (USER 11-05 "fix invariant"): an ORDER-INDEPENDENT sum of doubles.
//!
//! A plain floating-point sum depends on the order of its terms, so every global or
//! shell sum taken in (MeshBlock, k, j, i) order changes in the last bits with the
//! MeshBlock layout and the rank count.  ReproSum is an exact binned accumulator
//! (ReproBLAS idea, truncation variant): the terms are split on a GLOBAL exponent grid
//! (multiples of RS_W bits) into RS_K integer folds; fold f of the accumulator with top
//! grid exponent g counts units of 2^(g - f RS_W).  A term x (|x| < 2^(g + RS_W)) is
//! deposited by truncation: y = x 2^-g, q0 = trunc(y), y = (y - q0) 2^RS_W, q1 = ...;
//! every step is exact (power-of-two scaling, the fraction of a double is exact), the
//! split of x depends only on (x, g), and a larger g only adds zero folds on top and
//! drops the lowest (whose parts the larger g truncates away anyway).  The int64 fold
//! sums are therefore the same for ANY order of the terms and any way of merging
//! partial accumulators: the converted value depends only on the multiset of terms.
//! Bits below 2^(g - RS_K RS_W) of each term are dropped (~2^-95 of the largest term);
//! the int64 folds take 2^31 deposits of 2^32 units.  Non-finite terms go to a plain
//! side sum (inf/nan propagate).  Usable on host and device (Kokkos reducer below).

#include "athena.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace repro {

//! prod-1108 <mesh>/layout_invariant_mode: true = "block" (DEFAULT): every sum is taken
//! per MeshBlock in a FIXED order set by the code (independent of how many MeshBlocks a
//! rank holds and of the launch configuration), and only the per-MeshBlock partials go
//! into the exact accumulator below (a restart changes the rank assignment of the
//! MeshBlocks, never the MeshBlocks themselves).  false = "exact": every cell term
//! deposited exactly (layinv-1105).  Set once by the Mesh constructor.
inline bool block_mode = true;

constexpr int RS_K = 3;      // folds
constexpr int RS_W = 32;     // bits per fold
// the lowest top grid exponent: every fold scale 2^(g - f RS_W) and every remainder
// stays a NORMAL double (2^(-960 - 64) > 2^-1022 + 53 bits margin for the remainders)
constexpr int RS_GMIN = -960;

struct ReproSum {
  int g;                      // top grid exponent (a multiple of RS_W)
  std::int64_t f[RS_K];       // fold f counts units of 2^(g - f RS_W)
  Real sp;                    // the non-finite terms (plain)

  KOKKOS_INLINE_FUNCTION void Init() {
    g = RS_GMIN;
    for (int q = 0; q < RS_K; ++q) {f[q] = 0;}
    sp = 0.0;
  }
  // shift the folds to a higher top grid exponent gn (> g)
  KOKKOS_INLINE_FUNCTION void Raise(const int gn) {
    const int s = (gn - g)/RS_W;
    for (int q = RS_K - 1; q >= 0; --q) {f[q] = (q >= s) ? f[q - s] : 0;}
    g = gn;
  }
  KOKKOS_INLINE_FUNCTION void Add(const Real x) {
    if (x == 0.0) {return;}
    if (!(fabs(x) <= 1.0e308)) {sp += x; return;}   // inf / nan
    int ex;
    (void) frexp(x, &ex);      // |x| < 2^ex
    // the smallest grid g >= ex - RS_W (ceil division towards +inf)
    int gn = ex - RS_W;
    gn = (gn >= 0) ? ((gn + RS_W - 1)/RS_W)*RS_W : -((-gn)/RS_W)*RS_W;
    if (gn < RS_GMIN) {gn = RS_GMIN;}
    if (gn > g) {Raise(gn);}
    Real y = ldexp(x, -g);    // exact, |y| < 2^RS_W
    for (int q = 0; q < RS_K; ++q) {
      const std::int64_t n = static_cast<std::int64_t>(y);   // truncation
      f[q] += n;
      y = (y - static_cast<Real>(n))*4294967296.0;   // exact: fraction x 2^32
    }
  }
  KOKKOS_INLINE_FUNCTION void Join(const ReproSum &o) {
    ReproSum b = o;
    if (b.g > g) {Raise(b.g);}
    if (g > b.g) {b.Raise(g);}
    for (int q = 0; q < RS_K; ++q) {f[q] += b.f[q];}
    sp += b.sp;
  }
  // the value: carries normalised from the lowest fold up (exact integer arithmetic,
  // a function of the exact sum only), then a fixed sequence of double operations
  KOKKOS_INLINE_FUNCTION Real Value() const {
    std::int64_t c[RS_K];
    for (int q = 0; q < RS_K; ++q) {c[q] = f[q];}
    for (int q = RS_K - 1; q > 0; --q) {
      const std::int64_t hi = c[q] >> RS_W;   // floor division by 2^32
      c[q] -= hi*(static_cast<std::int64_t>(1) << RS_W);
      c[q-1] += hi;
    }
    Real v = 0.0;
    for (int q = RS_K - 1; q >= 0; --q) {
      v += ldexp(static_cast<Real>(c[q]), g - q*RS_W);
    }
    return v + sp;
  }
};

//! N sums and one max (the shape of the M1 Krylov reductions)
template <int N>
struct ReproVal {
  ReproSum s[N];
  Real mx;
};

//! Kokkos reducer of ReproVal<N> (same shape as Kokkos::Sum); 0 is the identity of mx
//! (the reduced maxima are of non-negative quantities)
template <int N, class Space>
struct ReproRed {
 public:
  using reducer = ReproRed<N, Space>;
  using value_type = ReproVal<N>;
  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool refs;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit ReproRed(value_type &v) : value(&v), refs(true) {}
  KOKKOS_INLINE_FUNCTION
  void join(value_type &d, const value_type &s) const {
    for (int q = 0; q < N; ++q) {d.s[q].Join(s.s[q]);}
    d.mx = (s.mx > d.mx) ? s.mx : d.mx;
  }
  KOKKOS_INLINE_FUNCTION
  void init(value_type &v) const {
    for (int q = 0; q < N; ++q) {v.s[q].Init();}
    v.mx = 0.0;
  }
  KOKKOS_INLINE_FUNCTION
  value_type &reference() const {return *value.data();}
  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {return value;}
  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const {return refs;}
};

//! merge n accumulators over all MPI ranks in place (Allgather of the raw words, then
//! Join: exact, so the rank order does not matter); nothing to do on one rank
inline void AllMerge(ReproSum *v, const int n) {
#if MPI_PARALLEL_ENABLED
  int nr = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nr);
  if (nr <= 1 || n <= 0) {return;}
  static_assert(sizeof(ReproSum) % sizeof(std::int64_t) == 0, "ReproSum padding");
  const int nw = static_cast<int>(sizeof(ReproSum)/sizeof(std::int64_t));
  std::vector<std::int64_t> mine(static_cast<std::size_t>(n)*nw);
  std::memcpy(mine.data(), v, static_cast<std::size_t>(n)*sizeof(ReproSum));
  std::vector<std::int64_t> all(mine.size()*nr);
  MPI_Allgather(mine.data(), n*nw, MPI_INT64_T, all.data(), n*nw, MPI_INT64_T,
                MPI_COMM_WORLD);
  for (int q = 0; q < n; ++q) {v[q].Init();}
  for (int r = 0; r < nr; ++r) {
    for (int q = 0; q < n; ++q) {
      ReproSum o;
      std::memcpy(&o, all.data() + (static_cast<std::size_t>(r)*n + q)*nw,
                  sizeof(ReproSum));
      v[q].Join(o);
    }
  }
#else
  (void) v; (void) n;
#endif
}

//! AllMerge of n accumulators plus the global max of mx >= 0 in ONE MPI call (prod-1108:
//! the max rides in the gathered words as its bit pattern)
inline void AllMergeMax(ReproSum *v, const int n, Real &mx) {
#if MPI_PARALLEL_ENABLED
  int nr = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &nr);
  if (nr <= 1) {return;}
  const int nw = static_cast<int>(sizeof(ReproSum)/sizeof(std::int64_t));
  const int nt = n*nw + 1;
  std::vector<std::int64_t> mine(nt);
  std::memcpy(mine.data(), v, static_cast<std::size_t>(n)*sizeof(ReproSum));
  std::memcpy(mine.data() + n*nw, &mx, sizeof(Real));
  std::vector<std::int64_t> all(static_cast<std::size_t>(nt)*nr);
  MPI_Allgather(mine.data(), nt, MPI_INT64_T, all.data(), nt, MPI_INT64_T,
                MPI_COMM_WORLD);
  for (int q = 0; q < n; ++q) {v[q].Init();}
  for (int r = 0; r < nr; ++r) {
    const std::int64_t *b = all.data() + static_cast<std::size_t>(r)*nt;
    for (int q = 0; q < n; ++q) {
      ReproSum o;
      std::memcpy(&o, b + q*nw, sizeof(ReproSum));
      v[q].Join(o);
    }
    Real m;
    std::memcpy(&m, b + n*nw, sizeof(Real));
    mx = (m > mx) ? m : mx;
  }
#else
  (void) v; (void) n; (void) mx;
#endif
}

//! convenience: n plain values -> n order-independent sums over all ranks, where each
//! rank deposits its terms with Add; Values() writes the results
struct ReproVec {
  std::vector<ReproSum> a;
  explicit ReproVec(const int n) : a(n) {for (auto &x : a) {x.Init();}}
  void Add(const int q, const Real x) {a[q].Add(x);}
  void Merge() {AllMerge(a.data(), static_cast<int>(a.size()));}
  Real Value(const int q) const {return a[q].Value();}
};

//! per-shell order-independent sums over all ranks: out[q*c1 + i] = sum over the
//! cells (m, k, j) of this rank and every other of f(q, m, k, j, i), i in [ilo, ihi]
//! (0 elsewhere).  One device thread per (q, i), serial over (m, k, j) (any order would
//! do: the accumulator is exact), then the exact rank merge on the host.
//! f is a device functor Real f(int q, int m, int k, int j, int i).
//! block_mode (prod-1108): one thread per (q, m, i) sums its MeshBlock's (k, j) cells in
//! increasing (k, j) order (a plain double sum: a function of the MeshBlock alone); one
//! thread per (q, i) then deposits the nmb per-MeshBlock partials exactly, and the exact
//! rank merge follows.  The result depends only on the set of MeshBlock partials.
template <class F>
void ShellSums(const char *name, const int nq, const int c1, const int nmb,
               const int ks, const int ke, const int js, const int je, const int ilo,
               const int ihi, const F &f, std::vector<Real> &out) {
  Kokkos::View<ReproSum**, DevMemSpace> rs(std::string(name), nq, c1);
  if (block_mode) {
    Kokkos::View<Real***, DevMemSpace> bp(Kokkos::view_alloc(std::string(name) + "_blk",
                                          Kokkos::WithoutInitializing), nq, nmb, c1);
    const int nqm = nq*nmb;
    Kokkos::parallel_for(std::string(name) + "_blk",
                         Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nqm*c1),
    KOKKOS_LAMBDA(const int idx) {
      const int i = idx % c1, qm = idx/c1;
      const int q = qm/nmb, m = qm - q*nmb;
      Real a = 0.0;
      if (i >= ilo && i <= ihi) {
        for (int k = ks; k <= ke; ++k) {
          for (int j = js; j <= je; ++j) {a += f(q, m, k, j, i);}
        }
      }
      bp(q, m, i) = a;
    });
    Kokkos::parallel_for(std::string(name) + "_dep",
                         Kokkos::RangePolicy<DevExeSpace>(DevExeSpace(), 0, nq*c1),
    KOKKOS_LAMBDA(const int idx) {
      const int q = idx/c1, i = idx % c1;
      ReproSum acc;
      acc.Init();
      for (int m = 0; m < nmb; ++m) {acc.Add(bp(q, m, i));}
      rs(q, i) = acc;
    });
  } else {
  // one team per (q, i), its threads over the (m, k, j) cells: the exact accumulator makes
  // the team's combination order irrelevant (a serial loop per shell was latency-bound on
  // the GPU: +7 % per cycle on the BSG deck with the vet_gd twin / half-range means)
  const int nk = ke - ks + 1, nj = je - js + 1;
  const int ncell = nmb*nk*nj;
  using TP = Kokkos::TeamPolicy<DevExeSpace>;
  Kokkos::parallel_for(std::string(name), TP(DevExeSpace(), nq*c1, Kokkos::AUTO),
  KOKKOS_LAMBDA(const TP::member_type &tm) {
    const int q = tm.league_rank()/c1, i = tm.league_rank() % c1;
    ReproVal<1> acc;
    acc.s[0].Init();
    acc.mx = 0.0;
    if (i >= ilo && i <= ihi) {
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ncell),
      [&](const int n, ReproVal<1> &a) {
        const int m = n/(nk*nj), r = n - m*nk*nj;
        const int k = ks + r/nj, j = js + r % nj;
        a.s[0].Add(f(q, m, k, j, i));
      }, ReproRed<1, Kokkos::HostSpace>(acc));
    }
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {rs(q, i) = acc.s[0];});
  });
  }
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rs);
  std::vector<ReproSum> v(static_cast<std::size_t>(nq)*c1);
  for (int q = 0; q < nq; ++q) {
    for (int i = 0; i < c1; ++i) {v[static_cast<std::size_t>(q)*c1 + i] = h(q, i);}
  }
  AllMerge(v.data(), static_cast<int>(v.size()));
  out.resize(v.size());
  for (std::size_t n = 0; n < v.size(); ++n) {out[n] = v[n].Value();}
}

}  // namespace repro

#endif  // UTILS_REPRO_SUM_HPP_
