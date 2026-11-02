//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_c2p_idem.cpp
//! \brief ckbval-1102 <hydro>/c2p_idempotent (read only when named, default false).
//!
//! WHY.  A conversion of a state that has already been converted is not a no-op today:
//! under etotgrav every call does u0(IEN) = (E - rho Phi) + rho Phi in place, the T(d,e)
//! solve is warm-started from the previous T and stops at a tolerance, and a floor that
//! fired once can fire again (or not) on the rounded state.  The ConToPrimSplit tasks
//! repeat the conversion just made (after RTOpSplitBvals; in "before_timeintegrator"),
//! so they nudge the state by round-off -- and sometimes by a table clamp -- every cycle.
//!
//! WHAT.  Exact idempotence by memo, not by changing the inversion: after every
//! conversion the outputs (u0, w0, wder, wtemp, every cell) are kept; a ConToPrimSplit
//! call whose inputs AND outputs are all bitwise those kept (nothing touched the state
//! since) returns at once.  Anything else -- any cell, any slot changed -- runs the
//! conversion exactly as today.  So with the key on a repeated conversion of an unchanged
//! state is a bitwise no-op, and <hydro>/split_c2p_once on vs off is bitwise.  Making the
//! inversion itself a fixed point (write-back only where a floor fired, T iterated to a
//! bitwise fixed point) was not taken: in a floored cell (e_floor + rho Phi) - rho Phi
//! need not round back to e_floor, so the second call can fire the floor or not, and a
//! Newton iterate is not guaranteed to have a bitwise fixed point.
//!
//! Only ConToPrimSplit calls are candidates (the RK-stage and split-refresh conversions
//! always follow a change); the frozen first conversion of a restart is never skipped.

#include <type_traits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"

namespace hydro {

namespace {
template <typename V>
void C2PIdemCopy(V &dst, const V &src) {
  if (src.size() == 0) return;
  if (dst.size() != src.size() || dst.extent(0) != src.extent(0)) {
    dst = V(Kokkos::view_alloc(Kokkos::WithoutInitializing, "c2pi_memo"), src.layout());
  }
  Kokkos::deep_copy(DevExeSpace(), dst, src);
}

KOKKOS_INLINE_FUNCTION
bool C2PIdemSame(const Real a, const Real b) {
  using U = std::conditional<sizeof(Real) == 8, unsigned long long, unsigned int>::type;
  return Kokkos::bit_cast<U>(a) == Kokkos::bit_cast<U>(b);
}

template <typename V>
int C2PIdemDiff(const V &a, const V &b) {
  if (a.size() == 0) return 0;
  if (a.size() != b.size()) return 1;
  const Real *pa = a.data();
  const Real *pb = b.data();
  int nd = 0;
  Kokkos::parallel_reduce("c2pi_check", Kokkos::RangePolicy<DevExeSpace,
                          Kokkos::IndexType<int64_t>>(DevExeSpace(), 0,
                          static_cast<int64_t>(a.size())),
  KOKKOS_LAMBDA(const int64_t n, int &s) {
    if (!C2PIdemSame(pa[n], pb[n])) s = 1;
  }, Kokkos::Max<int>(nd));
  return nd;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn bool Hydro::C2PIdemSkip
//! \brief true when this conversion would only repeat the last one (see above)

bool Hydro::C2PIdemSkip() {
  bool skip = c2pi_valid && !c2p_freeze_derived && !c2p_eint_rst;
  if (skip) {
    skip = (C2PIdemDiff(u0, c2pi_u) == 0) && (C2PIdemDiff(w0, c2pi_w) == 0) &&
           (C2PIdemDiff(wder, c2pi_d) == 0) && (C2PIdemDiff(wtemp, c2pi_t) == 0);
  }
  if (skip) {++c2pi_nskip;} else {++c2pi_nrun;}
  return skip;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::C2PIdemSave
//! \brief keep the outputs of the conversion that just ran

void Hydro::C2PIdemSave() {
  C2PIdemCopy(c2pi_u, u0);
  C2PIdemCopy(c2pi_w, w0);
  C2PIdemCopy(c2pi_d, wder);
  C2PIdemCopy(c2pi_t, wtemp);
  c2pi_valid = true;
}

}  // namespace hydro
