//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_splitbval.cpp
//! \brief ckbval-1102 <hydro>/split_bval_fast: the check that lets the ghost refresh after
//! the operator-split source (Hydro::RTOpSplitBvals) exchange u0(IEN) alone.

#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "driver/driver.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace hydro {

//----------------------------------------------------------------------------------------
//! \fn bool Hydro::SplitBvalFastOK
//! \brief ckbval-1102: may the ghost refresh after the operator-split source exchange
//! u0(IEN) alone and still be BITWISE the full refresh?
//!
//! The full refresh overwrites every ghost of every variable with the exchange of the
//! current active cells, then applies the physical BCs and converts the whole mesh.  The
//! energy-only refresh skips the exchange of the other variables X (density, momenta,
//! scalars) and keeps what the ghosts hold.  Everything after the exchange is the same
//! code on the same inputs, so the two agree bit for bit iff every ghost X now equals what
//! the full exchange would put there.
//!
//! S = u0 at the input of the first ConToPrim after the last refresh (sbf_s): ghosts
//! there hold exactly the exchange (+ BCs) of the active cells.  Accept iff, in EVERY cell
//! (active and ghost, all x1 included):
//!   - momenta and scalars are bit-identical to S;
//!   - the density is T(S), with T(d) = (d < dfloor ? dfloor : d) -- the ConToPrim
//!     density floor, applied the same way in every cell;
//!   - within ng+1 cells of a cubed-sphere panel seam (active cells and ghosts on that
//!     side, so the cube-vertex corners too) the density is bit-identical to S and was
//!     not below dfloor (T = identity there).
//! Why that is enough: a plain-copy ghost g of an active cell a had S_g = S_a, so
//! X_g = T(S_g) = T(S_a) = X_a, which the full exchange would copy.  A seam ghost is a
//! resample (basis transform of the tangential momenta, 3-point stencil, monotone clamp,
//! the density guard) of source cells within ng of the seam on the other panel; none of
//! them, nor the ghost, changed since S, so the resample gives what the ghost holds.  The
//! physical-BC ghosts are rebuilt by the BCs in both paths from these same values.  IEN,
//! which the source changed, is exchanged as always (its seam resample reads the
//! density, unchanged since S near seams).  The test is on the bits (-0.0 != +0.0).
//!
//! The check runs on every rank and the answer is reduced (MPI MAX of "bad"), so all
//! ranks pick the same message layout.

namespace {
KOKKOS_INLINE_FUNCTION
bool SbfSameBits(const Real a, const Real b) {
  using U = std::conditional<sizeof(Real) == 8, unsigned long long, unsigned int>::type;
  return Kokkos::bit_cast<U>(a) == Kokkos::bit_cast<U>(b);
}
}  // namespace

bool Hydro::SplitBvalFastOK() {
  Mesh *pm = pmy_pack->pmesh;
  int bad = (sbf_valid && !pm->use_polar_boundary && !pm->multilevel) ? 0 : 1;
  if (bad == 0) {
    auto &indcs = pm->mb_indcs;
    const int ng = indcs.ng;
    const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
    const int nfp = ng + 1;
    const bool cs = pm->use_cubed_sphere;
    const bool td = pm->three_d;
    const Real x2mn = pm->mesh_size.x2min, x2mx = pm->mesh_size.x2max;
    const Real x3mn = pm->mesh_size.x3min, x3mx = pm->mesh_size.x3max;
    auto size = pmy_pack->pmb->mb_size;
    const Real dfl = peos->eos_data.dfloor;
    const int nvar = nhydro + nscalars;
    auto u_ = u0;
    auto s_ = sbf_s;
    const int nmb = pmy_pack->nmb_thispack;
    const int n3 = u0.extent_int(2), n2 = u0.extent_int(3), n1 = u0.extent_int(4);
    int nbad = 0;
    Kokkos::parallel_reduce("sbf_check", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
        DevExeSpace(), {0, 0, 0, 0}, {nmb, n3, n2, n1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &sum) {
      bool fp = false;
      if (cs) {
        fp = (size.d_view(m).x2min <= x2mn && j < js + nfp) ||
             (size.d_view(m).x2max >= x2mx && j > je - nfp);
        if (td) {
          fp = fp || (size.d_view(m).x3min <= x3mn && k < ks + nfp) ||
                     (size.d_view(m).x3max >= x3mx && k > ke - nfp);
        }
      }
      for (int v = 0; v < nvar; ++v) {
        if (v == IEN) continue;
        const Real b = s_(m,v,k,j,i);
        Real e = b;
        if (v == IDN && b < dfl) {
          if (fp) {sum++; continue;}
          e = dfl;
        }
        if (!SbfSameBits(u_(m,v,k,j,i), e)) sum++;
      }
    }, Kokkos::Sum<int>(nbad));
    bad = (nbad > 0) ? 1 : 0;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  return bad == 0;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::SplitBvalSnap
//! \brief u0 at the input of the first ConToPrim after a ghost refresh (called from
//! ConToPrim when sbf_want), i.e. right after the exchange and the physical BCs.  Kept
//! only for the conversions a split refresh can follow: the last RK stage, or the split
//! refresh itself (no driver); any other refresh leaves no valid snapshot.

void Hydro::SplitBvalSnap(Driver *pdrive, int stage) {
  sbf_want = false;
  if (pdrive == nullptr || (!c2p_split && stage == pdrive->nexp_stages)) {
    if (sbf_s.extent(0) != u0.extent(0) || sbf_s.extent(1) != u0.extent(1) ||
        sbf_s.extent(2) != u0.extent(2) || sbf_s.extent(3) != u0.extent(3) ||
        sbf_s.extent(4) != u0.extent(4)) {
      Kokkos::realloc(sbf_s, u0.extent(0), u0.extent(1), u0.extent(2), u0.extent(3),
                      u0.extent(4));
    }
    Kokkos::deep_copy(DevExeSpace(), sbf_s, u0);
    sbf_valid = true;
  }
}

}  // namespace hydro
