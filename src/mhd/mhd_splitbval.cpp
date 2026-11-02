//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_splitbval.cpp
//! \brief ckbval-1102 <mhd>/split_bval_fast: the MHD twin of hydro/hydro_splitbval.cpp.
//! MHD::RTOpSplitBvals exchanges u0(IEN) alone when MHD::SplitBvalFastOK proves the
//! ghosts of every other u0 slot (density, momenta, x_e / floor tracers) already hold
//! what the full exchange would deliver.  The argument and the test are those of
//! Hydro::SplitBvalFastOK (see there), with these MHD points:
//!   - the split source (ck) changes u0(IEN) only; b0 is not exchanged by the refresh in
//!     either path (unchanged), so the face fields need no test;
//!   - the snapshot is u0 at the input of the first ConToPrim after the LAST u0 refresh,
//!     whichever it was: an RK stage, an RKG super-stage (MHD::RecvU also serves the RKG
//!     list, full or rkg_halo_ien_only-masked: either leaves every ghost equal to the
//!     exchange of the active cells), or the split refresh itself.  Under rkg_lean the
//!     intermediate super-stages skip both the exchange and the conversion, so the
//!     snapshot is taken at the final super-stage's conversion, right after its
//!     exchange.  The fused-Strang order (rkg_strang_fuse: leading RKG after
//!     before_timeintegrator) changes only which refresh is the last one;
//!   - anything that writes a slot other than IEN after the snapshot -- the tracer pass
//!     (TracerSource), a ConToPrim tracer rescale, the Alfven-speed density floor (not a
//!     constant), the no-infall momentum removal -- fails the bitwise test, and the
//!     refresh is then the full one;
//!   - hesdirk2 (no split source): RTOpSplitBvals is never called.

#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace mhd {

namespace {
KOKKOS_INLINE_FUNCTION
bool SbfSameBits(const Real a, const Real b) {
  using U = std::conditional<sizeof(Real) == 8, unsigned long long, unsigned int>::type;
  return Kokkos::bit_cast<U>(a) == Kokkos::bit_cast<U>(b);
}
}  // namespace

bool MHD::SplitBvalFastOK() {
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
    const int nvar = nmhd + nscalars;
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
//! \fn void MHD::SplitBvalSnap
//! \brief u0 at the input of the first ConToPrim after a u0 refresh (see above).

void MHD::SplitBvalSnap() {
  sbf_want = false;
  if (sbf_s.extent(0) != u0.extent(0) || sbf_s.extent(1) != u0.extent(1) ||
      sbf_s.extent(2) != u0.extent(2) || sbf_s.extent(3) != u0.extent(3) ||
      sbf_s.extent(4) != u0.extent(4)) {
    Kokkos::realloc(sbf_s, u0.extent(0), u0.extent(1), u0.extent(2), u0.extent(3),
                    u0.extent(4));
  }
  Kokkos::deep_copy(DevExeSpace(), sbf_s, u0);
  sbf_valid = true;
}

}  // namespace mhd
