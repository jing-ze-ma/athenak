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
// share-1104 P6: the test and the snapshot copy are shared with Hydro
#define SBV_CLASS MHD
#define SBV_NVAR nmhd
#include "hydro/splitbval_impl.hpp"
#undef SBV_CLASS
#undef SBV_NVAR

//----------------------------------------------------------------------------------------
//! \fn void MHD::SplitBvalSnap
//! \brief u0 at the input of the first ConToPrim after a u0 refresh (see above).

void MHD::SplitBvalSnap() {
  sbf_want = false;
  SbfSnapCopy(sbf_s, u0, sbf_valid);
}

}  // namespace mhd
