#ifndef MHD_RSOLVERS_CS_LOWBETA_LLF_HPP_
#define MHD_RSOLVERS_CS_LOWBETA_LLF_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_lowbeta_llf.hpp
//! \brief CUBED SPHERE: replace the HLLD flux by an LLF one in cells whose plasma beta is
//! below a threshold.
//!
//! WHY.  On the cubed sphere a uniform medium at rest carrying a uniform (hence
//! force-free) field is linearly UNSTABLE below beta ~ 0.05: kinetic energy grows out of
//! truncation noise until dt collapses.  Measured on
//! inputs/tests/cubed_sphere_mhd_strat.athinput, iprob = 8, the arm that stops it is
//! DISSIPATION, not accuracy -- hlle (which for a state at rest is identical to LLF)
//! takes nx2 = 32, beta = 0.02 from "dt collapses, dead at t = 2.5" to saturated at
//! KE/V ~ 2.4e-6 with dt still at 81% of its starting value, while ppm4 changes nothing
//! and a well-balanced geometric source changes nothing on its own.  HLLD is the least
//! dissipative of the solvers on exactly the modes involved (its Alfven and slow-mode
//! structure degenerates as the NORMAL field weakens, which on a curvilinear grid with a
//! strong TANGENTIAL field is the common case).
//!
//! WHY NOT JUST RUN HLLE.  The deep atmosphere of the problems this grid exists for is
//! high beta, subsonic and long-lived, and that is where HLLD's smaller dissipation on
//! contacts and Alfven waves is worth having.  Switching per CELL keeps it there and
//! spends the dissipation only where the scheme is not trustworthy anyway.  This is the
//! same trade mhd_fluxes.cpp already makes at a POLAR boundary, where it falls back from
//! HLLD to HLLE for the whole row.
//!
//! The L/R states are the RECONSTRUCTED ones, so the switch costs dissipation, not order.

#include "athena.hpp"
#include "eos/eos.hpp"
#include "llf_mhd_singlestate.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void CSLowBetaLLF
//! \brief overwrite flx/ey/ez with the LLF flux wherever beta < beta_thresh.
//!
//! Must run AFTER the Riemann solver and BEFORE GnomonicEquiangleFlux*, which rotates the
//! momentum flux into the covariant slots: everything here is in the face frame the
//! solver worked in.  Slot mapping and the sign convention on ey/ez follow mhd_fofc.cpp,
//! which calls the same single-state solver.

KOKKOS_INLINE_FUNCTION
void CSLowBetaLLF(TeamMember_t const &member, const EOS_Data &eos, const bool gen,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
     const Real beta_thresh,
     const ScrArray2D<Real> &wl, const ScrArray2D<Real> &wr,
     const ScrArray2D<Real> &bl, const ScrArray2D<Real> &br,
     const ScrArray2D<Real> &dl, const ScrArray2D<Real> &dr,
     const DvceArray4D<Real> &bx,
     DvceArray5D<Real> flx, DvceArray4D<Real> ey, DvceArray4D<Real> ez) {
  const int ivy = IVX + ((ivx-IVX)+1)%3;
  const int ivz = IVX + ((ivx-IVX)+2)%3;
  const int iby = ((ivx-IVX) + 1)%3;
  const int ibz = ((ivx-IVX) + 2)%3;

  par_for_inner(member, il, iu, [&](const int i) {
    const Real bxi = bx(m,k,j,i);
    MHDPrim1D wli, wri;
    wli.d = wl(IDN,i); wli.vx = wl(ivx,i); wli.vy = wl(ivy,i); wli.vz = wl(ivz,i);
    wri.d = wr(IDN,i); wri.vx = wr(ivx,i); wri.vy = wr(ivy,i); wri.vz = wr(ivz,i);
    wli.by = bl(iby,i); wli.bz = bl(ibz,i);
    wri.by = br(iby,i); wri.bz = br(ibz,i);
    if (eos.is_ideal) { wli.e = wl(IEN,i); wri.e = wr(IEN,i); }

    // beta from the FACE state, averaged over the two sides.  For a general EOS the
    // pressure was evaluated once per cell in ConsToPrim and reconstructed into dl/dr.
    const Real pl = gen ? dl(IDPR,i) : eos.IdealGasPressure(wl(IEN,i));
    const Real pr = gen ? dr(IDPR,i) : eos.IdealGasPressure(wr(IEN,i));
    const Real bsql = SQR(bxi) + SQR(wli.by) + SQR(wli.bz);
    const Real bsqr = SQR(bxi) + SQR(wri.by) + SQR(wri.bz);
    const Real bsq = 0.5*(bsql + bsqr);
    // beta = 2p/B^2, written as a product so a vanishing field never divides
    if ((pl + pr) >= beta_thresh*bsq) return;

    MHDCons1D flux;
    if (gen) {
      SingleStateLLF_GenMHD(wli, wri, bxi, pl, pr, dl(IDG1,i), dr(IDG1,i), flux);
    } else {
      SingleStateLLF_MHD(wli, wri, bxi, eos, flux);
    }
    flx(m,IDN,k,j,i) = flux.d;
    flx(m,ivx,k,j,i) = flux.mx;
    flx(m,ivy,k,j,i) = flux.my;
    flx(m,ivz,k,j,i) = flux.mz;
    if (eos.is_ideal) { flx(m,IEN,k,j,i) = flux.e; }
    ey(m,k,j,i) = flux.by;
    ez(m,k,j,i) = flux.bz;
  });
  return;
}
}  // namespace mhd
#endif  // MHD_RSOLVERS_CS_LOWBETA_LLF_HPP_
