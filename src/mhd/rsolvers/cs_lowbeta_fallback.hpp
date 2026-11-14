#ifndef MHD_RSOLVERS_CS_LOWBETA_FALLBACK_HPP_
#define MHD_RSOLVERS_CS_LOWBETA_FALLBACK_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file cs_lowbeta_fallback.hpp
//! \brief Cubed sphere: replaces the HLLD flux by the more dissipative HLLE flux on faces
//! whose plasma beta (averaged over the two face states) is below
//! <mhd>/cs_lowbeta_fallback. Suppresses the low-beta instability of near force-free
//! fields on the cubed sphere.

#include "athena.hpp"
#include "eos/eos.hpp"
#include "hlle_mhd_singlestate.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn bool CSLowBetaFires
//! \brief THE face test of the fallback (mean face beta below
//! beta_thresh, written as a product so a vanishing field never divides), shared by
//! CSLowBetaFallback and UctFaceStore (emf_scheme = uct | uct_hlld), so the corner-EMF
//! coefficients switch to HLL on exactly the faces whose flux was replaced.

KOKKOS_INLINE_FUNCTION
bool CSLowBetaFires(const Real pl, const Real pr, const Real bxi, const Real byl,
                    const Real bzl, const Real byr, const Real bzr,
                    const Real beta_thresh) {
  const Real bsql = SQR(bxi) + SQR(byl) + SQR(bzl);
  const Real bsqr = SQR(bxi) + SQR(byr) + SQR(bzr);
  return ((pl + pr) < beta_thresh*(0.5*(bsql + bsqr)));
}

//----------------------------------------------------------------------------------------
//! \fn void CSLowBetaFallback
//! \brief overwrite flx/ey/ez with the HLLE flux wherever beta < beta_thresh.
//!
//! Must run AFTER the Riemann solver and BEFORE GnomonicEquiangleFlux*, which rotates the
//! momentum flux into the covariant slots: everything here is in the face frame the
//! solver worked in. Slot mapping and the sign convention on ey/ez follow mhd_fofc.cpp,
//! which calls the LLF twin of the same single-state solver.

KOKKOS_INLINE_FUNCTION
void CSLowBetaFallback(TeamMember_t const &member, const EOS_Data &eos, const bool gen,
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

    // beta from the FACE state, averaged over the two sides. For a general EOS the
    // pressure was evaluated once per cell in ConsToPrim and reconstructed into dl/dr.
    // An isothermal EOS has no energy slot (IEN is out of range): p = d c_s^2.
    const Real pl = gen ? dl(IDPR,i) : (eos.is_ideal ? eos.IdealGasPressure(wl(IEN,i))
                                                     : wl(IDN,i)*SQR(eos.iso_cs));
    const Real pr = gen ? dr(IDPR,i) : (eos.is_ideal ? eos.IdealGasPressure(wr(IEN,i))
                                                     : wr(IDN,i)*SQR(eos.iso_cs));
    // beta = 2p/B^2, written as a product so a vanishing field never divides
    const bool fire = CSLowBetaFires(pl, pr, bxi, wli.by, wli.bz, wri.by, wri.bz,
                                     beta_thresh);

    if (!fire) return;

    MHDCons1D flux;
    if (gen) {
      SingleStateHLLE_GenMHD(wli, wri, bxi, pl, pr, dl(IDG1,i), dr(IDG1,i), flux);
    } else {
      SingleStateHLLE_MHD(wli, wri, bxi, eos, flux);
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
#endif  // MHD_RSOLVERS_CS_LOWBETA_FALLBACK_HPP_
