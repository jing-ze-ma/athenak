#ifndef MHD_RSOLVERS_HLLE_MHD_SINGLESTATE_HPP_
#define MHD_RSOLVERS_HLLE_MHD_SINGLESTATE_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hlle_mhd_singlestate.hpp
//! \brief The HLLE Riemann solver for NR MHD on a SINGLE L/R state, the counterpart of
//! llf_mhd_singlestate.hpp.  Callers that need one face's flux -- FOFC, and the
//! cubed-sphere low-beta fallback -- otherwise have only LLF, which is markedly more
//! diffusive wherever the flow is fast: LLF spreads both sides at |v| + c_f, while HLLE
//! keeps the signed fan and reduces to full upwinding once the flow is supersonic.
//!
//! The wave speeds are the L/R estimate, S_L = min(vl - cl, vr - cr) and
//! S_R = max(vl + cl, vr + cr), NOT the Roe average that hlle_mhd.hpp adds as a third
//! candidate for an ideal gas.  That is deliberate and not an approximation of
//! convenience: hlle_mhd.hpp itself uses exactly this estimate on the general-EOS path,
//! and dropping the Roe candidate can only widen the fan, i.e. err toward MORE
//! dissipation -- the safe direction for a fallback.
//!
//! Sign convention on flux.by / flux.bz matches SingleStateLLF_MHD, so callers can assign
//! ey = flux.by and ez = flux.bz directly, as mhd_fofc.cpp does.

#include "athena.hpp"
#include "eos/eos.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void SingleStateHLLEFluxes
//! \brief the shared algebra, once the pressures and fast speeds are known.

KOKKOS_INLINE_FUNCTION
void SingleStateHLLEFluxes(const MHDPrim1D &wl, const MHDPrim1D &wr, const Real &bxi,
                           const Real pl, const Real pr, const Real cl, const Real cr,
                           const bool is_ideal, MHDCons1D &flux) {
  const Real pbl = 0.5*(bxi*bxi + SQR(wl.by) + SQR(wl.bz));
  const Real pbr = 0.5*(bxi*bxi + SQR(wr.by) + SQR(wr.bz));

  // wave speeds, and the guard hlle_mhd.hpp applies for converging supersonic flow
  const Real al = fmin((wl.vx - cl), (wr.vx - cr));
  const Real ar = fmax((wl.vx + cl), (wr.vx + cr));
  const Real bp = ar > 0.0 ? ar : 1.0e-20;
  const Real bm = al < 0.0 ? al : -1.0e-20;

  const Real vxl = wl.vx - bm;
  const Real vxr = wr.vx - bp;

  MHDCons1D fl, fr;
  fl.d  = wl.d*vxl;
  fr.d  = wr.d*vxr;
  fl.mx = wl.d*wl.vx*vxl + pbl - SQR(bxi);
  fr.mx = wr.d*wr.vx*vxr + pbr - SQR(bxi);
  fl.my = wl.d*wl.vy*vxl - bxi*wl.by;
  fr.my = wr.d*wr.vy*vxr - bxi*wr.by;
  fl.mz = wl.d*wl.vz*vxl - bxi*wl.bz;
  fr.mz = wr.d*wr.vz*vxr - bxi*wr.bz;
  // pl, pr are the GAS pressure for an ideal gas and iso_cs^2*d for an isothermal one,
  // so this term is common to both; only the energy flux is ideal-only.
  fl.mx += pl;
  fr.mx += pr;
  if (is_ideal) {
    const Real el = wl.e + 0.5*wl.d*(SQR(wl.vx)+SQR(wl.vy)+SQR(wl.vz)) + pbl;
    const Real er = wr.e + 0.5*wr.d*(SQR(wr.vx)+SQR(wr.vy)+SQR(wr.vz)) + pbr;
    fl.e = el*vxl + wl.vx*(pl + pbl - bxi*bxi) - bxi*(wl.by*wl.vy + wl.bz*wl.vz);
    fr.e = er*vxr + wr.vx*(pr + pbr - bxi*bxi) - bxi*(wr.by*wr.vy + wr.bz*wr.vz);
  }
  fl.by = wl.by*vxl - bxi*wl.vy;
  fr.by = wr.by*vxr - bxi*wr.vy;
  fl.bz = wl.bz*vxl - bxi*wl.vz;
  fr.bz = wr.bz*vxr - bxi*wr.vz;

  Real tmp = 0.0;
  if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);

  flux.d  = 0.5*(fl.d  + fr.d ) + (fl.d  - fr.d )*tmp;
  flux.mx = 0.5*(fl.mx + fr.mx) + (fl.mx - fr.mx)*tmp;
  flux.my = 0.5*(fl.my + fr.my) + (fl.my - fr.my)*tmp;
  flux.mz = 0.5*(fl.mz + fr.mz) + (fl.mz - fr.mz)*tmp;
  if (is_ideal) flux.e = 0.5*(fl.e + fr.e) + (fl.e - fr.e)*tmp;
  flux.by = -0.5*(fl.by + fr.by) - (fl.by - fr.by)*tmp;
  flux.bz =  0.5*(fl.bz + fr.bz) + (fl.bz - fr.bz)*tmp;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleStateHLLE_MHD
//! \brief ideal gas.

KOKKOS_INLINE_FUNCTION
void SingleStateHLLE_MHD(const MHDPrim1D &wl, const MHDPrim1D &wr, const Real &bxi,
                         const EOS_Data &eos, MHDCons1D &flux) {
  Real pl = 0.0, pr = 0.0, cl, cr;
  if (eos.is_ideal) {
    pl = eos.IdealGasPressure(wl.e);
    pr = eos.IdealGasPressure(wr.e);
    cl = eos.IdealMHDFastSpeed(wl.d, pl, bxi, wl.by, wl.bz);
    cr = eos.IdealMHDFastSpeed(wr.d, pr, bxi, wr.by, wr.bz);
  } else {
    cl = eos.IdealMHDFastSpeed(wl.d, bxi, wl.by, wl.bz);
    cr = eos.IdealMHDFastSpeed(wr.d, bxi, wr.by, wr.bz);
    pl = SQR(eos.iso_cs)*wl.d;
    pr = SQR(eos.iso_cs)*wr.d;
  }
  SingleStateHLLEFluxes(wl, wr, bxi, pl, pr, cl, cr, eos.is_ideal, flux);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SingleStateHLLE_GenMHD
//! \brief GENERAL EOS: the pressures and first adiabatic exponents were evaluated once
//! per cell in ConsToPrim and reconstructed to this interface, as in
//! SingleStateLLF_GenMHD.  The fast speed uses Gamma_1 p in place of gamma p.

KOKKOS_INLINE_FUNCTION
void SingleStateHLLE_GenMHD(const MHDPrim1D &wl, const MHDPrim1D &wr, const Real &bxi,
                            const Real pl, const Real pr, const Real g1l, const Real g1r,
                            MHDCons1D &flux) {
  Real asq = g1l*pl;
  Real ctsq = SQR(wl.by) + SQR(wl.bz);
  Real qsq = bxi*bxi + ctsq + asq;
  Real tmp = bxi*bxi + ctsq - asq;
  const Real cl = sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ctsq))/wl.d);
  asq = g1r*pr;
  ctsq = SQR(wr.by) + SQR(wr.bz);
  qsq = bxi*bxi + ctsq + asq;
  tmp = bxi*bxi + ctsq - asq;
  const Real cr = sqrt(0.5*(qsq + sqrt(tmp*tmp + 4.0*asq*ctsq))/wr.d);
  SingleStateHLLEFluxes(wl, wr, bxi, pl, pr, cl, cr, true, flux);
  return;
}
}  // namespace mhd
#endif  // MHD_RSOLVERS_HLLE_MHD_SINGLESTATE_HPP_
