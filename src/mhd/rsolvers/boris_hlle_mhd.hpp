#ifndef MHD_RSOLVERS_BORIS_HLLE_MHD_HPP_
#define MHD_RSOLVERS_BORIS_HLLE_MHD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file boris_hlle_mhd.hpp
//! \brief <mhd>/boris_c (boris-1102): the face flux of the Alfven-speed (simplified
//! Boris) limiter.  See the note on EOS_Data::boris_c2inv for the formulation.
//!
//! On every face where (v_A/c_eff)^2 > boris_sw on either side, the flux written by the
//! Riemann solver (and by CSLowBetaFallback) is replaced by an HLLE flux of the BORIS
//! system:
//!  - outer speeds S_L/S_R from the fast speed with B -> sqrt(f) B, an upper bound of the
//!    fast speed of MHD with perpendicular inertia rho/f (same determinant, larger trace
//!    of the 2x2 fast/slow block);
//!  - physical fluxes F_L, F_R are the standard MHD ones (the simplified Boris system
//!    keeps the MHD fluxes and changes only the conserved momentum);
//!  - the HLL dissipation S_L S_R (U_R - U_L) uses the Boris momentum
//!    M_B = rho (v_par + v_perp/f) instead of rho v, so that after MHD::BorisCorrect
//!    multiplies the perpendicular momentum change by f the effective dissipation of
//!    v_perp is (S/2) dv_perp, i.e. the scheme is the plain HLLE of the limited system
//!    (with rho v it would be f S/2: an undamped perpendicular-velocity checkerboard).
//! The total energy dissipation uses the standard total energy (the code's E excludes the
//! electric-field energy rho v_perp^2 (1/f - 1)/2, see MHD::BorisCorrect).
//! One flux per face from the two reconstructed states, symmetric in L <-> R: the update
//! telescopes (mass, energy and the RIEMANN momentum flux are conservative; the momentum
//! non-conservation of the limiter is only the cell-wise BorisCorrect, measured there).
//! Runs in the face frame of the solver, BEFORE the gnomonic rotation (as
//! CSLowBetaFallback).  Ideal or general EOS (isothermal refused at start-up).

#include "athena.hpp"
#include "eos/eos.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn bool BorisFaceTest
//! \brief uctaudit-1103: THE test of a Boris face (fmax of the two sides' (v_A/c_eff)^2
//! > boris_sw) and its per-side sqrt(f), shared by BorisHLLEFace (the face flux) and
//! UctFaceStore (the emf_scheme = uct | uct_hlld corner-EMF coefficients), so both take
//! the same decision and the same limited fast speeds bit for bit.
//! b2 = |B|^2 = bx*bx + by*by + bz*bz, v2 = |v|^2 of each side.

KOKKOS_INLINE_FUNCTION
bool BorisFaceTest(const EOS_Data &eos, const Real dl, const Real b2l, const Real pl,
                   const Real v2l, const Real dr, const Real b2r, const Real pr,
                   const Real v2r, Real &fl, Real &fr) {
  const Real xl = eos.BorisX(dl, b2l, pl, v2l);
  const Real xr = eos.BorisX(dr, b2r, pr, v2r);
  if (!(fmax(xl, xr) > eos.boris_sw)) return false;
  fl = eos.BorisF(xl); fr = eos.BorisF(xr);
  return true;
}

KOKKOS_INLINE_FUNCTION
void BorisHLLEFace(TeamMember_t const &member, const EOS_Data &eos, const bool gen,
     const int m, const int k, const int j, const int il, const int iu, const int ivx,
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
    const Real dlv = wl(IDN,i), drv = wr(IDN,i);
    const Real byl = bl(iby,i), bzl = bl(ibz,i), byr = br(iby,i), bzr = br(ibz,i);
    Real pl, pr, g1l, g1r;
    if (gen) {
      pl = dl(IDPR,i); pr = dr(IDPR,i); g1l = dl(IDG1,i); g1r = dr(IDG1,i);
    } else {
      pl = eos.IdealGasPressure(wl(IEN,i)); pr = eos.IdealGasPressure(wr(IEN,i));
      g1l = eos.gamma; g1r = eos.gamma;
    }
    const Real b2l = bxi*bxi + byl*byl + bzl*bzl;
    const Real b2r = bxi*bxi + byr*byr + bzr*bzr;
    const Real v2l = SQR(wl(IVX,i)) + SQR(wl(IVY,i)) + SQR(wl(IVZ,i));
    const Real v2r = SQR(wr(IVX,i)) + SQR(wr(IVY,i)) + SQR(wr(IVZ,i));
    Real fl_, fr_;
    if (!BorisFaceTest(eos, dlv, b2l, pl, v2l, drv, b2r, pr, v2r, fl_, fr_)) return;
    const Real sl = sqrt(fl_), sr = sqrt(fr_);
    const Real cl = eos.FastSpeedFromP(dlv, pl, g1l, sl*bxi, sl*byl, sl*bzl);
    const Real cr = eos.FastSpeedFromP(drv, pr, g1r, sr*bxi, sr*byr, sr*bzr);

    const Real vxl = wl(ivx,i), vyl = wl(ivy,i), vzl = wl(ivz,i);
    const Real vxr = wr(ivx,i), vyr = wr(ivy,i), vzr = wr(ivz,i);
    const Real pbl = 0.5*b2l, pbr = 0.5*b2r;

    const Real al = fmin((vxl - cl), (vxr - cr));
    const Real ar = fmax((vxl + cl), (vxr + cr));
    const Real bp = ar > 0.0 ? ar : 1.0e-20;
    const Real bm = al < 0.0 ? al : -1.0e-20;

    // <mhd>/boris_face_fmax: one face factor for the Boris momentum (the wave speeds
    // above keep each side's own f)
    if (eos.boris_face_fmax) { fl_ = fmax(fl_, fr_); fr_ = fl_; }
    // extra Boris momentum rho (1/f - 1) v_perp of each side
    Real exl[3] = {0.0, 0.0, 0.0}, exr[3] = {0.0, 0.0, 0.0};
    if (b2l > 0.0 && !eos.boris_diss_plain) {
      const Real vb = (vxl*bxi + vyl*byl + vzl*bzl)/b2l;
      const Real w = dlv*(1.0/fl_ - 1.0);
      exl[0] = w*(vxl - vb*bxi); exl[1] = w*(vyl - vb*byl); exl[2] = w*(vzl - vb*bzl);
    }
    if (b2r > 0.0 && !eos.boris_diss_plain) {
      const Real vb = (vxr*bxi + vyr*byr + vzr*bzr)/b2r;
      const Real w = drv*(1.0/fr_ - 1.0);
      exr[0] = w*(vxr - vb*bxi); exr[1] = w*(vyr - vb*byr); exr[2] = w*(vzr - vb*bzr);
    }

    const Real uxl = vxl - bm;
    const Real uxr = vxr - bp;
    MHDCons1D fl, fr, flux;
    fl.d  = dlv*uxl;
    fr.d  = drv*uxr;
    fl.mx = dlv*vxl*uxl + pbl - bxi*bxi + pl - bm*exl[0];
    fr.mx = drv*vxr*uxr + pbr - bxi*bxi + pr - bp*exr[0];
    fl.my = dlv*vyl*uxl - bxi*byl - bm*exl[1];
    fr.my = drv*vyr*uxr - bxi*byr - bp*exr[1];
    fl.mz = dlv*vzl*uxl - bxi*bzl - bm*exl[2];
    fr.mz = drv*vzr*uxr - bxi*bzr - bp*exr[2];
    const Real eil = wl(IEN,i);   // internal energy density (ideal and general EOS)
    const Real eir = wr(IEN,i);
    const Real el = eil + 0.5*dlv*(vxl*vxl + vyl*vyl + vzl*vzl) + pbl;
    const Real er = eir + 0.5*drv*(vxr*vxr + vyr*vyr + vzr*vzr) + pbr;
    fl.e = el*uxl + vxl*(pl + pbl - bxi*bxi) - bxi*(byl*vyl + bzl*vzl);
    fr.e = er*uxr + vxr*(pr + pbr - bxi*bxi) - bxi*(byr*vyr + bzr*vzr);
    fl.by = byl*uxl - bxi*vyl;
    fr.by = byr*uxr - bxi*vyr;
    fl.bz = bzl*uxl - bxi*vzl;
    fr.bz = bzr*uxr - bxi*vzr;

    Real tmp = 0.0;
    if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);
    flux.d  = 0.5*(fl.d  + fr.d ) + (fl.d  - fr.d )*tmp;
    flux.mx = 0.5*(fl.mx + fr.mx) + (fl.mx - fr.mx)*tmp;
    flux.my = 0.5*(fl.my + fr.my) + (fl.my - fr.my)*tmp;
    flux.mz = 0.5*(fl.mz + fr.mz) + (fl.mz - fr.mz)*tmp;
    flux.e  = 0.5*(fl.e  + fr.e ) + (fl.e  - fr.e )*tmp;
    flux.by = -0.5*(fl.by + fr.by) - (fl.by - fr.by)*tmp;
    flux.bz =  0.5*(fl.bz + fr.bz) + (fl.bz - fr.bz)*tmp;

    flx(m,IDN,k,j,i) = flux.d;
    flx(m,ivx,k,j,i) = flux.mx;
    flx(m,ivy,k,j,i) = flux.my;
    flx(m,ivz,k,j,i) = flux.mz;
    flx(m,IEN,k,j,i) = flux.e;
    ey(m,k,j,i) = flux.by;
    ez(m,k,j,i) = flux.bz;
  });
  return;
}
}  // namespace mhd
#endif  // MHD_RSOLVERS_BORIS_HLLE_MHD_HPP_
