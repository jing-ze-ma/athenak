//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_tracers.cpp
//! \brief the MHD passive tracers (tracer-1101): the advected electron fraction
//! <mhd>/xe_tracer (xetracer-1101; model in xe_tracer.hpp) and the floor-mass tracer
//! <mhd>/floor_tracer (floortrc-1101).  Slot layout in mhd_tracers.hpp.
//!
//! One task, TracerSource, runs once per cycle in "before_timeintegrator" (after the RT
//! half step, before SaveMHDState), over EVERY cell including ghosts (pointwise, so a
//! ghost that mirrors an interior cell gets the identical value):
//!   1. FloorTracerDrain (floor tracer on): remove floor-added mass where the floor
//!      does not bind;
//!   2. XeTracerUpdate (x_e on): recombination relaxation toward Saha over dt;
//!   3. TracerX1Ghosts: the radial ghosts at physical x1 boundaries take the edge cell's
//!      specific tracer value (also done on u0 after every physical BC call);
//!   4. eta (Ohmic/ambipolar) from the relaxed z (x_e on).
//! The floor tracer's SOURCE is not here: Ideal/GeneralMHD::ConsToPrim add every density
//! floor increment (dfloor and the Alfven floor) to it, at every conversion.

#include <iostream>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/resistivity.hpp"
#include "pgen/pgen.hpp"
#include "mhd/mhd.hpp"
#include "mhd/mhd_tracers.hpp"
#include "mhd/xe_tracer.hpp"
#include "utils/c2p_track.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn void TracerPrimFromCons (file scope)
//! \brief trcrst-1103: w0 = u0/u0(IDN) in the tracer slots, every cell incl. ghosts (the
//! expression ConsToPrim and ConsToPrimFrozen use for the scalars)

static void TracerPrimFromCons(MHD *pm) {
  if (pm->tracers.NTracer() == 0) return;
  // u0 and w0 are (nmb_thispack, nvar, n3, n2, n1), ghosts included
  const int nmb = pm->u0.extent_int(0);
  const int n3m1 = pm->u0.extent_int(2) - 1;
  const int n2m1 = pm->u0.extent_int(3) - 1;
  const int n1m1 = pm->u0.extent_int(4) - 1;
  const int t0 = (pm->ixe >= 0) ? pm->ixe : pm->ifltr;
  const int t1 = (pm->ifltr >= 0) ? pm->ifltr : pm->ixe;
  auto u0_ = pm->u0;
  auto w0_ = pm->w0;
  par_for("mhd_tracer_w0", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real ud = u0_(m,IDN,k,j,i);
    for (int n=t0; n<=t1; ++n) {
      w0_(m,n,k,j,i) = u0_(m,n,k,j,i)/ud;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MHD::TracerSource
//! \brief the per-cycle tracer pass: floor drain, x_e relaxation over the cycle dt,
//! x1 ghosts, then eta from the new z

TaskStatus MHD::TracerSource(Driver *pdrive, int stage) {
  if (trc_strang && stage == 0) return TaskStatus::complete;   // Driver: TracerSourceHalf
  // rkgfuse-1102 tracer_strang: stage = -1|-2 = first|second half (TracerSourceHalf)
  if (floor_tracer && stage >= -1) FloorTracerDrain();
  if (xe_tracer) XeTracerUpdate(pmy_pack->pmesh->dt*((stage < 0) ? 0.5 : 1.0), false);
  TracerX1Ghosts(true);
  // trcrst-1103: the second tracer half closes the cycle, so the state it leaves is the
  // one a restart file stores -- and the file holds u0 only: a restart rebuilds the
  // tracer primitives as u0(n)/u0(IDN) (ConsToPrim[Frozen]).  XeTracerUpdate and
  // TracerX1Ghosts set w0(n) = y and u0(n) = rho y, and (rho y)/rho != y in the last
  // bit, which went into eta (here, and the RefreshEta of the next leading RKG step) and
  // made a restart differ from the straight run.  Re-derive w0 of the tracer slots from
  // u0 exactly as the restart does, every cell incl. ghosts (half 2 only: the first
  // half is followed by the RK stages, whose ConToPrim overwrites w0 anyway).
  if (stage == -2) TracerPrimFromCons(this);
  if (xe_tracer && presist != nullptr) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int ng = indcs.ng;
    const int n1m1 = indcs.nx1 + 2*ng - 1;
    const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    const Real rgas = pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas;
    if (presist->OhmicNeedsUpdate()) {
      presist->SetResistivity(w0, peos->eos_data, rgas, presist->eta_b,
                              0, n1m1, 0, n2m1, 0, n3m1);
    }
    if (presist->AmbipolarNeedsUpdate()) {
      presist->SetAmbipolar(w0, peos->eos_data, rgas, 0, n1m1, 0, n2m1, 0, n3m1);
    }
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::TracerSourceHalf
//! \brief rkgfuse-1102 <mhd>/tracer_strang: one half (1|2) of the Strang-split tracer pass

void MHD::TracerSourceHalf(int half) {
  (void) TracerSource(nullptr, -half);
}


//----------------------------------------------------------------------------------------
//! \fn void MHD::XeTracerUpdate
//! \brief the tracer is y = n_e/n_H (electrons per H nucleus, n_H = xh rho/(a_H m_u)).
//! init = true: y = local Saha x_eq n_tot/n_H.  Otherwise one backward-Euler
//! recombination step of n_e = y n_H over dt (code units) toward n_eq = x_eq n_tot,
//! capped at max(n_eq, gas-phase donors).  Writes w0 and u0 (= rho y) of the tracer in
//! every cell; the x1 boundary ghosts are the caller's (TracerX1Ghosts).  A density
//! floor that raises rho in ConsToPrim scales u0 = rho y with rho (EOS_Data::ixe), so the
//! floor mass inherits the cell's y; the floor-tracer drain also keeps y (it scales s).

void MHD::XeTracerUpdate(Real dt, bool init) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb = pmy_pack->nmb_thispack;
  const EOS_Data eos = peos->eos_data;
  const bool gen = eos.IsGeneral();
  const Real gm1 = eos.gamma - 1.0;
  const Real rgas = pmy_pack->pmesh->pgen->hot_jupiter_param.Rgas;
  if (!gen && !(rgas > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<mhd>/xe_tracer on an ideal gas needs <problem>/Rgas"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const xetracer::XeTracerParams par = xe_par;
  const Real dt_s = dt*par.time_cgs;
  const Real nh_fac = par.xh/(eos_cgs::a_hyd*eos_cgs::m_u);
  const int iz = ixe;
  auto u0_ = u0;
  auto w0_ = w0;
  auto &wtemp_ = wtemp;
  par_for("xe_tracer", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real kb = 1.380649e-16;
    const Real rho = w0_(m,IDN,k,j,i);
    Real t, ntot, rhoc;
    if (gen) {
      const Real tcode = wtemp_(m,k,j,i);
      t = tcode*eos.temp_cgs;
      const Real mu = eos.MeanMolecularWeight(rho, w0_(m,IEN,k,j,i), tcode);
      ntot = rho*eos.pres_cgs/(mu*eos.temp_cgs*kb);
      rhoc = rho*eos.dens_cgs;
    } else {
      t = gm1*w0_(m,IEN,k,j,i)/rgas/rho;
      ntot = rho*rgas/kb;
      rhoc = rho;
    }
    // a cell without a valid state (an unfilled cs corner ghost: T = 0) keeps its y
    const Real nh = nh_fac*rhoc;
    if (!(t > 0.0) || !(ntot > 0.0) || !(nh > 0.0)) {
      if (init) {
        w0_(m,iz,k,j,i) = 0.0;
        u0_(m,iz,k,j,i) = 0.0;
      }
      return;
    }
    // the tracer is y = n_e/n_H (electrons per H nucleus): n_H is proportional to rho, so
    // rho y is the electron number density times a constant and advection conserves
    // electrons.  (z = n_e/n_tot is NOT per-mass: n_tot/rho changes with mu, ~1.3 -> 2.3
    // between dissociated day and molecular night.)  Saha gives x_eq = n_e/n_tot.
    const Real xeq = eos.ElectronFractionKelvin(rho, t);
    const Real neq = xeq*ntot;
    Real z;   // (= y)
    if (init) {
      z = neq/nh;
    } else {
      const Real a = par.alpha_fac*xetracer::AlphaRR(t)*dt_s;
      Real n = xetracer::Relax(w0_(m,iz,k,j,i)*nh, neq, a);
      const Real ncap = fmax(neq, xetracer::DonorCap(par, t, ntot*kb*t, nh));
      n = fmin(n, ncap);
      z = n/nh;
    }
    w0_(m,iz,k,j,i) = z;
    u0_(m,iz,k,j,i) = u0_(m,IDN,k,j,i)*z;
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::TracerX1Ghosts
//! \brief zero gradient of the SPECIFIC tracer value s/rho (x_e: y = n_e/n_H; floor: the
//! floor mass fraction) in the radial ghosts at physical x1 boundaries (the dhj inner
//! wall and top, hse or outflow1).  The user/physical BCs set the MHD variables of those
//! ghosts but not the scalars, so without this the ghost scalar would keep a stale value
//! (for the floor tracer: every floor increment ConsToPrim ever made in that ghost) and
//! feed it into the reconstruction of the edge cells.  Only the tracer slots are touched;
//! user scalars keep their old treatment.  u0 from the edge cell's u0 ratio; with set_w
//! also w0.

void MHD::TracerX1Ghosts(bool set_w) {
  if (tracers.NTracer() == 0) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb = pmy_pack->nmb_thispack;
  const int is = indcs.is, ie = indcs.ie;
  const int t0 = (ixe >= 0) ? ixe : ifltr;
  const int t1 = (ifltr >= 0) ? ifltr : ixe;
  auto u0_ = u0;
  auto w0_ = w0;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  par_for("mhd_tracer_x1bc", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(int m, int k, int j) {
    const BoundaryFlag fi = mb_bcs.d_view(m, BoundaryFace::inner_x1);
    const BoundaryFlag fo = mb_bcs.d_view(m, BoundaryFace::outer_x1);
    const bool do_in = (fi != BoundaryFlag::block && fi != BoundaryFlag::periodic &&
                        fi != BoundaryFlag::panel);
    const bool do_out = (fo != BoundaryFlag::block && fo != BoundaryFlag::periodic &&
                         fo != BoundaryFlag::panel);
    for (int n=t0; n<=t1; ++n) {
      if (do_in) {
        const Real q = u0_(m,n,k,j,is)/u0_(m,IDN,k,j,is);
        for (int g=1; g<=ng; ++g) {
          u0_(m,n,k,j,is-g) = u0_(m,IDN,k,j,is-g)*q;
          if (set_w) w0_(m,n,k,j,is-g) = q;
        }
      }
      if (do_out) {
        const Real q = u0_(m,n,k,j,ie)/u0_(m,IDN,k,j,ie);
        for (int g=1; g<=ng; ++g) {
          u0_(m,n,k,j,ie+g) = u0_(m,IDN,k,j,ie+g)*q;
          if (set_w) w0_(m,n,k,j,ie+g) = q;
        }
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::FloorTracerDrain
//! \brief <mhd>/floor_tracer: remove floor-added mass (scalar ifltr, s) once per cycle
//! where the density floor does not bind.  With rho_fl = max(dfloor, B^2/v_A,max^2)
//! (the two floors ConsToPrim applies):
//!     dm = min(s, fltr_fmax rho, rho - fltr_margin rho_fl)      (no drain if dm <= 0)
//! so (a) only tagged floor mass is removed (dm <= s), (b) a drained cell stays at least
//! fltr_margin x above the floor (it cannot be pushed back onto the floor, which would
//! re-add the mass the next conversion: a drain/floor pump), (c) at most fltr_fmax of
//! the cell mass goes per cycle, so the pressure change per cycle is <= fltr_fmax.
//! The cell keeps its velocity, its specific internal energy and its specific values
//! of every other scalar: with f = dm/rho, rho, the momenta, e_int and the other
//! scalars' conserved values scale by (1 - f); the total energy loses
//! de = f (e_int + KE) (+ dm phi with etotgrav, whose u0(IEN) contains rho phi); the
//! magnetic energy is untouched.  The kinetic energy is m_i v^i/2 (covariant u0
//! momentum, contravariant w0 velocity on the cubed sphere; orthonormal elsewhere).  w0
//! is updated to the new state; the general-EOS pressure cache scales by (1 - f)
//! (constant T; exact for an ideal gas, O(f) otherwise until the next conversion
//! re-solves it); wtemp is kept.
//! Every cell including ghosts is drained (pointwise; ghosts that mirror an interior
//! cell get the identical result).  The ACTIVE-cell removal is booked twice: into
//! fltr_dm/fltr_de (history Mdrain/Edrain) and, when the floor bookkeeping (c2p_track)
//! is on, into c2p_acc rows 0/1 with weight 1, so Mfloor/Efloor and the energy ledgers
//! see the NET floor change (added - drained).

void MHD::FloorTracerDrain() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int nc1 = indcs.nx1 + 2*ng;
  const int nc2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  const int nc3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmy_pack->nmb_thispack;
  if (c2p_track) c2p_track::Before(pmy_pack, u0, c2p_ref_e, c2p_ref_d);
  auto u0_ = u0;
  auto w0_ = w0;
  auto bcc_ = bcc0;
  auto wder_ = wder;
  auto vol = pmy_pack->pcoord->volume;
  const bool hasder = (wder.size() > 0);
  const bool etg = use_etotgrav;
  DvceArray4D<Real> phi_;
  if (etg) phi_ = phicc0;
  const int isl = ifltr;
  const int ns0 = nmhd, ns1 = nmhd + nscalars;
  const Real dfl = peos->eos_data.dfloor;
  const Real va2i = peos->eos_data.va2max_inv;
  const Real marg = fltr_margin, fmx = fltr_fmax;
  const int ncell = nmb*nc3*nc2*nc1;
  Real sum_dm = 0.0, sum_de = 0.0;
  Kokkos::parallel_reduce("mhd_fltr_drain",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, ncell),
  KOKKOS_LAMBDA(const int idx, Real &adm, Real &ade) {
    const int m = idx/(nc3*nc2*nc1);
    int r = idx - m*(nc3*nc2*nc1);
    const int k = r/(nc2*nc1);
    r -= k*(nc2*nc1);
    const int j = r/nc1;
    const int i = r - j*nc1;
    const Real rho = u0_(m,IDN,k,j,i);
    const Real s = fmin(fmax(u0_(m,isl,k,j,i), 0.0), rho);
    const Real b2 = SQR(bcc_(m,IBX,k,j,i)) + SQR(bcc_(m,IBY,k,j,i))
                    + SQR(bcc_(m,IBZ,k,j,i));
    const Real rfl = fmax(dfl, b2*va2i);
    const Real dm = fmin(s, fmin(fmx*rho, rho - marg*rfl));
    // (no early return: the reducers must see every index)
    if (rho > 0.0 && dm > 0.0) {
      const Real f = dm/rho, g = 1.0 - f;
      // (not `ke`: that name is the last active x3 index, used below)
      const Real ekin = 0.5*(u0_(m,IM1,k,j,i)*w0_(m,IVX,k,j,i)
                             + u0_(m,IM2,k,j,i)*w0_(m,IVY,k,j,i)
                             + u0_(m,IM3,k,j,i)*w0_(m,IVZ,k,j,i));
      Real de = f*(w0_(m,IEN,k,j,i) + ekin);
      if (etg) de += dm*phi_(m,k,j,i);
      u0_(m,IDN,k,j,i) = rho - dm;
      u0_(m,IM1,k,j,i) *= g;
      u0_(m,IM2,k,j,i) *= g;
      u0_(m,IM3,k,j,i) *= g;
      u0_(m,IEN,k,j,i) -= de;
      for (int n=ns0; n<ns1; ++n) {
        if (n != isl) u0_(m,n,k,j,i) *= g;
      }
      const Real snew = s - dm;
      u0_(m,isl,k,j,i) = snew;
      w0_(m,IDN,k,j,i) = rho - dm;
      w0_(m,IEN,k,j,i) *= g;
      w0_(m,isl,k,j,i) = snew/(rho - dm);
      if (hasder) wder_(m,IDPR,k,j,i) *= g;
      if (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke) {
        const Real v = vol(m,k,j,i);
        adm += dm*v;
        ade += de*v;
      }
    }
  }, Kokkos::Sum<Real>(sum_dm), Kokkos::Sum<Real>(sum_de));
  fltr_dm += sum_dm;
  fltr_de += sum_de;
  if (c2p_track) {
    c2p_track::After(pmy_pack, u0, c2p_ref_e, c2p_ref_d, c2p_acc, 1.0, false);
  }
}

}  // namespace mhd
