//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_steady.cpp
//! \brief the steady start (clean_init = steady), the first face fluxes of a run and of a
//! restart from the old scheme.
//!
//!  called from: rad_vet_tasks.cpp, driver/driver.cpp
//!  calls into: rad_vet.cpp, rad_vet_bvals.cpp, rad_vet_gas.cpp, rad_vet_geometry.cpp,
//!    rad_vet_init.cpp, rad_vet_tasks.cpp

#include <math.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "driver/driver.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanSteadyStart
//! \brief <rad_vet>/clean_init = steady (fresh start, at the first dt): relaxes the
//! initial radiation to the initial matter and sources.
//! The gas is frozen (no hydro, no feedback) and backward-Euler radiation steps of
//! clean_init_dtfac x dt are taken until max(|dE|, |dF|/c) < clean_init_tol max E; as
//! dt -> infinity the fixed point is the scheme's own steady solution (closure included).
//! A failed pass (or E < 0) is undone and the step shrinks by 1000, down to dt.  The
//! passes then go on at the run's dt, to the tolerance or until the change contracts by
//! less than 1/2 per pass.  Only E and F change, except the gas internal energy with
//! force_reference = clean or icr_fmeas; the closure history is dropped afterwards.

void RadiationVet::CleanSteadyStart(Driver *pdrive) {
  if (!cln_steady) return;
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  const Real dt0 = pm->dt, cl = c_light;
  const bool sfb = gas_feedback;
  const int st2 = t2_solve;
  const Real nc0 = vet_ncall;
  // gas frozen, backward-Euler steps without substeps; big: the large step.  With
  // force_reference = clean or icr_fmeas the gas internal energy relaxes with E (rho,
  // momentum frozen) to radiative equilibrium; icr_fmeas also keeps each pass's x1 force
  // less the reference, dm - h rho a_ref, in icr_dm (icr_dm_dt = h, CleanWriteBack)
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool relax = (force_ref == VET_FREF_CLEAN || icr_fmeas) && sfb && fl.on;
  if (relax && icr_fmeas) {
    Kokkos::realloc(icr_dm, nmb1 + 1, u0.extent(2), u0.extent(3), u0.extent(4));
  }
  auto uh = fl.u0;
  gas_feedback = relax;
  cln_vzero = relax;
  cln_mdcap = (force_ref == VET_FREF_WB_MD);   // every pass's write-back: the last wins
  t2_solve = VET_T2S_NONE;
  nsub = 1;
  bool big = true, ok = false;
  int n2 = 0;
  // icr_fmeas: reach the large step by x10 per pass from dt (the gas solves fail when
  // started at the top)
  const bool lad = relax && icr_fmeas;
  Real dp = 0.0, fac = lad ? 1.0 : cln_steady_dtfac, emax = 0.0;
  cln_insteady = true;
  auto u0_ = u0;
  DvceArray5D<Real> uo("cln_steady_u", nmb1 + 1, 4, u0.extent(2), u0.extent(3),
                       u0.extent(4));
  DvceArray5D<Real> ug("cln_steady_g", relax ? nmb1 + 1 : 1, 4, u0.extent(2),
                       u0.extent(3), u0.extent(4));
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int ntot = (nmb1 + 1)*nk*nj*ni;
  const Real dv = MeshBlock0Dx(0)*MeshBlock0Dx(1)*MeshBlock0Dx(2);
  const bool sph = pm->use_spherical_polar || pm->use_cubed_sphere;   // volume sums
  auto vol = pmy_pack->pcoord->volume;
  // totals of E and F (code units x volume)
  auto sums = [&](Real t[4]) {
    for (int c = 0; c < 4; ++c) {
      Real a = 0.0;
      Kokkos::parallel_reduce("cln_steady_sum", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                              ntot),
      KOKKOS_LAMBDA(const int q, Real &sl) {
        const int m = q/(ni*nj*nk), k = ks + (q/(ni*nj))%nk;
        const int j = js + (q/ni)%nj, i = is + q%ni;
        sl += (sph ? vol(m,k,j,i) : dv)*u0_(m,VET_E+c,k,j,i);
      }, Kokkos::Sum<Real>(a));
      t[c] = CleanReduce(a, false);
    }
  };
  // u0 <-> uo with ghosts, plus face F and the write-back captures (cln_fref, icr_dm)
  DvceArray4D<Real> fo[4];
  DvceArray5D<Real> ro;
  auto copy = [&](const bool save) {
    DvceArray4D<Real> fv[4] = {f0x1, f0x2, f0x3, icr_dm};
    for (int d = 0; d < 4; ++d) {
      if (fo[d].size() != fv[d].size()) fo[d] = DvceArray4D<Real>("cln_steady_f",
                                                                  fv[d].layout());
      if (save) Kokkos::deep_copy(fo[d], fv[d]);
      if (!save) Kokkos::deep_copy(fv[d], fo[d]);
    }
    if (cln_mdcap && ro.size() != cln_fref.size()) {
      ro = DvceArray5D<Real>("cln_steady_r", cln_fref.layout());
    }
    if (cln_mdcap) Kokkos::deep_copy(save ? ro : cln_fref, save ? cln_fref : ro);
    Kokkos::parallel_for("cln_steady_copy", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                         uo.size()),
    KOKKOS_LAMBDA(const int q) {
      const int n4 = uo.extent(4), n3 = uo.extent(3), n2 = uo.extent(2);
      const int i = q%n4, j = (q/n4)%n3, k = (q/(n4*n3))%n2, c = (q/(n4*n3*n2))%4;
      const int m = q/(n4*n3*n2*4);
      const int g = (c == 0) ? IEN : IM1 + c - 1;
      if (save) {
        uo(m,c,k,j,i) = u0_(m,VET_E+c,k,j,i);
        if (relax) ug(m,c,k,j,i) = uh(m,g,k,j,i);
      } else {
        u0_(m,VET_E+c,k,j,i) = uo(m,c,k,j,i);
        if (relax) uh(m,g,k,j,i) = ug(m,c,k,j,i);
      }
    });
    if (relax && !save) FluidRefresh(pdrive);
  };
  Real t0[4], t1[4];
  sums(t0);
  int n = 0;
  Real d = 0.0;
  for (; n < cln_steady_n; ++n) {
    pm->dt = big ? fac*dt0 : dt0;
    dt_sub = pm->dt;
    if (icr_fmeas) icr_dm_dt = pm->dt;
    copy(true);
    cln_tc = -1.0e300;   // predict this pass's closure, not the last dt's
    const Real nf = cln_nfail + cln_nneg;
    // one radiation step: the ordinary m1 task lists (closure, solve, ghost exchange)
    for (int stage = 1; stage <= nstage; ++stage) {
      RunStage(pdrive, stage);
    }
    (void) ApplyClosureLimits(pdrive, 1);
    if (cln_nfail + cln_nneg > nf) {
      copy(false);
      dtr_fail = false;
      cln_have_n = cln_have_p = cln_have_c = cln_have_1 = false;
      if (global_variable::my_rank == 0) {
        std::printf("<rad_vet> clean_init steady pass %d: solve failed or E < 0 at %.0e "
                    "dt, undone\n", n, big ? fac : 1.0);
      }
      if (!big || fac <= 1.0) break;
      fac = fmax(1.0e-3*fac, 1.0);
      continue;
    }
    // its largest change of E and F/c over max E
    Real dm = 0.0, em = 0.0;
    Kokkos::parallel_reduce("cln_steady_chg", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            ntot),
    KOKKOS_LAMBDA(const int q, Real &dl, Real &el) {
      const int m = q/(ni*nj*nk), k = ks + (q/(ni*nj))%nk;
      const int j = js + (q/ni)%nj, i = is + q%ni;
      Real x = fabs(u0_(m,VET_E,k,j,i) - uo(m,0,k,j,i));
      for (int c = 1; c < 4; ++c) {
        x = fmax(x, fabs(u0_(m,VET_E+c,k,j,i) - uo(m,c,k,j,i))/cl);
      }
      dl = fmax(dl, x);
      el = fmax(el, fmax(u0_(m,VET_E,k,j,i), uo(m,0,k,j,i)));
    }, Kokkos::Max<Real>(dm), Kokkos::Max<Real>(em));
    emax = fmax(emax, CleanReduce(em, true));
    d = CleanReduce(dm, true)/fmax(emax, 1.0e-300);
    if (global_variable::my_rank == 0) {
      std::printf("<rad_vet> clean_init steady pass %d: max(|dE|, |dF|/c)/max E = %.3e\n",
                  n, d);
    }
    if (!big && (d < cln_steady_tol || (++n2 > 2 && d > 0.5*dp))) break;
    if (d < cln_steady_tol && (!lad || fac >= cln_steady_dtfac)) big = false;
    if (lad) fac = fmin(10.0*fac, cln_steady_dtfac);
    ok = ok || !big;
    dp = d;
  }
  cln_insteady = false;
  // restore the run's setup; the first step starts from this state without history
  sums(t1);
  if (cln_mdcap) {   // wb_arad_md: the reference in the ghosts, read by the hydro walks
    CleanGhostFill(cln_fref, 0, 9, true, 0);
    if (pm->use_spherical_polar && pm->multi_d) {   // across a pole: no theta halves
      auto fr = cln_fref;
      auto bc = pmy_pack->pmb->mb_bcs.d_view;
      const int n1 = fr.extent_int(4), n3 = fr.extent_int(2), ng = indcs.ng;
      par_for("cln_mdpole", DevExeSpace(), 0, nmb1, 0, n3-1, 0, ng-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int g, const int i) {
        for (int s = 0; s < 2; ++s) {
          if (bc(m,s ? BoundaryFace::outer_x2 : BoundaryFace::inner_x2) !=
              BoundaryFlag::polar) continue;
          const int j = s ? je + 1 + g : js - 1 - g;
          fr(m,5,k,j,i) = fr(m,6,k,j,i) = 0.0;
        }
      });
    }
    if (pm->use_cubed_sphere) {
      // cubed sphere: the lateral walks run along (e_xi, e_eta) and a seam ghost holds
      // the neighbour's axes, so each lateral half is a.(X_face - X_centre)
      auto fr = cln_fref;
      auto pnl = pmy_pack->pmb->mb_panel.d_view;
      auto x1v = pmy_pack->pcoord->x1v, x2v = pmy_pack->pcoord->x2v;
      auto x3v = pmy_pack->pcoord->x3v;
      auto x2f = pmy_pack->pcoord->xx2f, x3f = pmy_pack->pcoord->xx3f;
      const int n1 = fr.extent_int(4), n2 = fr.extent_int(3), n3 = fr.extent_int(2);
      par_for("cln_mdcs", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real c[3], f[3];
        const Real q = 0.25*M_PI;
        cubed_sphere::PanelToCart(pnl(m), q*x2v(m,j), q*x3v(m,k), c);
        for (int sd = 0; sd < 4; ++sd) {   // xi -, +, eta -, +
          const Real a2 = (sd < 2) ? x2f(m,j+sd) : x2v(m,j);
          const Real a3 = (sd < 2) ? x3v(m,k) : x3f(m,k+sd-2);
          cubed_sphere::PanelToCart(pnl(m), q*a2, q*a3, f);
          Real p = 0.0;
          for (int b = 0; b < 3; ++b) p += fr(m,b,k,j,i)*x1v(m,i)*(f[b] - c[b]);
          fr(m,5+sd,k,j,i) = (sd % 2) ? p : -p;
        }
      });
    }
    cln_fref_set = true;
    cln_mdcap = false;
  }
  gas_feedback = sfb;
  cln_vzero = false;
  t2_solve = st2;
  t2_ok = false;
  dtr_fail = false;
  cln_have_n = cln_have_p = cln_have_c = cln_have_1 = false;
  cln_tc = -1.0e300;
  pm->dt = dt0;
  (void) SetSubsteps(dt0);
  if (global_variable::my_rank == 0) {
    std::printf("<rad_vet> clean_init steady: %s after %d passes, residual %.3e "
                "(tol %.1e; %.0f formal solutions); radiation E %.10e -> %.10e, F/c "
                "%.4e %.4e %.4e -> %.4e %.4e %.4e (volume sums)\n",
                ok ? "converged" : "NOT converged",
                std::min(n + 1, cln_steady_n), d, cln_steady_tol, vet_ncall - nc0, t0[0],
                t1[0], t0[1]/cl, t0[2]/cl, t0[3]/cl, t1[1]/cl, t1[2]/cl, t1[3]/cl);
    if (!ok) std::printf("CLEAN WARNING clean_init steady: not converged (residual %.3e, "
                         "tol %.1e)\n", d, cln_steady_tol);
  }
  // a last change that is not small is a diverged start, not an unconverged one
  if (!(d < 1.0e-2)) {
    CleanFatal("clean_init steady: the last pass changed E or F/c by " +
               std::to_string(d) + " of max E (diverged: clean_init_dtfac, _n)");
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanInitFaceFluxes
//! \brief face fluxes as the mean of the two adjacent cell fluxes

void RadiationVet::CleanInitFaceFluxes() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u0_ = u0;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    CleanFor("cln_f0", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      f0(m,k,j,i) = 0.5*(u0_(m,VET_F1+dir,k,j,i) + u0_(m,VET_F1+dir,k-dk,j-dj,i-di));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CleanRstFaces
//! \brief staggered restart from a file without the clean state (the old scheme): its
//! faces hold the comoving normal flux F0 and its cells F = mean(F0) + a E, a = v + v.D,
//! while the staggered faces hold the lab-frame flux.  Each interior face adds the old
//! scheme's enthalpy flux a_f E_f (implicit_enthalpy: plm, central, or upwind donor) with
//! a = (F - mean F0)/E per cell (Cartesian on the cubed sphere, exchanged across seams);
//! physical-boundary faces keep theirs.

void RadiationVet::CleanRstFaces() {
  if (!cln_rstold) return;
  if (!cln_init) CleanInit();
  CleanWoBind();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto u0_ = u0;
  auto w = cln_w;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  auto cw = pmy_pack->pcoord->cos_cell, sw = pmy_pack->pcoord->sin_cell;
  const Real efl = e_floor;
  const bool fre = cln_freal;
  DvceArray1D<Real> rc("cln_rstr", 2);
  // a of every cell (CLW_DM, 3) and its E (CLW_B), as the old scheme's write-back
  par_for("cln_rsta", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real e = fmax(u0_(m,VET_E,k,j,i), efl);
    Real g[3] = {0.5*(f1(m,k,j,i) + f1(m,k,j,i+1)), 0.0, 0.0};
    if (x.nf > 1) g[1] = 0.5*(f2(m,k,j,i) + f2(m,k,j+1,i));
    if (x.nf > 2) g[2] = 0.5*(f3(m,k,j,i) + f3(m,k+1,j,i));
    w(m,CLW_B,k,j,i) = e;
    for (int d = 0; d < 3; ++d) w(m,CLW_DM+d,k,j,i) = (u0_(m,VET_F1+d,k,j,i) - g[d])/e;
    if (!x.cub) return;
    // covariant e_xi, e_eta components (face normals g = sin (a^xi, a^eta))
    const Real c = cw(m,k,j), sn = sw(m,k,j), s2 = 1.0 - c*c;
    const Real a2 = (u0_(m,VET_F2,k,j,i) - c*u0_(m,VET_F3,k,j,i))/s2 - g[1]/sn;
    const Real a3 = (u0_(m,VET_F3,k,j,i) - c*u0_(m,VET_F2,k,j,i))/s2 - g[2]/sn;
    const Real a1 = u0_(m,VET_F1,k,j,i) - g[0];
    for (int b = 0; b < 3; ++b) {
      w(m,CLW_DM+b,k,j,i) = (a1*x.g.At<CLG_R>(m,b,k,j,i) + a2*x.g.At<CLG_E>(m,b,k,j,i) +
                             a3*x.g.At<CLG_E>(m,3+b,k,j,i))/e;
    }
  });
  CleanGhostFill(cln_w, CLW_DM, 4);
  for (int dir = 0; dir < x.nf; ++dir) {
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("cln_rstf", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (CleanBoundaryFace(x, dir, m, k, j, i)) return;
      // cell o = -2, -1 | 0, 1 of the face: its a along the face normal, and E
      auto a = [&](const int o) {
        return CleanStagN(x, CLW_DM, dir, m, k, j, i, k+o*dk, j+o*dj, i+o*di);
      };
      auto e = [&](const int o) {return w(m,CLW_B,k+o*dk,j+o*dj,i+o*di);};
      Real af = 0.5*(a(-1) + a(0)), ef = 0.5*(e(-1) + e(0)), dum;
      if (!CleanBoundaryFace(x, dir, m, k-dk, j-dj, i-di) &&
          !CleanBoundaryFace(x, dir, m, k+dk, j+dj, i+di)) {
        Real afl, afr;
        PLM(a(-2), a(-1), a(0), afl, dum);
        PLM(a(-1), a(0), a(1), dum, afr);
        af = 0.5*(afl + afr);
        if (af > 0.0) PLM(e(-2), e(-1), e(0), ef, dum);
        if (af < 0.0) PLM(e(-1), e(0), e(1), dum, ef);
      }
      f0(m,k,j,i) += af*ef;
      // clean_stag_freal: realisable, |F| <= c E of its donor.  In an old E-floor hole
      // the file's F, F0 are not those of its E (the old scheme would have clipped them),
      // so the excess is dropped: counted and its momentum sum |dF| V/c^2 reported
      const Real fc = fmin(fmax(f0(m,k,j,i), -x.c*e(0)), x.c*e(-1));
      if (fre && fc != f0(m,k,j,i)) {
        const Real v = CleanVol(x, m, k, j, i);
        Kokkos::atomic_add(&rc(0), 1.0);
        Kokkos::atomic_add(&rc(1), fabs(f0(m,k,j,i) - fc)*v/(x.c*x.c));
        f0(m,k,j,i) = fc;
      }
    });
  }
  if (fre) {
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rc);
    const Real n = CleanReduce(h(0), false), p = CleanReduce(h(1), false);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> clean_stag_freal RESTART CONVERSION (one time, not a "
                << "deposit): faces of the old scheme's file above c E_donor taken at it="
                << n << ", radiation momentum removed sum |dF| V/c^2=" << p
                << " (its energy equivalent c sum |dF| V/c^2=" << p*x.c
                << "; E unchanged)" << std::endl;
    }
  }
  // the cell flux of these faces, as CleanWriteBack stores it (cubed sphere: the frame of
  // y, not the old covariant panel components)
  par_for("cln_rstc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), (x.nf > 1) ? f2(m,k,j,i) : 0.0,
                       (x.nf > 1) ? f2(m,k,j+1,i) : 0.0, (x.nf > 2) ? f3(m,k,j,i) : 0.0,
                       (x.nf > 2) ? f3(m,k+1,j,i) : 0.0};
    Real fv[3];
    CleanStagCellVector(x, m, k, j, i, c, fv);
    CleanStagHydro(x, m, k, j, i, fv);
    for (int d = 0; d < 3; ++d) u0_(m,VET_F1+d,k,j,i) = fv[d];
  });
  Kokkos::deep_copy(cln_w.s, 0.0);
  Kokkos::deep_copy(cln_w.o, 0.0);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> clean restart from a file without the clean state: the faces "
              << "take its enthalpy flux a E (lab-frame F)" << std::endl;
  }
}

}  // namespace radvet
