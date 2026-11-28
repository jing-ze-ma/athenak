//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dfloor_flux.cpp
//! \brief <hydro|mhd>/dfloor_flux_limit: the outflow limiter that keeps rho >= dfloor by
//! construction, and its momentum (DflMomentum), templated over the fluid (Hydro, MHD).
//! Under <mhd> the mass flux it scales carries the gas's (e + p + m.v/2)/rho only: the
//! field's energy flux (Poynting) and the EMFs are not mass-carried and stay; the
//! kinetic energy DflMomentum takes becomes heat at fixed E (|B|^2/2 untouched).

#include <cfloat>
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/gnomonic_kernels.hpp"
#include "coordinates/gnomonic_raisevel.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "bvals/bvals.hpp"
#include "utils/p_ledger.hpp"

namespace hydro {

namespace {
int NFluid(const Hydro &h) {return h.nhydro;}
int NFluid(const mhd::MHD &h) {return h.nmhd;}

//! the outflow limiter's face correction: the mass flux fd at a face whose donor has the
//! factor t < 1 is scaled by t, and the momentum (mr: the donor's m/rho, covariant on the
//! cubed sphere), energy (h: (e + p + m.v/2)/rho, + the face potential under etotgrav)
//! and scalars (sr) that mass carries leave with it; the pressure part of the flux is
//! left alone (well-balanced).  A face the limiter makes massless (t = 0) is CLOSED as
//! StarWallFlux closes a wall: the rest of its flux but the normal momentum (dr: IM1/2/3
//! of an x1/2/3 face), i.e. FOFC's LLF / HLLC dissipative part, is set to 0, and on the
//! cubed sphere the normal covariant flux loses the cos (cc) x tangential part (IM3 / IM2
//! of an x2 / x3 face) it holds.  Both sides of the face hold the same t: exact.
KOKKOS_INLINE_FUNCTION
void DflFace(const DvceArray5D<Real> &f, const int m, const int k, const int j,
             const int i, const Real t, const Real m1, const Real m2, const Real m3,
             const Real h, const DvceArray5D<Real> &w, const int nvar, const int kd,
             const int jd, const int id, const int dr, const Real cc) {
  const Real df = (1.0 - t)*f(m,IDN,k,j,i);
  f(m,IDN,k,j,i) -= df;
  f(m,IM1,k,j,i) -= df*m1;
  f(m,IM2,k,j,i) -= df*m2;
  f(m,IM3,k,j,i) -= df*m3;
  f(m,IEN,k,j,i) -= df*h;
  for (int n = IEN + 1; n < nvar; ++n) f(m,n,k,j,i) -= df*w(m,n,kd,jd,id);
  if (t > 0.0) return;
  f(m,dr,k,j,i) -= cc*f(m,(dr == 2) ? IM3 : IM2,k,j,i);
  for (int n = 0; n < nvar; ++n) {
    if (n != dr) f(m,n,k,j,i) = 0.0;
  }
}
//! the donor's specific values: 0 th, 1-3 m/rho, 4 (e + p + m.v/2)/rho
KOKKOS_INLINE_FUNCTION
Real DflSpec(const int q, const DvceArray5D<Real> &u, const DvceArray5D<Real> &w,
             const DvceArray5D<Real> &wd, const DvceArray5D<Real> &th, const bool gen,
             const Real gm1, const int m, const int k, const int j, const int i) {
  if (q == 0) return th(m,0,k,j,i);
  const Real d = u(m,IDN,k,j,i);
  if (q < 4) return u(m,q,k,j,i)/d;
  const Real e = w(m,IEN,k,j,i), p = gen ? wd(m,IDPR,k,j,i) : gm1*e;
  return (e + p + 0.5*(u(m,IM1,k,j,i)*w(m,IVX,k,j,i) + u(m,IM2,k,j,i)*w(m,IVY,k,j,i) +
                       u(m,IM3,k,j,i)*w(m,IVZ,k,j,i)))/d;
}


//----------------------------------------------------------------------------------------
//! \fn void Hydro::DfloorFluxInit
//! \brief the limiter's arrays and exchanges (first use)

template <class F>
void DflInit(F &flu, MeshBlockPack *pp) {
  auto *pm = pp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng, nmb = pp->nmb_thispack;
  const bool md = pm->multi_d, td = pm->three_d, cs = pm->use_cubed_sphere;
  const int n1 = indcs.nx1 + 2*ng, n2 = md ? indcs.nx2 + 2*ng : 1;
  const int n3 = td ? indcs.nx3 + 2*ng : 1;
  if (flu.dfl_th.extent_int(0) == nmb) return;
  // (a closed face (th = 0) zeroes its whole energy / tangential flux: no diffusive part)
  if (pm->multilevel || NFluid(flu) != 5 || flu.pcond != nullptr ||
      flu.pvisc != nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << "<hydro>/"
              << "dfloor_flux_limit: uniform mesh, adiabatic hydro without conduction "
              << "or viscosity only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Kokkos::realloc(flu.dfl_th, nmb, 1, n3, n2, n1);
  Kokkos::realloc(flu.dfl_thc, nmb, 1, 1, 1, 1);
  flu.pdfl_f0 = new DvceFaceFld4D<Real>("dfl_f0", 1, 1, 1, 1);   // (uflx's shape)
  for (int d = 0; d < 3; ++d) {
    auto &f0 = *flu.pdfl_f0;
    auto &q = (d == 0) ? f0.x1f : ((d == 1) ? f0.x2f : f0.x3f);
    const auto &f = (d == 0) ? flu.uflx.x1f : ((d == 1) ? flu.uflx.x2f : flu.uflx.x3f);
    Kokkos::realloc(q, nmb, f.extent_int(2), f.extent_int(3), f.extent_int(4));
  }
  for (int c = 0; c < 4; ++c) {
    if (c == 3 || c < ((cs || pm->use_spherical_polar) ? 1 : 3)) {
      Kokkos::realloc(flu.dfl_g[c], nmb, 1, n3, n2, n1);
    }
  }
  flu.pbval_dfl->InitializeBuffers(1);
  if (cs) {   // the panel-seam faces (bits: x2 lo / hi, x3 lo / hi); 5, 6: DflMomentum's
    flu.pbval_dfs->InitializeBuffers(7);   // g and donor rho
    flu.pdfl_fs = new DvceFaceFld5D<Real>("dfl_fs", nmb, 7, n3, n2, n1);
    Kokkos::realloc(flu.pdfl_fs->x1f, nmb, 7, 1, 1, 1);
    Kokkos::realloc(flu.dfl_sm, nmb);
    auto hm = Kokkos::create_mirror_view(flu.dfl_sm);
    const int nb_[4] = {8, 12, 24, 28};
    for (int m = 0; m < nmb; ++m) {
      hm(m) = 0;
      for (int q = 0; q < 4; ++q) {
        const auto nb = pp->pmb->nghbr.h_view(m, nb_[q]);
        if (nb.gid >= 0 && nb.panel != pp->pmb->mb_panel.h_view(m)) {
          hm(m) |= (1 << q);
        }
      }
    }
    Kokkos::deep_copy(flu.dfl_sm, hm);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::DfloorFluxLimit
//! \brief <hydro>/dfloor_flux_limit: the stage's outflows of every cell are scaled so
//! that the RK update cannot take rho below dfloor (inflows ignored, so the bound holds
//! for any inflow; convex RK combinations keep it).  Per cell,
//! th = min(1, (gam0 rho0 + gam1 rho1 - dfloor)/(beta dt sum_out A F/V)); the factors
//! go to the ghosts by one exchange, and every face takes its donor cell's factor and
//! specific values, so the two blocks of a block face hold the same flux and mass,
//! momentum and energy stay exact (no mass is created).  On the cubed sphere the donor
//! across a panel seam (resampled ghosts) comes from the seam record pdfl_fs (seam_raw,
//! the m/rho pair re-expressed on the receiver's basis).  Called by RKUpdate, on the
//! final fluxes (the seam fluxes already averaged).

template <class F>
void DflLimit(F &flu, MeshBlockPack *pp, Driver *pdrive, int stage) {
  auto *pm = pp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pp->nmb_thispack;
  const bool md = pm->multi_d, td = pm->three_d, cs = pm->use_cubed_sphere;
  DflInit(flu, pp);
  const Real g0 = pdrive->gam0[stage-1], g1 = pdrive->gam1[stage-1];
  const Real bdt = pdrive->beta[stage-1]*pm->dt;
  const Real dfl = flu.peos->eos_data.dfloor;
  const bool curv = cs || pm->use_spherical_polar;
  auto &vol = pp->pcoord->volume;
  auto &a1 = pp->pcoord->area.x1f;
  auto &a2 = pp->pcoord->area.x2f;
  auto &a3 = pp->pcoord->area.x3f;
  auto &size = pp->pmb->mb_size;
  auto f1 = flu.uflx.x1f, f2 = flu.uflx.x2f, f3 = flu.uflx.x3f;
  auto u0_ = flu.u0, u1_ = flu.u1, th = flu.dfl_th;
  Kokkos::deep_copy(th, 1.0);
  Real nlim = 0.0;
  Kokkos::parallel_reduce("dfl_th", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(),
                          {0, ks, js, is}, {nmb, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sn) {
    const auto &sz = size.d_view(m);
    const Real w1 = curv ? 1.0 : 1.0/sz.dx1, w2 = curv ? 1.0 : 1.0/sz.dx2;
    const Real w3 = curv ? 1.0 : 1.0/sz.dx3;
    Real q = (fmax(f1(m,IDN,k,j,i+1), 0.0)*(curv ? a1(m,k,j,i+1) : w1) +
              fmax(-f1(m,IDN,k,j,i), 0.0)*(curv ? a1(m,k,j,i) : w1));
    if (md) q += (fmax(f2(m,IDN,k,j+1,i), 0.0)*(curv ? a2(m,k,j+1,i) : w2) +
                  fmax(-f2(m,IDN,k,j,i), 0.0)*(curv ? a2(m,k,j,i) : w2));
    if (td) q += (fmax(f3(m,IDN,k+1,j,i), 0.0)*(curv ? a3(m,k+1,j,i) : w3) +
                  fmax(-f3(m,IDN,k,j,i), 0.0)*(curv ? a3(m,k,j,i) : w3));
    q *= bdt/(curv ? vol(m,k,j,i) : 1.0);
    // (less the round-off of the update itself, so the stored rho is >= dfloor)
    const Real tot = g0*u0_(m,IDN,k,j,i) + g1*u1_(m,IDN,k,j,i);
    const Real room = tot - dfl - 8.0*DBL_EPSILON*(tot + q);
    if (q > room && q > 0.0) {   // (a cell without outflow keeps th = 1)
      th(m,0,k,j,i) = (room > 0.0) ? room/q : 0.0;
      sn += 1.0;
    }
  }, Kokkos::Sum<Real>(nlim));
  flu.dfl_n += nlim;
  while (flu.pbval_dfl->InitRecv(1) == TaskStatus::incomplete) {}
  while (flu.pbval_dfl->PackAndSendCC(th, flu.dfl_thc) == TaskStatus::incomplete) {}
  while (flu.pbval_dfl->RecvAndUnpackCC(th, flu.dfl_thc) == TaskStatus::incomplete) {}
  while (flu.pbval_dfl->ClearSend() == TaskStatus::incomplete) {}
  while (flu.pbval_dfl->ClearRecv() == TaskStatus::incomplete) {}
  auto w0_ = flu.w0, wd = flu.wder;
  const int nvar = NFluid(flu) + flu.nscalars;
  const bool gen = flu.peos->eos_data.IsGeneral(), eg = flu.use_etotgrav;
  const Real gm1 = flu.peos->eos_data.gamma - 1.0;
  auto p1 = flu.phi0.x1f, p2 = flu.phi0.x2f, p3 = flu.phi0.x3f;
  auto cxi = cs ? pp->pcoord->cos_face_xi : DvceArray3D<Real>("dfl_c", 1, 1, 1);
  auto cet = cs ? pp->pcoord->cos_face_eta : cxi;
  DvceArray5D<Real> r2 = f2, r3 = f3;   // (the seam record's faces on the cubed sphere)
  DvceArray1D<int> sm;
  if (cs) {
    r2 = flu.pdfl_fs->x2f;
    r3 = flu.pdfl_fs->x3f;
    sm = flu.dfl_sm;
    par_for("dfl_rec2", DevExeSpace(), 0, nmb - 1, ks, ke, 0, 1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
      const int jf = sd ? je + 1 : js, jc = sd ? je : js;
      for (int q = 0; q < 5; ++q) {
        r2(m,q,k,jf,i) = DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, k, jc, i);
      }
    });
    par_for("dfl_rec3", DevExeSpace(), 0, nmb - 1, 0, 1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
      const int kf = sd ? ke + 1 : ks, kc = sd ? ke : ks;
      for (int q = 0; q < 5; ++q) {
        r3(m,q,kf,j,i) = DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, kc, j, i);
      }
    });
    flu.pbval_dfs->InitFluxSeamRecv(7);
    flu.pbval_dfs->PackAndSendFluxSeamCC(*flu.pdfl_fs);
    auto &fs = *flu.pdfl_fs;
    while (flu.pbval_dfs->RecvAndUnpackFluxSeamCC(fs) != TaskStatus::complete) {}
    flu.pbval_dfs->ClearFluxSend();
    flu.pbval_dfs->ClearFluxRecv();
  }
  // (the mass fluxes before the limiter: DflMomentum's outflow weights)
  const int idn = IDN;
  auto all = Kokkos::ALL;
  auto &f0 = *flu.pdfl_f0;
  Kokkos::deep_copy(DevExeSpace(), f0.x1f, Kokkos::subview(f1, all, idn, all, all,
                                                                 all));
  if (md) Kokkos::deep_copy(DevExeSpace(), f0.x2f, Kokkos::subview(f2, all, idn, all,
                                                                         all, all));
  if (td) Kokkos::deep_copy(DevExeSpace(), f0.x3f, Kokkos::subview(f3, all, idn, all,
                                                                         all, all));
  // dr = 1, 2, 3: the face's direction; the donor is (kd, jd, id); a seam face's outer
  // donor (the ghost side) comes from the record
  auto face = [=] KOKKOS_FUNCTION (const DvceArray5D<Real> &f, const DvceArray5D<Real> &r,
                                   const int dr, const int m, const int k, const int j,
                                   const int i, const Real phif, const Real cc) {
    const bool up = (f(m,IDN,k,j,i) > 0.0);   // the donor is the lower cell
    const int kd = (dr == 3 && up) ? k - 1 : k, jd = (dr == 2 && up) ? j - 1 : j;
    const int id = (dr == 1 && up) ? i - 1 : i;
    const int bit = (dr == 2) ? ((j == js) ? 1 : ((j == je + 1) ? 2 : 0))
                              : ((dr == 3) ? ((k == ks) ? 4 : ((k == ke + 1) ? 8 : 0)) : 0);
    const bool srec = cs && bit != 0 && (sm(m) & bit) &&
                      (((bit & 5) && up) || ((bit & 10) && !up));
    Real v[5];
    for (int q = 0; q < 5; ++q) {
      v[q] = srec ? r(m,q,k,j,i) : DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, kd, jd, id);
    }
    if (!(v[0] < 1.0)) return;
    DflFace(f, m, k, j, i, v[0], v[1], v[2], v[3], v[4] + phif, w0_, nvar, kd, jd, id, dr,
            cc);
  };
  par_for("dfl_f1", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    face(f1, f1, 1, m, k, j, i, eg ? p1(m,k,j,i) : 0.0, 0.0);
  });
  if (md) {
    par_for("dfl_f2", DevExeSpace(), 0, nmb - 1, ks, ke, js, je + 1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      face(f2, r2, 2, m, k, j, i, eg ? p2(m,k,j,i) : 0.0, cs ? cxi(m,k,j) : 0.0);
    });
  }
  if (td) {
    par_for("dfl_f3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      face(f3, r3, 3, m, k, j, i, eg ? p3(m,k,j,i) : 0.0, cs ? cet(m,k,j) : 0.0);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::DflMomentum
//! \brief <hydro>/dfloor_flux_limit, the momentum: a cell whose outflow was scaled by
//! th < 1 kept that mass, so its velocity is scaled to the mass flux it really had: the
//! momentum but its L_z part pf phi_hat (EkinScale's, floor_from_ekin) times th, E fixed
//! (the kinetic energy taken becomes heat).  Otherwise gravity and pressure keep
//! accelerating held gas whose work (the potential and the enthalpy its mass flux
//! carries, etotgrav) never arrives: KE > E, floors, infall to the velocity ceiling.
//! The momentum taken leaves by a face floor flux: its part along P (Cartesian: all;
//! spherical polar, cubed sphere: m_r) goes to the cells its outflow went to, shared by
//! the outflows A F before the limiter (pdfl_f0; g = P taken / sum, one exchange and on
//! the cubed sphere the seam record, so both sides of a face move the same P), E fixed;
//! a physical boundary face's share leaves the mesh (p_ledger: hflx).  After every
//! momentum source of the stage; the x1 momentum and the heat are booked.

template <class F>
void DflMomentumT(F &flu, MeshBlockPack *pp) {
  auto *pm = pp->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks, nx1 = indcs.nx1;
  const int ie = indcs.ie, je = indcs.je, ke = indcs.ke;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3, nmb = pp->nmb_thispack;
  const bool cs = pm->use_cubed_sphere, sp = pm->use_spherical_polar, curv = cs || sp;
  const bool md = pm->multi_d, td = pm->three_d;
  const int nc = curv ? 1 : 3;
  auto &vol = pp->pcoord->volume;
  auto a1 = pp->pcoord->area.x1f, a2 = pp->pcoord->area.x2f;
  auto a3 = pp->pcoord->area.x3f;
  auto cosc = cs ? pp->pcoord->cos_cell : DvceArray3D<Real>("dfl_c", 1, 1, 1);
  auto &size = pp->pmb->mb_size;
  auto pnl = pp->pmb->mb_panel.d_view;
  auto bcs = pp->pmb->mb_bcs;
  auto u0_ = flu.u0, th = flu.dfl_th, w0_ = flu.w0;
  auto q1 = flu.pdfl_f0->x1f, q2 = flu.pdfl_f0->x2f, q3 = flu.pdfl_f0->x3f;
  for (int c = 0; c < 4; ++c) {
    if (c == 3 || c < nc) Kokkos::deep_copy(DevExeSpace(), flu.dfl_g[c], 0.0);
  }
  auto *dg = flu.dfl_g;
  auto g0 = dg[0], g1 = curv ? g0 : dg[1], g2 = curv ? g0 : dg[2], g3 = dg[3];
  // the outflow weight of the face of cell (k,j,i) along d (0..2), side s (0 lo, 1 hi):
  // A F before the limiter out of the cell (Cartesian per volume); bd: physical face
  auto wout = [=] KOKKOS_FUNCTION (const int m, const int d, const int s, const int k,
                                   const int j, const int i, bool &bd) {
    const int kf = k + s*(d == 2), jf = j + s*(d == 1), iff = i + s*(d == 0);
    const auto &sz = size.d_view(m);
    const Real f = (d == 0) ? q1(m,kf,jf,iff) : ((d == 1) ? q2(m,kf,jf,iff)
                                                          : q3(m,kf,jf,iff));
    const Real a = curv ? ((d == 0) ? a1(m,kf,jf,iff) : ((d == 1) ? a2(m,kf,jf,iff)
                                                                 : a3(m,kf,jf,iff)))
                        : 1.0/((d == 0) ? sz.dx1 : ((d == 1) ? sz.dx2 : sz.dx3));
    const int e = (d == 0) ? (s ? ie : is) : ((d == 1) ? (s ? je : js) : (s ? ke : ks));
    const int x = (d == 0) ? i : ((d == 1) ? j : k);
    const BoundaryFlag fl = bcs.d_view(m, 2*d + s);
    bd = x == e && fl != BoundaryFlag::block && fl != BoundaryFlag::periodic &&
         fl != BoundaryFlag::panel && fl != BoundaryFlag::polar;
    return fmax(s ? f : -f, 0.0)*a;
  };
  const int nd = td ? 3 : (md ? 2 : 1);
  Real dp = 0.0, dk = 0.0, ob[3] = {0.0, 0.0, 0.0};
  Kokkos::parallel_reduce("dfl_mom", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nx3*nx2*nx1),
  KOKKOS_LAMBDA(const int n, Real &sdp, Real &sk, Real &b0, Real &b1, Real &b2) {
    const int m = n/(nx3*nx2*nx1), k = (n/(nx2*nx1))%nx3 + ks;
    const int j = (n/nx1)%nx2 + js, i = n%nx1 + is;
    const Real t = th(m,0,k,j,i);
    if (!(t < 1.0)) return;
    const auto &sz = size.d_view(m);
    const Real x2 = CellCenterX(j - js, nx2, sz.x2min, sz.x2max);
    const Real x3 = CellCenterX(k - ks, nx3, sz.x3min, sz.x3max);
    // phi_hat: spherical polar e3, the cubed sphere its covariant (0, f2, f3),
    // Cartesian (-y, x)/R of the cell
    Real f[3] = {0.0, 0.0, sp ? 1.0 : 0.0};
    if (cs) {
      GnomonicPhiHat(pnl(m), 0.25*M_PI*x2, 0.25*M_PI*x3, f[1], f[2]);
    } else if (!sp) {
      const Real x1 = CellCenterX(i - is, nx1, sz.x1min, sz.x1max);
      const Real r = sqrt(x1*x1 + x2*x2);
      if (r > 0.0) {f[0] = -x2/r; f[1] = x1/r;}
    }
    const Real c = cs ? cosc(m,k,j) : 0.0, det = 1.0 - c*c;
    const Real m1 = u0_(m,IM1,k,j,i), m2 = u0_(m,IM2,k,j,i), m3 = u0_(m,IM3,k,j,i);
    const Real msq = cs ? m1*m1 + (m2*m2 - 2.0*c*m2*m3 + m3*m3)/det
                        : m1*m1 + m2*m2 + m3*m3;
    const Real pf = cs ? (m2*(f[1] - c*f[2]) + m3*(f[2] - c*f[1]))/det
                       : m1*f[0] + m2*f[1] + m3*f[2];
    u0_(m,IM1,k,j,i) = t*(m1 - pf*f[0]) + pf*f[0];
    u0_(m,IM2,k,j,i) = t*(m2 - pf*f[1]) + pf*f[1];
    u0_(m,IM3,k,j,i) = t*(m3 - pf*f[2]) + pf*f[2];
    const Real v = curv ? vol(m,k,j,i) : 1.0;
    sdp += (m1 - u0_(m,IM1,k,j,i))*v*(curv ? 1.0 : sz.dx1*sz.dx2*sz.dx3);
    sk += 0.5*(1.0 - t*t)*fmax(msq - pf*pf, 0.0)/u0_(m,IDN,k,j,i)*
          (curv ? v : sz.dx1*sz.dx2*sz.dx3);
    // (the P taken per outflow weight; the physical faces' share leaves)
    Real qa = 0.0, qb = 0.0;
    for (int d = 0; d < nd; ++d) {
      for (int s = 0; s < 2; ++s) {
        bool bd;
        const Real w = wout(m, d, s, k, j, i, bd);
        qa += w;
        if (bd) qb += w;
      }
    }
    if (!(qa > 0.0)) return;
    const Real dm[3] = {m1 - u0_(m,IM1,k,j,i), m2 - u0_(m,IM2,k,j,i),
                        m3 - u0_(m,IM3,k,j,i)};
    g0(m,0,k,j,i) = dm[0]*v/qa;
    g3(m,0,k,j,i) = u0_(m,IDN,k,j,i);
    if (!curv) {g1(m,0,k,j,i) = dm[1]/qa; g2(m,0,k,j,i) = dm[2]/qa;}
    const Real bw = qb/qa*(curv ? v : sz.dx1*sz.dx2*sz.dx3);
    b0 += dm[0]*bw;
    if (!curv) {b1 += dm[1]*bw; b2 += dm[2]*bw;}
  }, Kokkos::Sum<Real>(dp), Kokkos::Sum<Real>(dk), Kokkos::Sum<Real>(ob[0]),
     Kokkos::Sum<Real>(ob[1]), Kokkos::Sum<Real>(ob[2]));
  flu.dfl_dp[0] += dp;
  flu.dfl_dp[1] += dk;
  if (!flu.dfl_pass) return;   // (dfloor_flux_pass = false: the momentum taken vanishes)
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, ob, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  p_ledger::Outflow(ob);
  for (int c = 0; c < 4; ++c) {   // (g and rho of the donors into the ghosts)
    if (c != 3 && c >= nc) continue;
    while (flu.pbval_dfl->InitRecv(1) == TaskStatus::incomplete) {}
    auto &gc = flu.dfl_g[c];
    while (flu.pbval_dfl->PackAndSendCC(gc, flu.dfl_thc) == TaskStatus::incomplete) {}
    while (flu.pbval_dfl->RecvAndUnpackCC(gc, flu.dfl_thc) == TaskStatus::incomplete) {}
    while (flu.pbval_dfl->ClearSend() == TaskStatus::incomplete) {}
    while (flu.pbval_dfl->ClearRecv() == TaskStatus::incomplete) {}
  }
  DvceArray5D<Real> r2 = g0, r3 = g0;
  DvceArray1D<int> sm;
  if (cs) {   // (a seam face's donor across the panel seam: its g by the seam record)
    r2 = flu.pdfl_fs->x2f;
    r3 = flu.pdfl_fs->x3f;
    sm = flu.dfl_sm;
    par_for("dfl_grec2", DevExeSpace(), 0, nmb - 1, ks, ke, 0, 1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
      r2(m,5,k,sd ? je + 1 : js,i) = g0(m,0,k,sd ? je : js,i);
      r2(m,6,k,sd ? je + 1 : js,i) = g3(m,0,k,sd ? je : js,i);
    });
    par_for("dfl_grec3", DevExeSpace(), 0, nmb - 1, 0, 1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
      r3(m,5,sd ? ke + 1 : ks,j,i) = g0(m,0,sd ? ke : ks,j,i);
      r3(m,6,sd ? ke + 1 : ks,j,i) = g3(m,0,sd ? ke : ks,j,i);
    });
    flu.pbval_dfs->InitFluxSeamRecv(7);
    flu.pbval_dfs->PackAndSendFluxSeamCC(*flu.pdfl_fs);
    auto &fs = *flu.pdfl_fs;
    while (flu.pbval_dfs->RecvAndUnpackFluxSeamCC(fs) != TaskStatus::complete) {}
    flu.pbval_dfs->ClearFluxSend();
    flu.pbval_dfs->ClearFluxRecv();
  }
  // every cell at least as dense as the donor takes g of each donor whose outflow came
  // in, times that weight (a lighter cell does not: that P is not moved, counted), so no
  // floor gas is driven by the held gas's momentum
  Real drop = 0.0;
  const int nji = nx2*nx1;
  Kokkos::parallel_reduce("dfl_gin", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nx3*nji),
  KOKKOS_LAMBDA(const int idx, Real &sd) {
    const int m = idx/(nx3*nji), k = (idx/nji)%nx3 + ks, j = (idx/nx1)%nx2 + js;
    const int i = idx%nx1 + is;
    const auto &sz = size.d_view(m);
    const Real vol0 = sz.dx1*sz.dx2*sz.dx3;   // (Cartesian: the weights are per volume)
    Real in[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < nd; ++d) {
      for (int s = 0; s < 2; ++s) {
        const int o = 2*s - 1;   // the neighbour across the face
        const int kn = k + o*(d == 2), jn = j + o*(d == 1), in_ = i + o*(d == 0);
        bool bd, bn;   // (its outflow into this cell: the face seen from the neighbour)
        (void) wout(m, d, s, k, j, i, bd);
        const Real w = wout(m, d, 1 - s, kn, jn, in_, bn);
        if (bd || !(w > 0.0)) continue;
        const int bit = (d == 1) ? ((j == js && !s) ? 1 : ((j == je && s) ? 2 : 0)) :
                        ((d == 2) ? ((k == ks && !s) ? 4 : ((k == ke && s) ? 8 : 0)) : 0);
        const bool srec = cs && bit != 0 && (sm(m) & bit);
        const int kf = k + s*(d == 2), jf = j + s*(d == 1), iff = i + s*(d == 0);
        const Real rd = srec ? ((d == 1) ? r2(m,6,kf,jf,iff) : r3(m,6,kf,jf,iff))
                             : g3(m,0,kn,jn,in_);
        const Real ga = srec ? ((d == 1) ? r2(m,5,kf,jf,iff) : r3(m,5,kf,jf,iff))
                             : g0(m,0,kn,jn,in_);
        if (u0_(m,IDN,k,j,i) < rd) {
          sd += w*fabs(ga)*(curv ? 1.0 : vol0);
          continue;
        }
        in[0] += w*ga;
        if (!curv) {in[1] += w*g1(m,0,kn,jn,in_); in[2] += w*g2(m,0,kn,jn,in_);}
      }
    }
    const Real iv = curv ? 1.0/vol(m,k,j,i) : 1.0;
    // (the kinetic energy it adds at most 1e-2 of the cell's eint: else scaled, the rest
    // not moved, counted)
    const Real d = u0_(m,IDN,k,j,i), dq[3] = {in[0]*iv, in[1], in[2]};
    Real a = 0.0, b = 0.0;
    for (int n = 0; n < nc; ++n) {
      a += 0.5*dq[n]*dq[n]/d;
      b += u0_(m,IM1+n,k,j,i)*dq[n]/d;
    }
    const Real cap = 1.0e-2*fmax(w0_(m,IEN,k,j,i), 0.0);
    const Real s = (a + b <= cap) ? 1.0 : ((a > 0.0) ? fmin(fmax((-b + sqrt(fmax(b*b +
                   4.0*a*cap, 0.0)))/(2.0*a), 0.0), 1.0) : 0.0);
    if (s < 1.0) sd += (1.0 - s)*sqrt(2.0*a*d)/iv*(curv ? 1.0 : vol0);
    u0_(m,IM1,k,j,i) += s*dq[0];
    if (!curv) {u0_(m,IM2,k,j,i) += s*dq[1]; u0_(m,IM3,k,j,i) += s*dq[2];}
  }, Kokkos::Sum<Real>(drop));
  flu.dfl_dp[2] += drop;
}

}  // namespace

void Hydro::DfloorFluxInit() {DflInit(*this, pmy_pack);}
void Hydro::DfloorFluxLimit(Driver *d, int stage) {
  DflLimit(*this, pmy_pack, d, stage);
}
void Hydro::DflMomentum() {DflMomentumT(*this, pmy_pack);}

}  // namespace hydro

namespace mhd {
void MHD::DfloorFluxInit() {hydro::DflInit(*this, pmy_pack);}
void MHD::DfloorFluxLimit(Driver *d, int stage) {
  hydro::DflLimit(*this, pmy_pack, d, stage);
}
void MHD::DflMomentum() {hydro::DflMomentumT(*this, pmy_pack);}
}  // namespace mhd
