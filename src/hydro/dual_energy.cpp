//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dual_energy.cpp
//! \brief see hydro/dual_energy.hpp.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "mesh/meshblock_pack.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_hyd.hpp"
#include "hydro/dual_energy.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace fluid {

namespace {
bool Named(ParameterInput *pin, const std::string &b, const std::string &k) {
  return pin->DoesParameterExist(b, k);
}
void Fatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Read

void DualEnergy::Read(ParameterInput *pin, const std::string &block, const EOS_Data &eos,
                      int nfluid, int du_off) {
  blk = block;
  on = Named(pin, blk, "dual_energy") && pin->GetBoolean(blk, "dual_energy");
  // floor_pay: default rad_column for a fresh run with eos = general and <rad_vet>
  // (recorded; a restart file without the key keeps none), else none
  const bool rad = eos.IsGeneral() && pin->DoesBlockExist("rad_vet") &&
                   !global_variable::restart_run;
  if (rad || Named(pin, blk, "floor_pay")) {
    const std::string s = pin->GetOrAddString(blk, "floor_pay", "rad_column");
    if (s == "rad" || s == "rad_column") {
      pay = true;
      pcol = (s == "rad_column");
    } else if (s != "none") {
      Fatal("<" + blk + ">/floor_pay = " + s + ": use rad, rad_column or none");
    }
  }
  // floor_pay_retry (read only when named; default false): energy the radiation cannot
  // pay is not created, the cycle is taken again with a smaller dt (time/dt_retry)
  if (Named(pin, blk, "floor_pay_retry")) {
    pretry = pin->GetBoolean(blk, "floor_pay_retry");
    if (pretry && !pay) Fatal("<" + blk + ">/floor_pay_retry needs floor_pay = rad*");
  }
  // floor_pay_skip (read only when named; default false): the dfloor / vceil cells too
  if (Named(pin, blk, "floor_pay_skip")) {
    pskip = pin->GetBoolean(blk, "floor_pay_skip");
    if (pskip && !pay) Fatal("<" + blk + ">/floor_pay_skip needs floor_pay = rad*");
  }
  retest = Named(pin, blk, "fofc_retest") && pin->GetBoolean(blk, "fofc_retest");
  if (Named(pin, blk, "floor_audit")) report = pin->GetInteger(blk, "floor_audit");
  if (Named(pin, blk, "dual_energy_eta")) eta = pin->GetReal(blk, "dual_energy_eta");
  any = on || pay || retest || (report > 0);
  if (!any) return;
  if (!eos.IsGeneral()) {
    Fatal("<" + blk + ">/dual_energy, floor_pay, fofc_retest, floor_audit need eos = "
          "general (the gamma-law mode of the general EOS reproduces the ideal gas)");
  }
  if (on && (eos.efloor_from_ekin || eos.ekin_heat)) {
    Fatal("<" + blk + ">/dual_energy with efloor_from_ekin / floor_from_ekin is not "
          "supported");
  }
  if (on && !(eta > 0.0 && eta < 1.0)) {
    Fatal("<" + blk + ">/dual_energy_eta must be in (0, 1)");
  }
  if (on) {
    if (du_off < 0) Fatal("dual_energy: no slots in the tracer layout");
    iden = nfluid + du_off;
    ivol = iden + 1;
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<" << blk << "> dual_energy " << on << " (eta " << eta
              << ", slots " << iden << "," << ivol << ") floor_pay "
              << (pay ? (pcol ? "rad_column" : "rad") : "none")
              << " fofc_retest " << retest << " floor_audit " << report
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Allocate

void DualEnergy::Allocate(MeshBlockPack *pp, const DvceArray5D<Real> &u0) {
  const int nmb = u0.extent_int(0);
  const int n3 = u0.extent_int(2), n2 = u0.extent_int(3), n1 = u0.extent_int(4);
  if (init && pre.extent_int(0) == nmb && pre.extent_int(1) == 5 &&
      pre.extent_int(2) == n3 &&
      pre.extent_int(3) == n2 && pre.extent_int(4) == n1) return;
  if (pay && pp->pradvet == nullptr) {
    Fatal("<" + blk + ">/floor_pay = rad needs <rad_vet>");
  }
  Kokkos::realloc(pre, nmb, 5, n3, n2, n1);
  Kokkos::realloc(pst, nmb, n3, n2, n1);
  if (pay) {
    Kokkos::realloc(payv, nmb, n3, n2, n1);
    Kokkos::deep_copy(payv, 0.0);
  }
  init = true;
}

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Pre

void DualEnergy::Pre(MeshBlockPack *pp, const DvceArray5D<Real> &u0) {
  if (!any) return;
  Allocate(pp, u0);
  const int nmb = u0.extent_int(0);
  const int n3m1 = u0.extent_int(2) - 1, n2m1 = u0.extent_int(3) - 1;
  const int n1m1 = u0.extent_int(4) - 1;
  auto pre_ = pre;
  par_for("duale_pre", DevExeSpace(), 0, nmb-1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    pre_(m,0,k,j,i) = u0(m,IEN,k,j,i);
    pre_(m,1,k,j,i) = u0(m,IDN,k,j,i);
    // the momentum components themselves: a recomputed |m|^2 can differ in the last bit
    // between two kernels (FMA contraction) and would mark untouched cells as changed
    pre_(m,2,k,j,i) = u0(m,IM1,k,j,i);
    pre_(m,3,k,j,i) = u0(m,IM2,k,j,i);
    pre_(m,4,k,j,i) = u0(m,IM3,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Post
//! \brief per cell: the internal energy the ConToPrim started from, e_c = E_pre - (E - e)
//! (E - e = KE (+ EB) exactly as the inversion formed it, metric-correct on the cubed
//! sphere); the dual selection; the audit; the payment record; the slot resync.

void DualEnergy::Post(MeshBlockPack *pp, DvceArray5D<Real> &u0, DvceArray5D<Real> &w0,
                      DvceArray5D<Real> &wder, DvceArray4D<Real> &wtemp,
                      const EOS_Data &eos) {
  if (!any) return;
  const int nmb = u0.extent_int(0);
  const int n3 = u0.extent_int(2), n2 = u0.extent_int(3), n1 = u0.extent_int(4);
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const bool curv = pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar;
  auto vol_ = curv ? pp->pcoord->volume : DvceArray4D<Real>("duale_v", 1, 1, 1, 1);
  auto &size = pp->pmb->mb_size;
  const bool sel_ = on && armed && valid;
  const bool rs_ = on, pay_ = pay, psk_ = pay && pskip, edg_ = (report > 0);
  const Real eta_ = eta;
  const Real wst_ = wst;   // energies are counted and paid with their end-of-step weight
  const int iden_ = iden, ivol_ = ivol;
  auto pre_ = pre;
  auto pst_ = pst;
  auto payv_ = payv;
  const bool grv_ = (phi.extent_int(0) > 0);
  auto phi_ = phi;
  auto eos_ = eos;
  const Real tfl_ = eos.tfloor;
  const Real tmin_ = eos.TableTempMin();
  const Real tmax_ = eos.TableTempMax();
  const int nji = n2*n1, nkji = n3*nji;
  Real s[16] = {};
  auto x1v_ = pp->pcoord->x1v;
  Kokkos::parallel_reduce("duale_post", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &s0, Real &s1, Real &s2, Real &s3, Real &s4,
                Real &s5, Real &s6, Real &s7, Real &s8, Real &s9, Real &s10,
                Real &s11, Real &s12, Real &s13, Real &s14, Real &s15) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji;
    const int j = (idx - m*nkji - k*nji)/n1;
    const int i = idx - m*nkji - k*nji - j*n1;
    const Real d = u0(m,IDN,k,j,i);
    if (!(d > 0.0)) return;
    const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
    const bool same = (d == pre_(m,1,k,j,i)) && (u0(m,IM1,k,j,i) == pre_(m,2,k,j,i)) &&
                      (u0(m,IM2,k,j,i) == pre_(m,3,k,j,i)) &&
                      (u0(m,IM3,k,j,i) == pre_(m,4,k,j,i));
    const Real epost = w0(m,IEN,k,j,i);
    const Real etpost = u0(m,IEN,k,j,i);
    const Real etpre = pre_(m,0,k,j,i);
    const Real mech = etpost - epost;
    const Real econs = (etpost == etpre) ? epost : (etpre - mech);
    Real efin = epost, esel = econs;
    bool dsel = false;
    if (same) {
      const bool kedom = !(econs > 0.0) || (econs < eta_*etpre);
      if (sel_ && kedom) {
        const Real se = u0(m,iden_,k,j,i), sv = u0(m,ivol_,k,j,i);
        if (se > 0.0 && sv > 0.0) {
          const Real ed = se*exp(pst_(m,k,j,i)*(sv - 1.0)/se);
          if (ed > 0.0 && Kokkos::isfinite(ed)) {
            dsel = true;
            esel = ed;
          }
        }
      }
      if (dsel) {
        // the floors of the inversion, on the selected energy (own: E rebuilt below)
        Real eint = esel, ekin = mech, etot = esel + mech, fvk = 1.0;
        Real temp = -1.0, pg = 0.0, g1 = 0.0, dde = 0.0;
        bool efl = false, tfl = false, tcl = false, tset = false, msc = false;
        bool ffd = false, stale = false;
        GeneralHydFloorsPost(eos_, d, 1.0/d, wtemp(m,k,j,i), true, false, eint, ekin,
                             etot, fvk, temp, pg, g1, efl, tfl, tcl, tset, msc, ffd,
                             stale, dde);
        w0(m,IEN,k,j,i) = eint;
        u0(m,IEN,k,j,i) = eint + mech;
        wder(m,IDPR,k,j,i) = pg;
        wder(m,IDG1,k,j,i) = g1;
        wtemp(m,k,j,i) = temp;
        efin = eint;
      }
      if (act) {
        const Real vol = curv ? vol_(m,k,j,i)
                              : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        if (kedom) s0 += 1.0;
        if (!(econs > 0.0)) s1 += 1.0;
        if (dsel) {
          s2 += 1.0;
          s3 += wst_*(esel - econs)*vol;
        }
        if (efin != esel) {
          s4 += 1.0;
          s5 += wst_*(efin - esel)*vol;
          if (kedom) s10 += 1.0;
          s12 += x1v_(m,i);   // (where: the floored cells' mean r, and energy weighted)
          s13 += wst_*(efin - esel)*vol*x1v_(m,i);
        }
        if (edg_) {
          if (tfl_ > 0.0 && econs < eos_.EnergyFromTemperature(d, tfl_)) s6 += 1.0;
          if (tmin_ > 0.0 && econs < eos_.EnergyFromTemperature(d, tmin_)) s7 += 1.0;
          if (tmax_ > 0.0 && efin > eos_.EnergyFromTemperature(d, tmax_)) s11 += 1.0;
          if (dsel && tfl_ > 0.0 && esel < eos_.EnergyFromTemperature(d, tfl_)) s9 += 1.0;
        }
      }
      // every cell, ghosts included: a ghost holds its neighbour's state and makes the
      // same selection and floors, so its payment mirrors the neighbour's (the radiation
      // ghosts then stay consistent without another exchange)
      if (pay_) payv_(m,k,j,i) += wst_*(efin - econs);
    } else {
      // a density floor / velocity ceiling cell: its total energy (with the floor mass's
      // potential energy under etotgrav) and mass change
      const Real dd = d - pre_(m,1,k,j,i);
      const Real de = etpost - etpre + (grv_ ? dd*phi_(m,k,j,i) : 0.0);
      if (act) {
        const Real vol = curv ? vol_(m,k,j,i)
                              : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
        s8 += 1.0;
        s14 += wst_*de*vol;
        s15 += wst_*dd*vol;
      }
      // (the floor mass's Phi is counted, not paid: it depends on the potential's zero)
      if (psk_) payv_(m,k,j,i) += wst_*(etpost - etpre);
    }
    if (rs_) {
      u0(m,iden_,k,j,i) = efin;
      u0(m,ivol_,k,j,i) = 1.0;
      w0(m,iden_,k,j,i) = efin/d;
      w0(m,ivol_,k,j,i) = 1.0/d;
      pst_(m,k,j,i) = wder(m,IDPR,k,j,i);
    }
  }, Kokkos::Sum<Real>(s[0]), Kokkos::Sum<Real>(s[1]), Kokkos::Sum<Real>(s[2]),
     Kokkos::Sum<Real>(s[3]), Kokkos::Sum<Real>(s[4]), Kokkos::Sum<Real>(s[5]),
     Kokkos::Sum<Real>(s[6]), Kokkos::Sum<Real>(s[7]), Kokkos::Sum<Real>(s[8]),
     Kokkos::Sum<Real>(s[9]), Kokkos::Sum<Real>(s[10]), Kokkos::Sum<Real>(s[11]),
     Kokkos::Sum<Real>(s[12]), Kokkos::Sum<Real>(s[13]), Kokkos::Sum<Real>(s[14]),
     Kokkos::Sum<Real>(s[15]));
  const int b = armed ? 0 : DC_REFRESH;
  hc[b + DC_N] += static_cast<Real>(pp->nmb_thispack)*indcs.nx1*indcs.nx2*indcs.nx3;
  hc[b + DC_KEDOM] += s[0];
  hc[b + DC_LE0] += s[1];
  hc[b + DC_SEL] += s[2];
  hc[b + DC_SEL_DE] += s[3];
  hc[b + DC_FLR] += s[4];
  hc[b + DC_FLR_DE] += s[5];
  hc[b + DC_SUBT] += s[6];
  hc[b + DC_SUBEDGE] += s[7];
  hc[b + DC_SKIP] += s[8];
  hc[b + DC_SUBT_SEL] += s[9];
  hc[b + DC_FLR_KED] += s[10];
  hc[armed ? DC_TOPE_U : DC_TOPE_R] += s[11];
  hc[b + DC_FLR_R] += s[12];
  hc[b + DC_FLR_RDE] += s[13];
  hc[b + DC_SKIP_DE] += s[14];
  hc[b + DC_SKIP_DM] += s[15];
  armed = false;
  if (on) valid = true;
}

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Report

void DualEnergy::Report(MeshBlockPack *pp) {
  if (report <= 0) return;
  const int cyc = pp->pmesh->ncycle;
  if ((cyc % report) != 0 || cyc == last_report) return;
  last_report = cyc;
  Real c[DC_NCNT];
  for (int n = 0; n < DC_NCNT; ++n) c[n] = hc[n];
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, c, DC_NCNT, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank != 0) return;
  std::ostringstream os;   // (its own precision: std::cout's is set elsewhere)
  os.precision(10);
  for (int h = 0; h < 2; ++h) {
    const int b = h ? DC_REFRESH : 0;
    os << "<" << blk << "> floor_audit cycle " << cyc << (h ? " [refresh]" :
              " [update]") << " cell-c2p " << c[b+DC_N] << " e_c<=0 " << c[b+DC_LE0]
              << " kedom " << c[b+DC_KEDOM] << " sub_tfloor " << c[b+DC_SUBT]
              << " sub_edge " << c[b+DC_SUBEDGE] << " dual " << c[b+DC_SEL]
              << " dual_de " << c[b+DC_SEL_DE] << " dual_sub_tfloor " << c[b+DC_SUBT_SEL]
              << " floored " << c[b+DC_FLR] << " (kedom " << c[b+DC_FLR_KED]
              << ") floor_de " << c[b+DC_FLR_DE] << " <r> " << c[b+DC_FLR_R]/
                 fmax(c[b+DC_FLR], 1.0) << " <r>_de " << c[b+DC_FLR_RDE]/
                 (c[b+DC_FLR_DE] != 0.0 ? c[b+DC_FLR_DE] : 1.0) << " skip(dfloor/vceil) "
              << c[b+DC_SKIP] << " skip_de " << c[b+DC_SKIP_DE] << " skip_dm "
              << c[b+DC_SKIP_DM] << " above_table_top " << c[h ? DC_TOPE_R : DC_TOPE_U]
              << std::endl;
  }
  os << "<" << blk << "> floor_audit cycle " << cyc << " pay: paid " << c[DC_PAID]
            << " unpaid " << c[DC_UNPAID] << " given_to_rad " << c[DC_GIVEN]
            << " cells " << c[DC_NPAY] << " | fofc_retest: retested " << c[DC_RETEST]
            << " unflagged " << c[DC_UNFLAG] << std::endl;
  std::cout << os.str() << std::flush;
}

//----------------------------------------------------------------------------------------
//! \fn void DualPayFromRad

Real DualPayFromRad(MeshBlockPack *pp, DualEnergy &de, DvceArray5D<Real> &urad,
                    const int ie, const Real ch_cl, const Real efl) {
  if (!de.pay || !de.init) return 0.0;
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, ie_ = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ni = de.payv.extent_int(3), nj = de.payv.extent_int(2);
  const int nk = de.payv.extent_int(1);
  const int nmb = pp->nmb_thispack;
  const bool curv = pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar;
  auto vol_ = curv ? pp->pcoord->volume : DvceArray4D<Real>("duale_v", 1, 1, 1, 1);
  auto &size = pp->pmb->mb_size;
  auto payv_ = de.payv;
  const int nji = nj*ni, nkji = nk*nji;
  const bool col = de.pcol;
  Real paid = 0.0, unpaid = 0.0, npay = 0.0, given = 0.0;
  Kokkos::parallel_reduce("duale_pay", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &sp, Real &su, Real &sn, Real &sg) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji;
    const int j = (idx - m*nkji - k*nji)/ni;
    const int i = idx - m*nkji - k*nji - j*ni;
    const Real need = payv_(m,k,j,i);
    if (need == 0.0) return;
    payv_(m,k,j,i) = 0.0;
    const bool act = (i >= is && i <= ie_ && j >= js && j <= je && k >= ks && k <= ke);
    const Real vol = curv ? vol_(m,k,j,i)
                          : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    const Real er = urad(m,ie,k,j,i);
    const Real fa = act ? 1.0 : 0.0;
    if (act) sn += 1.0;
    if (need > 0.0) {
      const Real want = ch_cl*need;                   // radiation units
      const Real take = fmin(want, fmax(er - efl, 0.0));
      urad(m,ie,k,j,i) = er - take;
      sp += fa*(take/ch_cl)*vol;
      if (col && act) {   // (rad_column: the rest from the column, below)
        payv_(m,k,j,i) = (want - take)*vol;
      } else {
        su += fa*((want - take)/ch_cl)*vol;
      }
    } else {
      urad(m,ie,k,j,i) = er - ch_cl*need;
      sg += fa*(-need)*vol;
    }
  }, Kokkos::Sum<Real>(paid), Kokkos::Sum<Real>(unpaid), Kokkos::Sum<Real>(npay),
     Kokkos::Sum<Real>(given));
  // floor_pay = rad_column: what the cell's E_rad could not cover is taken from the E_rad
  // above e_floor of the nearest cells of its radial column, i -+ d, the pair split by
  // their energy above the floor (implicit_pos_gas_cons's column pass, CleanWriteBack)
  if (col) {
    Real cp = 0.0, cu = 0.0;
    Kokkos::parallel_reduce("duale_pay_col", Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                            DevExeSpace(), {0, ks, js}, {nmb, ke + 1, je + 1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, Real &sp, Real &su) {
      auto vol = [&](const int q) {
        return curv ? vol_(m,k,j,q)
                    : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
      };
      for (int i = is; i <= ie_; ++i) {
        Real def = payv_(m,k,j,i);
        if (!(def > 0.0)) continue;
        payv_(m,k,j,i) = 0.0;
        for (int d = 1; d <= ie_ - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(urad(m,ie,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie_) ? fmax(urad(m,ie,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          if (aa > 0.0) {
            urad(m,ie,k,j,qa) = fmax(urad(m,ie,k,j,qa) - (tk*aa/at)/vol(qa), efl);
          }
          if (ab > 0.0) {
            urad(m,ie,k,j,qb) = fmax(urad(m,ie,k,j,qb) - (tk*ab/at)/vol(qb), efl);
          }
          def -= tk;
          sp += tk/ch_cl;
        }
        su += def/ch_cl;
      }
    }, Kokkos::Sum<Real>(cp), Kokkos::Sum<Real>(cu));
    paid += cp;
    unpaid += cu;
  }
  // (no silent fallback: energy the column could not give is created, reported; with
  // floor_pay_retry the Driver takes the cycle again instead)
  Real ug = 0.0;
  if (col || de.pretry) {
    ug = unpaid;
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &ug, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (ug > 0.0 && global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d: floor_pay = rad_column left %.6e of the fluid "
                  "floors' energy unpaid (%s)\n", pp->pmesh->ncycle, ug,
                  de.pretry ? "the cycle is taken again" : "created");
    }
  }
  de.hc[DC_PAID] += paid;
  de.hc[DC_UNPAID] += unpaid;
  de.hc[DC_NPAY] += npay;
  de.hc[DC_GIVEN] += given;
  return ug;
}

}  // namespace fluid
