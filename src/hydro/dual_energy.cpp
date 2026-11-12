//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dual_energy.cpp
//! \brief duale-1112: see hydro/dual_energy.hpp.

#include <cmath>
#include <iostream>
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
  if (Named(pin, blk, "floor_pay")) {
    const std::string s = pin->GetString(blk, "floor_pay");
    if (s == "rad") {
      pay = true;
    } else if (s != "none") {
      Fatal("<" + blk + ">/floor_pay = " + s + ": use rad or none");
    }
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
  if (on && eos.efloor_from_ekin) {
    Fatal("<" + blk + ">/dual_energy with efloor_from_ekin is not supported");
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
    std::cout << "<" << blk << "> duale-1112: dual_energy " << on << " (eta " << eta
              << ", slots " << iden << "," << ivol << ") floor_pay " << (pay ? "rad" :
              "none") << " fofc_retest " << retest << " floor_audit " << report
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void DualEnergy::Allocate

void DualEnergy::Allocate(MeshBlockPack *pp, const DvceArray5D<Real> &u0) {
  const int nmb = u0.extent_int(0);
  const int n3 = u0.extent_int(2), n2 = u0.extent_int(3), n1 = u0.extent_int(4);
  if (init && pre.extent_int(0) == nmb && pre.extent_int(2) == n3 &&
      pre.extent_int(3) == n2 && pre.extent_int(4) == n1) return;
  if (pay && pp->pradm1 == nullptr) {
    Fatal("<" + blk + ">/floor_pay = rad needs <rad_m1>");
  }
  Kokkos::realloc(pre, nmb, 3, n3, n2, n1);
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
    pre_(m,2,k,j,i) = SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) +
                      SQR(u0(m,IM3,k,j,i));
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
  const bool rs_ = on, pay_ = pay, edg_ = (report > 0);
  const Real eta_ = eta;
  const int iden_ = iden, ivol_ = ivol;
  auto pre_ = pre;
  auto pst_ = pst;
  auto payv_ = payv;
  auto eos_ = eos;
  const Real tfl_ = eos.tfloor;
  const Real tmin_ = eos.TableTempMin();
  const int nji = n2*n1, nkji = n3*nji;
  Real s[11] = {};
  Kokkos::parallel_reduce("duale_post", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &s0, Real &s1, Real &s2, Real &s3, Real &s4,
                Real &s5, Real &s6, Real &s7, Real &s8, Real &s9, Real &s10) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji;
    const int j = (idx - m*nkji - k*nji)/n1;
    const int i = idx - m*nkji - k*nji - j*n1;
    const Real d = u0(m,IDN,k,j,i);
    if (!(d > 0.0)) return;
    const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
    const Real msq = SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i)) + SQR(u0(m,IM3,k,j,i));
    const bool same = (d == pre_(m,1,k,j,i)) && (msq == pre_(m,2,k,j,i));
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
          s3 += (esel - econs)*vol;
        }
        if (efin != esel) {
          s4 += 1.0;
          s5 += (efin - esel)*vol;
          if (kedom) s10 += 1.0;
        }
        if (edg_) {
          if (tfl_ > 0.0 && econs < eos_.EnergyFromTemperature(d, tfl_)) s6 += 1.0;
          if (tmin_ > 0.0 && econs < eos_.EnergyFromTemperature(d, tmin_)) s7 += 1.0;
          if (dsel && tfl_ > 0.0 && esel < eos_.EnergyFromTemperature(d, tfl_)) s9 += 1.0;
        }
        if (pay_) payv_(m,k,j,i) += efin - econs;
      }
    } else if (act) {
      s8 += 1.0;
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
     Kokkos::Sum<Real>(s[9]), Kokkos::Sum<Real>(s[10]));
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
  const bool was_armed = armed;
  armed = false;
  if (on) valid = true;
  if (was_armed) Report(pp);
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
  for (int h = 0; h < 2; ++h) {
    const int b = h ? DC_REFRESH : 0;
    std::cout << "<" << blk << "> floor_audit cycle " << cyc << (h ? " [refresh]" :
              " [update]") << " cell-c2p " << c[b+DC_N] << " e_c<=0 " << c[b+DC_LE0]
              << " kedom " << c[b+DC_KEDOM] << " sub_tfloor " << c[b+DC_SUBT]
              << " sub_edge " << c[b+DC_SUBEDGE] << " dual " << c[b+DC_SEL]
              << " dual_de " << c[b+DC_SEL_DE] << " dual_sub_tfloor " << c[b+DC_SUBT_SEL]
              << " floored " << c[b+DC_FLR] << " (kedom " << c[b+DC_FLR_KED]
              << ") floor_de " << c[b+DC_FLR_DE] << " skip(dfloor/vceil) "
              << c[b+DC_SKIP] << std::endl;
  }
  std::cout << "<" << blk << "> floor_audit cycle " << cyc << " pay: paid " << c[DC_PAID]
            << " unpaid " << c[DC_UNPAID] << " given_to_rad " << c[DC_GIVEN]
            << " cells " << c[DC_NPAY] << " | fofc_retest: retested " << c[DC_RETEST]
            << " unflagged " << c[DC_UNFLAG] << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn void DualPayFromRad

void DualPayFromRad(MeshBlockPack *pp, DualEnergy &de, DvceArray5D<Real> &urad,
                    const int ie, const Real ch_cl, const Real efl) {
  if (!de.pay || !de.init) return;
  auto &indcs = pp->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.nx1, nj = indcs.nx2, nk = indcs.nx3;
  const int nmb = pp->nmb_thispack;
  const bool curv = pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar;
  auto vol_ = curv ? pp->pcoord->volume : DvceArray4D<Real>("duale_v", 1, 1, 1, 1);
  auto &size = pp->pmb->mb_size;
  auto payv_ = de.payv;
  const int nji = nj*ni, nkji = nk*nji;
  Real paid = 0.0, unpaid = 0.0, npay = 0.0, given = 0.0;
  Kokkos::parallel_reduce("duale_pay", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int &idx, Real &sp, Real &su, Real &sn, Real &sg) {
    const int m = idx/nkji;
    const int k = (idx - m*nkji)/nji + ks;
    const int j = (idx - m*nkji - (k - ks)*nji)/ni + js;
    const int i = (idx - m*nkji - (k - ks)*nji - (j - js)*ni) + is;
    const Real need = payv_(m,k,j,i);
    if (need == 0.0) return;
    payv_(m,k,j,i) = 0.0;
    const Real vol = curv ? vol_(m,k,j,i)
                          : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    sn += 1.0;
    const Real er = urad(m,ie,k,j,i);
    if (need > 0.0) {
      const Real want = ch_cl*need;                   // radiation units
      const Real take = fmin(want, fmax(er - efl, 0.0));
      urad(m,ie,k,j,i) = er - take;
      sp += (take/ch_cl)*vol;
      su += ((want - take)/ch_cl)*vol;
    } else {
      urad(m,ie,k,j,i) = er - ch_cl*need;
      sg += (-need)*vol;
    }
  }, Kokkos::Sum<Real>(paid), Kokkos::Sum<Real>(unpaid), Kokkos::Sum<Real>(npay),
     Kokkos::Sum<Real>(given));
  de.hc[DC_PAID] += paid;
  de.hc[DC_UNPAID] += unpaid;
  de.hc[DC_NPAY] += npay;
  de.hc[DC_GIVEN] += given;
}

}  // namespace fluid
