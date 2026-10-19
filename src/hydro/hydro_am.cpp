//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_am.cpp
//! \brief cs-open-1018: the angular-momentum module AmModule (hydro_am.hpp), shared by
//! <hydro> and <mhd>: am_transport, am_vector (+ spin reservoir, am_spin_return),
//! frame_omega (+ frame_exact).  BASES: on the cubed sphere the stored momentum u0 IM2/IM3
//! is COVARIANT on the cell's unit tangents (m_a = p . e_a, g_ab = [[1, c], [c, 1]],
//! c = cos_cell, det g = 1 - c^2), the primitive velocity w0 IVY/IVZ CONTRAVARIANT
//! (v = v^a e_a), the MHD cell field bcc the ORTHONORMAL triple (B.rhat, B.e_xi,
//! B.(e_eta - c e_xi)/s) (GnomonicEquiangleRaiseVelMHD); the face fluxes of IM2/IM3 are
//! covariant on the FACE's basis (raised with the face's cos).  Every torque r x p is formed
//! from the CARTESIAN vector; a Cartesian change goes back to the covariant momentum by
//! p . e_a.  Spherical polar: orthonormal (e_r, e_theta, e_phi) throughout.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "rad_m1/rad_m1.hpp"
#include "hydro/hydro_am.hpp"

namespace hydro {

bool AmModule::Named(ParameterInput *pin, const std::string &b) {
  return pin->DoesParameterExist(b, "am_transport") ||
         pin->DoesParameterExist(b, "am_vector") ||
         pin->DoesParameterExist(b, "frame_omega") ||
         pin->DoesParameterExist(b, "am_basis_check");
}

bool AmModule::DefaultOn(MeshBlockPack *ppack, ParameterInput *pin, const std::string &b) {
  if (pin->DoesParameterExist(b, "am_transport")) return false;
  Mesh *pm = ppack->pmesh;
  if (!(pm->use_cubed_sphere || pm->use_spherical_polar)) return false;   // silent
  const bool root = (global_variable::my_rank == 0);
  const std::string ev = pin->DoesParameterExist("time","evolution")
                         ? pin->GetString("time","evolution") : "";
  std::string why = "";
  if (ev != "dynamic") why = "<time>/evolution is not dynamic";
  if (pm->use_spherical_polar && pm->mb_indcs.nx2 < 2) why = "1-D spherical polar";
  if (pm->use_spherical_polar && pm->mb_indcs.nx2 >= 2 && !pm->three_d &&
      ppack->pcoord->sp_cart_polar_momentum) {
    why = "2-D spherical polar with sp_cart_polar_momentum";
  }
  if (!why.empty()) {
    if (root) {
      std::cout << "<" << b << ">/am_transport not named: default off (" << why << ")"
                << std::endl;
    }
    return false;
  }
  if (global_variable::restart_run) {
    if (root) {
      std::cout << "<" << b << ">/am_transport: restart input has no <" << b
                << ">/am_transport: keeping the old default (off); set it true to switch"
                << std::endl;
    }
    return false;
  }
  pin->GetOrAddBoolean(b, "am_transport", true);
  if (root) {
    std::cout << "<" << b << ">/am_transport = true (default on the "
              << (pm->use_cubed_sphere ? "cubed sphere" : "spherical-polar mesh")
              << "; set it false to switch the angular-momentum transport correction off)"
              << std::endl;
  }
  return true;
}

void AmModule::YieldToPgen(ParameterInput *pin, const std::string &who) {
  if (!am_transport) return;
  if (!am_default || am_vector) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<" << blk << ">/am_transport (or am_vector) = true together "
              << "with " << who << ": the transport correction would be applied twice; "
              << "set one of them false" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  am_transport = false;
  pin->SetBoolean(blk, "am_transport", false);
  if (global_variable::my_rank == 0) {
    std::cout << "<" << blk << ">/am_transport: default off (" << who << " is on)"
              << std::endl;
  }
}

AmModule *AmModuleOf(MeshBlockPack *pp) {
  if (pp->phydro != nullptr && pp->phydro->pam != nullptr) return pp->phydro->pam;
  if (pp->pmhd != nullptr && pp->pmhd->pam != nullptr) return pp->pmhd->pam;
  return nullptr;
}

AmModule::AmModule(MeshBlockPack *ppack, ParameterInput *pin, const std::string &b) :
    pmy_pack(ppack), blk(b) {
  if (pin->DoesParameterExist(blk,"am_transport")) {
    am_transport = pin->GetBoolean(blk,"am_transport");
  }
  if (pin->DoesParameterExist(blk,"frame_omega")) {
    std::string so = pin->GetString(blk,"frame_omega");
    for (auto &ch : so) {if (ch == ',') ch = ' ';}
    std::istringstream iso(so);
    iso >> frame_om[0] >> frame_om[1] >> frame_om[2];
    frame_on = true;
    if (pin->DoesParameterExist(blk,"frame_exact")) {
      frame_exact = pin->GetBoolean(blk,"frame_exact");
    }
  }
  if (pin->DoesParameterExist(blk,"am_vector")) {
    am_vector = pin->GetBoolean(blk,"am_vector");
    if (am_vector) {
      am_transport = true;
      auto &ic = ppack->pmesh->mb_indcs;
      const int nmbs = std::max(ppack->nmb_thispack, ppack->pmesh->nmb_maxperrank);
      const int n3 = ic.nx3 + 2*ic.ng, n2 = ic.nx2 + 2*ic.ng, n1 = ic.nx1 + 2*ic.ng;
      Kokkos::realloc(am_spin, nmbs, n3, n2, n1);
      Kokkos::realloc(am_spin1, nmbs, n3, n2, n1);
      Kokkos::realloc(am_spins, nmbs, n3, n2, n1);
      Kokkos::deep_copy(am_spin, 0.0);
      Kokkos::deep_copy(am_spin1, 0.0);
      Kokkos::deep_copy(am_spins, 0.0);
      if (pin->DoesParameterExist(blk,"am_spin_return")) {
        am_eta = pin->GetReal(blk,"am_spin_return");
      }
      if (am_eta > 0.0) {
        Kokkos::realloc(spw, nmbs, 1, n3, n2, n1);
        Kokkos::deep_copy(spw, 0.0);
        Kokkos::realloc(spw_c, nmbs, 1, 1, 1, 1);
        pbval_spin = new MeshBoundaryValuesCC(ppack, pin, false);
        pbval_spin->InitializeBuffers(1);
        pbval_spin->SetVectorPairs(1, {});
        pbval_spin->cs_noresample = true;
      }
    }
  }
  if (pin->DoesParameterExist(blk,"am_basis_check") &&
      pin->GetBoolean(blk,"am_basis_check") && ppack->pmesh->use_cubed_sphere) {
    BasisCheck();
  }
}

AmModule::~AmModule() {
  if (pbval_spin != nullptr) {delete pbval_spin;}
}

void AmModule::Apply(Driver *pdrive, int stage, Real bdt, const AmFluid &fl) {
  if (am_transport) {AmTransport(bdt, pdrive, stage, fl);}
  if (frame_on) {FrameRotation(bdt, pdrive, stage, fl);}
  if (am_vector && am_eta > 0.0 && stage == pdrive->nexp_stages) {AmSpinReturn(fl);}
}

//----------------------------------------------------------------------------------------
//! \fn void AmModule::AmTransport
//! \brief cs-open-1018 <hydro>/am_transport: add T_cons - T_code (hydro_am.hpp) of every
//! active cell to its zonal momentum.  The source variants it does not replicate are a
//! FATAL (as the deep_hot_jupiter_rt problem/am_transport).

void AmModule::AmTransport(Real bdt, Driver *pdrive, int stage, const AmFluid &fl) {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  if (!am_checked) {
    // amdef-1019: T_code re-evaluates every source variant (wellbalance_static[_reconst],
    // cs/sp_wellbalanced_src, sp_cart_polar/all_momentum; sp_face_avg only rewrites the
    // fluxes the update and T_cons both read); what is left is the mesh
    std::string bad = "";
    if (!(pm->use_cubed_sphere || pm->use_spherical_polar)) {
      bad = "a mesh that is neither cubed sphere nor spherical polar";
    }
    if (pm->use_spherical_polar && pm->mb_indcs.nx2 < 2) bad = "a 1-D spherical-polar mesh";
    if (pm->use_spherical_polar && !pm->three_d && pco->sp_cart_polar_momentum) {
      bad = "sp_cart_polar_momentum in 2-D";
    }
    if (!bad.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<" << blk << ">/am_transport does not support " << bad
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    am_checked = true;
  }
  HydroAmV v;
  AmFillV(v, fl);
  auto &ic = pm->mb_indcs;
  auto u0_ = fl.u0;
  if (am_vector) {
    // the spin reservoir follows the RK register logic of u0 (u1 = u0 at stage 1)
    if (pdrive->integrator == "rk4") {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro>/am_vector: integrator rk4 not supported"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (stage == 1) {Kokkos::deep_copy(DevExeSpace(), am_spin1, am_spin);}
    const Real g0 = pdrive->gam0[stage-1], g1 = pdrive->gam1[stage-1];
    auto sp_ = am_spin;
    auto sp1_ = am_spin1;
    auto sps_ = am_spins;
    par_for("hydro_am_vector", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
            ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dT[3], rc[3], b1[3], b2[3], cc;
      HydroAmVecTorque(v, m, k, j, i, dT, rc, b1, b2, cc);
      const Real r2 = rc[0]*rc[0] + rc[1]*rc[1] + rc[2]*rc[2];
      const Real ir = 1.0/sqrt(r2);
      const Real sr = (dT[0]*rc[0] + dT[1]*rc[1] + dT[2]*rc[2])*ir;   // rhat . dT
      // the tangential momentum whose torque is the perpendicular part of dT:
      // dp = (dT x rc)/(|rc|^2 V)
      Real dp[3];
      HydroAmCross(dT, rc, dp);
      const Real iv = 1.0/v.vol(m,k,j,i);
      for (int q = 0; q < 3; ++q) {dp[q] *= iv/r2;}
      const Real p1 = dp[0]*b1[0] + dp[1]*b1[1] + dp[2]*b1[2];
      const Real p2 = dp[0]*b2[0] + dp[1]*b2[1] + dp[2]*b2[2];
      // covariant on the cell basis (cs: m_a = p.e_a; sp: orthonormal)
      u0_(m,IM2,k,j,i) += p1*bdt;
      u0_(m,IM3,k,j,i) += p2*bdt;
      sps_(m,k,j,i) = sp_(m,k,j,i);   // the stage's own state (frame_exact)
      sp_(m,k,j,i) = g0*sp_(m,k,j,i) + g1*sp1_(m,k,j,i) + sr*bdt;
      (void) cc;
    });
    return;
  }
  par_for("hydro_am_transport", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
          ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtq, c2, c3;
    if (!HydroAmTorque(v, m, k, j, i, dtq, c2, c3)) return;
    u0_(m,IM2,k,j,i) += dtq*c2*bdt;
    u0_(m,IM3,k,j,i) += dtq*c3*bdt;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmModule::AmSpinReturn
//! \brief cs-open-1018 <hydro>/am_spin_return (see hydro.hpp).  Face f between cell c
//! and its lateral neighbour n (directions rhat_c, rhat_n, d = |rhat_c - rhat_n|,
//! s = (rhat_c - rhat_n)/d): c takes dL = lam s, n takes -lam s, lam = -eta (S_c + S_n)
//! d / (16 (1 - rhat_c.rhat_n)) (S the spin), so that each face lowers both spins by
//! ~eta (S_c + S_n)/16.  A cell carries dL exactly: rhat.dL to its spin, the rest as the
//! tangential momentum (dL x rc)/(|rc|^2 V).  Faces at a physical or periodic lateral
//! boundary (where the two sides do not share a geometry) take nothing.

void AmModule::AmSpinReturn(const AmFluid &fl) {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  auto &ic = pm->mb_indcs;
  const int is = ic.is, ie = ic.ie, js = ic.js, je = ic.je, ks = ic.ks, ke = ic.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  // the spins in the ghosts (mirror cells at a seam)
  auto sp_ = am_spin;
  auto sw_ = spw;
  par_for("hydro_am_spin_cp", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    sw_(m,0,k,j,i) = sp_(m,k,j,i);
  });
  while (pbval_spin->InitRecv(1) == TaskStatus::incomplete) {}
  while (pbval_spin->PackAndSendCC(spw, spw_c) == TaskStatus::incomplete) {}
  while (pbval_spin->RecvAndUnpackCC(spw, spw_c) == TaskStatus::incomplete) {}
  while (pbval_spin->ClearSend() == TaskStatus::incomplete) {}
  while (pbval_spin->ClearRecv() == TaskStatus::incomplete) {}
  const bool cs = pm->use_cubed_sphere;
  const Real eta = am_eta;
  auto u0_ = fl.u0;
  auto vol = pco->volume;
  auto x1v = pco->x1v;
  auto x2v = pco->x2v;
  auto x3v = pco->x3v;
  auto &size = pmy_pack->pmb->mb_size;
  auto &pan = pmy_pack->pmb->mb_panel;
  auto &bcs = pmy_pack->pmb->mb_bcs;
  const int nx2 = ic.nx2, nx3 = ic.nx3;
  par_for("hydro_am_spin_ret", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the direction of cell (k', j') of this block (extended panel / sp coordinates)
    auto dirn = [&](const int kq, const int jq, Real *q, Real *b1, Real *b2) {
      if (cs) {
        const Real xi = 0.25*M_PI*CellCenterX(jq-js, nx2, size.d_view(m).x2min,
                                              size.d_view(m).x2max);
        const Real et = 0.25*M_PI*CellCenterX(kq-ks, nx3, size.d_view(m).x3min,
                                              size.d_view(m).x3max);
        cubed_sphere::PanelToCart(pan.d_view(m), xi, et, q);
        cubed_sphere::PanelTangents(pan.d_view(m), xi, et, b1, b2);
      } else {
        const Real th = x2v(m,jq), ph = x3v(m,kq);
        q[0] = sin(th)*cos(ph); q[1] = sin(th)*sin(ph); q[2] = cos(th);
        b1[0] = cos(th)*cos(ph); b1[1] = cos(th)*sin(ph); b1[2] = -sin(th);
        b2[0] = -sin(ph); b2[1] = cos(ph); b2[2] = 0.0;
      }
    };
    Real rh[3], b1[3], b2[3];
    dirn(k, j, rh, b1, b2);
    const Real r = x1v(m,i);
    Real dL[3] = {0.0, 0.0, 0.0};
    const Real sc = sw_(m,0,k,j,i);
    for (int f = 0; f < 4; ++f) {
      const int dj = (f == 0) ? -1 : ((f == 1) ? 1 : 0);
      const int dk = (f == 2) ? -1 : ((f == 3) ? 1 : 0);
      const bool edge = (f == 0 && j == js) || (f == 1 && j == je) ||
                        (f == 2 && k == ks) || (f == 3 && k == ke);
      bool seam = false;
      if (edge) {
        const BoundaryFace bf = (f == 0) ? BoundaryFace::inner_x2
                              : (f == 1) ? BoundaryFace::outer_x2
                              : (f == 2) ? BoundaryFace::inner_x3 : BoundaryFace::outer_x3;
        const BoundaryFlag fl = bcs.d_view(m,bf);
        if (fl == BoundaryFlag::panel) {
          seam = true;
        } else if (fl != BoundaryFlag::block) {
          continue;   // physical / periodic: nothing crosses
        }
      }
      Real rn[3], t1[3], t2[3];
      if (seam) {
        // the mirror cell: rhat_c reflected across the seam plane (normal = gradient of
        // the panel coordinate at the face, which holds the origin)
        const int jf = (f == 1) ? (j + 1) : j, kf = (f == 3) ? (k + 1) : k;
        Real a[3], b[3], n0[3];
        cubed_sphere::PanelFrame(pan.d_view(m), a, b, n0);
        Real nf[3];
        if (f < 2) {
          const Real x = tan(0.25*M_PI*LeftEdgeX(jf-js, nx2, size.d_view(m).x2min,
                                                 size.d_view(m).x2max));
          for (int q = 0; q < 3; ++q) {nf[q] = a[q] - x*n0[q];}
        } else {
          const Real y = tan(0.25*M_PI*LeftEdgeX(kf-ks, nx3, size.d_view(m).x3min,
                                                 size.d_view(m).x3max));
          for (int q = 0; q < 3; ++q) {nf[q] = b[q] - y*n0[q];}
        }
        const Real nn = 1.0/sqrt(nf[0]*nf[0] + nf[1]*nf[1] + nf[2]*nf[2]);
        for (int q = 0; q < 3; ++q) {nf[q] *= nn;}
        const Real pr = rh[0]*nf[0] + rh[1]*nf[1] + rh[2]*nf[2];
        for (int q = 0; q < 3; ++q) {rn[q] = rh[q] - 2.0*pr*nf[q];}
      } else {
        dirn(k + dk, j + dj, rn, t1, t2);
      }
      Real sv[3];
      for (int q = 0; q < 3; ++q) {sv[q] = rh[q] - rn[q];}
      const Real d = sqrt(sv[0]*sv[0] + sv[1]*sv[1] + sv[2]*sv[2]);
      const Real omc = 1.0 - (rh[0]*rn[0] + rh[1]*rn[1] + rh[2]*rn[2]);
      if (!(d > 0.0) || !(omc > 0.0)) continue;
      const Real sn = sw_(m,0,k+dk,j+dj,i);
      const Real lam = -eta*(sc + sn)*d/(16.0*omc);
      for (int q = 0; q < 3; ++q) {dL[q] += lam*sv[q]/d;}
    }
    // this cell's share: radial part to the spin, the rest to the tangential momentum
    const Real srad = dL[0]*rh[0] + dL[1]*rh[1] + dL[2]*rh[2];
    Real rc[3] = {r*rh[0], r*rh[1], r*rh[2]};
    Real dp[3];
    HydroAmCross(dL, rc, dp);
    const Real iv = 1.0/(vol(m,k,j,i)*r*r);
    Real p1 = 0.0, p2 = 0.0;
    for (int q = 0; q < 3; ++q) {p1 += dp[q]*iv*b1[q]; p2 += dp[q]*iv*b2[q];}
    u0_(m,IM2,k,j,i) += p1;
    u0_(m,IM3,k,j,i) += p2;
    sp_(m,k,j,i) += srad;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void AmModule::FrameRotation
//! \brief cs-open-1018 <hydro>/frame_omega: the rotating-frame force (HydroAmFrame) on the
//! cell momentum and its work (the centrifugal part) on the total energy; under
//! <hydro>/am_vector also the L-conserving torque correction (perpendicular part in the
//! tangential momentum, radial part in the spin; the energy is not touched by it).

void AmModule::FrameRotation(Real bdt, Driver *pdrive, int stage, const AmFluid &fl) {
  Mesh *pm = pmy_pack->pmesh;
  if (!(pm->use_cubed_sphere || pm->use_spherical_polar)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/frame_omega needs the cubed sphere or spherical "
              << "polar" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int nst = pdrive->nexp_stages;
  const bool fex = frame_exact && am_vector;
  if (frame_exact && (!am_vector || nst > 2)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<hydro>/frame_exact needs <hydro>/am_vector and integrator "
              << "rk1 | rk2" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  HydroAmV v;
  AmFillV(v, fl);
  auto &ic = pm->mb_indcs;
  auto u0_ = fl.u0;
  auto u1_ = fl.u1;
  const Real o0 = frame_om[0], o1 = frame_om[1], o2 = frame_om[2];
  const bool amc = am_vector;
  auto sp_ = am_spin;
  auto sp1_ = am_spin1;
  auto sps_ = am_spins;
  // the radiation's L (F / (c chat)) of the cell is part of L_abs: the frame terms book
  // its rotation on the gas (the M1 equations carry no frame terms); u0 F covariant on
  // the cs panel basis (raised with cos_cell), orthonormal on sp
  // (only when the radiation exchanges momentum with this gas: coupling + gas_feedback)
  const bool hrad = amc && (pmy_pack->pradm1 != nullptr) &&
                    pmy_pack->pradm1->coupling && pmy_pack->pradm1->gas_feedback;
  auto ur_ = hrad ? pmy_pack->pradm1->u0 : DvceArray5D<Real>();
  const Real icc = hrad ? 1.0/(pmy_pack->pradm1->c_light*pmy_pack->pradm1->chat) : 0.0;
  if (fex && frame_ln.extent_int(0) == 0) {
    Kokkos::realloc(frame_ln, pmy_pack->nmb_thispack, 3, ic.nx3 + 2*ic.ng,
                    ic.nx2 + 2*ic.ng, ic.nx1 + 2*ic.ng);
  }
  auto ln_ = frame_ln;
  const bool rec = fex && (stage == 1);       // record L^n
  const bool fin = fex && (stage == nst);     // add the exact-rotation remainder
  const Real dt = pm->dt;
  // the exact rotation R(-Omega dt) (Rodrigues) and A = -Omega x
  const Real om = sqrt(o0*o0 + o1*o1 + o2*o2);
  const Real kx = (om > 0.0) ? o0/om : 0.0, ky = (om > 0.0) ? o1/om : 0.0;
  const Real kz = (om > 0.0) ? o2/om : 0.0;
  const Real cth = cos(-om*dt), sth = sin(-om*dt);
  const bool one = (nst == 1);
  par_for("hydro_frame_rot", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1,
          ic.ks, ic.ke, ic.js, ic.je, ic.is, ic.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real omv[3] = {o0, o1, o2};
    Real f[3], dT[3], rc[3], rh[3], b1[3], b2[3], cc;
    HydroAmFrame(v, omv, amc, m, k, j, i, f, dT, rc, rh, b1, b2, cc);
    const Real v1 = v.w0(m,IVX,k,j,i), v2 = v.w0(m,IVY,k,j,i), v3 = v.w0(m,IVZ,k,j,i);
    Real pw = 0.0;
    for (int q = 0; q < 3; ++q) {pw += f[q]*(v1*rh[q] + v2*b1[q] + v3*b2[q]);}
    Real fr = 0.0, f1 = 0.0, f2 = 0.0;
    for (int q = 0; q < 3; ++q) {fr += f[q]*rh[q]; f1 += f[q]*b1[q]; f2 += f[q]*b2[q];}
    u0_(m,IM1,k,j,i) += fr*bdt;
    u0_(m,IM2,k,j,i) += f1*bdt;
    u0_(m,IM3,k,j,i) += f2*bdt;
    u0_(m,IEN,k,j,i) += pw*bdt;
    if (amc) {
      const Real dv = v.vol(m,k,j,i);
      const Real s2 = 1.0 - cc*cc;
      // a covariant cell vector (radial, x2, x3) -> Cartesian
      auto cart = [&](const Real a1, const Real b2c, const Real b3c, Real *o) {
        const Real a2 = (b2c - cc*b3c)/s2, a3 = (b3c - cc*b2c)/s2;
        for (int q = 0; q < 3; ++q) {o[q] = a1*rh[q] + a2*b1[q] + a3*b2[q];}
      };
      // the spin and radiation parts of L_abs (current): -Omega x (S rhat + L_rad)
      Real lx[3];
      // (the stage's state: the spin before this stage's RK combination)
      for (int q = 0; q < 3; ++q) {lx[q] = sps_(m,k,j,i)*rh[q];}
      if (hrad) {
        Real fc[3], lr[3];
        cart(ur_(m,radm1::M1_F1,k,j,i), ur_(m,radm1::M1_F2,k,j,i),
             ur_(m,radm1::M1_F3,k,j,i), fc);
        for (int q = 0; q < 3; ++q) {fc[q] *= dv*icc;}
        HydroAmCross(rc, fc, lr);
        for (int q = 0; q < 3; ++q) {lx[q] += lr[q];}
      }
      Real ol[3];
      HydroAmCross(omv, lx, ol);
      for (int q = 0; q < 3; ++q) {dT[q] -= ol[q];}
      Real dl[3];
      for (int q = 0; q < 3; ++q) {dl[q] = dT[q]*bdt;}
      if (rec) {
        // L_abs^n of the cell: gas V rc x (m + rho Omega x rc) from the step-start
        // conserved state u1 (covariant m), spin^n (am_spin1), radiation (u0 before the
        // step's solves)
        Real mc[3], orc[3], pa[3], la[3];
        cart(u1_(m,IM1,k,j,i), u1_(m,IM2,k,j,i), u1_(m,IM3,k,j,i), mc);
        HydroAmCross(omv, rc, orc);
        const Real rho = u1_(m,IDN,k,j,i);
        for (int q = 0; q < 3; ++q) {pa[q] = (mc[q] + rho*orc[q])*dv;}
        HydroAmCross(rc, pa, la);
        for (int q = 0; q < 3; ++q) {
          la[q] += sp1_(m,k,j,i)*rh[q] + (lx[q] - sps_(m,k,j,i)*rh[q]);
          ln_(m,q,k,j,i) = la[q];
        }
      }
      if (fin) {
        Real l[3] = {ln_(m,0,k,j,i), ln_(m,1,k,j,i), ln_(m,2,k,j,i)};
        // R l (Rodrigues about k by -|Omega| dt)
        const Real kl = kx*l[0] + ky*l[1] + kz*l[2];
        const Real kxl[3] = {ky*l[2] - kz*l[1], kz*l[0] - kx*l[2], kx*l[1] - ky*l[0]};
        Real al[3], aal[3];
        {Real t[3]; HydroAmCross(omv, l, t);
         for (int q = 0; q < 3; ++q) {al[q] = -t[q];}}
        {Real t[3]; HydroAmCross(omv, al, t);
         for (int q = 0; q < 3; ++q) {aal[q] = -t[q];}}
        for (int q = 0; q < 3; ++q) {
          const Real kq = (q == 0) ? kx : ((q == 1) ? ky : kz);
          const Real rl = l[q]*cth + kxl[q]*sth + kq*kl*(1.0 - cth);
          dl[q] += rl - l[q] - dt*al[q] - (one ? 0.0 : 0.5*dt*dt*aal[q]);
        }
      }
      const Real r2 = rc[0]*rc[0] + rc[1]*rc[1] + rc[2]*rc[2];
      Real dp[3];
      HydroAmCross(dl, rc, dp);
      const Real iv = 1.0/(dv*r2);
      Real p1 = 0.0, p2 = 0.0;
      for (int q = 0; q < 3; ++q) {p1 += dp[q]*iv*b1[q]; p2 += dp[q]*iv*b2[q];}
      u0_(m,IM2,k,j,i) += p1;
      u0_(m,IM3,k,j,i) += p2;
      sp_(m,k,j,i) += dl[0]*rh[0] + dl[1]*rh[1] + dl[2]*rh[2];
    }
  });
}

//! \brief the views of the angular-momentum kernels (hydro_am.hpp)
void AmModule::AmFillV(HydroAmV &v, const AmFluid &fl) {
  Mesh *pm = pmy_pack->pmesh;
  Coordinates *pco = pmy_pack->pcoord;
  auto &ic = pm->mb_indcs;
  v.w0 = fl.w0;
  v.eos = fl.eos;
  v.gen = v.eos.IsGeneral();
  v.wder = fl.wder;
  v.f1 = fl.uflx.x1f;  v.f2 = fl.uflx.x2f;  v.f3 = fl.uflx.x3f;
  v.mhd = fl.mhd;
  v.bcc = fl.bcc;
  v.a1 = pco->area.x1f;  v.a2 = pco->area.x2f;  v.a3 = pco->area.x3f;
  v.vol = pco->volume;
  v.xovr = pco->x_ov_rD;
  v.yovr = pco->y_ov_rC;
  v.ccell = pco->cos_cell;
  v.scell = pco->sin_cell;
  v.x1v = pco->x1v;  v.x1f = pco->xx1f;
  v.x2v = pco->x2v;  v.x2f = pco->xx2f;
  v.x3v = pco->x3v;  v.x3f = pco->xx3f;
  v.size = pmy_pack->pmb->mb_size;
  v.panel = pmy_pack->pmb->mb_panel;
  v.is = ic.is;  v.js = ic.js;  v.ks = ic.ks;  v.nx2 = ic.nx2;  v.nx3 = ic.nx3;
  v.cs = pm->use_cubed_sphere;
  // amdef-1019: the source variants
  v.je = ic.je;
  v.three_d = pm->three_d;
  v.wbs = fl.wbs;
  v.pwb = fl.pwb;
  v.wbsrc = v.cs ? pco->cs_wellbalanced_src : pco->sp_wellbalanced_src;
  v.spcart = (!v.cs) && pco->sp_cart_polar_momentum;
  v.spall = v.spcart && pco->sp_cart_all_momentum;
  if (v.wbsrc || v.spcart) {
    // (the coordinate source of this stage built it; AMR rebuilds it there every call)
    if (static_cast<int>(pco->wb_geom.extent(0)) != pmy_pack->nmb_thispack) {
      pco->BuildWBGeometry();
    }
    v.geom = pco->wb_geom;
  }
  v.bcs = pmy_pack->pmb->mb_bcs;
}

//----------------------------------------------------------------------------------------
//! \fn void AmModule::BasisCheck
//! \brief cs-open-1018 <hydro|mhd>/am_basis_check = true (test only, read only when named;
//! rank 0 prints at start-up): the vector bases the angular-momentum code uses, checked
//! pointwise on every panel over a grid that covers the interior, the seams, the cube
//! vertices and two ghost cells beyond each edge (N = 32), with deterministic pseudo-random
//! Cartesian vectors V (a velocity / momentum) and B (a field):
//!  (1) covariant on the unit tangents m_a = V.e_a -> contravariant g^ab m_b (g = [[1,c],
//!      [c,1]]) -> Cartesian V_r rhat + v^a e_a;
//!  (2) face-normal f_a = V.nhat_a (nhat_xi = (e_xi - c e_eta)/s) -> contravariant f_a/s;
//!  (3) the MHD cell triple (B.rhat, B.e_xi, B.(e_eta - c e_xi)/s) -> Cartesian, |B|^2 = the
//!      sum of its squares;
//!  (4) invariants |V|^2 = V_r^2 + m_a v^a, V.B = V_r B_r + m_a B^a, |B|^2 = B_r^2 + g_ab
//!      B^a B^b;
//!  (5) across a seam (a point of panel p that lies on panel q, q = FindPanel): the code's
//!      cubed_sphere::TransformMomentum (covariant) and TransformFieldToDstNormals
//!      (face-normal) against the direct projection on panel q's basis at that point.
//! Every error is relative to |V| (|B|); round-off (~1e-15) is the pass.

void AmModule::BasisCheck() {
  const int nn = 32;
  const Real h = 0.5*M_PI/nn;
  Real e[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nseam = 0;
  unsigned int seed = 12345u;
  auto rnd = [&]() {
    seed = 1103515245u*seed + 12345u;
    return static_cast<Real>((seed >> 8) & 0xffff)/32768.0 - 1.0;
  };
  for (int p = 0; p < 6; ++p) {
    for (int jj = -2; jj < nn + 2; ++jj) {
      for (int kk = -2; kk < nn + 2; ++kk) {
        const Real xi = -0.25*M_PI + (jj + 0.5)*h, et = -0.25*M_PI + (kk + 0.5)*h;
        Real rh[3], e1[3], e2[3];
        cubed_sphere::PanelToCart(p, xi, et, rh);
        cubed_sphere::PanelTangents(p, xi, et, e1, e2);
        const Real c = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
        const Real s2 = 1.0 - c*c, s = std::sqrt(s2);
        Real V[3], B[3];
        for (int q = 0; q < 3; ++q) {V[q] = rnd(); B[q] = rnd();}
        auto dot = [](const Real *a, const Real *b) {return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];};
        const Real nv = std::sqrt(dot(V, V)), nb = std::sqrt(dot(B, B));
        const Real vr = dot(V, rh), br = dot(B, rh);
        // (1)
        const Real m2 = dot(V, e1), m3 = dot(V, e2);
        const Real v2 = (m2 - c*m3)/s2, v3 = (m3 - c*m2)/s2;
        Real err = 0.0;
        for (int q = 0; q < 3; ++q) {
          err = std::max(err, std::fabs(vr*rh[q] + v2*e1[q] + v3*e2[q] - V[q]));
        }
        e[0] = std::max(e[0], err/nv);
        // (2)
        Real n1[3], n2[3];
        for (int q = 0; q < 3; ++q) {n1[q] = (e1[q] - c*e2[q])/s; n2[q] = (e2[q] - c*e1[q])/s;}
        const Real f2 = dot(V, n1), f3 = dot(V, n2);
        err = 0.0;
        for (int q = 0; q < 3; ++q) {
          err = std::max(err, std::fabs(vr*rh[q] + (f2/s)*e1[q] + (f3/s)*e2[q] - V[q]));
        }
        e[1] = std::max(e[1], err/nv);
        // (3)
        Real o3[3];
        for (int q = 0; q < 3; ++q) {o3[q] = (e2[q] - c*e1[q])/s;}
        const Real t1 = dot(B, e1), t2 = dot(B, o3);
        err = 0.0;
        for (int q = 0; q < 3; ++q) {err = std::max(err, std::fabs(br*rh[q] + t1*e1[q] + t2*o3[q] - B[q]));}
        e[2] = std::max(e[2], err/nb);
        e[3] = std::max(e[3], std::fabs(br*br + t1*t1 + t2*t2 - dot(B, B))/(nb*nb));
        // (4)
        const Real b2c = (dot(B, e1) - c*dot(B, e2))/s2, b3c = (dot(B, e2) - c*dot(B, e1))/s2;
        e[4] = std::max(e[4], std::fabs(vr*vr + m2*v2 + m3*v3 - dot(V, V))/(nv*nv));
        e[5] = std::max(e[5], std::fabs(vr*br + m2*b2c + m3*b3c - dot(V, B))/(nv*nb));
        const Real bm2 = dot(B, e1), bm3 = dot(B, e2);
        e[6] = std::max(e[6], std::fabs(br*br + bm2*b2c + bm3*b3c - dot(B, B))/(nb*nb));
        // (5) a point that belongs to another panel (the ghost region beyond a seam)
        const int q = cubed_sphere::FindPanel(rh);
        if (q != p) {
          nseam += 1;
          Real d1[3], d2[3], xq, eq;
          cubed_sphere::CartToPanel(q, rh, xq, eq);
          cubed_sphere::PanelTangents(q, xq, eq, d1, d2);
          const Real cq = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2], sq = std::sqrt(1.0 - cq*cq);
          Real mo2, mo3, fo2, fo3;
          cubed_sphere::TransformMomentum(p, q, xi, et, m2, m3, mo2, mo3);
          cubed_sphere::TransformFieldToDstNormals(p, q, xi, et, f2, f3, fo2, fo3);
          Real nq1[3], nq2[3];
          for (int r = 0; r < 3; ++r) {nq1[r] = (d1[r] - cq*d2[r])/sq; nq2[r] = (d2[r] - cq*d1[r])/sq;}
          e[7] = std::max(e[7], std::max(std::fabs(mo2 - dot(V, d1)), std::fabs(mo3 - dot(V, d2)))/nv);
          e[7] = std::max(e[7], std::max(std::fabs(fo2 - dot(V, nq1)), std::fabs(fo3 - dot(V, nq2)))/nv);
        }
      }
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<" << blk << ">/am_basis_check (6 panels, N = " << nn << " + 2 ghosts, "
              << nseam << " cross-seam points): covariant round trip " << e[0]
              << ", face-normal round trip " << e[1] << ", MHD triple " << e[2]
              << ", |B|^2 triple " << e[3] << ", |V|^2 " << e[4] << ", V.B " << e[5]
              << ", |B|^2 metric " << e[6] << ", seam transforms " << e[7] << std::endl;
  }
}

}  // namespace hydro
