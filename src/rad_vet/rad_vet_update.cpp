//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_update.cpp
//! \brief the write-back of a stage: E, cell and face fluxes, gas, stage slopes, L_z.
//!
//!  called from: rad_vet_tasks.cpp
//!  calls into: rad_vet.cpp, rad_vet_bvals.cpp, rad_vet_fct.cpp, rad_vet_flux.cpp,
//!    rad_vet_flux_seam.cpp, rad_vet_init.cpp, rad_vet_ledger.cpp

#include <cmath>
#include <cstdio>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "utils/scratch_view.hpp"

namespace radvet {

//! \brief staggered: the radial component of a vector v in the components of
//! FaceFrameGasMomentum (spherical polar: v[0]; the cubed sphere: r_hat . v)

KOKKOS_INLINE_FUNCTION
static Real RadialPart(const SolveArgs &x, const int m, const int k, const int j,
                       const int i, const Real v[3]) {
  if (!x.cub) return v[0];
  return x.g(m,CLG_R,k,j,i)*v[0] + x.g(m,CLG_R+1,k,j,i)*v[1] + x.g(m,CLG_R+2,k,j,i)*v[2];
}

//----------------------------------------------------------------------------------------
//! \fn Real CsTorque
//! \brief cubed sphere, clean_cs_lz_tel: the gas L_z change of cell (k,j,i) over the
//! stage in telescoping form; rl: the radiation's (halves of V F_f lambda_f/c^2 of the
//! lateral faces, lambda_f = (M_f x n_f)_z/A_f); xs: the stress part (exact less faces').
//! The gas takes the faces' absorption torque V cf lambda plus the exact stress torque of
//! P = D E: fluxes r_f A (n x P n)_z (radial) and (M_f x P n)_z (lateral), P the two-cell
//! mean (one-sided at a radial boundary).  Interior faces cancel; LzBook books the
//! boundary torque.  Radial boundary faces take no stress
//! torque (lever arm along n), so gas and radiation L_z are exact with no ledger.

KOKKOS_INLINE_FUNCTION
static Real CsTorque(const SolveArgs &x, const int e, const int m, const int k,
                     const int j, const int i, Real &rl, Real &xs) {
  const Real v = x.vol(m,k,j,i), ic2 = 1.0/(x.c*x.c);
  Real gl = 0.0, gs = 0.0, tc = 0.0, p0[6];
  rl = 0.0;
  for (int q = 1; q < 3; ++q) {   // lateral faces: absorption and momenta
    Real cf[2], wf[2], tq[2];
    FaceMomentumChange(x, e, q, m, k, j, i, nullptr, cf, wf, tq);
    for (int s = 0; s < 2; ++s) {
      const int kf = k + s*(q == 2), jf = j + s*(q == 1);
      const Real f = x.fs ? x.fc[q](m,CLF_N,kf,jf,i)
                          : FaceFlux(x, e, q, m, kf, jf, i, true);
      const Real df = (f - x.fc[q](m,CLF_G,kf,jf,i))*ic2;
      const Real lam = CsArm(x, q, m, kf, jf, i);
      gl += v*lam*cf[s];
      gs += 0.5*v*lam*(tq[s] - df);
      rl += 0.5*v*lam*df;
    }
  }
  for (int n = 0; n < 6; ++n) p0[n] = x.cs(m,CLN_D11+n,k,j,i)*x.w(m,e,k,j,i);
  for (int q = 0; q < 3; ++q) {   // stress torque through the six faces
    const int dk = (q == 2), dj = (q == 1), di = (q == 0);
    for (int s = 0; s < 2; ++s) {
      const int o = 2*s - 1, kf = k + s*dk, jf = j + s*dj, iff = i + s*di;
      const int ko = k + o*dk, jo = j + o*dj, io = i + o*di;
      const bool bnd = IsBoundaryFace(x, q, m, kf, jf, iff);
      if (bnd && q == 0) continue;
      Real pf[6], a[3], t[2];
      for (int n = 0; n < 6; ++n) {
        pf[n] = bnd ? p0[n] : 0.5*(p0[n] + x.cs(m,CLN_D11+n,ko,jo,io)*x.w(m,e,ko,jo,io));
      }
      for (int c = 0; c < 3; ++c) {
        a[c] = (q == 0) ? x.x1f(m,iff)*x.area[0](m,kf,jf,iff)*x.g.At<CLG_N>(m,c,kf,jf,iff)
                        : x.g.At<CLG_M>(m,3*(q-1)+c,kf,jf,iff);
      }
      for (int b = 0; b < 2; ++b) {   // P n, outward: x and y components
        t[b] = 0.0;
        for (int c = 0; c < 3; ++c) {
          t[b] += o*pf[(b == c) ? b : 2 + b + c]*x.g.At<CLG_N>(m,3*q+c,kf,jf,iff);
        }
      }
      tc += a[0]*t[1] - a[1]*t[0];
    }
  }
  xs = gs - x.h*tc;
  return gl + xs;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::WriteBack
//! \brief stores the solution: face fluxes, E, cell fluxes, gas momentum and energy, and
//! the stage slopes K = (U - U_old)/h.  The gas takes exactly the momentum and energy the
//! radiation gives up (less the reference force with force_reference), so both totals
//! change only by face fluxes; with wb_arad (force_reference_work = split) E pays the
//! reference work h d arad_ref.v; esrc is the radiation's alone.  Cubed sphere: the gas
//! takes CellMomentumChange's vector in the hydro's components, the output cell flux is
//! the face momenta's vector (frame of y): energy and L_z conserved to round-off.

void RadiationVet::WriteBack() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool slopes = (t2_solve != VET_T2S_NONE);
  const bool feedback = cln_feedback;
  const Real h = dt_sub;
  GhostFill(cln_w, CLW_X, 1);
  CsLateralStress(CLW_X);
  SolveArgs x = MakeArgs(), xh = x;
  xh.hi = cln_split;   // clean_stag_split: the faces store the momentum state
  auto kk = (t2_solve == VET_T2S_STAGE1) ? t2k2 : t2k1;
  if (marena_k2 && t2_solve == VET_T2S_STAGE1) {   // (as its own allocation: 0 elsewhere)
    Kokkos::deep_copy(DevExeSpace(), t2k2, 0.0);
    Kokkos::deep_copy(DevExeSpace(), t2pf2, 0.0);
  }
  // h div Phi of the solution (CLW_Q) before the F solve below moves the memory d
  // (clean_stag_fkin): E takes the residual of the operator it was solved with; with
  // clean_stag_split also the cell flux of Phi (the faces keep F^H)
  auto u0_ = u0;
  // face tie: each shared face's CLF_PT the lower gid's (CsFaceTie), so energy
  // telescopes across MeshBlock faces and the periodic wraps
  CsFaceTie(1, CLF_PT);
  // clean_stag_xdonor: the cell flux (radiation momentum) of the central cross flux and
  // face sources, as the gas momentum's stress: Phi's donor form equals it only at the
  // exact fixed point of q
  SolveArgs xc = x;
  xc.xq = -1;
  const SolveArgs *px = DevArgs(x, 0), *pxh = DevArgs(xh, 1);
  const SolveArgs *pxc = DevArgs(xc, 2);
  // time2_face_fct: the energy flux Phi of every face, which moved E, with its slope
  const bool pfo = t2_pfst && slopes && t2pf1.extent_int(0) > nmb1;
  auto pf = (t2_solve == VET_T2S_STAGE1) ? t2pf2 : t2pf1;
  // xdonor: each face's central cross flux once, for both its cells, in cln_xf's scratch
  // (its users write it before they read it)
  auto phq = cln_phic;   // the cell flux of the faces' Phi (output vet_phi)
  DvceArray5D<Real> xcv;
  if (xh.hi && x.xq >= 0) {
    xcv = ScratchView<DvceArray5D<Real>>("cln_xf", nmb1 + 1, 3, cln_w.s.extent_int(2),
                                         cln_w.s.extent_int(3), cln_w.s.extent_int(4));
    for (int n = 0; n < x.nf; ++n) {
      SolveFor("cln_wbxc", nmb1, ks, ke+(n == 2), js, je+(n == 1), is, ie+(n == 0),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        DirCase(n, [&](auto d) {
          constexpr int nd = decltype(d)::value;
          xcv(m,nd,k,j,i) = FaceFlux<nd>(*pxc, CLW_X, nd, m, k, j, i, true);
        });
      });
    }
  }
  // clean_wb_fct: each face's Phi for WbFct, components 0..2 of cln_wbfp
  const bool wf = feedback && !x.vzero;
  const DvceArray5D<Real> wfp = !wf ? DvceArray5D<Real>() :
      ScratchView<DvceArray5D<Real>>("cln_wbfp", nmb1 + 1, 4, u0.extent_int(2),
                                     u0.extent_int(3), u0.extent_int(4));
  SolveFor("cln_wbdiv", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const SolveArgs &x = *px, &xh = *pxh;   // device memory: no copy
    // h div Phi with the faces' Phi of CountUnrealizable (CLF_PT), FaceFlux on
    // physical faces
    Real dv = 0.0, pc[6];   // cubed sphere: each face's (CellVector)
    for (int n = 0; n < 3; ++n) {
      const int dk = (n == 2), dj = (n == 1), di = (n == 0);
      if (n >= x.nf) {
        if (xh.hi) phq(m,n,k,j,i) = 0.0;
        continue;
      }
      const Real fl = IsBoundaryFace(x, n, m, k, j, i) ?
          FaceFlux(x, CLW_X, n, m, k, j, i, true) : x.fc[n](m,CLF_PT,k,j,i);
      const bool bh = IsBoundaryFace(x, n, m, k+dk, j+dj, i+di);
      const Real fu = bh ? FaceFlux(x, CLW_X, n, m, k+dk, j+dj, i+di, true) :
                           x.fc[n](m,CLF_PT,k+dk,j+dj,i+di);
      dv += DivF(x, n, m, k, j, i, fl, fu);
      if (wf) {
        wfp(m,n,k,j,i) = fl;
        if (n == 0 ? i == x.ie : (n == 1 ? j == x.je : k == x.ke)) {
          wfp(m,n,k+dk,j+dj,i+di) = fu;
        }
      }
      if (xh.hi && x.xq >= 0) {
        phq(m,n,k,j,i) = 0.5*(xcv(m,n,k,j,i) + xcv(m,n,k+dk,j+dj,i+di));
      } else if (xh.hi) {
        phq(m,n,k,j,i) = 0.5*(fl + fu);
      }
      if (xh.hi && x.cub) {
        pc[2*n] = (x.xq >= 0) ? xcv(m,n,k,j,i) : fl;
        pc[2*n+1] = (x.xq >= 0) ? xcv(m,n,k+dk,j+dj,i+di) : fu;
      }
      if (pfo) pf(m,n,k,j,i) = fl;
      // the MeshBlock's high face: the Phi that moved E, one flux on both sides for
      // time2_face_fct
      if (pfo && (n == 0 ? i == x.ie : (n == 1 ? j == x.je : k == x.ke))) {
        pf(m,n,k+dk,j+dj,i+di) = fu;
      }
    }
    if (xh.hi && x.cub) {   // cubed sphere: the faces' Phi's vector, frame of y
      Real fv[3];
      CellVector(x, m, k, j, i, pc, fv);
      ToHydroFrame(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) phq(m,d,k,j,i) = fv[d];
    }
    x.w(m,CLW_Q,k,j,i) = x.h*dv;
  });
  if (pfo) {
    if (t2_solve == VET_T2S_STAGE1) {t2pf2_ok = true;} else {t2pf1_ok = true;}
  }
  // the history's boundary energy: + h A Phi of the physical x1 faces (BoundaryEnergyAcc)
  if (cln_lacc.extent_int(0) > nmb1) {
    auto la = cln_lacc;
    auto ar = pmy_pack->pcoord->area.x1f;   // x1 face area of the clean fluxes
    const bool stg = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
    par_for("cln_wbbnd", DevExeSpace(), 0, nmb1, 0, 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int s, const int k, const int j) {
      const int i = s ? ie + 1 : is;
      la(m,s,k,j) = !IsBoundaryFace(x, 0, m, k, j, i) ? 0.0 :
                    la(m,stg ? 2+s : s,k,j) + h*ar(m,k,j,i)*FaceFlux(x, CLW_X, 0, m,
                                                                     k, j, i, true);
    });
  }
  // face F of the F solve and its dissipation (lateral sums of X: above)
  FluxSolve(CLW_X, true, false);
  xh.fs = true;
  xh.fa = true;   // the tied T in CLF_PT (FluxTie)

  // clean_cs_lz_tel: CsTorque's stress torque xs is zero-sum noise per cell; the gas
  // takes only its net X = sum xs, as the solid-body rotation dOmega = X/I of least
  // kinetic energy, I = sum rho V R^2.  tqv: the gas's L_z gl and the ledger's; lxs, lxi:
  // xs and rho V R^2, lzs their sums X and I
  const bool lzt = cln_lzt && feedback;
  DvceArray5D<Real> tqv;
  Real lzs[2] = {0.0, 0.0};
  // (clean_lz_spin: X, the L_z of fcone's cut of the faces with an arm (FluxCone)
  // and of clean_wb_mcap's cap (a pass of the arm families' terms), the gas's as the
  // rigid rotation dOmega = X/I, I = sum rho V R^2; the cubed sphere: none, its gas phi
  // momentum is lz_tel's torque of the uncut, uncapped face terms)
  const bool spn = cln_lzsp && feedback && !x.vzero && !x.cub;
  Real lzx = 0.0;
  if (spn) {   // (per cell, then summed in a fixed order: each cell's own halves)
    SolveArgs xm = xh;
    xm.mcl = cln_wbmc;
    xm.mclc = DvceArray1D<Real>("cln_lzmc", 1);
    auto fcx = cln_fcx;
    auto lxx = ScratchView<DvceArray4D<Real>>("cln_lxx", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    CLEAN_HIPX(const SolveArgs *pxm = DevArgs(xm, 0), *pxh = DevArgs(xh, 1);)
    SolveFor("cln_lzmc", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const SolveArgs &xm = *pxm, &xh = *pxh;)
      Real s = 0.0;
      for (int d = 0; d < xh.nf; ++d) {
        if (!LzFamily(xh.sph, xh.cub, d)) continue;
        Real cf[2], wf[2], a = 0.0;
        const Real u = FaceMomentumChange(xh, CLW_X, d, m, k, j, i, nullptr, cf, wf);
        for (int sd = 0; sd < 2; ++sd) {   // (fcone's cut, as the gas's share of it)
          a += wf[sd]*fcx(m,d,k+sd*(d == 2),j+sd*(d == 1),i+sd*(d == 0));
        }
        if (xm.mcl) a += u - FaceMomentumChange(xm, CLW_X, d, m, k, j, i);
        s += a*DirLzArm(xh, d, m, k, j, i);
      }
      lxx(m,k,j,i) = s*CellVolume(xh, m, k, j, i);
    });
    lzx = SumCells(lxx);
    lzs[0] = lzx;
    if (!lzt) {
      auto ug = FluidRef::Get(pmy_pack).u0;
      auto lxi = ScratchView<DvceArray4D<Real>>("cln_lxi", nmb1 + 1, u0.extent_int(2),
                                                u0.extent_int(3), u0.extent_int(4));
      SolveFor("cln_lzsi", nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        lxi(m,k,j,i) = ug(m,IDN,k,j,i)*CellVolume(x, m, k, j, i)*
                       SQR(LzArm(x, m, k, j, i));
      });
      lzs[1] = SumCells(lxi);
    }
  }
  if (lzt) {
    tqv = ScratchView<DvceArray5D<Real>>("cln_tqv", nmb1 + 1, 3, u0.extent_int(2),
                                         u0.extent_int(3), u0.extent_int(4));
    auto lxs = ScratchView<DvceArray4D<Real>>("cln_lxs", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    auto lxi = ScratchView<DvceArray4D<Real>>("cln_lxi", nmb1 + 1, u0.extent_int(2),
                                              u0.extent_int(3), u0.extent_int(4));
    px = DevArgs(x, 0);   // xh changed: the F solve's F
    pxh = DevArgs(xh, 1);
    auto ug = FluidRef::Get(pmy_pack).u0;
    SolveFor("cln_wbtq", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const SolveArgs &x = *px, &xh = *pxh;)
      Real rl, xs, q[3][3];
      const Real gl = CsTorque(xh, CLW_X, m, k, j, i, rl, xs);
      const Real rc = CsFrame(x, m, k, j, i, q);
      lxs(m,k,j,i) = xs;
      lxi(m,k,j,i) = ug(m,IDN,k,j,i)*x.vol(m,k,j,i)*rc*rc;
      tqv(m,1,k,j,i) = gl - xs;   // the absorption torque
      tqv(m,2,k,j,i) = -(gl + rl);
    });
    lzs[0] += SumCells(lxs);
    lzs[1] = SumCells(lxi);
  }
  const Real lzw_ = (lzs[1] > 0.0) ? lzs[0]/lzs[1] : 0.0;   // (dOmega)
  if (spn && global_variable::my_rank == 0 && lzx != 0.0) {
    std::printf("CLEAN LZSPIN cycle %d t2s %d: L_z of the cut and capped momentum X "
                "%.6e, the gas's as dOmega %.6e (I %.6e)\n", pmy_pack->pmesh->ncycle,
                t2_solve, lzx, lzw_, lzs[1]);
  }

  // face fluxes of the solution
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int tf = VET_T2_F1 + dir;
    SolveFor("cln_wbface", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real f = fc(m,CLF_N,k,j,i);
      if (slopes) kk(m,tf,k,j,i) = (f - fc(m,CLF_G,k,j,i))/h;
      f0(m,k,j,i) = f;
    });
  }

  // E, cell fluxes, gas
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto f1 = f0x1, f2 = f0x2, f3 = f0x3;
  const int fref = force_ref;
  const bool wsp = (fref == VET_FREF_WB_ARAD || fref == VET_FREF_WB_MD);
  const bool fcap = (fref == VET_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  const bool wpa = wsp && FrefWaccOn();   // fref_wacc: RightHandSide
  auto wc = wpa ? FrefWacc() : DvceArray4D<Real>();
  const Real fra = FrefWaccStage() ? 1.0 : h/pmy_pack->pmesh->dt;
  const bool wex = wsp &&
                   time_scheme == VET_TIME_HESDIRK2;   // no slope carries a payment
  const bool fm = icr_fmeas && cln_vzero;   // SteadyStart: icr_dm
  auto fmd = icr_dm;
  const bool eso = esrc_on;
  auto es = esrc;
  const Real arad_ = arad;
  // wb_arad_md, clean_init = steady: the steady state's reference; cln_fref 0..2 its
  // acceleration (components of dm), 3 + 2 dir + side the potential the hydro walks of
  // <hydro>/wb_dirs = all take along dir: the line integral of 2 cf/(h rho) over each
  // half cell (side 0 face to centre, 1 centre to face).  None outside the Roche lobe
  // (problem/wb_lobe_cut)
  const bool mdc = cln_mdcap;
  const Real phx = fl.ph ? fl.ph->wb_phimax : (fl.pm ? fl.pm->wb_phimax : 1.0e300);
  auto phc = fl.phicc_wb;
  auto cx1 = pmy_pack->pcoord->x1v, cx2 = pmy_pack->pcoord->x2v;
  auto cx3 = pmy_pack->pcoord->x3v, cf1 = pmy_pack->pcoord->xx1f;
  auto cf2 = pmy_pack->pcoord->xx2f, cf3 = pmy_pack->pcoord->xx3f;
  const Real tgf = eos.tfloor;
  const bool lex = cln_wbxl && feedback;
  if (lex && !(tgf > 0.0)) RadVetFatal("<rad_vet>/clean_wb_exch_local needs a tfloor > "
                                       "0");
  const Real pgf = 1.0e-3, efl = e_floor;
  // clean_wb_floor_diag N: the first N floored cells of each write-back, their terms
  const int fdn = cln_wbfdg;
  DvceArray1D<Real> fdc;
  if (fdn > 0) fdc = DvceArray1D<Real>("cln_wbfdc", 1);
  const int fcy = pmy_pack->pmesh->ncycle, fts = t2_solve;
  auto pc_ = pos_cnt_d;
  const bool pflr = impl_pos_floor, pfc = pflr;
  // the E floor's raise neither gas nor column covered (created, VET_POS_FLR_UN): its
  // change in this write-back, CLEAN WARNING below
  Real un0 = 0.0;
  if (pflr) Kokkos::deep_copy(un0, Kokkos::subview(pc_, VET_POS_FLR_UN));
  if (pfc && pfc_def.extent_int(0) != nmb1 + 1) {
    auto &icp = pmy_pack->pmesh->mb_indcs;
    Kokkos::realloc(pfc_def, nmb1 + 1, icp.nx3 + 2*icp.ng, icp.nx2 + 2*icp.ng,
                    icp.nx1 + 2*icp.ng);
  }
  if (pfc) Kokkos::deep_copy(pfc_def, 0.0);
  auto pfd = pfc_def;
  const bool pgc = pfc;
  if (pgc && pgc_def.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(pgc_def, nmb1 + 1, pfc_def.extent_int(1), pfc_def.extent_int(2),
                    pfc_def.extent_int(3));
  }
  if (pgc) Kokkos::deep_copy(pgc_def, 0.0);
  auto pgd = pgc_def;
  const bool xwb = feedback && !cln_vzero;
  const bool cls = xwb;   // the closure in CLW_PW, printed below
  // (clean_wb_mcap: the gas's momentum at the passes' cap, that of the E row's work; the
  // terms capped and their momentum withheld per face family, CLEAN MCAP below)
  const bool wmc = cln_wbmc && feedback && !x.vzero;
  if (wmc) {
    xh.mcl = true;
    xh.mclc = DvceArray1D<Real>("cln_wbmcc", 4);
  }
  px = DevArgs(x, 0);   // the F solve took the slots; xh changed
  pxh = DevArgs(xh, 1);
  const bool lzb = cln_lzbook && feedback;
  DvceArray4D<Real> lzw;
  if (lzb) lzw = ScratchView<DvceArray4D<Real>>("cln_lzw", nmb1 + 1, u0.extent_int(2),
                                                 u0.extent_int(3), u0.extent_int(4));
  // (prb, clean_wb_fct's probe: no side effects, each cell's E without its face
  // transport, e + qf, into cln_wbfp 3 for WbFct, before the cells' write-back)
  // (clean_dbg_inf: once, in the first stage-1 write-back of that cycle)
  const bool dinf = (cln_dbginf == fcy) && t2_solve == VET_T2S_STAGE1 && !t2_inprobe &&
                    global_variable::my_rank == 0;
  if (dinf) cln_dbginf = -1;
  auto wbc = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                           const bool prb) CLEAN_HIP_INLINE {
    CLEAN_HIPX(const SolveArgs &x = *px, &xh = *pxh;)   // in device memory
    Real e = x.w(m,CLW_X,k,j,i);
    if (mdc) {
      Real dr[3], mg[3], cf[2], wf[2];
      CellMomentumChange(xh, CLW_X, m, k, j, i, dr, mg);
      const Real d = uh(m,IDN,k,j,i), s = (phc(m,k,j,i) > phx) ? 0.0 : 1.0/(h*d);
      for (int n = 0; n < 3; ++n) {
        fr(m,n,k,j,i) = s*dr[n];
        fr(m,3+2*n,k,j,i) = fr(m,4+2*n,k,j,i) = 0.0;
        if (n >= x.nf) continue;
        FaceMomentumChange(xh, CLW_X, n, m, k, j, i, nullptr, cf, wf);
        const int dk = (n == 2), dj = (n == 1), di = (n == 0);
        for (int sd = 0; sd < 2; ++sd) {
          const int kf = k + sd*dk, jf = j + sd*dj, iff = i + sd*di;
          const Real l = !(x.sph || x.cub) ? 0.5/x.idx[n] : ((n == 0) ?
              fabs(cf1(m,iff) - cx1(m,i)) : ((n == 1) ? cx1(m,i)*fabs(cf2(m,jf) -
              cx2(m,j)) : cx1(m,i)*sin(cx2(m,j))*fabs(cf3(m,kf) - cx3(m,k))));
          fr(m,3+2*n+sd,k,j,i) = 2.0*s*cf[sd]*l;
        }
      }
    }
    Real dm[3] = {0.0, 0.0, 0.0}, dh[3] = {0.0, 0.0, 0.0};
    Real den = 0.0, wsl = 0.0;
    Real dg[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // clean_wb_floor_diag: eg wk wp dw rd
    if (feedback) {
      Real mg[3], pg[3], rf[3], dlz = 0.0;
      CellMomentumChange(xh, CLW_X, m, k, j, i, dm, mg, &dlz);
      if (spn) {   // (clean_lz_spin: dOmega rho z x r; the ledger's withheld L_z back)
        const Real w = lzw_*uh(m,IDN,k,j,i), r = LzArm(x, m, k, j, i);
        if (x.sph) {
          dm[2] += w*r;
        } else {
          dm[0] -= w*x.x2v(m,j);
          dm[1] += w*x.x1v(m,i);
        }
        dlz -= w*r*r*CellVolume(x, m, k, j, i);
      }
      if (lzt) {   // (clean_cs_lz_tel: the gas phi momentum of CsTorque, cln_wbtq)
        Real q[3][3];
        const Real rc = CsFrame(x, m, k, j, i, q);
        const Real dp = tqv(m,1,k,j,i)/(x.vol(m,k,j,i)*rc) + lzw_*uh(m,IDN,k,j,i)*rc -
                        (q[2][0]*dm[0] + q[2][1]*dm[1] + q[2][2]*dm[2]);
        // force_reference = clean: the reference withholds this phi row too (cln_fref 3),
        // its torque booked as that of cln_fref 0..2
        const Real dw = (fref != VET_FREF_CLEAN) ? 0.0 : (fcap ? dp : h*fr(m,3,k,j,i));
        if (fcap) fr(m,3,k,j,i) = dp/h;
        for (int n = 0; n < 3; ++n) dm[n] += (dp - dw)*q[2][n];
        dlz = tqv(m,2,k,j,i) + x.vol(m,k,j,i)*rc*dw;
      }
      FaceFrameGasMomentum(x, uh, m, k, j, i, pg);
      RefForce(x, fref, fr, aref, uh(m,IDN,k,j,i), m, k, j, i, rf);
      for (int n = 0; n < 3; ++n) {   // force_reference
        if (fref != VET_FREF_NONE) dm[n] -= h*rf[n];
      }
      if (lzb && fref == VET_FREF_CLEAN) {   // the reference torque withheld: booked
        Real q[3][3];
        const Real r = x.cub ? CsFrame(x, m, k, j, i, q) : x.x1v(m,i)*x.tg(m,2,j);
        dlz += h*x.vol(m,k,j,i)*r*(x.cub ? q[2][0]*rf[0] + q[2][1]*rf[1] + q[2][2]*rf[2]
                                         : rf[2]);
      }
      if (fm) fmd(m,k,j,i) = RadialPart(x, m, k, j, i, dm);
      if (x.vzero) dm[0] = dm[1] = dm[2] = dlz = 0.0;   // frozen gas momentum
      if (lzb) lzw(m,k,j,i) = dlz;   // clean_cs_lz_abs = 2: LzBook
      // the kinetic energy the gas gains (the E row's W at convergence); with wb_arad wr,
      // the reference work at the velocity after the kick, which E pays
      Real wk = 0.0, wr = 0.0;
      for (int n = 0; n < 3; ++n) {
        wk += (pg[n] + 0.5*dm[n])*dm[n]/uh(m,IDN,k,j,i);
        if (wsp && !x.vzero) wr += h*rf[n]*(pg[n] + dm[n])/uh(m,IDN,k,j,i);
        dh[n] = dm[n];
      }
      if (x.cub) {   // the hydro's components of dm, from the frame of y
        Real dy[3] = {dm[0], dm[1], dm[2]};
        ToHydroFrame(x, m, k, j, i, dy);
        GasChange(x, m, k, j, i, dy, dh);
      }
      for (int n = 0; n < 3 && !prb; ++n) uh(m,IM1+n,k,j,i) += dh[n];
      // gas energy: exactly what the radiation gives up besides its flux divergence, so
      // total energy changes only by the face fluxes
      Real div = x.idx[0]*(f1(m,k,j,i+1) - f1(m,k,j,i));
      if (multi_d) div += x.idx[1]*(f2(m,k,j+1,i) - f2(m,k,j,i));
      if (three_d) div += x.idx[2]*(f3(m,k+1,j,i) - f3(m,k,j,i));
      // clean_stag_split: the faces store F^H, E moved by the energy flux; esrc is the
      // radiation's external source (the gas takes the exchange only)
      den = -(e - x.w(m,CLW_EO,k,j,i)) - ((xh.hi || x.sph || x.cub) ?
            x.w(m,CLW_Q,k,j,i) : h*div) + (eso ? h*es(m,k,j,i) : 0.0);
      // steady start: the eint of the iterate's T, not the residual of a step h -> inf
      if (x.vzero) {
        Real ei, p, chi_rho, chi_t, cv;
        eos.ThermoAt(uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), ei, p, chi_rho, chi_t, cv);
        den = ei - x.w(m,CLW_EI,k,j,i);
      } else {
        // the gas takes the exchange of the solution at the iterate's T, the beta terms
        // and W; E the rest (its fluxes): the linear solve's residual goes to E, not eint
        const Real t = x.w(m,CLW_T,k,j,i);
        const Real ex = h*x.c*(x.op(m,VET_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_*
                               SQR(SQR(t)) -
                               x.op(m,VET_OP_E,k,j,i)*fmax(e, 0.0));
        // the residual over the row's local diagonal (all of it in a scattering cell,
        // little where the exchange is stiff).  wpa: E pays fref_wacc; the slope keeps no
        // payment (wex: FrefWaccStage) or the payment wr, re-applied by the next stages
        const Real wp = wpa ? fra*wc(m,k,j,i) : wr;
        // clean_gas_tfloor: the gas eint changes by that of its T root (>= tfloor),
        // e(rho, T) - eint, not by the exchange (two terms ~E >> eint where r ~ 1):
        // positive by construction; E takes all the rest.  The floor only here (inside
        // the Newton passes it caused 2-cycles); a root below it raises the gas to
        // e(rho, tfloor) with E's energy, counted as implicit_pos_gas
        Real eg = x.w(m,CLW_EP,k,j,i) - ex;
        // clean_wb_exch_local: the cell's exchange solved exactly from the transported
        // E* = X + den - wk - wp at the solve's opacities, G and s = EP: the root
        // T >= tfloor of (1 + b)(e(rho, T) - eint) + c4 T^4 - b E* - s (b = h c kE,
        // c4 = h c kP G a), increasing in T, so unique; E = (E* - s + c4 T^4)/(1 + b) > 0
        // for E* - s > 0 and the gas takes E* - E; the face FCT's probe keeps the
        // iterate's T
        const Real es = x.w(m,CLW_X,k,j,i) + den - wk - wp;
        const Real sp = x.w(m,CLW_EP,k,j,i);
        const Real lb = x.h*x.c*x.op(m,VET_OP_E,k,j,i);
        if (lex && !prb) {
          const Real d = uh(m,IDN,k,j,i), e0 = x.w(m,CLW_EI,k,j,i);
          const Real c4 = x.h*x.c*x.op(m,VET_OP_P,k,j,i)*x.w(m,CLW_GM,k,j,i)*arad_;
          auto fx = [&](const Real tt, Real &df) {
            Real ee, p, chi_rho, chi_t, cv;
            eos.ThermoAt(d, tt, ee, p, chi_rho, chi_t, cv);
            df = (1.0 + lb)*d*cv + 4.0*c4*tt*tt*tt;
            return (1.0 + lb)*(ee - e0) + c4*SQR(SQR(tt)) - lb*es - sp;
          };
          Real df, tl = tgf, th = fmax(t, tgf), tr = tgf;
          const bool flo = !(fx(tl, df) < 0.0);
          if (!flo) {
            for (int n = 0; n < 2100 && fx(th, df) < 0.0; ++n) {tl = th; th *= 2.0;}
            tr = th;
            Real fr = fx(tr, df);
            for (int n = 0; n < 200 && fr != 0.0 && th - tl > 4.0e-16*th; ++n) {
              if (fr > 0.0) {th = tr;} else {tl = tr;}
              const Real tn = tr - fr/df;
              tr = (tn > tl && tn < th) ? tn : 0.5*(tl + th);
              fr = fx(tr, df);
            }
          }
          // a root at tfloor: the exchange at T = tfloor; the gas's rest to e(rho,
          // tfloor) is implicit_pos_gas's raise, below
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(d, tr, ei, p, chi_rho, chi_t, cv);
          eg = flo ? (lb*es + sp - c4*SQR(SQR(tgf)))/(1.0 + lb) : ei - e0;
        } else {
          Real ei, p, chi_rho, chi_t, cv;
          eos.ThermoAt(uh(m,IDN,k,j,i), fmax(t, tgf), ei, p, chi_rho, chi_t, cv);
          eg = ei - x.w(m,CLW_EI,k,j,i);
          if (t < tgf && !prb) {
            Real e0;
            eos.ThermoAt(uh(m,IDN,k,j,i), t, e0, p, chi_rho, chi_t, cv);
            const Real vol = CellVolume(x, m, k, j, i);
            Kokkos::atomic_add(&pc_(VET_POS_GAS), 1.0);
            Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), (ei - e0)*vol);
          }
        }
        // the work the row took (and q X) but the gas's kinetic energy did not receive
        // goes to E whole, eint the heating only
        const Real dw = xwb ? x.w(m,CLW_XW,k,j,i) - wk - wp +
                              x.w(m,CLW_PD,k,j,i)*x.w(m,CLW_X,k,j,i) : 0.0;
        const Real rd = den - eg - wk - wp - dw + dw;
        if (cls) x.w(m,CLW_PW,k,j,i) = fabs(den - eg - wk - wp)/fmax(fabs(e) +
            fabs(x.w(m,CLW_EO,k,j,i)) + fabs(x.w(m,CLW_EI,k,j,i)), 1.0e-300);
        e += rd;
        if (fdn > 0) {dg[0] = eg; dg[1] = wk; dg[2] = wp; dg[3] = dw; dg[4] = rd;}
        if (prb) {   // less E's round-off and e_floor: a limited cell stays above it
          const Real q = x.w(m,CLW_Q,k,j,i);
          wfp(m,3,k,j,i) = e + q - efl -
              1.0e-12*(fabs(e) + fabs(q) + fabs(eg) + fabs(wk) +
              fabs(wp) + fabs(x.w(m,CLW_X,k,j,i)) + fabs(x.w(m,CLW_EO,k,j,i)));
          return;
        }
        den -= rd + wp;
        wsl = wex ? wp : wp - wr;
        // a gas eint below e_min = max(e(rho, tfloor), 1e-3 x its old eint) is raised to
        // it with energy from the cell's E above e_floor (gas + E kept; counted)
        Real ef = 0.0, pr, chi_rho, chi_t, cv;
        if (eos.tfloor > 0.0) {
          eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
        }
        // from the gas's true eint, the stage start's deficit below 0 included
        const Real ei = x.w(m,CLW_EI,k,j,i), emin = fmax(ef, pgf*fmax(ei, 0.0));
        const Real need = emin - ei - (den - wk);
        const Real take = fmin(need, e - efl);
        const Real vol = CellVolume(x, m, k, j, i);
        const bool tk = take > 1.0e-12*emin;   // beyond eint's round-off at e_min
        if (tk) {
          den += take;
          e -= take;
          Kokkos::atomic_add(&pc_(VET_POS_GAS), 1.0);
          Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), take*vol);
        }
        // implicit_pos_gas_cons: what E could not cover, from the column, below
        const Real un = need - (tk ? take : 0.0);
        if (pgc && un > 1.0e-12*emin) pgd(m,k,j,i) = un*vol;
      }
      uh(m,IEN,k,j,i) += den;
      x.w(m,CLW_JV,k,j,i) = den - wk;   // its change of the gas eint (cln_neint)
    }
    // (implicit_pos_floor: an E at or below e_floor is
    // raised to it with the energy of the cell's gas above e_min; what the gas cannot
    // cover is taken from the radiation of the cell's column (implicit_pos_floor_cons,
    // below), else counted as created; counted).  A non-finite E (a diverged iterate) is
    // never raised: it stays non-finite, and the check below rejects the solve
    if (dinf && !prb && m == 0 && k == ks && j == js && i == is) e = -INFINITY;
    if (pflr && e <= efl && Kokkos::isfinite(e)) {
      const Real vol = CellVolume(x, m, k, j, i);
      const bool gas = feedback && !x.vzero;
      Real take = 0.0;
      if (gas) {
        Real ef = 0.0, pr, chi_rho, chi_t, cv;
        if (eos.tfloor > 0.0) {
          eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
        }
        const Real ei = x.w(m,CLW_EI,k,j,i), en = ei + x.w(m,CLW_JV,k,j,i);
        // (the gas keeps its eint: the raise is the column's)
        const Real emin = fmax(fmax(ef, pgf*fmax(ei, 0.0)), fmax(en, 0.0));
        take = fmin(efl - e, fmax(en - emin, 0.0));
        den -= take;
        uh(m,IEN,k,j,i) -= take;
        x.w(m,CLW_JV,k,j,i) -= take;
        Kokkos::atomic_add(&pc_(VET_POS_FLR_DE), take*vol);
      }
      if (pfc) {
        pfd(m,k,j,i) = (efl - e - take)*vol;
      } else {
        Kokkos::atomic_add(&pc_(VET_POS_FLR_UN), (efl - e - take)*vol);
      }
      if (fdn > 0 && Kokkos::atomic_fetch_add(&fdc(0), 1.0) < fdn) {
        printf("CLEAN WBFDIAG cycle %d t2s %d m %d kji %d %d %d rho %.4e T %.4e EO %.4e "
               "X %.4e E %.4e Q %.4e EI %.4e eg %.4e wk %.4e wp %.4e dw %.4e "
               "rd %.4e den %.4e take %.4e dm %.4e %.4e %.4e vol %.4e\n", fcy, fts, m, k,
               j, i, uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), x.w(m,CLW_EO,k,j,i),
               x.w(m,CLW_X,k,j,i), e, x.w(m,CLW_Q,k,j,i), x.w(m,CLW_EI,k,j,i), dg[0],
               dg[1], dg[2], dg[3], dg[4], den, take, dm[0], dm[1], dm[2], vol);
      }
      e = efl;
      Kokkos::atomic_add(&pc_(VET_POS_FLR), 1.0);
    }
    if (slopes) {
      kk(m,VET_T2_E,k,j,i) = (e + wsl - x.w(m,CLW_EO,k,j,i))/h;
      for (int n = 0; n < 3; ++n) kk(m,VET_T2_M1+n,k,j,i) = dh[n]/h;
      kk(m,VET_T2_EN,k,j,i) = den/h;
    }
    u0_(m,VET_E,k,j,i) = e;
    // the cell flux: the face momenta's cell vector (frame of y on the cubed sphere);
    // Phi's in cln_phic, above
    {
      const Real c[6] = {f1(m,k,j,i), f1(m,k,j,i+1), multi_d ? f2(m,k,j,i) : 0.0,
                         multi_d ? f2(m,k,j+1,i) : 0.0, three_d ? f3(m,k,j,i) : 0.0,
                         three_d ? f3(m,k+1,j,i) : 0.0};
      Real fv[3];
      CellVector(x, m, k, j, i, c, fv);
      ToHydroFrame(x, m, k, j, i, fv);
      for (int d = 0; d < 3; ++d) {
        u0_(m,VET_F1+d,k,j,i) = fv[d];
        if (!xh.hi) phq(m,d,k,j,i) = fv[d];   // no split: Phi is the face momentum
      }
    }
  };
  if (wf) {
    SolveFor("cln_wbprb", nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      wbc(m, k, j, i, true);
    });
    WbFct(x, px, wfp, pfo, pf);
  }
  if (wmc) Kokkos::deep_copy(xh.mclc, 0.0);   // (the write-back's own, not the probe's)
  SolveFor("cln_wbcell", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    wbc(m, k, j, i, false);
  });
  if (wmc) {
    auto hm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), xh.mclc);
    Real r[4];
    for (int n = 0; n < 4; ++n) r[n] = GlobalReduce(hm(n), false);
    if (global_variable::my_rank == 0 && r[0] > 0.0) {
      std::printf("CLEAN MCAP cycle %d t2s %d: %.0f face momentum terms at the cap in the "
                  "write-back, gas momentum withheld sum dp V (x1 x2 x3) %.6e %.6e %.6e\n",
                  pmy_pack->pmesh->ncycle, t2_solve, r[0], r[1], r[2], r[3]);
    }
  }
  if (lzb) {   // (the L_z ledger: the stage's L_z as its slope, then booked)
    const Real s = SumCells(lzw);
    if (slopes) (t2_solve == VET_T2S_STAGE1 ? cln_lzs2 : cln_lzs1) = s/h;
    LzBook(s, t2_solve == VET_T2S_STAGE2);
  }
  // xt_exact: the stage's closure |dE + de + dK + h div Phi| over E + E_old + eint, max
  // over the cells
  if (cls) {
    auto w = cln_w;
    const int n1 = ie - is + 1, nji = (je - js + 1)*n1, nkji = (ke - ks + 1)*nji;
    Real cm = 0.0;
    Kokkos::parallel_reduce("cln_wbclose", SolveRp256(DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx, Real &c) {
      const int m = idx/nkji, k = ks + (idx%nkji)/nji, j = js + (idx%nji)/n1;
      const Real v = w(m,CLW_PW,k,j,is + idx%n1);
      c = fmax(c, Kokkos::isfinite(v) ? v : INFINITY);   // Max drops a NaN
    }, Kokkos::Max<Real>(cm));
    cm = GlobalReduce(cm, true);
    if (global_variable::my_rank == 0 && (cln_verbose > 0 || !std::isfinite(cm))) {
      std::printf("CLEAN %sXT cycle %d t2s %d write-back closure %.3e\n",
                  std::isfinite(cm) ? "" : "WARNING ",
                  pmy_pack->pmesh->ncycle, t2_solve, cm);
    }
  }
  // implicit_pos_floor_cons: the uncovered raise of a cell from the radiation above
  // e_floor of the nearest cells of its radial column, i -+ d, split by their energy
  // above the floor, with their stage slopes.
  if (pfc) {
    par_for("cln_pos_flr_cons", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      CLEAN_HIPX(const SolveArgs &x = *px;)
      auto vol = [&](const int q) {
        return CellVolume(x, m, k, j, q);
      };
      Real mv = 0.0, un = 0.0;
      for (int i = is; i <= ie; ++i) {
        Real def = pfd(m,k,j,i);
        for (int d = 1; d <= ie - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(u0_(m,VET_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,VET_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,VET_E,k,j,q) = fmax(u0_(m,VET_E,k,j,q) - de, efl);
            if (slopes) kk(m,VET_T2_E,k,j,q) -= de/h;
          }
          def -= tk;
          mv += tk;
        }
        un += fmax(def, 0.0);
      }
      if (mv > 0.0) Kokkos::atomic_add(&pc_(VET_POS_FLR_MV), mv);
      if (un > 0.0) Kokkos::atomic_add(&pc_(VET_POS_FLR_UN), un);
      if (!pgc) return;
      // implicit_pos_gas_cons: the gas raise E could not cover, the same way, into the
      // cell's gas; what the column cannot give stays below e_min (the EOS floor ledger)
      Real gm = 0.0;
      for (int i = is; i <= ie; ++i) {
        Real def = pgd(m,k,j,i);
        for (int d = 1; d <= ie - is && def > 0.0; ++d) {
          const int qa = i - d, qb = i + d;
          const Real aa = (qa >= is) ? fmax(u0_(m,VET_E,k,j,qa) - efl, 0.0)*vol(qa) : 0.0;
          const Real ab = (qb <= ie) ? fmax(u0_(m,VET_E,k,j,qb) - efl, 0.0)*vol(qb) : 0.0;
          const Real at = aa + ab;
          if (!(at > 0.0)) continue;
          const Real tk = fmin(def, at);
          for (int sd = 0; sd < 2; ++sd) {
            const int q = sd ? qb : qa;
            const Real a = sd ? ab : aa;
            if (!(a > 0.0)) continue;
            const Real de = (tk*a/at)/vol(q);
            u0_(m,VET_E,k,j,q) = fmax(u0_(m,VET_E,k,j,q) - de, efl);
            if (slopes) kk(m,VET_T2_E,k,j,q) -= de/h;
          }
          const Real dg = tk/vol(i);
          uh(m,IEN,k,j,i) += dg;
          x.w(m,CLW_JV,k,j,i) += dg;
          if (slopes) kk(m,VET_T2_EN,k,j,i) += dg/h;
          def -= tk;
          gm += tk;
        }
      }
      if (gm > 0.0) Kokkos::atomic_add(&pc_(VET_POS_GAS_DE), gm);
    });
  }
  if (pflr) {
    Real un = 0.0;
    Kokkos::deep_copy(un, Kokkos::subview(pc_, VET_POS_FLR_UN));
    un = GlobalReduce(un - un0, false);
    cln_flrun = un;
    cln_flrid = 4*pmy_pack->pmesh->ncycle + t2_solve;
    if (un > 0.0 && global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d t2s %d: implicit_pos_floor created %.6e of E "
                  "(neither the gas nor the column covered the raise)\n",
                  pmy_pack->pmesh->ncycle, t2_solve, un);
    }
  }
  // cells the write-back left with gas eint < 0 or below the EOS floor (the solve is then
  // not accepted; an eint the hydro left there is not the solve's)
  cln_neint = 0;
  if (!feedback) return;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nx1 = id.nx1, nx2 = id.nx2, nx3 = id.nx3, nji = nx2*nx1, nkji = nx3*nji;
  Real nb = 0.0, nh = 0.0;
  Kokkos::parallel_reduce("cln_neint", SolveRp256(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n, Real &nf) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji + ks;
    const int j = (idx%nji)/nx1 + js, i = idx%nx1 + is;
    Real p[3];
    GasMomentum(x, uh, m, k, j, i, p);
    const Real ke = 0.5*(SQR(p[0]) + SQR(p[1]) + SQR(p[2]))/uh(m,IDN,k,j,i);
    // and below the eint of the EOS temperature floor (the hydro floor would add energy
    // the radiation's budget does not see)
    Real ef = 0.0;
    if (eos.tfloor > 0.0) {
      Real pr, chi_rho, chi_t, cv;
      eos.ThermoAt(uh(m,IDN,k,j,i), eos.tfloor, ef, pr, chi_rho, chi_t, cv);
    }
    const Real tz = 1.0e-12*(fabs(uh(m,IEN,k,j,i)) + ke + ef);   // eint's round-off
    // and any non-finite E, F or gas state, also at the floor: rejected and retried,
    // never accepted
    const bool fin = Kokkos::isfinite(u0_(m,VET_E,k,j,i) + u0_(m,VET_F1,k,j,i) +
                                      u0_(m,VET_F1+1,k,j,i) + u0_(m,VET_F1+2,k,j,i) +
                                      uh(m,IDN,k,j,i) + uh(m,IEN,k,j,i) + ke);
    const bool low = !(GasEint(x, uh, ke, m, k, j, i) > ef - tz);
    if ((low && !(x.w(m,CLW_JV,k,j,i) >= 0.0) && x.w(m,CLW_EI,k,j,i) > ef + tz) ||
        (low && x.w(m,CLW_EI,k,j,i) < 0.0) ||
        u0_(m,VET_E,k,j,i) < 0.0 || !fin) n += 1.0;
  }, Kokkos::Sum<Real>(nb), Kokkos::Sum<Real>(nh));
  cln_neint = static_cast<int>(GlobalReduce(nb, false));
  nh = GlobalReduce(nh, false);
  if (nh > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: %.0f cells left with gas eint below the "
                "EOS floor (implicit_pos_gas off: the hydro's floor raises them, energy "
                "created)\n", pmy_pack->pmesh->ncycle, t2_solve, nh);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::UpdateState
//! \brief the new state (WriteBack)

void RadiationVet::UpdateState(const SolvePass &p) {
  // a failed solve or a probe does not change the energy flux's memory
  cln_dok = (p.nit >= 0) && !t2_inprobe;
  WriteBack();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VerdictState
//! \brief what the new state means: gas eint below its floor (a failed stage, dt_retry or
//! a warning); K1 and the reference force of an accepted solve

void RadiationVet::VerdictState(const SolvePass &p) {
  // a write-back that leaves gas eint < 0 is rejected like a failed solve (stage redone
  // by backward Euler, a backward-Euler step by dt_retry)
  if (cln_feedback && cln_neint > 0 && !p.t2stage && dtr_on) {
    if (!dtr_fail) {
      dtr_why = 4;   // reason eint<=0
      dtr_resid = cln_neint;
    }
    dtr_fail = true;
  }
  if (cln_feedback && cln_neint > 0 && !p.t2stage && !dtr_on && !t2_inprobe &&
      global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d t2s %d: backward-Euler step accepted with %d "
                "cells of gas eint below 0 / its floor, E < 0 or not finite (no "
                "dt_retry)\n", pmy_pack->pmesh->ncycle, t2_solve, cln_neint);
  }
  if (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2) {
    // implicit_pos_floor(_cons): an E < 0 raised by the conservative write-back floor
    // does not fail the stage (counted)
    const bool flr = impl_pos_floor;
    t2_fail = (p.nit < 0) || (p.nneg > 0 && !flr) || (cln_feedback && cln_neint > 0) ||
              (t2_solve == VET_T2S_STAGE1 && pmy_pack->pmesh->ncycle == t2_dbg_fail);
    if (cln_lzbook && cln_feedback && t2_solve == VET_T2S_STAGE2 && !t2_fail) {
      LzBook(0.0, false);   // the accepted stage 2's L_z, pending until here
    }
  }
  // K1 = f(U^n) for the next predictor and time error: from an accepted backward-Euler
  // or stage-2 solve; a failed one leaves none (BE redo and auto-BE steps keep the last)
  if (t2_solve == VET_T2S_BESTORE || t2_solve == VET_T2S_STAGE2) {
    cln_k1ok = (p.nit >= 0) && !dtr_fail && !(t2_solve == VET_T2S_STAGE2 && t2_fail);
  }
  // force_reference = clean: the reference is the force of the converged, accepted first
  // coupled solve (each iterate recaptures it; a failed one is redone)
  if (force_ref == VET_FREF_CLEAN && cln_feedback && !cln_vzero && p.nit >= 0 &&
      p.nneg == 0 && !(p.t2stage && t2_fail) && !dtr_fail && !t2_inprobe) {
    cln_fref_set = true;
  }
}

}  // namespace radvet
