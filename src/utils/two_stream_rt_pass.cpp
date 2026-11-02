//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt_pass.cpp
//! \brief picket_fence_two_stream_RT_pass, the two-stream column pass (declared in
//! two_stream_rt.hpp).  tusplit-1102: moved out of two_stream_rt.hpp verbatim so that it
//! is compiled once, in its own translation unit, instead of in every problem generator
//! that includes the header.

#include <cstdio>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>

#include "utils/two_stream_rt.hpp"
#include "utils/ck_2nd.hpp"     // ck2nd-1102 ck_vef_tint (read only when named)

namespace two_stream_rt {
// ck2nd-1102 (utils/ck_2nd.hpp (4b)): the VEF refresh test, in simulated time when named
inline bool CkVefDue2(Mesh *pm, const bool ckfop, const bool ckfst) {
  if (!(ck_vef_dt > 0.0)) return CkVefDue(pm, ckfop, ckfst);
  const bool stp = !ck_vef_sync || !ckfop || ckfst;
  if (ck_cad_dt > 0.0) {   // ck2ndt-1102: stateless, on the multiples of ck_vef_dt
    return stp && (ck_vef_last < 0 || ck_vef_s_ptr == nullptr ||
                   CkSchedCross(pm->time, pm->dt, ck_vef_dt));
  }
  return stp && (ck_vef_last < 0 || ck_vef_s_ptr == nullptr ||
                 pm->time + pm->dt >= ck_vef_tlast + ck_vef_dt);
}
}  // namespace two_stream_rt

namespace two_stream_rt {

void picket_fence_two_stream_RT_pass(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    // Noti+2023; Lee+2021

    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;

    // --- PLANE-PARALLEL geometry (problem/rt_plane_parallel).  See the flag's note: the
    // three radial geometry Views above are placeholders on a Cartesian mesh, so in this
    // mode every read of them goes through X1V / X1F / DX1 instead, which rebuild the
    // same quantities from the MeshBlock's RegionSize.  With the flag off each helper
    // returns the View element it replaced, so the radial path is unchanged bit for bit.
    const bool pp_ = rt_plane_parallel;
    if (pp_ && !(rt_grey && rt_split)) {
      std::cout << "### FATAL ERROR in two_stream_rt: rt_plane_parallel is implemented "
                << "for the GREY split sweep only (rt_grey = true, which implies "
                << "rt_split); the picket-fence and correlated-k kernels carry radial "
                << "stellar-beam geometry." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (rt_layer_legacy && !rt_plane_parallel && rt_grey) {
      std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_layer_legacy (the old "
                << "staggered whole-cell layers) was NOT converted to the spherical "
                << "form of the grey sweep and still carries the J = A I substitution, "
                << "which is wrong in the diffusion limit by 1 - H_T/(2r).  Use the "
                << "centre-to-centre layers (rt_layer_legacy = false) on a curvilinear "
                << "mesh." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (pp_ && (use_cubed_sphere_ || use_spherical_polar)) {
      std::cout << "### FATAL ERROR in two_stream_rt: rt_plane_parallel is for a "
                << "Cartesian mesh, but this mesh is spherical-polar or cubed-sphere."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/ck_beam_sph (see the flag's note): the pseudo-spherical direct beam.  Pure
    // geometry, so meaningless on a Cartesian mesh and refused there; belongs to the ck
    // kernel; and the legacy staggered ck layers carry their own copy of the beam, which
    // was not converted.
    if (ck_beam_sph) {
      if (!(use_cubed_sphere_ || use_spherical_polar)) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_beam_sph is the "
                  << "PSEUDO-SPHERICAL stellar beam and has no meaning on a Cartesian "
                  << "mesh.  Leave it false there." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!rt_ck) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_beam_sph converts the "
                  << "CORRELATED-K beam, but problem/rt_ck is false." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rt_layer_legacy) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_beam_sph with "
                  << "problem/rt_layer_legacy: the staggered whole-cell correlated-k "
                  << "layers carry their own copy of the beam and were not converted."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      static bool ckb_announced = false;
      if (!ckb_announced && global_variable::my_rank == 0) {
        ckb_announced = true;
        std::cout << "### two_stream_rt: PSEUDO-SPHERICAL stellar beam ON "
                  << "(problem/ck_beam_sph): straight rays through the shells, no secant "
                  << "clamp, twilight (mu0 < 0) lit down to the grazing radius."
                  << std::endl;
      }
    }
    // problem/ck_spherical (see the flag's note): the spherical form of the correlated-k
    // thermal two-stream.  It is geometry, so it is meaningless on a Cartesian mesh and
    // REFUSED there rather than silently inert; it belongs to the ck kernel alone; and
    // the legacy staggered ck layers were never converted.
    if (ck_spherical) {
      if (!(use_cubed_sphere_ || use_spherical_polar)) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_spherical is the "
                  << "SPHERICAL form of the correlated-k two-stream and has no meaning "
                  << "on a Cartesian mesh.  Leave it false there." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!rt_ck) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_spherical converts "
                  << "the CORRELATED-K kernel, but problem/rt_ck is false.  The grey "
                  << "sweep is spherical already (8bca3dfa) and needs no switch."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rt_layer_legacy) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_spherical with "
                  << "problem/rt_layer_legacy: the staggered whole-cell correlated-k "
                  << "layers were NOT converted to the spherical form.  Use the "
                  << "centre-to-centre layers (rt_layer_legacy = false)." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      static bool cksph_announced = false;
      if (!cksph_announced && global_variable::my_rank == 0) {
        cksph_announced = true;
        std::cout << "### two_stream_rt: correlated-k SPHERICAL form ON "
                  << "(problem/ck_spherical): face mixing on both rays, Phi = A F "
                  << "single-valued for the THERMAL two-stream; the stellar beam is "
                  << "unchanged (a parallel pencil is not confined to the column)."
                  << std::endl;
      }
    }
    // problem/ck_sweep_form (see the flag's note): the probe-free two-pass forms.  They
    // ARE the spherical face coupling solved exactly, so they have nothing to do without
    // it; they are instantiated only at ck_sweep_cache = 2; and only tm carries the
    // implicit column solve's Jacobian (the JAC instantiation) and the frozen operator
    // (the FOP instantiation, DESIGN_tm.md phase T2).
    if (ck_sweep_form != 0) {
      if (!ck_spherical) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_sweep_form solves "
                  << "the SPHERICAL face coupling exactly, and problem/ck_spherical is "
                  << "false.  Without the coupling the plane-parallel sweep is already "
                  << "the exact two-pass form; leave ck_sweep_form = 0 there."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ck_sweep_cache != 2) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_sweep_form is "
                  << "instantiated at problem/ck_sweep_cache = 2 only (the default "
                  << "wherever it can run), to keep the kernel's instantiation count "
                  << "down.  Got ck_sweep_cache = " << ck_sweep_cache << "."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (ck_implicit && ck_sweep_form != 1) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_sweep_form = 2 (sd) "
                  << "with problem/ck_implicit: the column solve's tridiagonal has been "
                  << "derived for the four-pass and the tm (ck_sweep_form = 1) "
                  << "recurrences only." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // problem/ck_impl_lin: the linear re-apply kernel exists for the tm sweep with the
    // stored operator (the triple included) and nothing else.
    if (ck_impl_lin) {
      if (!ck_implicit || ck_sweep_form != 1 || !ck_impl_frozen_op ||
          !ck_impl_frozen_cof || ck_impl_refresh_kappa || RT_FP32 != 0) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_impl_lin needs "
                  << "ck_implicit, ck_sweep_form = 1, ck_impl_frozen_op, "
                  << "ck_impl_frozen_cof, ck_impl_refresh_kappa = false and a "
                  << "double-precision chain (RT_FP32 = 0)." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    // problem/ck_impl_glob lives in the fused step only
    if (ck_impl_glob > 0 && (!ck_impl_fuse || ck_impl_debug > 0)) {
      std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_impl_glob needs "
                << "problem/ck_impl_fuse (and ck_impl_debug <= 0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/ck_impl_esc refines ls_sub; problem/ck_impl_aa lives in the fused step
    if ((ck_impl_esc > 0 && ck_impl_glob != 2)
        || (ck_impl_aa > 0 && (!ck_impl_fuse || ck_impl_debug > 0))
        || ck_impl_aa < 0 || ck_impl_aa > CK_AA_MAX) {
      std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_impl_esc needs "
                << "ck_impl_glob = ls_sub, and problem/ck_impl_aa (0.." << CK_AA_MAX
                << ") needs ck_impl_fuse (and ck_impl_debug <= 0)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/ck_impl_cvsec lives in the fused step only
    if (ck_impl_cvsec && !ck_impl_fuse) {
      std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_impl_cvsec needs "
                << "problem/ck_impl_fuse." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // problem/ck_impl_jac_lin: the Jacobian kernel reads the factorisation the linear
    // kernel stores and writes the one-chain-per-thread partials
    if (ck_impl_jac_lin) {
      if (!ck_impl_lin || ck_impl_lin_thr != 1) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_impl_jac_lin needs "
                  << "ck_impl_lin with ck_impl_lin_thr = 1; got ck_impl_lin="
                  << ck_impl_lin << " ck_impl_lin_thr=" << ck_impl_lin_thr << "."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (ck_sweep_form != 0) {
      static bool ckfrm_announced = false;
      if (!ckfrm_announced && global_variable::my_rank == 0) {
        ckfrm_announced = true;
        std::cout << "### two_stream_rt: correlated-k PROBE-FREE sweep ON "
                  << "(problem/ck_sweep_form = " << ck_sweep_form << ", "
                  << ((ck_sweep_form == 1) ? "tm: transfer-matrix/adding in (u,d)"
                                           : "sd: S = (u+d)/2, D = (u-d)/2")
                  << "): two column passes, the face coupling exact." << std::endl;
      }
    }
    const int is_pp = is;
    const int nx1_pp = indcs.nx1;
    auto X1V = KOKKOS_LAMBDA (const int m, const int i) {
      return pp_ ? CellCenterX(i-is_pp, nx1_pp, size.d_view(m).x1min,
                               size.d_view(m).x1max)
                 : x1v_(m,i);
    };
    auto X1F = TsrtX1F{pp_, is_pp, nx1_pp, size, x1f_};
    auto dx1_ = pmbp->pcoord->dx1;
    auto DX1 = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return pp_ ? size.d_view(m).dx1 : dx1_(m,k,j,i);
    };
    // ---- THE TRANSVERSE ARC LENGTHS.  On a CURVILINEAR mesh size.dx2/dx3 are the
    // COORDINATE spacings -- on the cubed sphere the gnomonic angles xi, eta in [-1,1],
    // i.e. dimensionless -- while pcoord->dx2/dx3 are the per-cell physical arc lengths
    // (~ r dxi), the same distinction conduction.cpp:1069 makes with its `curv` flag.
    // Every transverse DERIVATIVE in this file must divide by the arc length; the
    // coordinate spacing is off by ~r (2e11 on the He star).  On a Cartesian mesh
    // pcoord->dx2/dx3 are 1x1x1x1 dummies, so the flag also decides which is readable.
    const bool curv_ = pmbp->pmesh->use_cubed_sphere || pmbp->pmesh->use_spherical_polar;
    auto dx2_c_ = pmbp->pcoord->dx2;
    auto dx3_c_ = pmbp->pcoord->dx3;
    auto DX2 = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return curv_ ? dx2_c_(m,k,j,i) : size.d_view(m).dx2;
    };
    auto DX3 = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return curv_ ? dx3_c_(m,k,j,i) : size.d_view(m).dx3;
    };

    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    // problem/rt_use_cons: where the solver reads the cell's thermodynamic state from.
    // eiN/rhoN are the ONLY way this routine touches e and rho below, so the switch
    // reaches the precompute, the sweeps, the equilibrium, the Newton, the positivity
    // guard, the rescue floor and LimitRTSource together -- there is no half-converted
    // path.  With the switch off they return exactly w0(...), so the default is
    // bit-identical by construction.
    const bool usecons_ = rt_use_cons;
    const bool bface_on = rt_bface;
    const bool cs_uc_ = pm->use_cubed_sphere;
    auto cosc_uc_ = pmbp->pcoord->cos_cell;
    const bool etg_uc_ = (pmbp->phydro != nullptr) ? pmbp->phydro->use_etotgrav
                       : ((pmbp->pmhd != nullptr) ? pmbp->pmhd->use_etotgrav : false);
    DvceArray4D<Real> phicc_uc_ = CkDum<DvceArray4D<Real>>("rt_phi_dummy");
    if (pmbp->phydro != nullptr) {
      phicc_uc_ = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      phicc_uc_ = pmbp->pmhd->phicc0;
    }
    auto u0_uc_ = u0;
    auto w0_uc_ = w0;
    // MHD: u0(IEN) also holds the magnetic energy, so the extraction has to subtract it
    // or every temperature this solver forms is the temperature of e + ME.  bcc0 is the
    // cell-centred form of the CURRENT b0 (see MagEnergyCC): the RT apply and its
    // precompute run from MHD::MHDSrcTerms, inside the stage but BEFORE MHD::CT, so b0
    // has not moved since the ConToPrim that filled bcc0 and re-averaging the faces here
    // would return the same numbers.  On the cubed sphere bcc0 is already in the
    // orthonormal frame, so no metric enters.
    const bool mhd_uc_ = (pmbp->pmhd != nullptr);
    DvceArray5D<Real> bcc_uc_ = CkDum<DvceArray5D<Real>>("rt_bcc_dummy");
    if (pmbp->pmhd != nullptr) {
      bcc_uc_ = pmbp->pmhd->bcc0;
    }
    // see the note on rt_eiclamp_cnt: e <= 0 (or NaN) out of the conserved state is
    // clamped to e(rho, tfloor) before any temperature or Planck function is formed
    if (rt_eiclamp_cnt == nullptr) {
      rt_eiclamp_cnt = new DvceArray1D<int>("rt_eiclamp", 1);
      Kokkos::deep_copy(*rt_eiclamp_cnt, 0);
    }
    auto eicl_g = *rt_eiclamp_cnt;
    auto eos_uc_ = eos;
    auto eiN = TsrtEiN{usecons_, w0_uc_, u0_uc_, cs_uc_, cosc_uc_, etg_uc_, phicc_uc_,
                       mhd_uc_, bcc_uc_, eicl_g, eos_uc_};
    // problem/ck_impl_nosync: these two once-only warnings are checked on the first
    // pass of a call only (one read-back each instead of one per pass)
    const bool clchk_ = !(ck_impl_nosync && ck_impl_pass > 0);
    if (usecons_ && !rt_eiclamp_warned && clchk_) {
      auto eicl_h = Kokkos::create_mirror_view(eicl_g);
      Kokkos::deep_copy(eicl_h, eicl_g);
      if (eicl_h(0) > 0) {
        rt_eiclamp_warned = true;
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING in two_stream_rt: rt_use_cons gave a non-positive "
                    << "internal energy in " << eicl_h(0) << " cell read(s); clamped to "
                    << "e(rho,tfloor) before the temperature/Planck evaluation. "
                    << "Reported once." << std::endl;
        }
      }
    }
    // see rt_stclamp_cnt: a non-positive or non-finite density, and a (p,T) the EOS
    // could not form from it, are clamped before any Planck function or optical depth
    if (rt_stclamp_cnt == nullptr) {
      rt_stclamp_cnt = new DvceArray1D<int>("rt_stclamp", 1);
      Kokkos::deep_copy(*rt_stclamp_cnt, 0);
    }
    auto stcl_g = *rt_stclamp_cnt;
    if (!rt_stclamp_warned && clchk_) {
      auto stcl_h = Kokkos::create_mirror_view(stcl_g);
      Kokkos::deep_copy(stcl_h, stcl_g);
      if (stcl_h(0) > 0) {
        rt_stclamp_warned = true;
        if (global_variable::my_rank == 0) {
          std::cout << "### WARNING in two_stream_rt: a non-finite or non-positive "
                    << "density/(p,T) state was read in " << stcl_h(0) << " cell(s); "
                    << "clamped to the floor state before the Planck function and the "
                    << "optical depth, so the column below it is not poisoned. "
                    << "Reported once." << std::endl;
        }
      }
    }
    const Real dfl_uc_ = eos.dfloor;
    auto rhoN = TsrtRhoN{usecons_, u0_uc_, w0_uc_, stcl_g, dfl_uc_};
    {
      static bool uc_announced = false;
      if (!uc_announced && global_variable::my_rank == 0) {
        uc_announced = true;
        std::cout << "### two_stream_rt: thermodynamic state read from "
                  << (usecons_ ? "u0 (CONSERVED, problem/rt_use_cons = true)"
                               : "w0 (previous ConToPrim; rt_use_cons = false)")
                  << std::endl;
      }
    }

    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;

    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;

    // ================================================================================
    // SPHERICAL DILUTION: WHY THE DEPOSIT *AND* THE PROPAGATION CARRY THE AREA
    //
    // SUPERSEDED FOR THE GREY CENTRE-TO-CENTRE SWEEP AND FOR MODE 3.  The J = A I form
    // derived below is exact in the TRANSPARENT limit and wrong in the DIFFUSION limit
    // by F_2s/F_exact = 1 - H_T/(2r) (tests_r2/thick measures it to 1e-4), because it
    // puts the area on BOTH two-stream moments while the exact pair puts it on the
    // antisymmetric one alone.  The grey centre-to-centre sweep and the mode-3 column
    // solve now carry PLAIN intensities with the face mixing of BTF/MIXF instead -- see
    // the SPHERICAL FORM note on RTCol3::Bt in two_stream_column_implicit.hpp, which is
    // the derivation that applies.  What follows still describes the LEGACY whole-cell
    // layers (refused on a curvilinear mesh), the picket-fence chain kernel and the
    // correlated-k kernel, which were not converted.
    //
    // The sweep integrates, per ray direction mu, the plane-parallel transfer equation
    //     mu dI/dr = -kappa rho (I - B),
    // and used to deposit -(F_top - F_bot)/dx1 with F = 2 pi w mu I.  On a RADIAL mesh
    // that is wrong twice over, and the two errors have to be fixed together:
    //   (a) the deposit is a plane-parallel divergence.  What the shell actually loses
    //       is (A_top F_top - A_bot F_bot)/V, and over r_out/r_in = 2.5 the area ratio
    //       is 6.25 -- the emergent luminosity of a global star was wrong by that much;
    //   (b) the PROPAGATION does not dilute.  A transparent shell returns I_out = I_in,
    //       i.e. F_out = F_in, so L = A F grows like r^2 -- and had only (a) been fixed
    //       that same transparent shell would be handed a spurious cooling
    //       -(A_top - A_bot) F/V = -2 F/r.  NEITHER FORM ALONE IS CONSERVATIVE.
    //
    // The exact statement is the zeroth moment, div F = kappa rho (4 pi B - c E), i.e.
    // (1/r^2) d(r^2 F)/dr on the left.  Multiply the transfer equation by r^2 and write
    // J = r^2 I (equivalently J = A I; the column's solid angle cancels):
    //     mu dJ/dr = -kappa rho (J - A B),
    // which is the SAME equation the layer solve already integrates exactly, with the
    // source function A(r) B(r) in place of B(r).  So the sweep keeps its form and
    // carries the AREA-WEIGHTED intensity J instead of I:
    //   * every layer source endpoint is multiplied by the area AT THAT ENDPOINT -- a
    //     face endpoint by area.x1f (AFC), a cell-centre endpoint by the cell's mean
    //     area V/dx1 (ACC);
    //   * the boundary data are scaled the same way (the top incoming intensity by
    //     A(ie+1), the internal flux entering the bottom face by A(is));
    //   * the face flux handed to the rest of the code stays a flux PER UNIT AREA,
    //     recovered as F = 2 pi w mu J / A;
    //   * the deposit is (J_in - J_out), summed over rays, divided by the cell VOLUME.
    //     It telescopes to exactly -(A_top F_top - A_bot F_bot)/V, so the column
    //     conserves energy cell by cell and L = A F is constant wherever the local
    //     balance vanishes -- which is the 1-D radiative-equilibrium gate.
    // The ray-curvature term of the true spherical transfer equation, (1-mu^2)/r dI/dmu,
    // is NOT represented: this is the standard radial-ray approximation, exact in the
    // two limits that matter here (diffusion, and a radially streaming flux) and
    // conservative everywhere.
    //
    // PLANE-PARALLEL INERTNESS.  Under problem/rt_plane_parallel the helpers return
    // A = 1 exactly and V = the same dx1 the expression used before, so every new
    // factor is a multiplication or a division by 1.0 -- exact in IEEE arithmetic --
    // and the Cartesian box is bitwise what it was.  The helpers are also the only
    // reads of the Coordinates area/volume Views, which are 1x1x1x1 placeholders on a
    // Cartesian mesh (coordinates.cpp:63-94) exactly like x1v/dx1, so the ternary is
    // what keeps them unread there.
    //
    // WHAT IS NOT CONVERTED: the monolithic (rt_split = false) sweep, the correlated-k
    // kernel and the mode-1/2 nearest-neighbour Jacobian still carry the plane-parallel
    // divergence.  They serve the box and the r_out/r_in < 1.3 hot-Jupiter/red-giant
    // domains, where the error is at the per-cent level.  The GREY SPLIT sweep
    // (rt_grey + rt_split), the apply kernel and the mode-3 column solve -- the path a
    // global star runs on -- are area/volume correct.
    // ================================================================================
    // the face area, the cell's MEAN area V/dx1 (the weight of a cell-centre source
    // endpoint), and the cell volume.  THE TWO VOLUME HELPERS DIFFER ONLY IN THEIR
    // PLANE-PARALLEL BRANCH, and each kernel must use the one matching the radial width
    // it already spells: VLS for the picket-fence chain kernel, which indexes the raw
    // dx1 View, and VLA wherever the code goes through DX1 (the grey sweep and the
    // apply kernel).  That is what makes the substitution bitwise inert on the box --
    // getting it the wrong way round reads the 1x1x1x1 dx1 placeholder there.
    auto AFC = TsrtAFC{pp_, area1};
    auto ACC = TsrtACC{pp_, volume, dx1};
    auto VLS = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return pp_ ? dx1(m,k,j,i) : volume(m,k,j,i);
    };
    auto VLA = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return pp_ ? DX1(m,k,j,i) : volume(m,k,j,i);
    };
    // ---- THE SPHERICAL FORM of the grey centre-to-centre sweep.  The J = A I
    // substitution above is kept ONLY by the legacy layers, the picket-fence chain and
    // the correlated-k kernel; the grey c2c sweep -- what a global star actually runs --
    // carries PLAIN intensities and puts the geometry in the FACE conditions instead.
    // See the derivation on RTCol3::Bt in two_stream_column_implicit.hpp: the area
    // belongs to the ANTISYMMETRIC moment alone, so at every face S = (I+ + I-)/2 and
    // A D = A (I+ - I-)/2 are passed continuously, which is the single correction
    //     c = beta (d_above - u_below),   beta = (A_above - A_below)/(A_above + A_below)
    // added to BOTH sides.  Inside a cell the frame area is constant (A_i = V_i/dx_i),
    // so the layer solve and the deposit are the plane-parallel ones they always were.
    // Both helpers are identically inert under rt_plane_parallel: beta is exactly 0 and
    // AUN/ACN are the very area factors the old expressions read.
    auto BTF = TsrtBTF{pp_, area1, volume, dx1, ie};
    auto MIXF = TsrtMIXF{};
    // the area factors the OLD expression still spells out and the spherical form no
    // longer carries: exactly 1.0, read from the Views under rt_plane_parallel so that
    // the box's expression DAG -- and hipcc's FMA contraction -- is untouched
    auto AUN = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return pp_ ? AFC(m,k,j,i) : 1.0;
    };
    auto ACN = KOKKOS_LAMBDA (const int m, const int k, const int j, const int i) {
      return pp_ ? ACC(m,k,j,i) : 1.0;
    };

//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real met = 0.0;

//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    const bool tide = pm->pgen->hot_jupiter_param.stellar_tide;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    Real met = pm->pgen->hot_jupiter_param.met;

    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;

    Real Tirr = Teq*sqrt(2);
    Real Tirr4 = SQR(SQR(Tirr));
    Real Fstar = boltz_sigma*Tirr4;
    Real Tint;
    if (rt_tint_override > 0.0) {
      Tint = rt_tint_override;
    } else {
      get_Tint(Teq, Tint);
    }
    Real Tint4 = SQR(SQR(Tint));
    bool int_at_cut = rt_int_at_cut;
    const Real bot_flux = rt_bot_flux;          // see rt_bottom_flux
    const bool layer_legacy = rt_layer_legacy;  // see rt_layer_legacy
    const bool cksph_ = ck_spherical;           // see ck_spherical
    const bool ckbsph_ = ck_beam_sph;           // see ck_beam_sph
    const int ckcache_ = ck_sweep_cache;         // see ck_sweep_cache
    const int ckform_ = ck_sweep_form;           // see ck_sweep_form
    const bool top_re = rt_top_re;
    const bool top_vac = rt_top_vacuum;
    Real Iint = boltz_sigma/M_PI*Tint4;

    const int nchain_rt = rt_nchain;

    Real mug[2];
    Real wg[2];
    mug[0] = 0.21132487;
    mug[1] = 0.78867513;
    wg[0] = 0.5;
    wg[1] = 0.5;

    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);
    if (rt_albedo_set >= 0.0) albedo = rt_albedo_set;   // problem/albedo
    rt_hist_fstar = Fstar;                               // for problem/flux_hst
    rt_hist_albedo = albedo;


    // ================================================================================
    // CHAIN-PARALLEL SPLIT PATH  (problem/rt_split)
    // --------------------------------------------------------------------------------
    // The monolithic kernel below is parallel over (m,k,j) only -- 8192 columns here,
    // 128 wavefronts on 1216 SIMDs -- so it runs at MeanOccupancyPerCU ~0.5 and
    // VALUBusy 4-6 %, and the frequency chains are a SERIAL loop inside each thread.
    // That is why the cost is exactly linear in the chain count: it is exp() latency on
    // a dependency chain with no other wave to hide it, not a throughput limit.
    //
    // This path splits the same arithmetic into three kernels so the chain-block index
    // becomes a parallel dimension:
    //   A "rt_pre"    (m,k,j)       chain-independent: tau, B, Q_v and the picket-fence
    //                               coefficients, written to global arrays.
    //   B "rt_chain"  (m,blk,k,j)   one thread per column per block of RT_NB chains.
    //                               nblk x more threads than the monolithic kernel.
    //   C "rt_apply"  (m,k,j,i)     sum the per-block fluxes in block order and apply.
    //
    // Kernel C sums the blocks in the SAME order the serial loop accumulated them,
    // starting from the same 0.0, so with one block (production, 4 chains) the result is
    // bitwise identical to the monolithic path. With the harness active the extra blocks
    // contribute exactly 0.0, so it stays bitwise identical at any chain count -- which
    // is what makes the speed-up measurable against an unchanged answer.
    if (rt_split) {
      constexpr int NC = RT_NB;
      // grey runs one chain: one band, one column sweep, whatever the angular quadrature
      const int nblk = rt_grey ? 1 : (nchain_rt + NC - 1)/NC;
      if (rt_tau_ptr == nullptr) {
        const int nmb = pmbp->nmb_thispack;
        // Deliberately leaked, like the k-table: a namespace-scope View would outlive
        // Kokkos::finalize(). rt_Fb is the only large one, nmb*nblk*n3*n2*n1 Reals.
        rt_tau_ptr = new DvceArray4D<Real>("rt_tau", nmb, n3, n2, n1);
        rt_B_ptr   = new DvceArray4D<Real>("rt_B",   nmb, n3, n2, n1);
        rt_Qv_ptr  = new DvceArray4D<Real>("rt_Qv",  nmb, n3, n2, n1);
        // slot 4 (roche-dhj-1031): the Roche factor |grad Phi|/g_pm of the column top for
        // the tau = kappa p/g closure (EffGravAt x cf(4)); 1 unless problem/roche
        rt_cf_ptr  = new DvceArray4D<Real>("rt_cf",  nmb, n3, n2, 5);
        Kokkos::deep_copy(Kokkos::subview(*rt_cf_ptr, Kokkos::ALL, Kokkos::ALL,
                                          Kokkos::ALL, 4), 1.0);
        // (m, slot, i, k, j). Threads in a wave vary in j -- see the flattening in
        // par_for -- so j has to be the FASTEST array index. With j second-slowest, as
        // (m,slot,k,j,i) had it, adjacent lanes were 544 bytes apart and each one pulled
        // its own cache line: 8 useful bytes out of every 64 fetched.
        rt_Fb_ptr  = new DvceArray5D<Real>("rt_Fb",  nmb, nblk, n1, n3, n2);
        rt_Em_ptr  = new DvceArray5D<Real>("rt_Em",  nmb, nblk, n1, n3, n2);
        rt_Src_ptr = new DvceArray5D<Real>("rt_Src", nmb, nblk, n1, n3, n2);
        if (rt_ck || rt_grey) {
          const int nb_a = rt_ck ? CK_NB : 1;
          rt_kc_ptr = new DvceArray5D<Real>("rt_kc", nmb, nb_a, n1, n3, n2);
          rt_Bb_ptr = new DvceArray5D<Real>("rt_Bb", nmb, nb_a, n1, n3, n2);
          rt_T_ptr  = new DvceArray4D<Real>("rt_T",  nmb, n3, n2, n1);
          rt_pb_ptr = new DvceArray4D<Real>("rt_pb", nmb, n3, n2, n1);
          ck_pbar_ptr = rt_pb_ptr;     // ck_impl_tol_ptop reads it in the column step
          rt_xT_ptr = new DvceArray4D<Real>("rt_xT", nmb, n3, n2, n1);
          rt_xP_ptr = new DvceArray4D<Real>("rt_xP", nmb, n3, n2, n1);
          rt_icut_ptr = new DvceArray3D<int>("rt_icut", nmb, n3, n2);
          rt_Qb_ptr = new DvceArray5D<Real>("rt_Qb", nmb, nblk, n1, n3, n2);
          RtFhRstAfterAlloc(nmb);   // problem/flux_hst_rst: the restart's cache
          // problem/ck_implicit: the column solve's own scratch.  Nothing here exists
          // with the switch off (see utils/two_stream_column_ck.hpp).
          if (ck_implicit) {
            CkImplAlloc(nmb, nb_a, nblk*NC, n1, n2, n3, nblk);
            CkRstAfterAlloc();     // restart: ck_thk / ck_cv (two_stream_ck_rst.hpp)
          }
        }
        if (rt_diag) {
          rt_diag_ptr = new DvceArray5D<Real>("rt_diag", nmb, 7, n3, n2, n1);
        }
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: RT split path ON, " << nblk
                    << " chain block(s) of " << NC << " -> "
                    << nblk << "x the thread count of the monolithic kernel"
                    << std::endl;
        }
      }
      auto tau_g = *rt_tau_ptr;
      auto B_g   = *rt_B_ptr;
      auto Qv_g  = *rt_Qv_ptr;
      auto cf_g  = *rt_cf_ptr;
      auto Fb_g  = *rt_Fb_ptr;
      auto Em_g  = *rt_Em_ptr;
      auto Src_g = *rt_Src_ptr;
      // ---- problem/ck_implicit: the column solve's state.  Every one of these is a
      // 1-element dummy with the switch off, captured and never read.
      const bool ckimp_ = ck_implicit;
      const bool ckpass0_ = ckimp_ && (ck_impl_pass <= 0);
      // the opacity is FROZEN over the step by default: pass 0 builds it, later passes
      // rebuild only the Planck functions and their dB/dT.  See ck_impl_refresh_kappa.
      // (not const: problem/ck_impl_xstep may turn pass 0 into a frozen pass once
      // rt_pre_tp has measured the state; see the block after rt_pre_tp)
      bool ckfrz_ = ckimp_ && (ck_impl_pass > 0) && !ck_impl_refresh_kappa;
      // ---- problem/ck_impl_reuse_jac: does THIS pass build the tridiagonal?  With the
      // chord method only the first pass that takes a Newton step does (pass 1 when
      // pass 0 is a ck_impl_seed guess, pass 0 otherwise); every later pass leaves
      // ckjac_g and ckdb_g exactly as that pass wrote them, which removes the ~10
      // atomic_add per (cell, chain) of the assembly AND the CK_NB pair of
      // Planck-fraction look-ups per cell that dB_b/dT costs.  See
      // utils/two_stream_column_ck.hpp.
      // ck_impl_glob = ls_sub: also the pass after a sub-step restart / advance
      // problem/ck_impl_jreuse (lever 3): the chord with an adaptive rebuild -- the
      // driver sets ck_impl_jac_again when the residual stopped contracting
      // (not const: ck_impl_jreuse_xc may drop the build of a re-applying call's pass 0,
      // decided with ck_impl_xstep after rt_pre_tp)
      bool ckjacp_ = ckimp_ && ((ck_impl_reuse_jac == 0 && !(ck_impl_jreuse > 0.0))
                                      || (ck_impl_pass == CkImplJacPass())
                                      || (ck_impl_jac_again && ck_impl_pass > 0));
      // problem/ck_impl_jac_lin: the tridiagonal of a pass that builds it comes from the
      // stored factorisation of the linear kernel, never from the JAC chain kernel
      const bool ckjl_ = ckimp_ && ck_impl_jac_lin;
      auto ckjac_g = ckimp_ ? *ck_jac_ptr
                            : CkDum<DvceArray5D<Real>>("ck_jac_dummy");
      auto ckdb_g = ckimp_ ? *ck_dbdt_ptr
                           : CkDum<DvceArray5D<Real>>("ck_dbdt_dummy");
      auto cksrc_g = ckimp_ ? *ck_src_ptr : CkDum<DvceArray4D<Real>>("ck_src_dummy");
      auto ckei_g = ckimp_ ? *ck_ei_ptr : CkDum<DvceArray4D<Real>>("ck_ei_dummy");
      auto ckest_g = ckimp_ ? *ck_estar_ptr : CkDum<DvceArray4D<Real>>("ck_est_dummy");
      auto ckem_g = ckimp_ ? *ck_em_ptr : CkDum<DvceArray4D<Real>>("ck_em_dummy");
      auto ckthk_g = ckimp_ ? *ck_thk_ptr : CkDum<DvceArray4D<Real>>("ck_thk_dummy");
      auto ckthu_g = ckimp_ ? *ck_thu_ptr : CkDum<DvceArray4D<Real>>("ck_thu_dummy");
      auto ckdone_g = ckimp_ ? *ck_done_ptr : CkDum<DvceArray3D<Real>>("ck_done_dummy");
      auto cksd_g = ckimp_ ? *ck_seed_ptr : CkDum<DvceArray4D<Real>>("ck_seed_dummy");
      // ---- problem/ck_impl_frozen_op: STORE on the first pass, RE-APPLY on the rest.
      // The stored quantities are the bits pass 0 computed and the optical depths are
      // re-formed from them by the same expressions, so a frozen pass is bitwise the pass
      // it replaces.  Refused with ck_impl_refresh_kappa (which wants the opacity rebuilt
      // every pass) and with the legacy staggered layers (which ck_implicit already
      // refuses).
      // ckfop_ picks the kernel's FOP INSTANTIATION at launch (see launch_ck_cache);
      // ckfst_/ckfus_ are read inside the FOP = 1 kernel only.
      const bool ckfop_ = ckimp_ && ck_impl_frozen_op && !ck_impl_refresh_kappa
                          && !rt_layer_legacy;
      bool ckfst_ = ckfop_ && (ck_impl_pass <= 0);
      bool ckfus_ = ckfop_ && (ck_impl_pass > 0);
      const bool ckfcf_ = ck_impl_frozen_cof;
      auto ckkro_g = (ckfop_ && ck_kro_ptr != nullptr) ? *ck_kro_ptr
                   : CkDum<DvceArray5D<Real>>("ck_kro_dummy");
      auto ckc0_g = (ckfop_ && ck_c0_ptr != nullptr) ? *ck_c0_ptr
                  : CkDum<DvceArray5D<Real>>("ck_c0_dummy");
      auto ckci_g = (ckfop_ && ck_ci_ptr != nullptr) ? *ck_ci_ptr
                  : CkDum<DvceArray5D<Real>>("ck_ci_dummy");
      auto ckco_g = (ckfop_ && ck_co_ptr != nullptr) ? *ck_co_ptr
                  : CkDum<DvceArray5D<Real>>("ck_co_dummy");
      auto cktpf_g = (ckfop_ && ck_tpf_ptr != nullptr) ? *ck_tpf_ptr
                   : CkDum<DvceArray4D<Real>>("ck_tpf_dummy");
      // ck_impl_tau_min / ck_impl_colskip, as plain values for the device lambdas
      const Real cktaumin_ = ck_impl_tau_min;
      const Real ckarat_ = ck_impl_arat;
      const bool ckskip_ = ckimp_ && ck_impl_colskip;
      const bool ck_on = rt_ck;
      // the grey path shares the correlated-k scaffolding: the per-cell (T, p) and
      // opacity precompute, the cut, the tau blend and the semi-implicit application.
      // band_on says "one of the two band solvers is running", ck_on says which.
      const bool grey_on = rt_grey;
      const bool band_on = ck_on || grey_on;
      // problem/rt_top_clamp: top slot reads cell ie instead of the hydro ghost
      const bool topclamp = rt_top_clamp;
      // the optical-depth blend with the conduction module's radiative diffusion: its
      // x1-face weight w (rad_w) says how much of each face's longwave flux the
      // two-stream still owns (1 - w); the column's RT bottom is the first face with
      // w = 1, and no internal flux is injected there (the diffusion carries it)
      Conduction *pcond_rt = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->pcond
                                                            : pm->pmb_pack->phydro->pcond;
      const bool taublend = (band_on && pcond_rt != nullptr &&
                             pcond_rt->rad_tau_mode);
      auto w_g = taublend ? pcond_rt->rad_w : CkDum<DvceArray4D<Real>>("rt_w_dummy");
      auto tauf_g = taublend ? pcond_rt->rad_tauf
                             : CkDum<DvceArray4D<Real>>("rt_tau_dummy");
      if (taublend) int_at_cut = false;
      auto kc_g   = (band_on) ? *rt_kc_ptr : Fb_g;
      auto Bb_g   = (band_on) ? *rt_Bb_ptr : Fb_g;
      auto T_g    = (band_on) ? *rt_T_ptr  : tau_g;
      auto pb_g   = (band_on) ? *rt_pb_ptr : tau_g;
      auto xT_g   = (band_on) ? *rt_xT_ptr : tau_g;
      auto xP_g   = (band_on) ? *rt_xP_ptr : tau_g;
      auto icut_g = (band_on) ? *rt_icut_ptr : CkDum<DvceArray3D<int>>("dummy");
      // ---- ck-fast2 lever 1 (problem/ck_dif_dtau, see two_stream_column_ck.hpp) ----
      // icc_g is the CHAIN cut every ck chain / linear / jacobian kernel reads; the RT
      // domain cut icut_g (the apply, the column solve, the cadence) is unchanged.  Off,
      // icc_g IS icut_g and the datum array is a dummy that is never read.
      const bool ckdif_ = ck_implicit && rt_ck && (ck_dif_dtau > 0.0);
      // ck_beam_par: the storing pass's pseudo-spherical beam in its own kernels
      const bool cbt_ = ck_beam_par;
      auto icc_g = ckdif_ ? *ck_ich_ptr : icut_g;
      auto difg_g = ckdif_ ? *ck_difg_ptr : CkDum<DvceArray4D<Real>>("ck_difg_d");
      if (ckdif_ && taublend) {
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_dif_dtau is "
                  << "the route-B deep handover (<hydro|mhd>/isotropic_conduction = "
                  << "none); it cannot be combined with the tau blend." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      auto Qb_g   = (band_on) ? *rt_Qb_ptr : Fb_g;
      // see rt_cell_report: the per-face streams the report needs, q = 0 only
      const bool report_on = rt_cell_report && band_on;
      if (report_on && rt_idn_ptr == nullptr) {
        const int nmb_r = pmbp->nmb_thispack;
        rt_idn_ptr = new DvceArray4D<Real>("rt_idn", nmb_r, n3, n2, n1);
        rt_iup_ptr = new DvceArray4D<Real>("rt_iup", nmb_r, n3, n2, n1);
      }
      auto idn_g = report_on ? *rt_idn_ptr : CkDum<DvceArray4D<Real>>("rt_idn_d");
      auto iup_g = report_on ? *rt_iup_ptr : CkDum<DvceArray4D<Real>>("rt_iup_d");
      // the grey opacity: the conduction module's own table if it has one, else the
      // Freedman fit, which is what the old grey path used unconditionally
      // rt_implicit_column = 3: the exact block-tridiagonal column solve.  It runs as a
      // separate kernel right after the sweep and owns the energy update, so the apply
      // block below skips de entirely (skip_de) and keeps only its diagnostics and the
      // radiative momentum source.
      const bool mode3_ = (rt_implicit_column == 3) && grey_on && rt_split;
      // problem/rt_col3_skip_sweep: the entry sweep is dead weight under mode 3 when the
      // blend weight is 0 everywhere.  The REQUIREMENTS are checked in the problem
      // generator; the one thing it cannot check there is that w really is 0 on every
      // face (the weights are built per cycle), so that is checked here, once.
      const bool skipsweep_ = mode3_ && rt_col3_skip_sweep;
      if (skipsweep_) {
        static bool wchecked = false;
        if (!wchecked) {
          wchecked = true;
          Real wmax = 0.0;
          if (taublend) {
            Kokkos::parallel_reduce("rt_c3_wchk",
              Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*(ke-ks+1)*(je-js+1)
                                                      *(ie+2-is)),
              KOKKOS_LAMBDA(const int idx, Real &mx) {
                const int nn1 = ie + 2 - is, nn2 = je - js + 1, nn3 = ke - ks + 1;
                const int i = is + (idx % nn1);
                const int j = js + ((idx/nn1) % nn2);
                const int k = ks + ((idx/(nn1*nn2)) % nn3);
                const int m = idx/(nn1*nn2*nn3);
                const Real w = fabs(w_g(m,k,j,i));
                if (w > mx) mx = w;
              }, Kokkos::Max<Real>(wmax));
          }
          if (wmax > 0.0) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_col3_skip_sweep "
                      << "needs the tau-blend weight to be 0 on EVERY x1 face (the "
                      << "handover the skipped sweep would feed is then identically "
                      << "zero), but max|w| = " << wmax << ". Set <hydro>/rad_tau_lo "
                      << "deeper than the bottom of the box, or drop the switch."
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
        }
      }
      const bool grey_ktab = grey_on && pcond_rt != nullptr &&
                             pcond_rt->rad_kappa_tab && pcond_rt->rad_kr_nT > 0;
      const bool grey_krho = grey_ktab && pcond_rt->rad_kappa_rho;
      auto grey_kt  = grey_ktab ? pcond_rt->rad_kr_tab : CkDum<DvceArray2D<Real>>("d");
      auto grey_klT = grey_ktab ? pcond_rt->rad_kr_lT : CkDum<DvceArray1D<Real>>("d");
      auto grey_klP = grey_ktab ? pcond_rt->rad_kr_lP : CkDum<DvceArray1D<Real>>("d");
      const int grey_nT = grey_ktab ? pcond_rt->rad_kr_nT : 0;
      const int grey_nP = grey_ktab ? pcond_rt->rad_kr_nP : 0;
      const Real grey_kfac = (pcond_rt != nullptr) ? pcond_rt->rad_kappa_fac : 1.0;
      const Real grey_kabove = (pcond_rt != nullptr) ? pcond_rt->rad_kappa_above : 0.0;
      // ...and its DENSITY form (rad_gate_rho), which is the one to use when the
      // artificial medium and the star exchange gas: kappa_eff = G kappa + (1-G) kabove
      // with G = RadGate(rho).  Applied to kc_g, which is the ONE array every sweep
      // below reads -- the top-ghost dtau, the BFace emissivity weighting and the
      // up/down sweeps all go through it, so gating it here gates the whole solver.
      const Real gate_rho = (pcond_rt != nullptr) ? pcond_rt->rad_gate_rho : 0.0;
      const Real gate_dex = (pcond_rt != nullptr) ? pcond_rt->rad_gate_dex : 0.5;
      auto ckswf  = (ck_on) ? *ck_swf_ptr : CkDum<DvceArray1D<Real>>("d");
      auto cklk = (ck_on) ? *ck_lk_ptr : CkDum<DvceArray4D<Real>>("d");
      auto cklT = (ck_on) ? *ck_lT_ptr : CkDum<DvceArray1D<Real>>("d");
      auto cklP = (ck_on) ? *ck_lP_ptr : CkDum<DvceArray1D<Real>>("d");
      auto ckgw = (ck_on) ? *ck_gw_ptr : CkDum<DvceArray1D<Real>>("d");
      auto ckwl = (ck_on) ? *ck_wl_ptr : CkDum<DvceArray1D<Real>>("d");
      auto ckpf = (ck_on) ? *ck_pf_ptr : CkDum<DvceArray2D<Real>>("d");
      auto cece = (ck_on) ? *ce_ptr : CkDum<DvceArray3D<Real>>("d");
      auto celT = (ck_on) ? *ce_lT_ptr : CkDum<DvceArray1D<Real>>("d");
      auto celP = (ck_on) ? *ce_lP_ptr : CkDum<DvceArray1D<Real>>("d");
      auto cian = (ck_on) ? *cia_nT_ptr : CkDum<DvceArray1D<int>>("d");
      auto ciaT = (ck_on) ? *cia_T_ptr : CkDum<DvceArray2D<Real>>("d");
      auto ciak = (ck_on) ? *cia_k_ptr : CkDum<DvceArray3D<Real>>("d");
      auto rayx = (ck_on) ? *ray_x_ptr : CkDum<DvceArray2D<Real>>("d");
      const int ckNT = ck_nT;
      const int ckNP = ck_nP;
      const int ceNT = ce_nT;
      const int ceNP = ce_nP;
      const Real pfl0 = ck_pf_lTmin;
      const Real pfid = ck_pf_idlT;
      const Real pcut = rt_ck_pcut;
      const int ck_nq_ = ck_nq;

      if (!band_on && n1 > RT_NNC) {
        // the correlated-k and grey paths size their columns by n1 at run time; the
        // picket-fence split path below still uses the fixed RT_NNC, so it is checked
        std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/rt_split with grey "
                  << "RT needs RT_NNC >= n1, but RT_NNC = " << RT_NNC
                  << " and n1 = " << n1
                  << ". Rebuild with -DRT_NNC=" << n1 << ", or use the monolithic path "
                  << "(problem/rt_split=false), which sizes its own arrays." << std::endl;
        std::exit(EXIT_FAILURE);
      }

      // ---- A (correlated-k): the per-cell work has NO radial dependency, so it does
      // not belong in a per-column kernel. Splitting it out takes it from 8192 threads at
      // 0.52 waves/CU and 2.8 % VALUBusy -- the same starvation the monolithic kernel
      // had -- to one thread per cell. Only mu0 and the cut index are per column, and
      // with correlated-k on, the grey optical depth sweep, the grey Planck function and
      // the three-band Q_v that the old rt_pre computed are all dead: nothing reads them.
      if (band_on) {
        const Real tmu_ = rt_test_mu0;
        // roche-dhj-1031 (problem/roche): cf(4) = |grad Phi_Roche|/g_pm at the column
        // top, for the unresolved-column closure; untouched (1) otherwise
        const bool rch_ = pm->pgen->hot_jupiter_param.roche;
        const Real rgm1_ = grav*ap*ap, rom2_ = SQR(omega);
        const Real rgms_ = pm->pgen->hot_jupiter_param.roche_gms;
        const Real ra_ = pm->pgen->hot_jupiter_param.roche_a;
        // problem/ck_impl_nosync: mu0 is geometry; the later passes of a call keep it
        if (!(ck_impl_nosync && ck_impl_pass > 0)) {
        par_for("rt_pre_geom", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          // PLANE-PARALLEL: there is no substellar direction, and x2v/x3v are 1x1
          // placeholder Views on a Cartesian mesh, so they must not be read at all.
          // mu0 reaches only EffGravAt's tidal term, which is off here.
          if (pp_) {
            cf_g(m,k,j,3) = 0.0;
            return;
          }
          const Real x2v = x2v_(m,j);
          const Real x3v = x3v_(m,k);
          Real lam, phi, theta;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -theta+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            lam = x3v*iap;
            theta = -lam+M_PI/2.0;
            phi = x2v*iap;
          }
          Real mu0 = sin(theta)*cos(phi);
          if (rch_) {
            const Real rt_ = x1v_(m,ie+1);
            const Real st_ = sin(theta);
            cf_g(m,k,j,4) = roche::GeffPlanet(rgm1_, rom2_, rgms_, ra_, rt_*st_*cos(phi),
                                              rt_*st_*sin(phi), rt_*cos(theta))
                            /(rgm1_/(rt_*rt_));
          }
          if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
          if (tmu_ >= -1.0) mu0 = tmu_;
          cf_g(m,k,j,3) = mu0;
        });
        }
        // ONE ANGULAR GHOST EACH SIDE under rt_rad_force: the transverse Prad grad w term
        // differences T_g over j+-1 and k+-1, so at a MeshBlock edge it read a T_g that
        // was never filled (0 or stale).  With the temperature gate on, that put a large
        // spurious horizontal force on every block-edge cell of the shell where dw/dlnT
        // is non-zero (He4 wedge, 0.988 R: v_h 3e6 cm/s and 20 % density errors at every
        // 24th cell by 0.5 turnover, neighbouring shells quiet), on the cubed sphere too.
        const int jg_tp = (rt_rad_force && pm->multi_d) ? 1 : 0;
        const int kg_tp = (rt_rad_force && pm->three_d) ? 1 : 0;
        par_for("rt_pre_tp", DevExeSpace(), 0, nmb1, ks-kg_tp, ke+kg_tp, js-jg_tp,
                je+jg_tp, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          // THE TOP SLOT READS THE TOP ACTIVE CELL, NOT THE HYDRO GHOST.  i = ie+1 is the
          // unresolved column above the domain, and physically it IS the continuation of
          // cell ie -- taking its state from ie rather than from w0's ghost decouples the
          // whole band solver from whatever the boundary condition put there.  R9's death
          // (t = 1.99e5, hydro dt 37 s -> 7e-15 in one step, 1.64e6 non-finite cells over
          // i ~ 200..323 on every rank in 26 cycles) is what a single poisoned ghost does
          // once the down-sweep starts from it; nothing else spreads that fast.
          const int ii = (topclamp && i > ie) ? ie : i;
          Real pp, TT;
          const Real dd = rhoN(m,k,j,ii);
          PresTempFromEint(eos,gm1,Rgas,dd,eiN(m,k,j,ii),
                           TGuess(wtemp_, m, k, j, ii),pp,TT);
          if (RTBadState(pp, TT)) {          // see RTBadState: one NaN kills the column
            Kokkos::atomic_fetch_add(&stcl_g(0), 1);
            TT = 0.0;                        // the marker the opacity kernels read
            pp = 0.0;
          }
          T_g(m,k,j,i) = TT;
          pb_g(m,k,j,i) = pp*1.0e-6;
        });
        // ---- problem/ck_impl_xstep (lever 2): does THIS call re-apply the operator an
        // earlier call stored?  Decided on pass 0, on the state rt_pre_tp just measured.
        // A re-applying call runs its pass 0 exactly like a later frozen pass of the
        // storing call (no cut, no opacity, no chain sweep, no beam, no factorisation).
        if (ckfop_ && ck_impl_xstep > 0 && ck_impl_pass <= 0) {
          bool reu = (ck_impl_xs_cyc >= 0) && (ck_xsT_ptr != nullptr)
                     && ((ck_xstep_dt > 0.0)    // ck2ndt-1102 (ck_2nd.hpp (7))
                         ? !CkSchedCross(pm->time, pm->dt, ck_xstep_dt)
                         : (pm->ncycle - ck_impl_xs_cyc < ck_impl_xstep));
          const bool thr_on = (ck_impl_xstep_thr > 0.0);
          if (reu && thr_on) {
            auto t0_ = *ck_xsT_ptr;
            auto d0_ = *ck_xsD_ptr;
            Real dmx = 0.0;
            const int nk_ = ke - ks + 1, nj_ = je - js + 1, ni_ = ie + 2 - is;
            Kokkos::parallel_reduce("ck_xs_chk",
              Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nk_*nj_*ni_),
              KOKKOS_LAMBDA(const int idx, Real &mx) {
                const int i = is + (idx % ni_);
                const int j = js + ((idx/ni_) % nj_);
                const int k = ks + ((idx/(ni_*nj_)) % nk_);
                const int m = idx/(ni_*nj_*nk_);
                if (i < icut_g(m,k,j)) return;
                const Real t0 = t0_(m,k,j,i), d0 = d0_(m,k,j,i);
                const Real tt = T_g(m,k,j,i);
                const int ii = (topclamp && i > ie) ? ie : i;
                Real r = 0.0;
                if (t0 > 0.0 && tt > 0.0) r = fabs(tt/t0 - 1.0);
                else if (t0 > 0.0 || tt > 0.0) r = 1.0e30;   // a cell went bad or back
                if (d0 > 0.0) {
                  const Real rd = fabs(rhoN(m,k,j,ii)/d0 - 1.0);
                  if (rd > r) r = rd;
                }
                if (r > mx) mx = r;
              }, Kokkos::Max<Real>(dmx));
            ck_impl_xs_dmax = dmx;
            reu = (dmx <= ck_impl_xstep_thr);
          }
          // ck_impl_every_thr: a partial (guard) call always stores for its columns;
          // CkCadStep then drops the partial operator (xs_cyc = -1)
          if (ck_cad_partial) reu = false;
          ck_impl_reuse_op = reu;
          if (reu) {
            ckfrz_ = ckimp_ && !ck_impl_refresh_kappa;
            ckfst_ = false;
            ckfus_ = true;
            ++ck_impl_nreuse;
            // ck_impl_jreuse_xc: the chord Jacobian of the previous call carries over to
            // a call that re-applies the same operator (rebuilt as usual once the
            // residual stops contracting)
            if (ck_impl_jreuse > 0.0 && ck_impl_jreuse_xc && ck_impl_jac_built) {
              ckjacp_ = false;
            }
          } else {
            ck_impl_xs_cyc = pm->ncycle;
            ++ck_impl_nstore;
            CkXsSnap(pm, bdt);       // what this store read, for a restart's rebuild
            if (ck_xsT_ptr == nullptr) {
              ck_xsT_ptr = new DvceArray4D<Real>("ck_xsT", nmb1+1, n3, n2, n1);
              ck_xsD_ptr = new DvceArray4D<Real>("ck_xsD", nmb1+1, n3, n2, n1);
            }
            if (thr_on) {
              auto t0_ = *ck_xsT_ptr;
              auto d0_ = *ck_xsD_ptr;
              par_for("ck_xs_rec", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
              KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
                const int ii = (topclamp && i > ie) ? ie : i;
                t0_(m,k,j,i) = T_g(m,k,j,i);
                d0_(m,k,j,i) = rhoN(m,k,j,ii);
              });
            }
          }
        }
        // ck_sph_face = 5: the state the refresh of THIS pass will read, for a restart
        // (two_stream_ck_rst_state.hpp; same point in the pass as CkXsSnap); the VEF
        // rebuild of a restart re-takes it, so a second restart is bitwise too
        if ((ck_sph_face == 5 || ck_sph_face == 6) &&
            (ck_rst_rebuild_active ? ck_vef_force
                                                       : CkVefDue2(pm, ckfop_, ckfst_))) {
          CkVefSnap(pm, bdt);
        }
        if (ckjacp_) ck_impl_jac_built = true;
        // problem/ck_impl_frozen_op: the CUT is frozen over the Newton passes with the
        // opacity.  The stored operator covers the cells pass 0 swept, so a cut that
        // moved down by one cell at a later pass would read a kappa rho this call never
        // stored.  With the switch off the kernel runs every pass, as it always did.
        // problem/ck_impl_xstep: and over the calls that re-apply a stored operator.
        if (!(ckfop_ && (ck_impl_pass > 0 || (ck_impl_xstep > 0 && ck_impl_reuse_op)))) {
        par_for("rt_pre_cut", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          // i = is is the bottom, so pressure falls as i rises: the cut is the deepest
          // cell still shallower than pcut, i.e. the first one scanning up
          int icut = ie+1;
          if (taublend) {
            // the RT reaches down to the first face it still owns a share of
            for (int i=is; i<ie+2; ++i) {
              if (w_g(m,k,j,i) < 1.0) { icut = i; break; }
            }
          } else {
            for (int i=is; i<ie+2; ++i) {
              if (pb_g(m,k,j,i) < pcut) { icut = i; break; }
            }
          }
          icut_g(m,k,j) = icut;
        });
        }
        if (grey_on) {
          par_for("rt_pre_opac_grey", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            if (i < icut_g(m,k,j)) return;        // deeper than the cut: never read
            const int ii = (topclamp && i > ie) ? ie : i;  // top slot: see rt_pre_tp
            // T = 0 is rt_pre_tp's marker for a state the EOS could not form (see
            // RTBadState).  Such a cell takes NO part in this RT call: zero opacity and
            // zero emission make the layer transparent (e0 = 0, so alp = bet = 0, the
            // streams pass through untouched, Src and Em pick up nothing and rt_apply's
            // Newton branch is skipped).  Giving it a FLOOR state instead would hand a
            // 0.1 K blackbody to a cell with 8000 K neighbours, which collapses dt.
            if (!(T_g(m,k,j,i) > 0.0)) {
              kc_g(m,0,i,k,j) = 0.0;
              Bb_g(m,0,i,k,j) = 0.0;
              return;
            }
            const Real TT = T_g(m,k,j,i);
            const Real pcgs = pb_g(m,k,j,i)*1.0e6;
            const Real rho = rhoN(m,k,j,ii);
            Real kr = grey_ktab
                   ? RosselandTable(grey_kt, grey_klT, grey_klP, grey_nT, grey_nP, TT,
                                    grey_krho ? rho : pcgs)
                   : RosselandFreedman2014(TT, pcgs, met);
            if (gate_rho > 0.0) {
              const Real g = RadGate(rho, gate_rho, gate_dex);
              kr = g*kr + (1.0 - g)*grey_kabove;
            }
            kc_g(m,0,i,k,j) = grey_kfac*kr;
            Bb_g(m,0,i,k,j) = boltz_sigma/M_PI*SQR(SQR(TT));
          });
        } else {
        const Real kgr_ = ck_test_kgrey;
        par_for("rt_pre_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          if (i < icut_g(m,k,j)) return;          // deeper than the cut: never read
          // ck_impl_colskip: a converged column is not swept again, so nothing it feeds
          // is rebuilt either
          if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
          const int ii = (topclamp && i > ie) ? ie : i;  // top slot: see rt_pre_tp
          if (!(T_g(m,k,j,i) > 0.0)) {            // unusable state: inert, see the grey
            xT_g(m,k,j,i) = 0.0;                  // kernel's note
            xP_g(m,k,j,i) = 0.0;
            for (int b=0; b<CK_NB; ++b) {
              kc_g(m,b,i,k,j) = 0.0;
              Bb_g(m,b,i,k,j) = 0.0;
              if (ckjacp_) ckdb_g(m,b,i,k,j) = 0.0;
            }
            if (ckimp_) ckthu_g(m,k,j,i) = 0.0;
            return;
          }
          const Real TT = T_g(m,k,j,i);
          // ---- problem/ck_implicit: dB_b/dT, and the FROZEN-OPACITY pass -------------
          // The Jacobian needs d(sigma T^4/pi * f_b(T))/dT.  f_b is the tabulated Planck
          // fraction of the band, so the derivative is taken by one-sided difference on
          // the same interpolant the sweep reads -- exact consistency with B_b matters
          // more here than the O(h) truncation, which enters the Jacobian alone.
          // On a pass with the opacity frozen (every pass after the first) this rebuilds
          // the Planck functions and their derivative and NOTHING else: kc_g, xT_g and
          // xP_g are left as pass 0 wrote them, which is what "frozen opacity" means.
          // A pass that REUSES the Jacobian (ck_impl_reuse_jac) still rebuilds B_b, which
          // is the iterate itself, but not dB_b/dT, which is the matrix.
          if (ckimp_ && (ckjacp_ || ckfrz_)) {
            const Real sT4a = boltz_sigma/M_PI*SQR(SQR(TT));
            const Real TTd = TT*1.0001;
            const Real sT4d = boltz_sigma/M_PI*SQR(SQR(TTd));
            for (int b=0; b<CK_NB; ++b) {
              const Real bb = sT4a*ck_planck_frac(ckpf, pfl0, pfid, TT, b);
              if (ckjacp_) {
                const Real bd = sT4d*ck_planck_frac(ckpf, pfl0, pfid, TTd, b);
                ckdb_g(m,b,i,k,j) = (bd - bb)/(TTd - TT);
              }
              if (ckfrz_) Bb_g(m,b,i,k,j) = bb;
            }
          }
          if (!ckfrz_) {
            const Real pbar = pb_g(m,k,j,i);
            int iT, iP;
            Real fT, fP;
            ck_tp_index(cklT, ckNT, log10(TT), iT, fT);
            ck_tp_index(cklP, ckNP, log10(pbar), iP, fP);
            xT_g(m,k,j,i) = static_cast<Real>(iT) + fT;
            xP_g(m,k,j,i) = static_cast<Real>(iP) + fP;
            Real kcb[CK_NB];
            ck_continuum(cece, celT, celP, ceNT, ceNP, cian, ciaT, ciak, rayx, ckwl,
                         TT, pbar, rhoN(m,k,j,ii), kcb);
            if (kgr_ > 0.0) {
              for (int b=0; b<CK_NB; ++b) kcb[b] = kgr_;
            }
            // the density gate.  The correlated-k path never carried the rad_kappa_rmax
            // radius test -- the inert corona is a grey-path feature -- but the gate is
            // a property of the GAS, so it must reach every band here as well or a ck
            // run would keep heating the medium the grey run refuses to touch.
            const Real gk = (gate_rho > 0.0)
                ? RadGate(rhoN(m,k,j,ii), gate_rho, gate_dex) : 1.0;
            const Real sigT4_pi = boltz_sigma/M_PI*SQR(SQR(TT));
            for (int b=0; b<CK_NB; ++b) {
              kc_g(m,b,i,k,j) = (gate_rho > 0.0)
                  ? (gk*kcb[b] + (1.0 - gk)*grey_kabove) : kcb[b];
              Bb_g(m,b,i,k,j) = sigT4_pi*ck_planck_frac(ckpf, pfl0, pfid, TT, b);
            }
          }
          // ---- problem/ck_impl_tau_min: THICK OR THIN?  ----------------------------
          // The cell's own optical depth, Planck-weighted over the bands, out of the
          // SAME continuum cache the sweep builds its layers with.  The per-g-point LINE
          // opacity is not in it -- it is a function of the g index, not of the cell --
          // so this UNDER-estimates dtau and errs toward "thin".  That is the safe
          // direction: a thin cell takes the bracketed per-cell implicit step and still
          // converges to the same balance, only with its neighbours lagged; a cell
          // wrongly called thick is the one that makes the Newton diverge.
          // Zeroing dB/dT is how the cell leaves the tridiagonal: every off-diagonal
          // entry a neighbour carries for it is proportional to it (see the jck block in
          // rt_chain_ck), and its own diagonal is too.
          if (ckimp_) {
            // the A/E verdict the previous pass's apply left, AND -- if
            // ck_impl_tau_min is set -- the cell's own CONTINUUM optical depth, which
            // is all that is knowable here (the line opacity is a g-point quantity).
            bool thk = (ckthk_g(m,k,j,i) > 0.0);
            if (thk && cktaumin_ > 0.0) {
              Real kb = 0.0, bs = 0.0;
              for (int b=0; b<CK_NB; ++b) {
                kb += kc_g(m,b,i,k,j)*Bb_g(m,b,i,k,j);
                bs += Bb_g(m,b,i,k,j);
              }
              const Real kmn = (bs > 0.0) ? (kb/bs) : 0.0;
              if (kmn*rhoN(m,k,j,ii)*DX1(m,k,j,ii) < cktaumin_) thk = false;
            }
            ckthu_g(m,k,j,i) = thk ? 1.0 : 0.0;
            // with ck_impl_reuse_jac the matrix is not this pass's, so its dB/dT is not
            // rewritten either; a cell that turns thin after the build keeps the
            // identity row the solve gives it, and the neighbours' frozen off-diagonals
            // pointing at it stay non-positive, so the M-matrix property is unaffected.
            if (!thk && ckjacp_) {
              for (int b=0; b<CK_NB; ++b) ckdb_g(m,b,i,k,j) = 0.0;
            }
          }
        });
        // ---- ck-fast2 lever 1 (problem/ck_dif_dtau): the chain cut and the deep
        // diffusion operator, frozen with the opacity (built on the passes that build
        // it).  Band by band: the first face above the RT cut whose centre-to-centre
        // layer is thinner than thr in some g-point, and below it the face
        // conductances G_b and the emission weights.  Then per column the chain cut
        // ich (the lowest such face over the bands, minus one, minus the margin) and
        // per chain the flux datum g_c at face ich.  See two_stream_column_ck.hpp.
        if (ckdif_ && !ckfrz_) {
          auto difG_ = *ck_difG_ptr;
          auto difE_ = *ck_difE_ptr;
          auto diff_ = *ck_diff_ptr;
          const Real dthr_ = ck_dif_dtau;
          const int dmrg_ = ck_dif_margin;
          par_for("ck_dif_band", DevExeSpace(), 0, nmb1, 0, CK_NB-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int b, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int ic = icut_g(m,k,j);
            int fail = ie + 1;
            if (ic > ie) {
              diff_(m,b,k,j) = fail;
              return;
            }
            // per g: sum over the angles of the chains' 2 wfc mu (flux) and 2 wfc/mu
            // (emission), the expressions rt_chain_ck forms its wfc and Em with
            Real wsg[CK_NG], esg[CK_NG], krl[CK_NG];
            for (int g=0; g<CK_NG; ++g) {
              if (ck_nq_ == 1) {
                const Real mu = 1.0/CK_DIFFUSIVITY;
                const Real wf = M_PI*ckgw(g);
                wsg[g] = 2.0*wf*mu;
                esg[g] = 2.0*(wf/mu);
              } else {
                wsg[g] = 0.0;
                esg[g] = 0.0;
                for (int q=0; q<2; ++q) {
                  const Real wf = 2.0*M_PI*wg[q]*mug[q]*ckgw(g);
                  wsg[g] += 2.0*wf*mug[q];
                  esg[g] += 2.0*(wf/mug[q]);
                }
              }
            }
            // kappa rho of (b, g) at cell i, by the chain's own expression (kapof/krof)
            auto krb = [&](const int i, Real *kr) {
              const Real xTv = xT_g(m,k,j,i);
              const Real xPv = xP_g(m,k,j,i);
              const int iT = static_cast<int>(xTv);
              const int iP = static_cast<int>(xPv);
              const Real fT = xTv - static_cast<Real>(iT);
              const Real fP = xPv - static_cast<Real>(iP);
              const Real rho = rhoN(m,k,j,i);
              for (int g=0; g<CK_NG; ++g) {
                kr[g] = (ck_kappa(cklk, iT, fT, iP, fP, b, g) + kc_g(m,b,i,k,j))*rho;
              }
            };
            if (!(T_g(m,k,j,ic) > 0.0)) {
              diff_(m,b,k,j) = ic;
              return;
            }
            krb(ic, krl);
            {
              Real e = 0.0;
              for (int g=0; g<CK_NG; ++g) e += esg[g]*krl[g];
              difE_(m,b,ic,k,j) = e;
            }
            for (int f=ic+1; f<=ie; ++f) {
              if (!(T_g(m,k,j,f) > 0.0)) {
                fail = f;
                break;
              }
              Real kru[CK_NG];
              krb(f, kru);
              const Real dzl = 0.5*dx1(m,k,j,f-1);
              const Real dzu = 0.5*dx1(m,k,j,f);
              bool thick = true;
              Real gs = 0.0, e = 0.0;
              for (int g=0; g<CK_NG; ++g) {
                const Real dtc = krl[g]*dzl + kru[g]*dzu;
                if (!(dtc >= dthr_)) thick = false;
                const Real wbf = BFaceW(kru[g], krl[g], bface_on)
                               + BFaceW(krl[g], kru[g], bface_on) - 1.0;
                gs += wsg[g]*wbf/dtc;
                e += esg[g]*kru[g];
              }
              if (!thick) {
                fail = f;
                break;
              }
              difG_(m,b,f,k,j) = gs;
              difE_(m,b,f,k,j) = e;
              for (int g=0; g<CK_NG; ++g) krl[g] = kru[g];
            }
            diff_(m,b,k,j) = fail;
          });
          par_for("ck_dif_cut", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int ic = icut_g(m,k,j);
            int fmin = ie + 1;
            for (int b=0; b<CK_NB; ++b) {
              const int f = diff_(m,b,k,j);
              if (f < fmin) fmin = f;
            }
            // ich: the highest face that is thick in every chain with every face
            // below it thick too; at least one deep cell, at least one chain cell
            int ich = fmin - 1 - dmrg_;
            if (ic > ie || ich < ic + 1) ich = ic;
            if (ich > ie) ich = ie;
            icc_g(m,k,j) = ich;
          });
          par_for("ck_dif_datum", DevExeSpace(), 0, nmb1, 0, nblk*RT_NB-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int ic = icut_g(m,k,j);
            const int ich = icc_g(m,k,j);
            if (ich <= ic || ich > ie) {
              difg_g(m,c,k,j) = -1.0;
              return;
            }
            int b, g;
            Real mu;
            if (ck_nq_ == 1) {
              g = c % CK_NG;
              b = c/CK_NG;
              mu = 1.0/CK_DIFFUSIVITY;
            } else {
              g = (c/2) % CK_NG;
              b = c/(2*CK_NG);
              mu = mug[c % 2];
            }
            auto kr1 = [&](const int i) {
              const Real xTv = xT_g(m,k,j,i);
              const Real xPv = xP_g(m,k,j,i);
              const int iT = static_cast<int>(xTv);
              const int iP = static_cast<int>(xPv);
              const Real fT = xTv - static_cast<Real>(iT);
              const Real fP = xPv - static_cast<Real>(iP);
              return (ck_kappa(cklk, iT, fT, iP, fP, b, g) + kc_g(m,b,i,k,j))
                     *rhoN(m,k,j,i);
            };
            const Real kl = kr1(ich-1);
            const Real ku = kr1(ich);
            const Real dtc = kl*(0.5*dx1(m,k,j,ich-1)) + ku*(0.5*dx1(m,k,j,ich));
            const Real wbf = BFaceW(ku, kl, bface_on) + BFaceW(kl, ku, bface_on) - 1.0;
            difg_g(m,c,k,j) = 2.0*mu*wbf/dtc;
          });
          if (ck_impl_verbose) {
            Real ncol = 0.0, fsum = 0.0;
            const int nk_ = ke - ks + 1, nj_ = je - js + 1;
            Kokkos::parallel_reduce("ck_dif_stat",
              Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1+1)*nk_*nj_),
              KOKKOS_LAMBDA(const int idx, Real &a, Real &s) {
                const int j = js + (idx % nj_);
                const int k = ks + ((idx/nj_) % nk_);
                const int m = idx/(nj_*nk_);
                const int d = icc_g(m,k,j) - icut_g(m,k,j);
                if (d > 0) {
                  a += 1.0;
                  s += static_cast<Real>(d);
                }
              }, ncol, fsum);
            ck_dif_ncol = static_cast<int64_t>(ncol);
            ck_dif_fsum = fsum;
          }
        }
        }
      } else {
      // ---- A: chain-independent per-column precompute -----------------------------
      par_for("rt_pre", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        // Subviews so the body below is the monolithic kernel's code unchanged: these
        // stand in for the private tau_down_r_f[NN], B[NN], Q_v[NN] arrays.
        auto tau_down_r_f = Kokkos::subview(tau_g, m, k, j, Kokkos::ALL);
        auto B            = Kokkos::subview(B_g,   m, k, j, Kokkos::ALL);
        auto Q_v          = Kokkos::subview(Qv_g,  m, k, j, Kokkos::ALL);
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);

        Real rtop = x1v_(m,ie+1);
        Real rbot = x1v_(m,is);

        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);

        Real mus = (mu0 > 0.0) ? mu0 : 0.0;
        Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
        Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
        get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);

        // 3 V Bands
        // top
        // the top slot is the continuation of the top ACTIVE cell, not the hydro ghost
        // (see rt_pre_tp): the band solver must not be able to read a poisoned ghost
        const int itop = topclamp ? ie : (ie+1);
        Real rho = rhoN(m,k,j,itop);
        Real p = PresFromEint(eos,gm1,rho,eiN(m,k,j,itop));
        Real T = TempKelvin(eos,Rgas,rho,eiN(m,k,j,itop),p);
        bool badtop = RTBadState(p, T);      // see RTBadState: one NaN kills the column
        if (badtop) {
          Kokkos::atomic_fetch_add(&stcl_g(0), 1);
          T = 0.0;
          p = 0.0;
        }
        B[ie+1] = badtop ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr = 0.0;
        if (!badtop) get_kapr(T, p, met, kapr);
        // tau = kappa p / g for the UNRESOLVED column above the domain. g must be the
        // value at the top, not the surface value: at r/ap = 1.5 they differ by 2.3x.
        // The tidal term belongs here too: it is the EFFECTIVE gravity that sets how
        // much mass the unresolved column above the domain holds. See EffGravAt.
        Real tau_r_f = RTTopDtau(kapr, p,
                                 EffGravAt(grav, ap, rtop, grav_pmass, omega, mu0,
                                           tide));
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = (kapr*rho != 0.0) ? tau_r_f/(kapr*rho) : 0.0;
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
        // beam transmission at the face ABOVE the cell being filled, carried down the
        // sweep so that differencing the flux across a cell costs no extra exp
        Real trp1 = trans1;
        Real trp2 = trans2;
        Real trp3 = trans3;
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = rhoN(m,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,eiN(m,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          bool badcell = RTBadState(p, T);   // see RTBadState: one NaN kills the column
          if (badcell) {
            Kokkos::atomic_fetch_add(&stcl_g(0), 1);
            T = 0.0;
            p = 0.0;
          }
          B[i] = badcell ? 0.0 : boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr = 0.0;
          if (!badcell) get_kapr(T, p, met, kapr);
          Real dr = dx1(m,k,j,i);
          tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
          Real r = x1f_(m,i);
////          Real delta = (drtop+(rtop-r))/r;
//          Real delta = dr/r;
//          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
//          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
          Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
          Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
//          Real trans1 = exp(-gamv1*tausl);
//          Real trans2 = exp(-gamv2*tausl);
//          Real trans3 = exp(-gamv3*tausl);
//          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
          Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
          // Deposit the flux DIFFERENCE across the cell. The old form,
          // kappa rho F exp(-tau) with tau at the lower face, is right only for a thin
          // layer: it returns u e^-u / (1 - e^-u) of what the cell actually absorbs, with
          // u = dtau/mu, which is 0.95 at u = 0.1 but 0.58 at u = 1. On this grid, about
          // 0.46 scale heights per cell, that put dtau/mu near one wherever it mattered
          // and lost about 24 % of the incident stellar flux. Found by comparing the
          // correlated-k version of the same expression against Exo-FMS on an identical
          // column. Written this way the column integral telescopes to
          // mu F (1 - e^-tau_total) exactly, and it still reduces to the old expression
          // as dtau -> 0.
          Real Qv = (1.0-albedo)*Fstar*(1.0/3.0)
                  * ((trp1-trans1)+(trp2-trans2)+(trp3-trans3))/(fac*dr);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
          trp1 = trans1;
          trp2 = trans2;
          trp3 = trans3;
        }
        cf_g(m,k,j,0) = gamir1;
        cf_g(m,k,j,1) = gamir2;
        cf_g(m,k,j,2) = beta;
        cf_g(m,k,j,3) = mu0;
      });
      }

      // ---- B (correlated-k): one thread per (column, chain block) ------------------
      // Chain c decodes as c = ((band*CK_NG) + g)*2 + quad, so a block of RT_NB = 4
      // consecutive chains is two g-points x two angular quadrature points of ONE band,
      // which is what lets the block share a band's continuum and Planck function.
      //
      // The longwave needs no cumulative optical depth: dtau is kappa*rho*dr, purely
      // local, so going from a grey tau scaled by gamma to a per-chain kappa is a lookup
      // and nothing structural. Only the sweep's lower limit changes, from is to icut.
      if (grey_on) {
        // ---- B (grey): one thread per column.  One band, one opacity, one sweep.
        //
        // Structurally this is the correlated-k kernel with the chain loop collapsed to
        // the angular quadrature: the same layer-integrated coefficients (alpha, beta on
        // the linear-in-tau source), the same column-above-the-domain top boundary, the
        // same thermalised intensity at the cut, and the same flux accumulation.  What it
        // does NOT carry is the band structure, and that is the point: with one opacity
        // the cell's emission and its absorption of the field scale together, so the
        // local balance has a stable fixed point.
        //
        // Em is the cell's own emission per unit volume and it is written EXACTLY here,
        // 4 sigma kappa rho T^4, rather than through the hemispheric-mean quadrature the
        // band kernel uses (17 % low).  It only sets the implicit relaxation rate, but
        // the rate is what decides whether a thin cell can be integrated at the
        // hydrodynamic timestep, so it is worth having right.
        // ck-tiers: the two whole-column arrays, I_down and I_upb (2 angles each), live
        // in rt_ckscr_ptr (CkScrRows: rows 0-1 and 2-3), sized by n1 at run time.
        auto launch_grey_chain = [&]() {
          const int sg2 = je - js + 1, sg3 = ke - ks + 1;
          const size_t sgtsz = static_cast<size_t>(4)*n1*CKS_W;
          const size_t sgntl = (static_cast<size_t>(nmb1 + 1)*sg3*sg2 + CKS_W - 1)/CKS_W;
          CkScrEnsure(sgntl*sgtsz);
          Real *sgbuf = rt_ckscr_ptr->data();
          par_for("rt_chain_grey", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            const int sgt = (m*sg3 + (k - ks))*sg2 + (j - js);
            for (int i=is; i<ie+2; ++i) {
              Fb_g(m,0,i,k,j) = 0.0;
              Qb_g(m,0,i,k,j) = 0.0;
              Em_g(m,0,i,k,j) = 0.0;
              Src_g(m,0,i,k,j) = 0.0;
            }
            const int icut = icut_g(m,k,j);
            if (icut > ie) return;                  // whole column deeper than the cut
            // THE DEEP-LIMIT GRADIENT AT THE CUT.  dB/dtau with tau increasing
            // downward, from the two cell centres straddling the cut face and the
            // Rosseland tau between those centres.  Zero if the cut sits on the top
            // cell (no i+1 to difference against).
            Real dbdtau_cut = 0.0;
            if (icut + 1 <= ie) {
              const Real krb0 = kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut);
              const Real krb1 = kc_g(m,0,icut+1,k,j)*rhoN(m,k,j,icut+1);
              const Real dtc = 0.5*(krb0*DX1(m,k,j,icut) + krb1*DX1(m,k,j,icut+1));
              if (dtc > 0.0) {
                dbdtau_cut = (Bb_g(m,0,icut,k,j) - Bb_g(m,0,icut+1,k,j))/dtc;
              }
            }
            // problem/rt_bottom_flux: the cut IS the bottom wall and carries the
            // imposed internal flux, so its gradient is set by that flux, not by the
            // column's own two deepest cells.  See rt_bot_flux.
            if (bot_flux > 0.0) dbdtau_cut = 3.0*bot_flux/(4.0*M_PI);
            // ck_nquad = 1 is the hemispheric mean (mu = 1/1.66), 2 the two-point
            // Gauss-Legendre quadrature the band solver offers on the same switch
            const int nq = (ck_nq_ > 1) ? 2 : 1;
            // THE LAYER between the centres of cells iL and iL+1 (see rt_layer_legacy):
            // its two half thicknesses, the Planck source at each centre, and the source
            // at the face that separates them, taken at its own tau inside the layer.
            // BFace is applied to BOTH endpoints and symmetrically, so a radiatively
            // inert neighbour leaves a layer emitting with its own Planck function, as it
            // did before -- and with equal opacities on both sides this is just the
            // straight line from B(iL) to B(iL+1).
            auto rt_layer = [&](const int iL, Real &dt_l, Real &dt_u,
                                Real &s_l, Real &s_u, Real &s_f) {
              const Real kl = kc_g(m,0,iL,k,j)*rhoN(m,k,j,iL);
              const Real ku = kc_g(m,0,iL+1,k,j)*rhoN(m,k,j,iL+1);
              const Real bl = Bb_g(m,0,iL,k,j);
              const Real bu = Bb_g(m,0,iL+1,k,j);
              dt_l = 0.5*kl*DX1(m,k,j,iL);
              dt_u = 0.5*ku*DX1(m,k,j,iL+1);
              s_l = BFace(ku, kl, bu, bl, bface_on);
              s_u = BFace(kl, ku, bl, bu, bface_on);
              const Real dtc = dt_l + dt_u;
              s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dt_l/dtc)) : (0.5*(s_l + s_u));
            };
            Real muq[2], wfq[2];
            if (nq == 1) {
              muq[0] = 1.0/CK_DIFFUSIVITY;
              wfq[0] = M_PI;                        // F = pi I
            } else {
              for (int q=0; q<2; ++q) {
                muq[q] = mug[q];
                wfq[q] = 2.0*M_PI*wg[q]*mug[q];
              }
            }
            // GEOMETRY.  The LEGACY whole-cell layers below still carry the J = A I
            // substitution of 1159a8f3 (area-weighted intensities against a source A B)
            // and are refused on a curvilinear mesh; the CENTRE-TO-CENTRE layers -- the
            // default, and the only branch a global star runs -- carry PLAIN intensities
            // and the face mixing of BTF/MIXF instead.  Under rt_plane_parallel the two
            // are the same code: every factor is 1.0 and beta is 0.
            const Real aft_ = AFC(m,k,j,ie+1);
            const Real afc_ = AFC(m,k,j,icut);
            const Real aftn_ = AUN(m,k,j,ie+1);
            const Real afcn_ = AUN(m,k,j,icut);
            // the top face's incoming datum: area-weighted for the legacy branch, a
            // plain intensity in the face's own frame for the centre-to-centre one
            const Real atop_ = layer_legacy ? aft_ : aftn_;
            auto I_down = CkScrRows<Real>(sgbuf, 0, sgtsz, n1, sgt);
            // Top: the unresolved hydrostatic column above the domain, p/g of it, at the
            // top cell's opacity -- the same construction the band solver uses.
            {
              const Real kap = kc_g(m,0,ie+1,k,j);
              const Real mu0 = cf_g(m,k,j,3);
              // rt_top_vacuum: nothing above the domain, so the column has no optical
              // depth and hands back no intensity whatever its source function is.
              const Real dtau = top_vac ? 0.0
                  : RTTopDtau(kap, pb_g(m,k,j,ie+1)*1.0e6,
                              EffGravAt(grav, ap, X1V(m,ie+1), grav_pmass,
                                        omega, mu0, tide)*cf_g(m,k,j,4));
              // The source of that column: the ghost's own Planck function, or -- with
              // rt_top_re -- the radiative-equilibrium value I_up/2, which is what a slab
              // with nothing but space above it must emit downward.  See rt_top_re.
              //
              // I_up at the top face does not depend on I_down anywhere: the up-sweep
              // reads only the intensity it starts with at the cut and the Planck
              // function.  So the equilibrium source can be had EXACTLY, with a scalar
              // probe sweep here and no lag and no stored state -- one extra pass over
              // the column, run only when this boundary is selected.
              Real bsrc[2];
              for (int q=0; q<nq; ++q) {
                bsrc[q] = Bb_g(m,0,ie+1,k,j);
              }
              if (top_re && !layer_legacy) {
                // the same probe on the centre-to-centre layers (see rt_layer_legacy):
                // the cut half, then every layer, then the top half
                const Real dcut = 0.5*kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut)
                                * DX1(m,k,j,icut);
                const Real bcut = Bb_g(m,0,icut,k,j);
                const Real dtop = 0.5*kc_g(m,0,ie,k,j)*rhoN(m,k,j,ie)*DX1(m,k,j,ie);
                const Real btop = Bb_g(m,0,ie,k,j);
                // the probe is the up-sweep, so it dilutes exactly as the up-sweep
                // does: it carries J and is divided back by the top face's area
                // the probe is the up-sweep, so it carries what the up-sweep carries:
                // plain intensities.  It does NOT apply the face mixing -- the ratio
                // I_up(top)/2 it returns is a boundary MODEL, and rt_top_re is off in
                // every production configuration; documented as an approximation here.
                const Real acc_ = ACN(m,k,j,icut), act_ = ACN(m,k,j,ie);
                for (int q=0; q<nq; ++q) {
                  const Real bfc = bcut + dbdtau_cut*dcut;
                  Real ip = (bfc + (int_at_cut ? Iint : 0.0)
                             + muq[q]*dbdtau_cut)*afcn_;
                  Real ab, em;
                  RTLayer(dcut, muq[q], bfc*afcn_, bcut*acc_, ip, ab, em);
                  for (int i=icut; i<ie; ++i) {
                    Real dt_l, dt_u, s_l, s_u, s_f;
                    rt_layer(i, dt_l, dt_u, s_l, s_u, s_f);
                    const Real acl = ACN(m,k,j,i), afm = AUN(m,k,j,i+1);
                    const Real acu = ACN(m,k,j,i+1);
                    RTLayer(dt_l, muq[q], s_l*acl, s_f*afm, ip, ab, em);
                    RTLayer(dt_u, muq[q], s_f*afm, s_u*acu, ip, ab, em);
                  }
                  RTLayer(dtop, muq[q], btop*act_, btop*aftn_, ip, ab, em);
                  bsrc[q] = 0.5*ip/aftn_;
                }
              } else if (top_re) {
                for (int q=0; q<nq; ++q) {
                  Real ip = (Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0)
                          + muq[q]*dbdtau_cut)*afc_;
                  for (int i=icut+1; i<ie+2; ++i) {
                    const Real krb = kc_g(m,0,i-1,k,j)*rhoN(m,k,j,i-1);
                    const Real x = krb*DX1(m,k,j,i-1)/muq[q];
                    const Real e0 = -expm1(-x);
                    const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                    const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                                  : (x/2.0 - SQR(x)/3.0);
                    // emissivity-weighted far endpoint: see BFace
                    const int iir = (i > ie) ? ie : i;
                    const Real bfar = BFace(krb, kc_g(m,0,i,k,j)*rhoN(m,k,j,iir),
                                            Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
                    ip = (1.0-e0)*ip + bet*bfar*ACC(m,k,j,i)
                       + gm*Bb_g(m,0,i-1,k,j)*ACC(m,k,j,i-1);
                  }
                  bsrc[q] = 0.5*ip/aft_;
                }
              }
              for (int q=0; q<nq; ++q) {
                I_down[q][ie+1] = (1.0 - exp(-dtau/muq[q]))*(bsrc[q]*atop_);
              }
              if (report_on) idn_g(m,k,j,ie+1) = I_down[0][ie+1];
            }
            // THE LAYERS.  See rt_layer_legacy: the default branch runs layers
            // between cell CENTRES; the legacy one the old staggered whole-cell
            // layers, bit for bit.
            Real I_up[2];
            if (layer_legacy) {
              // down-sweep on the OLD staggered whole-cell layers
              for (int i=ie; i>icut-1; --i) {
                const Real krb_d = kc_g(m,0,i,k,j)*rhoN(m,k,j,i);
                const Real dtau_i = krb_d*DX1(m,k,j,i);
                // the far endpoint of the layer's source function, weighted by emitting
                // matter: above rad_kappa_rmax the neighbour has kappa = 0 and a Planck
                // function 1e9x this cell's, which used to drain it to the floor in one
                // step (see BFace).  Identical to Bb_g(i+1) at equal opacity.
                const int iip = (i+1 > ie) ? ie : i+1;
                const Real bfar_d = BFace(krb_d, kc_g(m,0,i+1,k,j)*rhoN(m,k,j,iip),
                                          Bb_g(m,0,i,k,j), Bb_g(m,0,i+1,k,j), bface_on);
                const Real bfa_d = bfar_d*ACC(m,k,j,i+1);
                const Real bwn_d = Bb_g(m,0,i,k,j)*ACC(m,k,j,i);
                for (int q=0; q<nq; ++q) {
                  const Real x = dtau_i/muq[q];
                  const Real e0 = -expm1(-x);
                  const Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                               : (x/2.0 - SQR(x)/3.0);
                  const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                  // direct source: what this stream leaves in cell i, absorbed minus
                  // emitted
                  Src_g(m,0,i,k,j) += wfq[q]/VLA(m,k,j,i)
                                    *(e0*I_down[q][i+1]
                                      - (alp*bfa_d + bet*bwn_d));
                  I_down[q][i] = (1.0-e0)*I_down[q][i+1]
                               + alp*bfa_d + bet*bwn_d;
                }
                if (report_on) idn_g(m,k,j,i) = I_down[0][i];
              }
              // Bottom of the RT domain: thermalised, plus the internal flux if the
              // layers below are not carrying it themselves (see rt_int_at_cut).
              for (int q=0; q<nq; ++q) {
                I_up[q] = (Bb_g(m,0,icut,k,j) + (int_at_cut ? Iint : 0.0)
                        + muq[q]*dbdtau_cut)*afc_;
                Fb_g(m,0,icut,k,j) += wfq[q]*(I_up[q] - I_down[q][icut])/afc_;
              }
              if (report_on) iup_g(m,k,j,icut) = I_up[0];
              // up-sweep
              for (int i=icut+1; i<ie+2; ++i) {
                const Real kap = kc_g(m,0,i-1,k,j);
                const Real rho = rhoN(m,k,j,i-1);
                const Real dtau_i = kap*rho*DX1(m,k,j,i-1);
                // same emissivity weighting for the upward stream and for Em: see BFace
                const int iiu = (i > ie) ? ie : i;
                const Real bfar_u = BFace(kap*rho, kc_g(m,0,i,k,j)*rhoN(m,k,j,iiu),
                                          Bb_g(m,0,i-1,k,j), Bb_g(m,0,i,k,j), bface_on);
                const Real bfa_u = bfar_u*ACC(m,k,j,i);
                const Real bwn_u = Bb_g(m,0,i-1,k,j)*ACC(m,k,j,i-1);
                const Real afi_ = AFC(m,k,j,i);
                for (int q=0; q<nq; ++q) {
                  const Real x = dtau_i/muq[q];
                  const Real e0 = -expm1(-x);
                  const Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
                  const Real gm  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x)
                                               : (x/2.0 - SQR(x)/3.0);
                  const Real Iup_in = I_up[q];
                  // the same for the upward stream through layer i-1
                  Src_g(m,0,i-1,k,j) += wfq[q]/VLA(m,k,j,i-1)
                                      *(e0*Iup_in
                                        - (bet*bfa_u
                                           + gm*bwn_u));
                  I_up[q] = (1.0-e0)*Iup_in
                          + bet*bfa_u + gm*bwn_u;
                  Fb_g(m,0,i,k,j) += wfq[q]*(I_up[q] - I_down[q][i])/afi_;
                }
                if (report_on) iup_g(m,k,j,i) = I_up[0];
                Em_g(m,0,i-1,k,j) = 4.0*M_PI*kap*rho
                                  * 0.5*(bfar_u + Bb_g(m,0,i-1,k,j));
              }
            } else {
              // ---- the centre-to-centre layers (see rt_layer_legacy) ----
              // The top half layer: the upper half of cell ie, entered at the top face,
              // with the source held at the cell's own centre value -- there is nothing
              // above it to interpolate towards.
              const Real dt_top = 0.5*kc_g(m,0,ie,k,j)*rhoN(m,k,j,ie)*DX1(m,k,j,ie);
              const Real b_top = Bb_g(m,0,ie,k,j);
              // The cut half layer: the lower half of cell icut.  B is continued from
              // the centre to the cut face at the deep-limit gradient, the same
              // linear-in-tau behaviour the cut boundary itself assumes.  A thick
              // half layer forgets its entry source anyway,
              // since that coefficient falls off as 1/x.
              const Real dt_cut = 0.5*kc_g(m,0,icut,k,j)*rhoN(m,k,j,icut)
                                * DX1(m,k,j,icut);
              const Real b_cut = Bb_g(m,0,icut,k,j);
              const Real b_cutf = b_cut + dbdtau_cut*dt_cut;
              Real ab, em;
              // down-sweep: centre to centre, recording the face intensity in between
              Real Idn[2];
              // THE PROBE UP-SWEEP (curvilinear only, see the SPHERICAL FORM note).
              // The down-sweep's face mixing c = beta (d_above - u_below) needs the
              // UPWARD intensity at each face, which the down-sweep does not have; the
              // up-sweep that follows it does have the downward one, so only this
              // direction is lagged.  A plane-parallel probe pass supplies u_below to
              // O(beta), leaving the down ray's mixing exact to O(beta^2) -- beta is
              // dr/2r, so ~2e-5 on the test grids, which is what the gates measure.
              // Mode 3 needs none of this: it solves the coupled system exactly.
              auto I_upb = CkScrRows<Real>(sgbuf, 2, sgtsz, n1, sgt);
              if (!pp_) {
                for (int q=0; q<nq; ++q) {
                  Real ipb = b_cutf + (int_at_cut ? Iint : 0.0) + muq[q]*dbdtau_cut;
                  I_upb[q][icut] = ipb;
                  RTLayer(dt_cut, muq[q], b_cutf, b_cut, ipb, ab, em);
                  for (int i=icut; i<ie; ++i) {
                    Real dt_l, dt_u, s_l, s_u, s_f;
                    rt_layer(i, dt_l, dt_u, s_l, s_u, s_f);
                    RTLayer(dt_l, muq[q], s_l, s_f, ipb, ab, em);
                    I_upb[q][i+1] = ipb;
                    RTLayer(dt_u, muq[q], s_f, s_u, ipb, ab, em);
                  }
                  RTLayer(dt_top, muq[q], b_top, b_top, ipb, ab, em);
                  I_upb[q][ie+1] = ipb;
                }
              }
              // the top half layer runs from the TOP FACE to the CENTRE of cell ie
              const Real act_ = ACN(m,k,j,ie), acc_ = ACN(m,k,j,icut);
              const Real bt_top = BTF(m,k,j,ie+1,icut);
              const Real bt_cut = BTF(m,k,j,icut,icut);
              for (int q=0; q<nq; ++q) {
                Idn[q] = I_down[q][ie+1];
                if (!pp_) {
                  Idn[q] += MIXF(bt_top, I_upb[q][ie+1], I_down[q][ie+1]);
                }
                RTLayer(dt_top, muq[q], b_top*aftn_, b_top*act_, Idn[q], ab, em);
                Src_g(m,0,ie,k,j) += wfq[q]/DX1(m,k,j,ie)*(ab - em);
              }
              for (int i=ie; i>icut; --i) {
                Real dt_l, dt_u, s_l, s_u, s_f;
                rt_layer(i-1, dt_l, dt_u, s_l, s_u, s_f);
                const Real acd_u = ACN(m,k,j,i), afd_m = AUN(m,k,j,i);
                const Real acd_l = ACN(m,k,j,i-1);
                const Real bt_d = BTF(m,k,j,i,icut);
                for (int q=0; q<nq; ++q) {
                  RTLayer(dt_u, muq[q], s_u*acd_u, s_f*afd_m, Idn[q], ab, em);
                  Src_g(m,0,i,k,j) += wfq[q]/DX1(m,k,j,i)*(ab - em);
                  // the value ABOVE the face, in cell i's own frame: this is what the
                  // up-sweep needs to form the mixing exactly
                  I_down[q][i] = Idn[q];
                  if (!pp_) Idn[q] += MIXF(bt_d, I_upb[q][i], Idn[q]);
                  RTLayer(dt_l, muq[q], s_f*afd_m, s_l*acd_l, Idn[q], ab, em);
                  Src_g(m,0,i-1,k,j) += wfq[q]/DX1(m,k,j,i-1)*(ab - em);
                }
                if (report_on) idn_g(m,k,j,i) = I_down[0][i];
              }
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_cut, muq[q], b_cut*acc_, b_cutf*afcn_, Idn[q], ab, em);
                Src_g(m,0,icut,k,j) += wfq[q]/DX1(m,k,j,icut)*(ab - em);
                I_down[q][icut] = Idn[q];
              }
              if (report_on) idn_g(m,k,j,icut) = I_down[0][icut];
              // Bottom of the RT domain, now AT the cut face: thermalised, plus the
              // internal flux if the layers below are not carrying it themselves (see
              // rt_int_at_cut), plus the deep-limit gradient.
              for (int q=0; q<nq; ++q) {
                I_up[q] = (b_cutf + (int_at_cut ? Iint : 0.0)
                           + muq[q]*dbdtau_cut)*afcn_;
                // the cut face: Phi = A(icut) (u_below - d_below), reported per unit
                // area, so the factor is A(icut)/A(icut) = 1 and the division stands
                const Real cc = MIXF(bt_cut, I_up[q], I_down[q][icut]);
                Fb_g(m,0,icut,k,j) += wfq[q]*(I_up[q]
                                              - (I_down[q][icut] + cc))/afcn_;
                I_up[q] += cc;                       // -> cell icut's own frame
              }
              if (report_on) iup_g(m,k,j,icut) = I_up[0];
              // up-sweep
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_cut, muq[q], b_cutf*afcn_, b_cut*acc_, I_up[q], ab, em);
                Src_g(m,0,icut,k,j) += wfq[q]/DX1(m,k,j,icut)*(ab - em);
              }
              for (int i=icut; i<ie; ++i) {
                Real dt_l, dt_u, s_l, s_u, s_f;
                rt_layer(i, dt_l, dt_u, s_l, s_u, s_f);
                const Real acu_l = ACN(m,k,j,i), afu_m = AUN(m,k,j,i+1);
                const Real acu_u = ACN(m,k,j,i+1);
                const Real bt_u = BTF(m,k,j,i+1,icut);
                // Phi_f / A(f) with Phi built on the LOWER side: A_cell(i)/A(i+1),
                // which is 1/AFC(i+1) = 1 on a Cartesian mesh
                const Real fsc_u = pp_ ? 1.0/afu_m
                                       : ACC(m,k,j,i)/AFC(m,k,j,i+1);
                for (int q=0; q<nq; ++q) {
                  RTLayer(dt_l, muq[q], s_l*acu_l, s_f*afu_m, I_up[q], ab, em);
                  Src_g(m,0,i,k,j) += wfq[q]/DX1(m,k,j,i)*(ab - em);
                  const Real cc = MIXF(bt_u, I_up[q], I_down[q][i+1]);
                  Fb_g(m,0,i+1,k,j) += wfq[q]*(I_up[q]
                                               - (I_down[q][i+1] + cc))*fsc_u;
                  I_up[q] += cc;                     // -> cell i+1's own frame
                  if (report_on && q == 0) iup_g(m,k,j,i+1) = I_up[0];
                  RTLayer(dt_u, muq[q], s_f*afu_m, s_u*acu_u, I_up[q], ab, em);
                  Src_g(m,0,i+1,k,j) += wfq[q]/DX1(m,k,j,i+1)*(ab - em);
                }
              }
              const Real fsc_t = pp_ ? 1.0/aft_ : ACC(m,k,j,ie)/AFC(m,k,j,ie+1);
              for (int q=0; q<nq; ++q) {
                RTLayer(dt_top, muq[q], b_top*act_, b_top*aftn_, I_up[q], ab, em);
                Src_g(m,0,ie,k,j) += wfq[q]/DX1(m,k,j,ie)*(ab - em);
                const Real cc = MIXF(bt_top, I_up[q], I_down[q][ie+1]);
                Fb_g(m,0,ie+1,k,j) += wfq[q]*(I_up[q]
                                              - (I_down[q][ie+1] + cc))*fsc_t;
                I_up[q] += cc;
              }
              if (report_on) iup_g(m,k,j,ie+1) = I_up[0];
              // the cell's own emission, exactly 4 sigma kappa rho T^4: each cell now
              // owns both of its half layers, so nothing here is a two-centre average
              for (int i=icut; i<ie+1; ++i) {
                Em_g(m,0,i,k,j) = 4.0*M_PI*kc_g(m,0,i,k,j)*rhoN(m,k,j,i)
                                * Bb_g(m,0,i,k,j);
              }
            }
          });
        };
        if (skipsweep_) {
          // NOTHING to launch: kc_g/Bb_g/icut_g come from the pre-kernels, Qb_g is
          // identically 0 on the grey path (only the correlated-k chain accumulates
          // into it) and stays at its allocation zero, Fb_g is written by the column
          // solve below, and Src_g/Em_g stay zero with no consumer left.
        } else {
          launch_grey_chain();
        }
        // ---- rt_implicit_column = 3: THE EXACT IMPLICIT COLUMN SOLVE ----------------
        // The sweep above has just filled kc_g/Bb_g (the FROZEN opacity and the entry
        // Planck function), Fb_g/Src_g (which the handover term is formed from) and
        // Qb_g.  Everything the block solve needs is therefore in place; it re-solves
        // the same column with the intensities as unknowns and applies the energy
        // itself.  Fb_g/Src_g/Em_g are deliberately LEFT at the entry-state values, so
        // rad_f2s, rt_rad_force and every diagnostic below see exactly the numbers mode
        // 0 hands them.
        if (mode3_) {
          const int nmb_c3 = pmbp->nmb_thispack;
          // THE PARTITIONED SOLVER (problem/rt_impl_solver = pcr).  A team per column
          // needs the spike H and the two sweep factors on top of the serial layout, and
          // a per-segment array for the reduced system.  On a host backend a Kokkos team
          // is one thread, so the partition is forced to a single segment there -- which
          // is the serial block Thomas, to round-off.
          const bool c3par = (rt_impl_solver == 1);
          const bool c3rp = c3par && rt_impl_redpar;
          // problem/rt_impl_reuse enlarges BOTH workspaces (the per-cell factors v, M
          // and the reduced solve's own stored factors); with the switch off not a word
          // is added.  The serial (thomas) path does not implement it.
          if (rt_impl_mixed > 0 && !c3par) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_impl_mixed "
                      << "requires problem/rt_impl_solver = pcr (the partitioned column "
                      << "solver); the thomas path is double only." << std::endl;
            std::exit(EXIT_FAILURE);
          }
          if (rt_impl_reuse > 0 && !c3par) {
            std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_impl_reuse "
                      << "requires problem/rt_impl_solver = pcr (the partitioned column "
                      << "solver); the thomas path does not store the factorisation."
                      << std::endl;
            std::exit(EXIT_FAILURE);
          }
          const int c3nrd = (rt_impl_reuse > 0) ? (c3rp ? RTCOL3_NRDPU : RTCOL3_NRDU)
                                                : (c3rp ? RTCOL3_NRDP : RTCOL3_NRD);
          int c3nseg = c3par ? rt_impl_nseg : 1;
          if (std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
            c3nseg = 1;
          }
          const int c3nw = c3par ? ((rt_impl_reuse > 0) ? RTCOL3_NWPU : RTCOL3_NWP)
                                 : RTCOL3_NW;
          if (rt_c3wk_ptr == nullptr) {
            // the partitioned path TRANSPOSES the workspace so that the fast index is
            // the one the threads of a team differ in; see RTCol3::Wk
            rt_c3wk_ptr = c3par
                ? new DvceArray5D<Real>("rt_c3wk", nmb_c3, c3nw, n3, n2, n1)
                : new DvceArray5D<Real>("rt_c3wk", nmb_c3, c3nw, n1, n3, n2);
            rt_c3top_ptr = new DvceArray4D<Real>("rt_c3top", nmb_c3, 2, n3, n2);
            rt_c3stat_ptr = new DvceArray1D<Real>("rt_c3stat", 26);
            // the warm-start history.  Zero-initialised, so "no history" is the state of
            // every cell on the first call and the fallback fires there by construction.
            const int wn1 = (rt_impl_warm > 0) ? n1 : 1;
            const int wn2 = (rt_impl_warm > 0) ? n2 : 1;
            const int wn3 = (rt_impl_warm > 0) ? n3 : 1;
            const int wnm = (rt_impl_warm > 0) ? nmb_c3 : 1;
            rt_c3bp_ptr = new DvceArray4D<Real>("rt_c3bp", wnm, wn3, wn2, wn1);
            rt_c3bp2_ptr = new DvceArray4D<Real>("rt_c3bp2",
                                                 (rt_impl_warm > 1) ? wnm : 1,
                                                 (rt_impl_warm > 1) ? wn3 : 1,
                                                 (rt_impl_warm > 1) ? wn2 : 1,
                                                 (rt_impl_warm > 1) ? wn1 : 1);
            // ... unless this run is a restart whose file carried the history, in which
            // case it is copied in here and the first solve of the restarted run starts
            // where the straight run's solve at this cycle started.
            RtWarmStageConsume();
            rt_c3rd_ptr = new DvceArray4D<Real>("rt_c3rd", nmb_c3, n3, n2,
                                                c3par ? c3nseg*c3nrd : 1);
            // problem/rt_impl_mixed = 2: the compact SINGLE-precision factor workspace.
            // Sized 1 when the switch is not 2, so nothing is added otherwise.
            const int c3nwf = (rt_impl_mixed == 2)
                ? ((rt_impl_reuse > 0) ? RTCOL3_NWFU : RTCOL3_NWF) : 1;
            rt_c3wkf_ptr = (rt_impl_mixed == 2)
                ? new DvceArray5D<float>("rt_c3wkf", nmb_c3, c3nwf, n3, n2, n1)
                : new DvceArray5D<float>("rt_c3wkf", 1, 1, 1, 1, 1);
            if (global_variable::my_rank == 0) {
              const double wmb = static_cast<double>(nmb_c3)*c3nw*n1*n2*n3
                                 *sizeof(Real)/1.0e6;
              const double rmb = c3par ? (static_cast<double>(nmb_c3)*c3nseg*c3nrd
                                          *n2*n3*sizeof(Real)/1.0e6) : 0.0;
              std::cout << "### two_stream_rt: rt_implicit_column = 3, the EXACT "
                        << "block-tridiagonal column solve (intensities as unknowns); "
                        << "solver " << (c3par ? "pcr" : "thomas")
                        << " nseg " << c3nseg
                        << (c3rp ? " redpar" : "")
                        << " hybrid_tau " << rt_col3_hybrid_tau
                        << (rt_impl_reuse > 0
                            ? ((rt_impl_reuse > 1) ? " reuse(always)" : " reuse") : "")
                        << (rt_impl_mixed > 0
                            ? ((rt_impl_mixed > 1) ? " mixed(2)" : " mixed(1)") : "")
                        << "; workspace " << wmb << " + " << rmb << " MB" << std::endl;
            }
          }
          auto c3wk = *rt_c3wk_ptr;
          auto c3wkf = *rt_c3wkf_ptr;
          auto c3rd = *rt_c3rd_ptr;
          auto c3top = *rt_c3top_ptr;
          auto c3stat = *rt_c3stat_ptr;
          auto c3bp = *rt_c3bp_ptr;
          auto c3bp2 = *rt_c3bp2_ptr;
          Kokkos::deep_copy(c3stat, 0.0);
          // the angular quadrature, the same one the sweep just used
          const int nq3 = (ck_nq_ > 1) ? 2 : 1;
          Real mu3[2], wf3[2];
          if (nq3 == 1) {
            mu3[0] = 1.0/CK_DIFFUSIVITY;
            wf3[0] = M_PI;
            mu3[1] = 1.0;
            wf3[1] = 0.0;
          } else {
            for (int q=0; q<2; ++q) {
              mu3[q] = mug[q];
              wf3[q] = 2.0*M_PI*wg[q]*mug[q];
            }
          }
          // the frozen top-face downward intensity: the sweep's unresolved-column model
          const Real mu3a = mu3[0], mu3b = mu3[1];
          const bool tv3 = top_vac;
          par_for("rt_c3_top", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            const Real kap = kc_g(m,0,ie+1,k,j);
            const Real mu0 = cf_g(m,k,j,3);
            const Real dtau = tv3 ? 0.0
                : RTTopDtau(kap, pb_g(m,k,j,ie+1)*1.0e6,
                            EffGravAt(grav, ap, X1V(m,ie+1), grav_pmass, omega, mu0,
                                      tide)*cf_g(m,k,j,4));
            const Real bsrc = Bb_g(m,0,ie+1,k,j);
            c3top(m,0,k,j) = (1.0 - exp(-dtau/mu3a))*bsrc;
            c3top(m,1,k,j) = (1.0 - exp(-dtau/mu3b))*bsrc;
          });
          RTCol3 c3;
          c3.u0 = u0;
          c3.bcc = bcc_uc_;
          c3.phicc = phicc_uc_;
          c3.cosc = cosc_uc_;
          c3.kc = kc_g;
          c3.Bb = Bb_g;
          c3.Fb = Fb_g;
          c3.Src = Src_g;
          c3.Qb = Qb_g;
          c3.Tg = T_g;
          c3.wblend = taublend ? w_g : CkDum<DvceArray4D<Real>>("rt_c3w_d");
          c3.icut = icut_g;
          c3.wk = c3wk;
          c3.wkf = c3wkf;
          c3.mixed = rt_impl_mixed;
          c3.rd = c3rd;
          c3.nseg = c3nseg;
          c3.redpar = c3rp;
          c3.dtop = c3top;
          c3.stat = c3stat;
          c3.size = size;
          c3.dx1 = dx1_;
          // the radial geometry of the SPHERICAL DILUTION (see the note above): the
          // column solve carries the same area-weighted intensities the sweep does
          c3.vol = volume;
          c3.area1 = area1;
          c3.eos = eos;
          c3.bdt = bdt;
          c3.sigma = boltz_sigma;
          c3.Iint = Iint;
          c3.bot_flux = bot_flux;
          c3.mu[0] = mu3[0];
          c3.mu[1] = mu3[1];
          c3.wf[0] = wf3[0];
          c3.wf[1] = wf3[1];
          c3.tol = rt_impl_tol;
          c3.norm = rt_impl_norm;
          c3.norm_eps = rt_impl_norm_eps;
          c3.exjac = rt_impl_exjac;
          c3.dstop = rt_impl_dstop;
          c3.rescheck = rt_impl_rescheck;
          c3.reuse = rt_impl_reuse;
          c3.reuse_rho = rt_impl_reuse_rho;
          c3.warm = rt_impl_warm;
          c3.bprev = c3bp;
          c3.bprev2 = c3bp2;
          // the extrapolation ratio: this call's bdt over the last one's.  The first
          // call of a run has no history at all, so the value is irrelevant there.
          c3.dtr = (rt_c3_bdt_prev > 0.0) ? (bdt/rt_c3_bdt_prev) : 1.0;
          rt_c3_bdt_prev = bdt;
          c3.dfloor = eos.dfloor;
          c3.rgas = Rgas;
          c3.gm1 = gm1;
          c3.nq = nq3;
          c3.maxit = (rt_impl_maxit > 0) ? rt_impl_maxit : 8;
          c3.is = is;
          c3.ie = ie;
          c3.is_pp = is_pp;
          c3.nx1_pp = nx1_pp;
          c3.pp = pp_;
          c3.mhd = mhd_uc_;
          c3.etg = etg_uc_;
          c3.cs = cs_uc_;
          c3.bface = bface_on;
          c3.taublend = taublend;
          c3.int_at_cut = int_at_cut;
          c3.direct = rt_src_direct;
          c3.ex_iter = rt_col3_ex_iter;
          c3.hyb_tau = rt_col3_hybrid_tau;
          c3.wrflux = skipsweep_;
          c3.dump = rt_outer_verbose && (pm->ncycle == 0);
          if (c3par) {
            RTCol3TeamLaunch(c3, nmb1, ks, ke, js, je);
          } else {
            RTCol3Launch(c3, nmb1, ks, ke, js, je);
          }
          if (rt_outer_verbose && global_variable::my_rank == 0 &&
              rt_c3_lines < 4000 &&
              (rt_report_every <= 0 || pm->ncycle % rt_report_every == 0)) {
            ++rt_c3_lines;
            auto hs = Kokkos::create_mirror_view(c3stat);
            Kokkos::deep_copy(hs, c3stat);
            const Real ncol = (hs(1) > 0.0) ? hs(1) : 1.0;
            std::cout << "### rt_col3 ncycle=" << pm->ncycle << " t=" << pm->time
                      << " dt=" << bdt << " ncol=" << static_cast<int>(hs(1))
                      << " it_mean=" << hs(0)/ncol
                      << " it_max=" << static_cast<int>(hs(2))
                      << " nclamp=" << static_cast<int>(hs(3))
                      << " max|db/b|_it=" << hs(4)
                      << " max|du/u|_tot=" << hs(18)
                      << " eos_roundtrip=" << hs(5)
                      << " nsing=" << static_cast<int>(hs(6))
                      << " nbadE=" << static_cast<int>(hs(7))
                      << " nfloor=" << static_cast<int>(hs(8))
                      << " sum_de_dx=" << hs(9)
                      << " sum|de|dx=" << hs(10)
                      << " ndiagviol=" << static_cast<int>(hs(11))
                      << " nwarmfb=" << static_cast<int>(hs(21))
                      << " nrefac=" << hs(22)/ncol
                      << " nreusepass=" << hs(23)/ncol
                      << " nmixfb=" << static_cast<int>(hs(24))
                      << " resid_max=" << hs(19)
                      << " resid_mean=" << hs(20)/ncol
                      << " budget_rel="
                      << ((hs(10) > 0.0) ? (hs(9) - hs(12))/hs(10) : 0.0)
                      << " telescope_rel="
                      << ((hs(15) > 0.0) ? (hs(13) - hs(14))/hs(15) : 0.0)
                      << " Ftop_impl/Ftop_sweep="
                      << ((hs(17) != 0.0) ? hs(16)/hs(17) : 0.0) << std::endl;
          }
        }
      } else if (ck_on) {
        // The private intensity column has to be sized at COMPILE time, but the radial
        // extent is only known at run time, and an oversized one is not free: at n1 = 68
        // the chain kernel costs 454 ms with a 72-deep column, 503 at 136 and 540 at 272,
        // because the scratch footprint grows with it. So instead of one generous ceiling
        // that taxes every run, the kernel is instantiated at a few sizes and the
        // smallest
        // one that fits is dispatched at run time.
        // The SPHERICAL form is a COMPILE-TIME parameter of this kernel, not a runtime
        // branch, for two reasons: the plane-parallel instantiation is then literally the
        // code it always was -- every ck_spherical test folds away, so nothing in its
        // expression DAG can shift -- and its per-thread face-mixing column (Cmx below)
        // shrinks to one element instead of NC x NN.  See ck_spherical.
        // problem/ck_implicit: the chain blocks accumulate into ONE tridiagonal per
        // column, so it is zeroed here, once, before any of them run.  ck-determinism:
        // the JAC chain kernel adds into per-block slots (ck_jacq_ptr; each thread zeroes
        // its own), and ck_jacq_sum adds them into ck_jac in block order after it.
        const bool ckjq_ = ckjacp_ && !ckjl_;
        if (ckjq_) {
          const auto &jx = *ck_jac_ptr;
          if (ck_jacq_ptr == nullptr
              || ck_jacq_ptr->extent(0) != jx.extent(0)
              || ck_jacq_ptr->extent(1) != static_cast<size_t>(3*nblk)) {
            if (ck_jacq_ptr != nullptr) delete ck_jacq_ptr;
            ck_jacq_ptr = new DvceArray5D<Real>("ck_jacq", jx.extent(0), 3*nblk,
                                                jx.extent(2), jx.extent(3), jx.extent(4));
          }
        }
        auto ckjq_g = ckjq_ ? *ck_jacq_ptr : CkDum<DvceArray5D<Real>>("ck_jacq_d");
        if (ckjacp_) {
          par_for("ck_jac_zero", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
          KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            ckjac_g(m,0,k,j,i) = 0.0;
            ckjac_g(m,1,k,j,i) = 0.0;
            ckjac_g(m,2,k,j,i) = 0.0;
          });
        }
        // nvcc: no extended (device) lambda may sit inside a generic lambda, so the
        // body of this launcher is the function template TsrtCkChain above
        // picket_fence_two_stream_RT_pass; ck_chain_ctx hands it, by reference, the
        // locals it reads (unpacked there in the same order).
        // ---- problem/ck_sph_face = 5: the formal solution for f = K/J (see
        // ck_sph_face).  Once per ck_vef_every cycles, BEFORE the chain kernel, on the
        // state this RT call starts from; every later pass of the cycle (the Newton
        // passes) reads the same factors.  Nodes q = 0 .. n+1 of a column: the cut face,
        // the frames of cells icut .. ie (radius sqrt(A/dOmega)), the top face.  Rays:
        // ck_vef_ncore through the wall, then one tangent at every node.
        // ck_vef_sync: on a storing pass only (a re-applying pass keeps the factors its
        // stored operator was formed with); the first call always forms them
        // a restart's rebuild passes (two_stream_ck_rst_state.hpp) form them only when
        // they rebuild the refresh itself (ck_vef_force) or nothing is there yet
        const bool vdue_ = ck_rst_rebuild_active ? ck_vef_force
                                                 : CkVefDue2(pm, ckfop_, ckfst_);
        // (ck_sph_face = 6 reads the same formal solution)
        // ck2nd-1102 ck_vef_tint: does this pass re-form the factors
        const bool vref_ = (ck_sph_face == 5 || ck_sph_face == 6) &&
                           (ck_vef_s_ptr == nullptr || vdue_);
        if ((ck_sph_face == 5 || ck_sph_face == 6) &&
            (ck_vef_s_ptr == nullptr || vdue_)) {
          // a refresh forced on a re-applying pass (nothing formed yet, e.g. the first
          // call after a restart without VEF state re-applies an xstep operator) is not
          // a scheduled one: the next storing pass refreshes (and snapshots) properly
          if (!ck_rst_rebuild_active && vdue_) {
            ck_vef_last = static_cast<int>(pm->ncycle);
            ck_vef_tlast = pm->time + pm->dt;   // ck2nd-1102 ck_vef_dt
          }
          const int nv = n1 + 3;
          const int nc_ = CK_NB*CK_NG;
          if (ck_vef_s_ptr == nullptr ||
              ck_vef_s_ptr->extent(0) != static_cast<size_t>(nmb1 + 1)) {
            if (ck_vef_s_ptr != nullptr) {
              delete ck_vef_s_ptr;
              delete ck_vef_g_ptr;
              delete ck_vef_w_ptr;
            }
            ck_vef_s_ptr = new DvceArray5D<Real>("ck_vef_s", nmb1 + 1, nc_, nv, ke + 1,
                                                 je + 1);
            ck_vef_g_ptr = new DvceArray5D<Real>("ck_vef_g", nmb1 + 1, nc_, nv, ke + 1,
                                                 je + 1);
            ck_vef_w_ptr = new DvceArray5D<Real>("ck_vef_w", nmb1 + 1, nc_, 14*nv,
                                                 ke + 1, je + 1);
          }
          auto vs_ = *ck_vef_s_ptr;
          auto vg_ = *ck_vef_g_ptr;
          auto vw_ = *ck_vef_w_ptr;
          // problem/ck_vef_report (diagnostic): keep the factors being replaced
          const bool vrep_ = (ck_vef_report > 0 && ck_vef_last >= 0 && vdue_);
          DvceArray5D<Real> vso_, vgo_;
          if (vrep_) {
            vso_ = DvceArray5D<Real>("ck_vef_so", vs_.extent(0), vs_.extent(1),
                                     vs_.extent(2), vs_.extent(3), vs_.extent(4));
            vgo_ = DvceArray5D<Real>("ck_vef_go", vg_.extent(0), vg_.extent(1),
                                     vg_.extent(2), vg_.extent(3), vg_.extent(4));
            Kokkos::deep_copy(vso_, vs_);
            Kokkos::deep_copy(vgo_, vg_);
          }
          const int ncore = ck_vef_ncore;
          const Real xlim_ = ck_vef_xlim;
          par_for("ck_vef_formal", DevExeSpace(), 0, nmb1, 0, nc_ - 1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            const int b = c/CK_NG;
            const int g = c % CK_NG;
            const int n = ie - icut + 1;
            const int nn = n + 2;
            // work slots: 0 r, 1 kappa rho, 2 B, 3 I_in, 4 J, 5 K, 6 mu_prev,
            // 7 (I+ + I-)_prev, 8 W0, 9 W2, 10 W1 (top), 11 H (top), 12 J^- and
            // 13 H^- at the wall (the incoming half moments; slot 3 of the previous
            // ray is kept in slot 12's partner below)
            auto W = [&](const int s, const int q) -> Real & {
              return vw_(m,c,s*nv + q,k,j);
            };
            for (int q=0; q<nn; ++q) {
              Real r;
              if (q == 0) {
                r = X1F(m,icut);
              } else if (q == nn - 1) {
                r = X1F(m,ie+1);
              } else {
                const Real rm = X1F(m,icut+q-1), rp = X1F(m,icut+q);
                r = sqrt((rp*rp + rp*rm + rm*rm)/3.0);
                const int i = icut + q - 1;
                const Real xTv = xT_g(m,k,j,i);
                const Real xPv = xP_g(m,k,j,i);
                const int iT = static_cast<int>(xTv);
                const int iP = static_cast<int>(xPv);
                W(1,q) = (ck_kappa(cklk, iT, xTv - static_cast<Real>(iT), iP,
                                   xPv - static_cast<Real>(iP), b, g)
                          + kc_g(m,b,i,k,j))*rhoN(m,k,j,i);
                W(2,q) = Bb_g(m,b,i,k,j);
              }
              W(0,q) = r;
              for (int s=3; s<14; ++s) W(s,q) = 0.0;
            }
            const Real Ib = Bb_g(m,b,icut,k,j)
                          + (int_at_cut ? boltz_sigma/M_PI*Tint4
                             * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
            const Real r0 = W(0,0);
            // chord of impact parameter p between radii ra < rb (both >= p)
            auto chord = [&](const Real p, const Real ra, const Real rb) {
              const Real za = (ra > p) ? sqrt(ra*ra - p*p) : 0.0;
              const Real zb = (rb > p) ? sqrt(rb*rb - p*p) : 0.0;
              return zb - za;
            };
            // one cell's layer of the segment between nodes q and q+1, crossed
            // inward (in = true) or outward, applied to I
            auto seg = [&](const int q, const Real p, Real &I, const bool in) {
              int cl[2];
              Real ds[2];
              int ns = 1;
              if (q == 0) {
                cl[0] = 1;
                ds[0] = chord(p, W(0,0), W(0,1));
              } else if (q == nn - 2) {
                cl[0] = q;
                ds[0] = chord(p, W(0,q), W(0,q+1));
              } else {
                const Real rf = X1F(m,icut+q);
                cl[0] = q;
                ds[0] = chord(p, W(0,q), rf);
                cl[1] = q + 1;
                ds[1] = chord(p, rf, W(0,q+1));
                ns = 2;
              }
              for (int t=0; t<ns; ++t) {
                const int tt = in ? (ns - 1 - t) : t;
                const Real x = W(1,cl[tt])*ds[tt];
                const Real e = exp(-x);
                I = I*e + W(2,cl[tt])*(1.0 - e);
              }
            };
            // l = ncore is the wall's edge ray p = r_cut, taken as hitting the wall;
            // the tangent ray of node 0 is the same p missing it, so the step of
            // the field at the edge of the wall's cone falls between two rays
            const int nray = ncore + 1 + nn;
            for (int l=0; l<nray; ++l) {
              Real p;
              int qlo;
              if (l <= ncore) {
                p = r0*sin(0.5*M_PI*static_cast<Real>(l)/static_cast<Real>(ncore));
                qlo = 0;
              } else {
                qlo = l - ncore - 1;
                p = W(0,qlo);
              }
              // inward from the top (vacuum)
              Real I = 0.0;
              for (int q=nn-1; q>=qlo; --q) {
                W(3,q) = I;
                if (q > qlo) seg(q-1, p, I, true);
              }
              I = (l <= ncore) ? Ib : W(3,qlo);
              for (int q=qlo; q<nn; ++q) {
                const Real r = W(0,q);
                const Real mu2 = (r > p) ? (1.0 - (p/r)*(p/r)) : 0.0;
                const Real mu = (l <= ncore || q > qlo) ? sqrt((mu2 > 0.0) ? mu2 : 0.0)
                                                        : 0.0;
                const Real gs = I + W(3,q);
                if (l > 0) {
                  const Real dm = W(6,q) - mu;
                  const Real mp = W(6,q), gp = W(7,q);
                  W(4,q) += 0.5*dm*(gp + gs);
                  W(5,q) += 0.5*dm*(mp*mp*gp + mu*mu*gs);
                  W(8,q) += dm;
                  W(9,q) += 0.5*dm*(mp*mp + mu*mu);
                  if (q == nn - 1) {
                    W(10,q) += 0.5*dm*(mp + mu);
                    W(11,q) += 0.5*dm*(mp*(gp - 2.0*W(3,q)) + mu*(gs - 2.0*W(3,q)));
                  }
                  if (q == 0) {
                    // I^- of the previous ray at the wall is parked in W(2,0)
                    // (node 0 has no cell of its own)
                    const Real ip = W(2,0), ic = W(3,0);
                    W(12,0) += 0.5*dm*(ip + ic);
                    W(13,0) += 0.5*dm*(mp*ip + mu*ic);
                  }
                }
                W(6,q) = mu;
                W(7,q) = gs;
                if (q == 0) W(2,0) = W(3,0);
                if (q < nn - 1) seg(q, p, I, false);
              }
            }
            // f, s = sqrt(3 f), ln(Phi f) and the face factors
            Real lpf = 0.0, fprev = 1.0/3.0, rprev = r0, lphi = 0.0;
            for (int q=0; q<nn; ++q) {
              const Real Jq = W(4,q), Kq = W(5,q);
              Real f = 1.0/3.0;
              if (Jq > 0.0 && W(9,q) > 0.0) {
                f = (Kq/Jq)*(W(8,q)/(3.0*W(9,q)));
                f = (f < 1.0e-3) ? 1.0e-3 : ((f > 1.0) ? 1.0 : f);
              }
              const Real r = W(0,q);
              if (q > 0) {
                lphi += 0.5*((3.0*fprev - 1.0)/fprev + (3.0*f - 1.0)/f)*log(r/rprev);
              }
              const Real lpn = lphi + log(f);
              if (q == 0) {
                // the wall's Marshak factor beta = J^-/H^- (2 if nothing comes in)
                const Real hm = W(13,0), jm = W(12,0);
                vg_(m,c,icut-1,k,j) = (hm > 1.0e-12*Jq && jm > 0.0) ? jm/hm : 2.0;
              }
              vs_(m,c,icut+q,k,j) = sqrt(3.0*f);
              if (q > 0) vg_(m,c,icut+q-1,k,j) = exp(lpn - lpf);
              lpf = lpn;
              fprev = f;
              rprev = r;
              if (q == nn - 1) {
                // Auer's outer condition H = h J, h normalised to 1/2 for a
                // hemispherically isotropic field; zeta = (s - 2h)/(s + 2h)
                Real h = 0.5;
                if (Jq > 0.0 && W(10,q) > 0.0) {
                  h = (W(11,q)/Jq)*(W(8,q)/(2.0*W(10,q)));
                }
                const Real st = sqrt(f);          // mu_eff of the pair
                vg_(m,c,ie+2,k,j) = (st - h)/(st + h);
              }
            }
          });
          // problem/ck_vef_xlim = X > 1: the face factor x = gamma rho' of an interior
          // face (q = 2..nn-2) kept in [1/X, X] by limiting gamma (rho', the flux
          // continuity, is untouched, so the deposit still telescopes).  A separate
          // kernel, so that ck_vef_formal (and every run with the key off) compiles to
          // the same code as without the limiter: it re-forms f from the moments the
          // formal solution left in the work slots (4 J, 5 K, 8 W0, 9 W2; slot 0 r)
          // with the same expressions, so the limited gamma is the in-kernel one.
          if (xlim_ > 1.0) {
            par_for("ck_vef_xlim", DevExeSpace(), 0, nmb1, 0, nc_ - 1, ks, ke, js, je,
            KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
              const int icut = icc_g(m,k,j);
              if (icut > ie) return;
              const int nn = ie - icut + 3;
              auto W = [&](const int s, const int q) -> Real {
                return vw_(m,c,s*nv + q,k,j);
              };
              auto fq = [&](const int q) -> Real {
                const Real Jq = W(4,q), Kq = W(5,q);
                Real f = 1.0/3.0;
                if (Jq > 0.0 && W(9,q) > 0.0) {
                  f = (Kq/Jq)*(W(8,q)/(3.0*W(9,q)));
                  f = (f < 1.0e-3) ? 1.0e-3 : ((f > 1.0) ? 1.0 : f);
                }
                return f;
              };
              Real fprev = fq(1);
              for (int q=2; q<=nn-2; ++q) {
                const Real f = fq(q);
                Real gq = vg_(m,c,icut+q-1,k,j);
                const Real ra = W(0,q), rb = W(0,q-1);
                const Real rhop = (rb*rb)/(ra*ra)*sqrt(fprev/f);
                const Real xq = gq*rhop;
                if (xq > xlim_) gq = xlim_/rhop;
                if (xq < 1.0/xlim_) gq = 1.0/(xlim_*rhop);
                vg_(m,c,icut+q-1,k,j) = gq;
                fprev = f;
              }
            });
          }
          if (vrep_) CkVefReport(pm, vs_, vg_, vso_, vgo_, *rt_pb_ptr, icc_g, X1F, nmb1,
                                 nc_, ks, ke, js, je, ie);
        }
        // ck2nd-1102 problem/ck_vef_tint (utils/ck_2nd.hpp (4)); off: nothing
        if (ck_vef_tint && (ck_sph_face == 5 || ck_sph_face == 6) &&
            !ck_rst_rebuild_active) {
          CkVefTint(pm, vref_, ckfst_);
        }
        // problem/ck_sph_face; = 6 runs the chain kernel (a storing pass without
        // ck_beam_par, ck_impl_lin_check) on the face-5 map: its own operator lives in
        // the linear kernels (rt_chain_ck_lin6 overwrites Src, Fb and Em)
        const int dil_ = (ck_sph_face == 6) ? 5 : ck_sph_face;
        const int tpm_ = ck_sph_top;    // problem/ck_sph_top
        const bool vefa_ = (ck_sph_face == 5 || ck_sph_face == 6);
        auto vfs_ = vefa_ ? *ck_vef_s_ptr : CkDum<DvceArray5D<Real>>("ck_vef_s_d");
        auto vfg_ = vefa_ ? *ck_vef_g_ptr : CkDum<DvceArray5D<Real>>("ck_vef_g_d");
        // the linear re-apply kernels' names for the same (ck_sph_face = 5); vef6_ the
        // conservative M-matrix form of ck_sph_face = 6 (its own build, linear and
        // Jacobian kernels)
        const bool vefl_ = vefa_;
        const bool vef6_ = (ck_sph_face == 6);
        const bool wex_ = ck_wall_flux_exact;   // see ck_wall_flux_exact
        auto vfs_g = vfs_;
        auto vfg_g = vfg_;
        auto ck_chain_ctx = std::forward_as_tuple(ACC, AFC, albedo, ap, Bb_g, bface_on,
                                                  boltz_sigma, BTF, cbt_, cf_g, ckc0_g,
                                                  ckci_g, ckco_g, ckdb_g, ckdif_,
                                                  ckdone_g, ckfcf_, ckfst_, ckfus_, ckgw,
                                                  ckjq_g, ckjacp_, ckkro_g, cklk, ck_nq_,
                                                  ckpf, ckskip_, ckswf, cktpf_g, difg_g,
                                                  dx1, Em_g, Fb_g, Fstar, grav,
                                                  grav_pmass, icc_g, icut_g, ie,
                                                  int_at_cut, is, je, js, kc_g, ke, ks,
                                                  layer_legacy, MIXF, mug, n1, nblk, nmb1,
                                                  omega, pb_g, pfid, pfl0, Qb_g, rhoN,
                                                  Src_g, T_g, tide, Tint, Tint4, wg, X1F,
                                                  x1v_, xP_g, xT_g, dil_, tpm_,
                                                  vfs_, vfg_, wex_);
        // tusplit-1102: TsrtCkChain is instantiated for this type (two_stream_rt.hpp)
        static_assert(std::is_same_v<decltype(ck_chain_ctx), TsrtCkChainCtx>,
                      "ck_chain_ctx != TsrtCkChainCtx: update two_stream_rt.hpp");
        auto launch_ck_chain = [&](auto nn_tag, auto sph_tag, auto bsp_tag,
                                   auto cch_tag, auto frm_tag, auto fop_tag,
                                   auto jac_tag) {
          TsrtCkChain(ck_chain_ctx, nn_tag, sph_tag, bsp_tag, cch_tag, frm_tag, fop_tag,
                      jac_tag);
        };
        auto launch_ck_cch = [&](auto nn_tag, auto sph_tag, auto bsp_tag,
                                 auto frm_tag, auto fop_tag) {
          if (ckcache_ >= 2) {       // see ck_sweep_cache
            launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                            std::integral_constant<int, 2>{}, frm_tag, fop_tag,
                            std::false_type{});
          } else if (ckcache_ == 1) {
            launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                            std::integral_constant<int, 1>{}, frm_tag, fop_tag,
                            std::false_type{});
          } else {
            launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                            std::integral_constant<int, 0>{}, frm_tag, fop_tag,
                            std::false_type{});
          }
        };
        // problem/ck_impl_frozen_op: the second compile-time tag on this kernel.  It is
        // dispatched here rather than at the top so that it multiplies the cached forms
        // only: the probe-free ck_sweep_form recurrences are refused under ck_implicit,
        // which is the only way the frozen operator can be on.
        auto launch_ck_cache = [&](auto nn_tag, auto sph_tag, auto bsp_tag,
                                   auto frm_tag) {
          if (ckfop_) {
            launch_ck_cch(nn_tag, sph_tag, bsp_tag, frm_tag, std::true_type{});
          } else {
            launch_ck_cch(nn_tag, sph_tag, bsp_tag, frm_tag, std::false_type{});
          }
        };
        // problem/ck_sweep_form: the probe-free forms exist only where there IS a face
        // coupling to solve, so they are instantiated under the SPHERICAL tag alone --
        // and only at ck_sweep_cache = 2, which is the default wherever they can run.
        // That keeps the instantiation count (and the compile time of this kernel) down:
        // 4 of the 28 per radial tier, instead of the 8 a full product would add.  The
        // startup guard refuses the other cases rather than silently changing them.
        // (The other 24 are 2 SPH x 2 BSP x 3 CCH x 2 FOP.)
        auto launch_ck_form = [&](auto nn_tag, auto sph_tag, auto bsp_tag) {
          constexpr bool SPHF = decltype(sph_tag)::value;
          if constexpr (SPHF) {
            if (ckform_ == 1) {
              // tm: the JAC instantiation only on a pass that assembles the Jacobian.
              // Its per-face window (3 Reals per face and chain) no longer sits on the
              // GPU stack (ck-scratch), so it exists at every n1.
              // problem/ck_impl_frozen_op (phase T2): the FOP instantiations, JAC or not,
              // exist only under ck_implicit; the production kernel is FOP = JAC = 0.
              if (ckjacp_ && !ckjl_) {
                if (ckfop_) {
                  launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                                  std::integral_constant<int, 2>{},
                                  std::integral_constant<int, 1>{}, std::true_type{},
                                  std::true_type{});
                } else {
                  launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                                  std::integral_constant<int, 2>{},
                                  std::integral_constant<int, 1>{}, std::false_type{},
                                  std::true_type{});
                }
              } else if (ckfop_) {
                launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                                std::integral_constant<int, 2>{},
                                std::integral_constant<int, 1>{}, std::true_type{},
                                std::false_type{});
              } else {
                launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                                std::integral_constant<int, 2>{},
                                std::integral_constant<int, 1>{}, std::false_type{},
                                std::false_type{});
              }
              return;
            } else if (ckform_ == 2) {
              launch_ck_chain(nn_tag, sph_tag, bsp_tag,
                              std::integral_constant<int, 2>{},
                              std::integral_constant<int, 2>{}, std::false_type{},
                              std::false_type{});
              return;
            }
          }
          launch_ck_cache(nn_tag, sph_tag, bsp_tag,
                          std::integral_constant<int, 0>{});
        };
        // ck-scratch: the chain kernel no longer has a column tier -- its column arrays
        // are sized by n1 at run time (see CkScrCol) -- so it is instantiated once, with
        // the tag 0, for every n1.
        auto launch_ck_tier = [&](auto sph_tag, auto bsp_tag) {
          launch_ck_form(std::integral_constant<int, 0>{}, sph_tag, bsp_tag);
        };
        // =============== problem/ck_impl_lin: THE LINEAR RE-APPLY OF THE tm SWEEP =====
        // DESIGN_tm.md phase T3.  With the opacity, the cut and the beam frozen for the
        // call (ck_impl_frozen_op), the tm sweep is affine in B_b and its reflectance R
        // does not depend on B: pass 1's R recurrence is the factorisation of the column.
        // ck_lin_build packs, on the storing pass, everything the recurrences need that
        // does not depend on B into ONE View (lP, 6 slots per chain since ck-store:
        // R and 1/(1 + R beta) at every face, the three BFace/interpolation weights of
        // every layer and the emission coefficient; the half-layer triple is read from
        // ck_c0/ci/co) and the column geometry
        // into another (lG: beta, the area ratios, 1/dz); rt_chain_ck_lin then re-runs
        // only the Sc forward substitution and the ray back substitution on the new B_b.
        // Same recurrences and accumulation order as the tm body of rt_chain_ck; the
        // divides by (1 + R beta) and the BFace sums become multiplies by stored
        // reciprocals and weights, so the answer is the chain kernel's to round-off.  It
        // is a separate kernel with few captured Views (the chain kernel's register
        // file is full): the production chain kernel is not touched by it.  Qb_g is left
        // as the storing pass wrote it.
        const bool cklon_ = ckfop_ && ck_impl_lin && (ckform_ == 1);
        const bool cklbuild_ = cklon_ && ckfst_;
        // under ck_impl_jac_lin the pass that builds the Jacobian is a linear pass too
        const bool cklin_ = cklon_ && ckfus_ && (!ckjacp_ || ckjl_);
        // ck_impl_jac_lin: this pass builds the Jacobian from the factorisation (a pass
        // that builds it, once the factorisation exists; not a ck_impl_seed pass, whose
        // rows are identity rows)
        const bool ckjlp_ = ckjl_ && cklon_ && ckjacp_
                             && !(ck_impl_seed > 0 && ck_impl_pass == 0 && !ck_impl_jac0);
        const bool cklw_ = ck_impl_lw;
        const bool cklchk_ = cklin_ && (ck_impl_lin_check > 0);
        // ck-store: THE STORING PASS WITHOUT THE CHAIN KERNEL.  The storing pass of the
        // linear re-apply needs from rt_chain_ck only what depends on the opacity alone
        // -- kappa rho, the half-layer triple and the top factor per (cell, chain) --
        // and the fluxes of the first pass.  ck_coef stores the former, by the chain
        // kernel's own expressions (bitwise the same numbers), a thread per (band, g)
        // pair and no serial recurrence; ck_lin_build then factorises as before, and
        // the linear kernels of every later pass give the fluxes of this one too, i.e.
        // the storing pass's Src/Fb/Em change by round-off (the lin kernels are the
        // chain kernel to round-off, see ck_lin_check); Em is bitwise.  Only where the
        // chain kernel's storing pass carries nothing else: the beam in its own kernels
        // (ck_beam_par), no four-pass Jacobian, no deep diffusion, FP64.
        const bool ckstsp_ = cklbuild_ && ckfcf_ && !ckfus_ && cbt_ && ckbsph_ && cksph_
                             && (ckcache_ >= 2) && !(ckjacp_ && !ckjl_) && !ckdif_
                             && (RT_FP32 == 0) && ck_store_split;
        // ck-lin2: the half-layer triple (e0, cin, cout) is formed by the paired linear
        // kernel, ck_lin_build and rt_chain_ck_jlin (CkE0, CkCinCout), so ck_coef
        // stores it only for a reader that still reads it: the
        // chain kernel re-applying a frozen operator (a Jacobian pass without
        // ck_impl_jac_lin, ck_impl_lin_check) and the unpaired linear tiers
        const bool cklp_ = (ck_impl_lin_thr == 1) && (ck_nq_ == 2);
        const bool ckcst_ = !(cklp_ && ckjl_ && (ck_impl_lin_check == 0));
        const int nch_ = nblk*RT_NB;
        auto lP_g = cklon_ ? *ck_linP_ptr : CkDum<DvceArray5D<Real>>("ck_lP_d");
        auto lG_g = cklon_ ? *ck_linG_ptr : CkDum<DvceArray5D<Real>>("ck_lG_d");
        auto lC_g = cklon_ ? *ck_linC_ptr : CkDum<DvceArray2D<Real>>("ck_lC_d");
        // ck-tiers: Sc, the one whole-column array of ck_lin / ck_lin1, lives in
        // rt_ckscr_ptr (CkScrRows), sized by the actual n1: RT_NB rows per ck_lin thread,
        // one per ck_lin1 thread; the thread count is the same (nch_ per column), so
        // both need nch_*n1 Reals per column, 1/7 of what the chain kernel's groups take.
        const int scl2 = je - js + 1, scl3 = ke - ks + 1;
        auto launch_ck_lin = [&]() {
          const size_t sltsz = static_cast<size_t>(RT_NB)*n1*CKS_W;
          const size_t slntl = (static_cast<size_t>(nmb1 + 1)*nblk*scl3*scl2 + CKS_W - 1)
                               /CKS_W;
          CkScrEnsure(slntl*sltsz);
          Real *slbuf = rt_ckscr_ptr->data();
          par_for("rt_chain_ck_lin", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
            constexpr int NC = RT_NB;
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            for (int i=is; i<ie+2; ++i) {
              Fb_g(m,blk,i,k,j) = 0.0;
              Em_g(m,blk,i,k,j) = 0.0;
              Src_g(m,blk,i,k,j) = 0.0;
            }
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            int bandc[NC];
            for (int cc=0; cc<NC; ++cc) {
              const int c = blk*NC + cc;
              bandc[cc] = (ck_nq_ == 1) ? (c/CK_NG) : (c/(2*CK_NG));
            }
            // PASS 1: Sc up the column on the stored R; Sc on the below side of every
            // face is kept for pass 2
            const int slt = ((m*nblk + blk)*scl3 + (k - ks))*scl2 + (j - js);
            auto Sc = CkScrRows<Real>(slbuf, 0, sltsz, n1, slt);
            Real ss[NC], sfc[NC], suc[NC];
            for (int cc=0; cc<NC; ++cc) {
              const Real bcut = Bb_g(m,bandc[cc],icut,k,j);
              sfc[cc] = bcut;
              suc[cc] = bcut;
              ss[cc] = bcut + lC_g(1,blk*NC+cc);
              if (ckdif_) {                  // ck_dif_dtau: the flux datum
                const Real gdf = difg_g(m,blk*NC+cc,k,j);
                if (gdf >= 0.0) ss[cc] = gdf*(Bb_g(m,bandc[cc],icut-1,k,j) - bcut);
              }
            }
            for (int i=icut; i<ie+1; ++i) {
              const Real bt = lG_g(m,0,i,k,j);
              for (int cc=0; cc<NC; ++cc) {
                const int c = blk*NC + cc;
                const int cw = (ck_nq_ == 2) ? (c & ~1) : c;   // ck-nq2: pair weights
                const int b = bandc[cc];
                const Real bown = Bb_g(m,b,i,k,j);
                Real slv = bown, sfv = bown, suu = bown;
                if (i < ie) {
                  const Real bfar = Bb_g(m,b,i+1,k,j);
                  const Real wl = lP_g(m,2*nch_+cw,i,k,j);
                  const Real wu = lP_g(m,3*nch_+cw,i,k,j);
                  slv = wl*bown + (1.0 - wl)*bfar;
                  suu = wu*bfar + (1.0 - wu)*bown;
                  sfv = slv + (suu - slv)*lP_g(m,4*nch_+cw,i,k,j);
                }
                const Real e0 = ckc0_g(m,c,i,k,j);
                const Real cin = ckci_g(m,c,i,k,j);
                const Real cout = ckco_g(m,c,i,k,j);
                const Real tr = 1.0 - e0;
                const Real pl = cin*sfc[cc] + cout*suc[cc];
                const Real ql = cin*suc[cc] + cout*sfc[cc];
                const Real pu = cin*slv + cout*sfv;
                const Real qu = cin*sfv + cout*slv;
                Sc[cc][i] = ss[cc];
                const Real idn = CkIdn(lP_g(m,0*nch_+c,i,k,j), lG_g(m,0,i,k,j));
                const Real rn = (lP_g(m,0*nch_+c,i,k,j) + bt)*idn;
                Real sv = (1.0 - bt)*ss[cc]*idn;
                sv = tr*(rn*ql + sv) + pl;
                const Real r2 = tr*tr*rn;
                ss[cc] = tr*(r2*qu + sv) + pu;
                sfc[cc] = sfv;
                suc[cc] = suu;
              }
            }
            for (int cc=0; cc<NC; ++cc) Sc[cc][ie+1] = ss[cc];
            // PASS 2: the down ray, the face solve on the stored R, the up ray
            Real dcu[NC], ubf[NC], slc[NC], sfu[NC];
            {
              const Real bt = lG_g(m,0,ie+1,k,j);
              const Real fsc = lG_g(m,2,ie+1,k,j);
              Real fb = 0.0;
              for (int cc=0; cc<NC; ++cc) {
                const int c = blk*NC + cc;
                const int b = bandc[cc];
                const Real da = cktpf_g(m,c,k,j)*Bb_g(m,b,ie+1,k,j);
                const Real rr = lP_g(m,0*nch_+c,ie+1,k,j);
                const Real ub = (rr*(1.0 + bt)*da + Sc[cc][ie+1])
                              * CkIdn(lP_g(m,0*nch_+c,ie+1,k,j), lG_g(m,0,ie+1,k,j));
                const Real db = da + bt*(da - ub);
                fb += lC_g(0,c)*(ub - db)*fsc;
                dcu[cc] = db;
                ubf[cc] = ub;
                const Real btop = Bb_g(m,b,ie,k,j);
                slc[cc] = btop;
                sfu[cc] = btop;
              }
              Fb_g(m,blk,ie+1,k,j) = fb;
            }
            for (int i=ie; i>icut-1; --i) {
              const Real bt = lG_g(m,0,i,k,j);
              const Real rat = lG_g(m,1,i,k,j);
              const Real fsc = lG_g(m,2,i,k,j);
              const Real idz = lG_g(m,3,i,k,j);
              Real src = 0.0, fb = 0.0, em = 0.0;
              for (int cc=0; cc<NC; ++cc) {
                const int c = blk*NC + cc;
                const int cw = (ck_nq_ == 2) ? (c & ~1) : c;   // ck-nq2: pair weights
                const int b = bandc[cc];
                const Real bown = Bb_g(m,b,i,k,j);
                // the layer joining cells i-1 and i: pass 1's slv/suu/sfv of cell i-1
                Real suv = bown, sfv = bown, snl = bown;
                if (i > icut) {
                  const Real bfar = Bb_g(m,b,i-1,k,j);
                  const Real wl = lP_g(m,2*nch_+cw,i-1,k,j);
                  const Real wu = lP_g(m,3*nch_+cw,i-1,k,j);
                  snl = wl*bfar + (1.0 - wl)*bown;
                  suv = wu*bown + (1.0 - wu)*bfar;
                  sfv = snl + (suv - snl)*lP_g(m,4*nch_+cw,i-1,k,j);
                }
                const Real e0 = ckc0_g(m,c,i,k,j);
                const Real cin = ckci_g(m,c,i,k,j);
                const Real cout = ckco_g(m,c,i,k,j);
                const Real tr = 1.0 - e0;
                const Real wfc = lC_g(0,c);
                const Real wz = wfc*idz;
                // the down ray: the upper half of cell i, then the lower half
                Real dI = dcu[cc];
                Real eh = cin*sfu[cc] + cout*slc[cc];
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                eh = cin*suv + cout*sfv;
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                // face i
                const Real rr = lP_g(m,0*nch_+c,i,k,j);
                const Real ub = (rr*(1.0 + bt)*dI + Sc[cc][i])*CkIdn(rr, bt);
                const Real db = dI + bt*(dI - ub);
                const Real dm = ub - db;
                fb += wfc*dm*fsc;
                Real ua = dI + dm*rat;
                dcu[cc] = db;
                // the up ray: the lower half, then the upper half, and the closure
                eh = cin*sfv + cout*suv;
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                eh = cin*slc[cc] + cout*sfu[cc];
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                src += wz*(ua - ubf[cc]);
                ubf[cc] = ub;
                em += CkEmW(lC_g(2,c), ckkro_g(m,c,i,k,j))*bown;
                slc[cc] = snl;
                sfu[cc] = sfv;
              }
              Src_g(m,blk,i,k,j) = src;
              Fb_g(m,blk,i,k,j) = fb;
              Em_g(m,blk,i,k,j) = em;
            }
          });
        };
        // ck_impl_lin_thr = 1: the same recurrences, ONE CHAIN PER THREAD: four times
        // the threads of the block kernel above (which runs at the chain kernel's 1.16
        // waves per SIMD) and a quarter of the state.  Each thread writes its chain's
        // Src and Fb into the partials, and rt_chain_ck_lin_sum adds a block's chains
        // in chain order (Fb and Em in exactly the block kernel's order, Src per chain
        // first, i.e. equal to round-off).
        auto lps_g = (cklon_ && ck_lps_ptr != nullptr) ? *ck_lps_ptr
                   : CkDum<DvceArray5D<Real>>("ck_lps_d");
        auto lpf_g = (cklon_ && ck_lpf_ptr != nullptr) ? *ck_lpf_ptr
                   : CkDum<DvceArray5D<Real>>("ck_lpf_d");
        // the per-block sums of the linear re-apply (shared by lin1 and lin1p)
        // ck-next: psum -- rt_chain_ck_lin1p wrote the first angle pair of every block
        // pre-added (chain 0's slot holds f0 + f1, which is bitwise the running sum
        // (0 + f0) + f1 below): the sum skips chain 1's slot, the same numbers
        auto launch_ck_lin_sum = [&](const bool psum) {
          // ck-next: a thread per (block, face) instead of per block: the same sums
          par_for("rt_chain_ck_lin_sum", DevExeSpace(), 0, nmb1, 0, nblk-1, is, ie+1,
                  ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int blk, const int i, const int k,
                        const int j) {
            constexpr int NC = RT_NB;
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int icut = icc_g(m,k,j);
            {
              Real src = 0.0, fb = 0.0, em = 0.0;
              if (i >= icut && icut <= ie) {
                for (int cc=0; cc<NC; ++cc) {
                  const int c = blk*NC + cc;
                  const bool skp = psum && (cc == 1);
                  if (!skp) fb += lpf_g(m,c,i,k,j);
                  if (i <= ie) {
                    const int b = (ck_nq_ == 1) ? (c/CK_NG) : (c/(2*CK_NG));
                    if (!skp) src += lps_g(m,c,i,k,j);
                    // ck-next: 2 (wfc/mu) kappa rho, formed (CkEmW); the two angles
                    // of a pair share kappa rho
                    const int ck = (ck_nq_ == 2) ? (c & ~1) : c;
                    em += CkEmW(lC_g(2,c), ckkro_g(m,ck,i,k,j))
                          *Bb_g(m,b,i,k,j);
                  }
                }
              }
              Src_g(m,blk,i,k,j) = src;
              Fb_g(m,blk,i,k,j) = fb;
              Em_g(m,blk,i,k,j) = em;
            }
          });
        };
        auto launch_ck_lin1 = [&]() {
          const size_t sltsz = static_cast<size_t>(n1)*CKS_W;
          const size_t slntl = (static_cast<size_t>(nmb1 + 1)*nch_*scl3*scl2 + CKS_W - 1)
                               /CKS_W;
          CkScrEnsure(slntl*sltsz);
          Real *slbuf = rt_ckscr_ptr->data();
          CkParFor4("rt_chain_ck_lin1", cklw_, 0, nmb1, 0, nch_-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            const int b = (ck_nq_ == 1) ? (c/CK_NG) : (c/(2*CK_NG));
            // ck-nq2: lP slots 2-4 (kappa rho weights) are stored once per angle pair
            const int cw = (ck_nq_ == 2) ? (c & ~1) : c;
            const Real wfc = lC_g(0,c);
            const int slt = ((m*nch_ + c)*scl3 + (k - ks))*scl2 + (j - js);
            auto Sc = CkScrRows<Real>(slbuf, 0, sltsz, n1, slt)[0];
            const Real bcut = Bb_g(m,b,icut,k,j);
            Real sfc = bcut, suc = bcut, ss = bcut + lC_g(1,c);
            if (ckdif_) {                    // ck_dif_dtau: the flux datum
              const Real gdf = difg_g(m,c,k,j);
              if (gdf >= 0.0) ss = gdf*(Bb_g(m,b,icut-1,k,j) - bcut);
            }
            Real bown = bcut;
            for (int i=icut; i<ie+1; ++i) {
              const Real bt = lG_g(m,0,i,k,j);
              Real slv = bown, sfv = bown, suu = bown;
              Real bnext = bown;
              if (i < ie) {
                bnext = Bb_g(m,b,i+1,k,j);
                const Real wl = lP_g(m,2*nch_+cw,i,k,j);
                const Real wu = lP_g(m,3*nch_+cw,i,k,j);
                slv = wl*bown + (1.0 - wl)*bnext;
                suu = wu*bnext + (1.0 - wu)*bown;
                sfv = slv + (suu - slv)*lP_g(m,4*nch_+cw,i,k,j);
              }
              const Real e0 = ckc0_g(m,c,i,k,j);
              const Real cin = ckci_g(m,c,i,k,j);
              const Real cout = ckco_g(m,c,i,k,j);
              const Real tr = 1.0 - e0;
              const Real pl = cin*sfc + cout*suc;
              const Real ql = cin*suc + cout*sfc;
              const Real pu = cin*slv + cout*sfv;
              const Real qu = cin*sfv + cout*slv;
              Sc[i] = ss;
              const Real idn = CkIdn(lP_g(m,0*nch_+c,i,k,j), lG_g(m,0,i,k,j));
              const Real rn = (lP_g(m,0*nch_+c,i,k,j) + bt)*idn;
              Real sv = (1.0 - bt)*ss*idn;
              sv = tr*(rn*ql + sv) + pl;
              const Real r2 = tr*tr*rn;
              ss = tr*(r2*qu + sv) + pu;
              sfc = sfv;
              suc = suu;
              bown = bnext;
            }
            Sc[ie+1] = ss;
            Real dcu, ubf, slc, sfu;
            {
              const Real bt = lG_g(m,0,ie+1,k,j);
              const Real da = cktpf_g(m,c,k,j)*Bb_g(m,b,ie+1,k,j);
              const Real rr = lP_g(m,0*nch_+c,ie+1,k,j);
              const Real ub = (rr*(1.0 + bt)*da + Sc[ie+1])*CkIdn(rr, bt);
              const Real db = da + bt*(da - ub);
              lpf_g(m,c,ie+1,k,j) = wfc*(ub - db)*lG_g(m,2,ie+1,k,j);
              dcu = db;
              ubf = ub;
              slc = Bb_g(m,b,ie,k,j);
              sfu = slc;
            }
            bown = sfu;
            for (int i=ie; i>icut-1; --i) {
              const Real bt = lG_g(m,0,i,k,j);
              Real suv = bown, sfv = bown, snl = bown;
              Real bnext = bown;
              if (i > icut) {
                bnext = Bb_g(m,b,i-1,k,j);
                const Real wl = lP_g(m,2*nch_+cw,i-1,k,j);
                const Real wu = lP_g(m,3*nch_+cw,i-1,k,j);
                snl = wl*bnext + (1.0 - wl)*bown;
                suv = wu*bown + (1.0 - wu)*bnext;
                sfv = snl + (suv - snl)*lP_g(m,4*nch_+cw,i-1,k,j);
              }
              const Real e0 = ckc0_g(m,c,i,k,j);
              const Real cin = ckci_g(m,c,i,k,j);
              const Real cout = ckco_g(m,c,i,k,j);
              const Real tr = 1.0 - e0;
              const Real wz = wfc*lG_g(m,3,i,k,j);
              Real src = 0.0;
              Real dI = dcu;
              Real eh = cin*sfu + cout*slc;
              src += wz*(e0*dI - eh);
              dI = tr*dI + eh;
              eh = cin*suv + cout*sfv;
              src += wz*(e0*dI - eh);
              dI = tr*dI + eh;
              const Real rr = lP_g(m,0*nch_+c,i,k,j);
              const Real ub = (rr*(1.0 + bt)*dI + Sc[i])*CkIdn(rr, bt);
              const Real db = dI + bt*(dI - ub);
              const Real dm = ub - db;
              lpf_g(m,c,i,k,j) = wfc*dm*lG_g(m,2,i,k,j);
              Real ua = dI + dm*lG_g(m,1,i,k,j);
              dcu = db;
              eh = cin*sfv + cout*suv;
              src += wz*(e0*ua - eh);
              ua = tr*ua + eh;
              eh = cin*slc + cout*sfu;
              src += wz*(e0*ua - eh);
              ua = tr*ua + eh;
              src += wz*(ua - ubf);
              ubf = ub;
              lps_g(m,c,i,k,j) = src;
              slc = snl;
              sfu = sfv;
              bown = bnext;
            }
          });
          launch_ck_lin_sum(false);
        };
        // ---- ck-nq2: rt_chain_ck_lin1 with BOTH angles of a (band, g) pair in one
        // thread.  With ck_nquad = 2 the chains c = 2p and 2p+1 share the band, the
        // Planck column B_b, the column geometry lG and the three BFace/interpolation
        // weights (lP slots 2-4, which depend on kappa rho alone), so one thread loads
        // those once and carries the two angles' recurrences side by side.  Each
        // angle's arithmetic is term for term the one of rt_chain_ck_lin1, and the
        // per-chain partials go to the same lps/lpf slots, so the sum kernel and the
        // answer are BITWISE those of the unpaired kernel.
        auto launch_ck_lin1p = [&]() {
          const int npr = nch_/2;
          const size_t sltsz = static_cast<size_t>(2)*n1*CKS_W;
          const size_t slntl = (static_cast<size_t>(nmb1 + 1)*npr*scl3*scl2 + CKS_W - 1)
                               /CKS_W;
          CkScrEnsure(slntl*sltsz);
          Real *slbuf = rt_ckscr_ptr->data();
          CkParFor4("rt_chain_ck_lin1p", cklw_, 0, nmb1, 0, npr-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            const int c0 = 2*pp;
            // ck-next: the first pair of a chain block stores its two partials added
            // (see launch_ck_lin_sum)
            const bool pfst = (RT_NB % 2 == 0) && (c0 % RT_NB == 0);
            const int b = c0/(2*CK_NG);
            Real wfc[2];
            wfc[0] = lC_g(0,c0);
            wfc[1] = lC_g(0,c0+1);
            const int slt = ((m*npr + pp)*scl3 + (k - ks))*scl2 + (j - js);
            auto Scr = CkScrRows<Real>(slbuf, 0, sltsz, n1, slt);
            const Real bcut = Bb_g(m,b,icut,k,j);
            Real sfc = bcut, suc = bcut, ss[2];
            ss[0] = bcut + lC_g(1,c0);
            ss[1] = bcut + lC_g(1,c0+1);
            if (ckdif_) {                    // ck_dif_dtau: the flux datum
              for (int q=0; q<2; ++q) {
                const Real gdf = difg_g(m,c0+q,k,j);
                if (gdf >= 0.0) ss[q] = gdf*(Bb_g(m,b,icut-1,k,j) - bcut);
              }
            }
            Real bown = bcut;
            // ck-next: every load of cell i+1 is issued before cell i is worked (the
            // stores keep the compiler from hoisting them): the same arithmetic
            Real bt_n = lG_g(m,0,icut,k,j);
            Real bnx_n = bown;
            if (icut < ie) bnx_n = Bb_g(m,b,icut+1,k,j);
            // ck-next: cin, cout and 1/(1 + R beta) are formed from e0, kappa rho
            // (one per pair), dz and R (CkCinCout, CkIdn), not read.  ck-lin2: so are
            // the layer weights, lP slots 2-4 (CkLayW, from the kappa rho and dz of the
            // two cells: those of cell i+1 are loaded two cells ahead)
            Real kro_n = ckkro_g(m,c0,icut,k,j), dz_n = lG_g(m,4,icut,k,j);
            Real kro_u = 0.0, dz_u = 0.0;
            if (icut < ie) {
              kro_u = ckkro_g(m,c0,icut+1,k,j);
              dz_u = lG_g(m,4,icut+1,k,j);
            }
            // ck-lin2: e0 formed too (CkE0), not read
            Real r0_n[2];
            for (int q=0; q<2; ++q) r0_n[q] = lP_g(m,0*nch_+c0+q,icut,k,j);
            for (int i=icut; i<ie+1; ++i) {
              const Real bt = bt_n;
              const Real bnx = bnx_n;
              const Real krov = kro_n, dzv = dz_n, kruv = kro_u, dzuv = dz_u;
              Real e0v[2], r0v[2];
              for (int q=0; q<2; ++q) {
                e0v[q] = CkE0(krov, dzv, mug[q]);
                r0v[q] = r0_n[q];
              }
              if (i < ie) {
                bt_n = lG_g(m,0,i+1,k,j);
                kro_n = kro_u;
                dz_n = dz_u;
                if (i+1 < ie) {
                  bnx_n = Bb_g(m,b,i+2,k,j);
                  kro_u = ckkro_g(m,c0,i+2,k,j);
                  dz_u = lG_g(m,4,i+2,k,j);
                }
                for (int q=0; q<2; ++q) r0_n[q] = lP_g(m,0*nch_+c0+q,i+1,k,j);
              }
              Real slv = bown, sfv = bown, suu = bown;
              Real bnext = bown;
              if (i < ie) {
                Real wl, wu, ff;
                CkLayW(krov, kruv, dzv, dzuv, bface_on, wl, wu, ff);
                bnext = bnx;
                slv = wl*bown + (1.0 - wl)*bnext;
                suu = wu*bnext + (1.0 - wu)*bown;
                sfv = slv + (suu - slv)*ff;
              }
              for (int q=0; q<2; ++q) {
                const Real e0 = e0v[q];
                Real cin, cout;
                CkCinCout(krov, dzv, mug[q], e0, cin, cout);
                const Real tr = 1.0 - e0;
                const Real pl = cin*sfc + cout*suc;
                const Real ql = cin*suc + cout*sfc;
                const Real pu = cin*slv + cout*sfv;
                const Real qu = cin*sfv + cout*slv;
                Scr[q][i] = ss[q];
                const Real idn = CkIdn(r0v[q], bt);
                const Real rn = (r0v[q] + bt)*idn;
                Real sv = (1.0 - bt)*ss[q]*idn;
                sv = tr*(rn*ql + sv) + pl;
                const Real r2 = tr*tr*rn;
                ss[q] = tr*(r2*qu + sv) + pu;
              }
              sfc = sfv;
              suc = suu;
              bown = bnext;
            }
            Real dcu[2], ubf[2];
            Real slc, sfu;
            {
              const Real bt = lG_g(m,0,ie+1,k,j);
              Real fq[2];
              for (int q=0; q<2; ++q) {
                const int c = c0 + q;
                Scr[q][ie+1] = ss[q];
                const Real da = cktpf_g(m,c,k,j)*Bb_g(m,b,ie+1,k,j);
                const Real rr = lP_g(m,0*nch_+c,ie+1,k,j);
                const Real ub = (rr*(1.0 + bt)*da + Scr[q][ie+1])
                              *CkIdn(rr, bt);
                const Real db = da + bt*(da - ub);
                fq[q] = CkMulRn(wfc[q]*(ub - db), lG_g(m,2,ie+1,k,j));
                dcu[q] = db;
                ubf[q] = ub;
              }
              if (pfst) {
                lpf_g(m,c0,ie+1,k,j) = fq[0] + fq[1];
              } else {
                lpf_g(m,c0,ie+1,k,j) = fq[0];
                lpf_g(m,c0+1,ie+1,k,j) = fq[1];
              }
              slc = Bb_g(m,b,ie,k,j);
              sfu = slc;
            }
            bown = sfu;
            // ck-next: cell i-1's loads before cell i is worked, as in pass 1
            Real bt_p = lG_g(m,0,ie,k,j), g1_p = lG_g(m,1,ie,k,j);
            Real g2_p = lG_g(m,2,ie,k,j), g3_p = lG_g(m,3,ie,k,j);
            Real bnx_p = bown;
            // ck-lin2: the layer (i-1, i) weights formed (CkLayW): kappa rho and dz of
            // cell i-1 are loaded two cells ahead
            Real kro_p = ckkro_g(m,c0,ie,k,j), dz_p = lG_g(m,4,ie,k,j);
            Real kro_l = 0.0, dz_l = 0.0;
            if (ie > icut) {
              bnx_p = Bb_g(m,b,ie-1,k,j);
              kro_l = ckkro_g(m,c0,ie-1,k,j);
              dz_l = lG_g(m,4,ie-1,k,j);
            }
            Real r0_p[2], sc_p[2];
            for (int q=0; q<2; ++q) {
              r0_p[q] = lP_g(m,0*nch_+c0+q,ie,k,j);
              sc_p[q] = Scr[q][ie];
            }
            for (int i=ie; i>icut-1; --i) {
              const Real bt = bt_p;
              const Real g1 = g1_p, g2 = g2_p, g3 = g3_p;
              const Real bnx = bnx_p;
              const Real krov = kro_p, dzv = dz_p, krlv = kro_l, dzlv = dz_l;
              Real e0v[2], r0v[2], scv[2];
              for (int q=0; q<2; ++q) {
                e0v[q] = CkE0(krov, dzv, mug[q]);
                r0v[q] = r0_p[q];
                scv[q] = sc_p[q];
              }
              if (i > icut) {
                bt_p = lG_g(m,0,i-1,k,j);
                kro_p = kro_l;
                dz_p = dz_l;
                g1_p = lG_g(m,1,i-1,k,j);
                g2_p = lG_g(m,2,i-1,k,j);
                g3_p = lG_g(m,3,i-1,k,j);
                if (i-1 > icut) {
                  bnx_p = Bb_g(m,b,i-2,k,j);
                  kro_l = ckkro_g(m,c0,i-2,k,j);
                  dz_l = lG_g(m,4,i-2,k,j);
                }
                for (int q=0; q<2; ++q) {
                  r0_p[q] = lP_g(m,0*nch_+c0+q,i-1,k,j);
                  sc_p[q] = Scr[q][i-1];
                }
              }
              Real suv = bown, sfv = bown, snl = bown;
              Real bnext = bown;
              if (i > icut) {
                Real wl, wu, ff;
                CkLayW(krlv, krov, dzlv, dzv, bface_on, wl, wu, ff);
                bnext = bnx;
                snl = wl*bnext + (1.0 - wl)*bown;
                suv = wu*bown + (1.0 - wu)*bnext;
                sfv = snl + (suv - snl)*ff;
              }
              Real fq[2], sq[2];
              for (int q=0; q<2; ++q) {
                const Real e0 = e0v[q];
                Real cin, cout;
                CkCinCout(krov, dzv, mug[q], e0, cin, cout);
                const Real tr = 1.0 - e0;
                const Real wz = wfc[q]*g3;
                Real src = 0.0;
                Real dI = dcu[q];
                Real eh = cin*sfu + cout*slc;
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                eh = cin*suv + cout*sfv;
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                const Real rr = r0v[q];
                const Real ub = (rr*(1.0 + bt)*dI + scv[q])*CkIdn(rr, bt);
                const Real db = dI + bt*(dI - ub);
                const Real dm = ub - db;
                fq[q] = CkMulRn(wfc[q]*dm, g2);   // not contracted into the add
                Real ua = dI + dm*g1;
                dcu[q] = db;
                eh = cin*sfv + cout*suv;
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                eh = cin*slc + cout*sfu;
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                src += wz*(ua - ubf[q]);
                ubf[q] = ub;
                sq[q] = src;
              }
              if (pfst) {
                lpf_g(m,c0,i,k,j) = fq[0] + fq[1];
                lps_g(m,c0,i,k,j) = sq[0] + sq[1];
              } else {
                lpf_g(m,c0,i,k,j) = fq[0];
                lpf_g(m,c0+1,i,k,j) = fq[1];
                lps_g(m,c0,i,k,j) = sq[0];
                lps_g(m,c0+1,i,k,j) = sq[1];
              }
              slc = snl;
              sfu = sfv;
              bown = bnext;
            }
          });
          launch_ck_lin_sum(RT_NB % 2 == 0);
        };
        // ck_sph_face = 5: rt_chain_ck_lin1p on the VEF face (the chain kernel's tm form
        // with vef), a kernel of its own, so that the default one is the code it always
        // was
        // ck_sph_face = 6: the re-apply of the ladder (see ck_sph_face = 6 and
        // ck_f6_build) on the new B_b: the Norton sources down (I_hi = I_up + g_p K_eq,
        // I_lo = a_s I_hi + g_p K_eq, I_up = a_c I_lo), the face potentials up (the cut
        // node from the wall, K_lo = a_c (K_hi below + R_c I_lo), K_hi = a_s K_lo +
        // b_s I_hi), then the fluxes accumulated from the top, G = y_t K_top, G += dep,
        // dep = g_p (K_lo + K_hi - 2 K_eq) the cell's absorption minus emission: the
        // pair's Src = 4 pi g_w dOmega dep/V and its flux Phi/A at every face, so that
        // V Src = (A F)_i - (A F)_{i+1} to round-off.  One thread per (band, g) pair,
        // two scratch rows (I, then K).
        auto launch_ck_lin6 = [&]() {
          const int npr = nch_/2;
          const size_t sltsz = static_cast<size_t>(2)*n1*CKS_W;
          const size_t slntl = (static_cast<size_t>(nmb1 + 1)*npr*scl3*scl2 + CKS_W - 1)
                               /CKS_W;
          CkScrEnsure(slntl*sltsz);
          Real *slbuf = rt_ckscr_ptr->data();
          CkParFor4("rt_chain_ck_lin6", cklw_, 0, nmb1, 0, npr-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            const int c0 = 2*pp, c1 = c0 + 1;
            const bool pfst = (RT_NB % 2 == 0) && (c0 % RT_NB == 0);
            const int b = c0/(2*CK_NG);
            const int slt = ((m*npr + pp)*scl3 + (k - ks))*scl2 + (j - js);
            auto Scr = CkScrRows<Real>(slbuf, 0, sltsz, n1, slt);
            auto P = [&](const int sl, const int c, const int i) -> Real {
              return lP_g(m,sl*nch_+c,i,k,j);
            };
            // down: the Norton sources of both face nodes of every cell
            Real I = 0.0;
            for (int i=ie; i>=icut; --i) {
              const Real src = P(0,c0,i)*P(1,c0,i)*Bb_g(m,b,i,k,j);
              const Real ihi = I + src;
              const Real ilo = P(2,c0,i)*ihi + src;
              Scr[1][i] = ihi;
              Scr[0][i] = ilo;
              I = P(4,c0,i)*ilo;
            }
            // up: the face potentials
            Real khi = 0.0;
            for (int i=icut; i<=ie; ++i) {
              Real klo;
              if (i == icut) {
                // ck_wall_flux_exact: I_int scaled by c6 = beta (y_lo + y_w)/((2 +
                // beta) y_lo), so the wall current's I_int part is r^2 I_int/4
                Real iw = lC_g(1,c0);
                if (wex_) {
                  const Real bwl = vfg_g(m,pp,icut-1,k,j);
                  iw *= bwl*(P(0,c1,icut) + P(4,c1,icut))/((2.0 + bwl)*P(0,c1,icut));
                }
                klo = (Scr[0][i] + P(3,c1,icut)*(Bb_g(m,b,icut,k,j) + iw))
                      /(P(0,c1,icut) + P(4,c1,icut));
              } else {
                klo = P(4,c0,i)*(khi + P(5,c0,i)*Scr[0][i]);
              }
              khi = P(2,c0,i)*klo + P(3,c0,i)*Scr[1][i];
              Scr[0][i] = klo;
              Scr[1][i] = khi;
            }
            // the fluxes from the top, the deposit per cell
            const Real wpr = lC_g(2,c0) + lC_g(2,c1);           // 4 pi g_w
            const Real cnv = wpr*AFC(m,k,j,ie+1)/SQR(X1F(m,ie+1));
            Real G = P(3,c1,ie+1)*khi;
            lpf_g(m,c0,ie+1,k,j) = cnv*G/AFC(m,k,j,ie+1);
            if (!pfst) lpf_g(m,c1,ie+1,k,j) = 0.0;
            for (int i=ie; i>=icut; --i) {
              const Real keqb = P(1,c0,i)*Bb_g(m,b,i,k,j);
              const Real d = P(0,c0,i)*((Scr[0][i] - keqb) + (Scr[1][i] - keqb));
              G += d;
              lps_g(m,c0,i,k,j) = cnv*d/VLS(m,k,j,i);
              lpf_g(m,c0,i,k,j) = cnv*G/AFC(m,k,j,i);
              if (!pfst) {
                lps_g(m,c1,i,k,j) = 0.0;
                lpf_g(m,c1,i,k,j) = 0.0;
              }
            }
          });
          launch_ck_lin_sum(RT_NB % 2 == 0);
        };
        auto launch_ck_lin1p_vef = [&]() {
          const int npr = nch_/2;
          const size_t sltsz = static_cast<size_t>(2)*n1*CKS_W;
          const size_t slntl = (static_cast<size_t>(nmb1 + 1)*npr*scl3*scl2 + CKS_W - 1)
                               /CKS_W;
          CkScrEnsure(slntl*sltsz);
          Real *slbuf = rt_ckscr_ptr->data();
          CkParFor4("rt_chain_ck_lin1p", cklw_, 0, nmb1, 0, npr-1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int icut = icc_g(m,k,j);
            if (icut > ie) return;
            const int c0 = 2*pp;
            // ck-next: the first pair of a chain block stores its two partials added
            // (see launch_ck_lin_sum)
            const bool pfst = (RT_NB % 2 == 0) && (c0 % RT_NB == 0);
            const int b = c0/(2*CK_NG);
            Real wfc[2];
            wfc[0] = lC_g(0,c0);
            wfc[1] = lC_g(0,c0+1);
            const int slt = ((m*npr + pp)*scl3 + (k - ks))*scl2 + (j - js);
            auto Scr = CkScrRows<Real>(slbuf, 0, sltsz, n1, slt);
            const Real bcut = Bb_g(m,b,icut,k,j);
            Real sfc = bcut, suc = bcut, ss[2];
            ss[0] = bcut + lC_g(1,c0);
            ss[1] = bcut + lC_g(1,c0+1);
            if (ckdif_) {                    // ck_dif_dtau: the flux datum
              for (int q=0; q<2; ++q) {
                const Real gdf = difg_g(m,c0+q,k,j);
                if (gdf >= 0.0) ss[q] = gdf*(Bb_g(m,b,icut-1,k,j) - bcut);
              }
            }
            // ck_sph_face = 5: the Marshak-type wall, Sc = (B + Iint) (1 + beta/2)/(1 +
            // beta mu_eff) (its R is the one ck_lin_build stored at the cut); the pair
            // (b, g) = pp has one s = sqrt(3 f) per frame, its two angles one mu_eff
            if (vefl_) {
              const Real bw = vfs_g(m,pp,icut,k,j)/SQRT3_VEF;
              const Real bwb = vfg_g(m,pp,icut-1,k,j);
              const Real vs0 = (1.0 + 0.5*bwb)/(1.0 + bwb*bw);
              // ck_wall_flux_exact: I_int scaled so the wall's net flux is pi I_int
              const Real cwx = wex_ ? (1.0 + bwb*bw)/(2.0*bw*(1.0 + 0.5*bwb)) : 1.0;
              for (int q=0; q<2; ++q) ss[q] = (bcut + lC_g(1,c0+q)*cwx)*vs0;
            }
            Real bown = bcut;
            // ck-next: every load of cell i+1 is issued before cell i is worked (the
            // stores keep the compiler from hoisting them): the same arithmetic
            Real bt_n = lG_g(m,0,icut,k,j);
            Real bnx_n = bown;
            if (icut < ie) bnx_n = Bb_g(m,b,icut+1,k,j);
            // ck-next: cin, cout and 1/(1 + R beta) are formed from e0, kappa rho
            // (one per pair), dz and R (CkCinCout, CkIdn), not read.  ck-lin2: so are
            // the layer weights, lP slots 2-4 (CkLayW, from the kappa rho and dz of the
            // two cells: those of cell i+1 are loaded two cells ahead)
            Real kro_n = ckkro_g(m,c0,icut,k,j), dz_n = lG_g(m,4,icut,k,j);
            Real kro_u = 0.0, dz_u = 0.0;
            if (icut < ie) {
              kro_u = ckkro_g(m,c0,icut+1,k,j);
              dz_u = lG_g(m,4,icut+1,k,j);
            }
            // ck-lin2: e0 formed too (CkE0), not read
            Real r0_n[2];
            for (int q=0; q<2; ++q) r0_n[q] = lP_g(m,0*nch_+c0+q,icut,k,j);
            for (int i=icut; i<ie+1; ++i) {
              const Real bt = bt_n;
              const Real bnx = bnx_n;
              const Real krov = kro_n, dzv = dz_n, kruv = kro_u, dzuv = dz_u;
              Real e0v[2], r0v[2], muv[2];
              for (int q=0; q<2; ++q) {
                muv[q] = mug[q];
                r0v[q] = r0_n[q];
              }
              // ck_sph_face = 5: both angles at mu_eff = s/sqrt 3 of cell i; rp and
              // gamma of face i
              Real vrp = 1.0, vgm = 1.0;
              if (vefl_) {
                const Real vsa = vfs_g(m,pp,i+1,k,j);
                muv[0] = vsa/SQRT3_VEF;
                muv[1] = muv[0];
                vrp = (1.0 - bt)/(1.0 + bt)*vfs_g(m,pp,i,k,j)/vsa;
                vgm = vfg_g(m,pp,i,k,j);
              }
              for (int q=0; q<2; ++q) e0v[q] = CkE0(krov, dzv, muv[q]);
              if (i < ie) {
                bt_n = lG_g(m,0,i+1,k,j);
                kro_n = kro_u;
                dz_n = dz_u;
                if (i+1 < ie) {
                  bnx_n = Bb_g(m,b,i+2,k,j);
                  kro_u = ckkro_g(m,c0,i+2,k,j);
                  dz_u = lG_g(m,4,i+2,k,j);
                }
                for (int q=0; q<2; ++q) r0_n[q] = lP_g(m,0*nch_+c0+q,i+1,k,j);
              }
              Real slv = bown, sfv = bown, suu = bown;
              Real bnext = bown;
              if (i < ie) {
                Real wl, wu, ff;
                CkLayW(krov, kruv, dzv, dzuv, bface_on, wl, wu, ff);
                bnext = bnx;
                slv = wl*bown + (1.0 - wl)*bnext;
                suu = wu*bnext + (1.0 - wu)*bown;
                sfv = slv + (suu - slv)*ff;
              }
              for (int q=0; q<2; ++q) {
                const Real e0 = e0v[q];
                Real cin, cout;
                CkCinCout(krov, dzv, muv[q], e0, cin, cout);
                const Real tr = 1.0 - e0;
                const Real pl = cin*sfc + cout*suc;
                const Real ql = cin*suc + cout*sfc;
                const Real pu = cin*slv + cout*sfv;
                const Real qu = cin*sfv + cout*slv;
                Scr[q][i] = ss[q];
                Real rn, sv;
                if (vefl_) {
                  Real c1;
                  CkVefMap(r0v[q], vrp, vgm, rn, c1);
                  sv = c1*ss[q];
                } else {
                  const Real idn = CkIdn(r0v[q], bt);
                  rn = (r0v[q] + bt)*idn;
                  sv = (1.0 - bt)*ss[q]*idn;
                }
                sv = tr*(rn*ql + sv) + pl;
                const Real r2 = tr*tr*rn;
                ss[q] = tr*(r2*qu + sv) + pu;
              }
              sfc = sfv;
              suc = suu;
              bown = bnext;
            }
            Real dcu[2], ubf[2];
            Real slc, sfu;
            {
              const Real bt = lG_g(m,0,ie+1,k,j);
              Real fq[2];
              for (int q=0; q<2; ++q) {
                const int c = c0 + q;
                Scr[q][ie+1] = ss[q];
                // problem/ck_sph_top = 1: nothing comes down from above
                const Real da = (tpm_ == 1) ? 0.0
                              : cktpf_g(m,c,k,j)*Bb_g(m,b,ie+1,k,j);
                const Real rr = lP_g(m,0*nch_+c,ie+1,k,j);
                if (vefl_) {
                  // ck_sph_face = 5: the down ray reflected at the top (Auer's
                  // H = h J), d_a = (zeta Sc' + datum)/(1 - zeta R'), then the face
                  // solve of the jump conditions (the chain kernel's vsolve)
                  const Real vsb = vfs_g(m,pp,ie+1,k,j);
                  const Real vrp = (1.0 - bt)/(1.0 + bt)*vsb/vfs_g(m,pp,ie+2,k,j);
                  const Real vgm = vfg_g(m,pp,ie+1,k,j);
                  const Real zt = vfg_g(m,pp,ie+2,k,j);
                  const Real ev = (1.0 + rr) + vgm*vrp*(1.0 - rr);
                  const Real rup = 1.0 - 2.0*vgm*vrp*(1.0 - rr)/ev;
                  const Real sup = 2.0*vrp*ss[q]/ev;
                  const Real dat = (zt*sup + da)/(1.0 - zt*rup);
                  const Real db = (2.0*vgm*dat - (1.0 - vgm*vrp)*ss[q])/ev;
                  const Real ub = rr*db + ss[q];
                  fq[q] = CkMulRn(wfc[q]*vsb/(SQRT3_VEF*mug[q])*(ub - db),
                                  lG_g(m,2,ie+1,k,j));
                  dcu[q] = db;
                  ubf[q] = ub;
                  continue;
                }
                const Real ub = (rr*(1.0 + bt)*da + Scr[q][ie+1])
                              *CkIdn(rr, bt);
                const Real db = da + bt*(da - ub);
                fq[q] = CkMulRn(wfc[q]*(ub - db), lG_g(m,2,ie+1,k,j));
                dcu[q] = db;
                ubf[q] = ub;
              }
              if (pfst) {
                lpf_g(m,c0,ie+1,k,j) = fq[0] + fq[1];
              } else {
                lpf_g(m,c0,ie+1,k,j) = fq[0];
                lpf_g(m,c0+1,ie+1,k,j) = fq[1];
              }
              slc = Bb_g(m,b,ie,k,j);
              sfu = slc;
            }
            bown = sfu;
            // ck-next: cell i-1's loads before cell i is worked, as in pass 1
            Real bt_p = lG_g(m,0,ie,k,j), g1_p = lG_g(m,1,ie,k,j);
            Real g2_p = lG_g(m,2,ie,k,j), g3_p = lG_g(m,3,ie,k,j);
            Real bnx_p = bown;
            // ck-lin2: the layer (i-1, i) weights formed (CkLayW): kappa rho and dz of
            // cell i-1 are loaded two cells ahead
            Real kro_p = ckkro_g(m,c0,ie,k,j), dz_p = lG_g(m,4,ie,k,j);
            Real kro_l = 0.0, dz_l = 0.0;
            if (ie > icut) {
              bnx_p = Bb_g(m,b,ie-1,k,j);
              kro_l = ckkro_g(m,c0,ie-1,k,j);
              dz_l = lG_g(m,4,ie-1,k,j);
            }
            Real r0_p[2], sc_p[2];
            for (int q=0; q<2; ++q) {
              r0_p[q] = lP_g(m,0*nch_+c0+q,ie,k,j);
              sc_p[q] = Scr[q][ie];
            }
            for (int i=ie; i>icut-1; --i) {
              const Real bt = bt_p;
              const Real g1 = g1_p, g2 = g2_p, g3 = g3_p;
              const Real bnx = bnx_p;
              const Real krov = kro_p, dzv = dz_p, krlv = kro_l, dzlv = dz_l;
              Real e0v[2], r0v[2], scv[2], muv[2];
              for (int q=0; q<2; ++q) {
                muv[q] = mug[q];
                r0v[q] = r0_p[q];
                scv[q] = sc_p[q];
              }
              // ck_sph_face = 5: s below face i (vsb) and in cell i (vsa), gamma
              Real vsb = 1.0, vsa = 1.0, vgm = 1.0;
              if (vefl_) {
                vsb = vfs_g(m,pp,i,k,j);
                vsa = vfs_g(m,pp,i+1,k,j);
                vgm = vfg_g(m,pp,i,k,j);
                muv[0] = vsa/SQRT3_VEF;
                muv[1] = muv[0];
              }
              for (int q=0; q<2; ++q) e0v[q] = CkE0(krov, dzv, muv[q]);
              if (i > icut) {
                bt_p = lG_g(m,0,i-1,k,j);
                kro_p = kro_l;
                dz_p = dz_l;
                g1_p = lG_g(m,1,i-1,k,j);
                g2_p = lG_g(m,2,i-1,k,j);
                g3_p = lG_g(m,3,i-1,k,j);
                if (i-1 > icut) {
                  bnx_p = Bb_g(m,b,i-2,k,j);
                  kro_l = ckkro_g(m,c0,i-2,k,j);
                  dz_l = lG_g(m,4,i-2,k,j);
                }
                for (int q=0; q<2; ++q) {
                  r0_p[q] = lP_g(m,0*nch_+c0+q,i-1,k,j);
                  sc_p[q] = Scr[q][i-1];
                }
              }
              Real suv = bown, sfv = bown, snl = bown;
              Real bnext = bown;
              if (i > icut) {
                Real wl, wu, ff;
                CkLayW(krlv, krov, dzlv, dzv, bface_on, wl, wu, ff);
                bnext = bnx;
                snl = wl*bnext + (1.0 - wl)*bown;
                suv = wu*bown + (1.0 - wu)*bnext;
                sfv = snl + (suv - snl)*ff;
              }
              Real fq[2], sq[2];
              for (int q=0; q<2; ++q) {
                const Real e0 = e0v[q];
                Real cin, cout;
                CkCinCout(krov, dzv, muv[q], e0, cin, cout);
                const Real tr = 1.0 - e0;
                const Real wz = vefl_ ? wfc[q]*vsa/(SQRT3_VEF*mug[q])*g3 : wfc[q]*g3;
                Real src = 0.0;
                Real dI = dcu[q];
                Real eh = cin*sfu + cout*slc;
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                eh = cin*suv + cout*sfv;
                src += wz*(e0*dI - eh);
                dI = tr*dI + eh;
                const Real rr = r0v[q];
                Real ub, db, dm, ua;
                if (vefl_) {
                  // ck_sph_face = 5: the face solve of the jump conditions
                  const Real ratv = g1*vsb/vsa;
                  const Real ev = (1.0 + rr) + vgm*ratv*(1.0 - rr);
                  db = (2.0*vgm*dI - (1.0 - vgm*ratv)*scv[q])/ev;
                  ub = rr*db + scv[q];
                  dm = ub - db;
                  fq[q] = CkMulRn(wfc[q]*vsb/(SQRT3_VEF*mug[q])*dm, g2);
                  ua = dI + dm*ratv;
                } else {
                  ub = (rr*(1.0 + bt)*dI + scv[q])*CkIdn(rr, bt);
                  db = dI + bt*(dI - ub);
                  dm = ub - db;
                  fq[q] = CkMulRn(wfc[q]*dm, g2);   // not contracted into the add
                  ua = dI + dm*g1;
                }
                dcu[q] = db;
                eh = cin*sfv + cout*suv;
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                eh = cin*slc + cout*sfu;
                src += wz*(e0*ua - eh);
                ua = tr*ua + eh;
                src += wz*(ua - ubf[q]);
                ubf[q] = ub;
                sq[q] = src;
              }
              if (pfst) {
                lpf_g(m,c0,i,k,j) = fq[0] + fq[1];
                lps_g(m,c0,i,k,j) = sq[0] + sq[1];
              } else {
                lpf_g(m,c0,i,k,j) = fq[0];
                lpf_g(m,c0+1,i,k,j) = fq[1];
                lps_g(m,c0,i,k,j) = sq[0];
                lps_g(m,c0+1,i,k,j) = sq[1];
              }
              slc = snl;
              sfu = sfv;
              bown = bnext;
            }
          });
          launch_ck_lin_sum(RT_NB % 2 == 0);
        };
        const bool cklthr1_ = (ck_impl_lin_thr == 1);
        auto launch_ck_lin_tier = [&]() {
          if (cklthr1_ && ck_nq_ == 2) {
            if (vef6_) {
              launch_ck_lin6();
            } else if (vefl_) {
              launch_ck_lin1p_vef();
            } else {
              launch_ck_lin1p();
            }
          } else if (cklthr1_) {
            launch_ck_lin1();
          } else {
            launch_ck_lin();
          }
        };
        // ck-store: ck_coef (see ckstsp_).  Every value is the one the tm body of
        // rt_chain_ck stores on this pass: kapof (the k-table plus the continuum; the
        // q = 1 chain of a pair takes the q = 0 chain's number, as there), krof, the
        // `cofs` triple of the half layer 0.5 kappa rho dz and the top layer's factor.
        auto launch_ck_coef = [&]() {
          const int nqc = ck_nq_;
          const int ngr = nch_/nqc;           // chains per thread: the nq angles
          // ck-next: the cell inputs, read (k, j, i) -- strided by n1 across a warp --
          // by every one of the ngr threads of a column, are copied once to a j-fastest
          // View; ck_coef reads the copies (the same numbers: bitwise)
          if (ck_cin_ptr == nullptr) {
            ck_cin_ptr = new DvceArray5D<Real>("ck_cin", nmb1+1, 4, n1, n3, n2);
          }
          auto cin_g = *ck_cin_ptr;
          par_for("ck_coef_in", DevExeSpace(), 0, nmb1, is, ie, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            if (i < icc_g(m,k,j)) return;
            cin_g(m,0,i,k,j) = xT_g(m,k,j,i);
            cin_g(m,1,i,k,j) = xP_g(m,k,j,i);
            cin_g(m,2,i,k,j) = rhoN(m,k,j,i);
            cin_g(m,3,i,k,j) = dx1(m,k,j,i);
          });
          // ck-next: a thread per (pair, cell): no cell depends on another, and a
          // column loop per pair left the kernel waiting on memory (long scoreboard)
          par_for("ck_coef", DevExeSpace(), 0, nmb1, 0, ngr-1, is, ie+1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int p, const int i, const int k,
                        const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int c0 = p*nqc;
            const int blk = c0/RT_NB;
            // the chain kernel zeroes the beam deposit its block's beam kernels add to
            if (c0 % RT_NB == 0) Qb_g(m,blk,i,k,j) = 0.0;
            const int icut = icc_g(m,k,j);
            if (icut > ie || i < icut) return;
            const int g = (nqc == 1) ? (c0 % CK_NG) : ((c0/2) % CK_NG);
            const int b = (nqc == 1) ? (c0/CK_NG) : (c0/(2*CK_NG));
            Real muq[2];
            muq[0] = (nqc == 1) ? 1.0/CK_DIFFUSIVITY : mug[0];
            muq[1] = mug[1];
            if (i == ie+1) {   // the top layer, (1 - e^-dtau) of the column above
              const Real mu0 = cf_g(m,k,j,3);
              const Real ptop = pb_g(m,k,j,ie+1);
              const Real xTv = xT_g(m,k,j,ie+1);
              const Real xPv = xP_g(m,k,j,ie+1);
              const int iT = static_cast<int>(xTv);
              const int iP = static_cast<int>(xPv);
              const Real fT = xTv - static_cast<Real>(iT);
              const Real fP = xPv - static_cast<Real>(iP);
              for (int q=0; q<nqc; ++q) {
                const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, g)
                               + kc_g(m,b,ie+1,k,j);
                const Real dtau = RTTopDtau(kap, ptop*1.0e6,
                                            EffGravAt(grav, ap, x1v_(m,ie+1),
                                                      grav_pmass, omega, mu0, tide)
                                            *cf_g(m,k,j,4));
                const RtF trans = RT_EXP(-static_cast<RtF>(dtau/muq[q]));
                cktpf_g(m,c0+q,k,j) = static_cast<Real>(static_cast<RtF>(1.0)-trans);
                // problem/ck_sph_top = 1: no datum, so the linear kernels' cktpf B = 0
                if (tpm_ == 1) cktpf_g(m,c0+q,k,j) = 0.0;
              }
              return;
            }
            {
              const Real rho = cin_g(m,2,i,k,j);
              const Real dz = cin_g(m,3,i,k,j);
              const Real xTv = cin_g(m,0,i,k,j);
              const Real xPv = cin_g(m,1,i,k,j);
              const int iT = static_cast<int>(xTv);
              const int iP = static_cast<int>(xPv);
              const Real fT = xTv - static_cast<Real>(iT);
              const Real fP = xPv - static_cast<Real>(iP);
              const Real kap = ck_kappa(cklk, iT, fT, iP, fP, b, g) + kc_g(m,b,i,k,j);
              const Real kro = kap*rho;
              for (int q=0; q<nqc; ++q) {
                const int c = c0 + q;
                ckkro_g(m,c,i,k,j) = kro;
                // ck-lin2: e0 by the function the linear kernels re-form it with,
                // and, with the pair (ck-next: CkCinCout), stored only when a reader
                // still reads it (ckcst_)
                if (ckcst_) {
                  // ck_sph_face = 5: the triple at mu_eff = s/sqrt 3 of cell i (the
                  // chain kernel's dthe), for a frozen chain kernel that reads it
                  const Real mqv = vefl_ ? vfs_g(m,c0/nqc,i+1,k,j)/SQRT3_VEF : muq[q];
                  const Real e0 = CkE0(kro, dz, mqv);
                  ckc0_g(m,c,i,k,j) = e0;
                  Real cin, cout;
                  CkCinCout(kro, dz, mqv, e0, cin, cout);
                  ckci_g(m,c,i,k,j) = cin;
                  ckco_g(m,c,i,k,j) = cout;
                }
              }
            }
          });
        };
        auto launch_ck_full = [&]() {
          if (cksph_) {
            if (ckbsph_) {
              launch_ck_tier(std::true_type{}, std::true_type{});
            } else {
              launch_ck_tier(std::true_type{}, std::false_type{});
            }
          } else {
            if (ckbsph_) {
              launch_ck_tier(std::false_type{}, std::true_type{});
            } else {
              launch_ck_tier(std::false_type{}, std::false_type{});
            }
          }
          // ck-determinism: the JAC pass's per-block tridiagonals, added in block order
          if (ckjq_ && cksph_ && ckform_ == 1) {
            par_for("ck_jacq_sum", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
            KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              Real s0 = 0.0, s1 = 0.0, s2 = 0.0;
              for (int b=0; b<nblk; ++b) {
                s0 += ckjq_g(m,3*b,k,j,i);
                s1 += ckjq_g(m,3*b+1,k,j,i);
                s2 += ckjq_g(m,3*b+2,k,j,i);
              }
              ckjac_g(m,0,k,j,i) = s0;
              ckjac_g(m,1,k,j,i) = s1;
              ckjac_g(m,2,k,j,i) = s2;
            });
          }
        };
        {   // ck-tiers: every column kernel is sized by n1 at run time; no n1 cap
          if (cklin_) {
            // ck_impl_lin_check: the chain kernel on the same B first, kept aside
            auto lchk_g = cklchk_ ? *ck_lchk_ptr : CkDum<DvceArray5D<Real>>("lchk_d");
            if (cklchk_) {
              launch_ck_full();
              par_for("ck_lin_chk_cp", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
                for (int i=is; i<ie+2; ++i) {
                  lchk_g(m,blk,i,k,j) = Src_g(m,blk,i,k,j);
                  lchk_g(m,nblk+blk,i,k,j) = Fb_g(m,blk,i,k,j);
                  lchk_g(m,2*nblk+blk,i,k,j) = Em_g(m,blk,i,k,j);
                }
              });
            }
            launch_ck_lin_tier();
            if (cklchk_) {
              DvceArray1D<Real> lmx("ck_lin_chk_mx", 6);
              Kokkos::deep_copy(lmx, 0.0);
              par_for("ck_lin_chk", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
                Real d[3] = {0.0, 0.0, 0.0}, a[3] = {0.0, 0.0, 0.0};
                for (int i=is; i<ie+2; ++i) {
                  const Real v[3] = {Src_g(m,blk,i,k,j), Fb_g(m,blk,i,k,j),
                                     Em_g(m,blk,i,k,j)};
                  for (int q=0; q<3; ++q) {
                    const Real r = lchk_g(m,q*nblk+blk,i,k,j);
                    const Real dd = fabs(v[q] - r);
                    if (dd > d[q] || !(dd == dd)) d[q] = (dd == dd) ? dd : 1.0e300;
                    if (fabs(r) > a[q]) a[q] = fabs(r);
                  }
                }
                for (int q=0; q<3; ++q) {
                  Kokkos::atomic_max(&lmx(q), d[q]);
                  Kokkos::atomic_max(&lmx(3+q), a[q]);
                }
              });
              auto hl = Kokkos::create_mirror_view_and_copy(HostMemSpace(), lmx);
              std::cout << "### ck_lin_check rank=" << global_variable::my_rank
                        << " ncycle=" << pm->ncycle << " pass=" << ck_impl_pass
                        << " src=" << ((hl(3) > 0.0) ? hl(0)/hl(3) : hl(0))
                        << " fb=" << ((hl(4) > 0.0) ? hl(1)/hl(4) : hl(1))
                        << " em=" << ((hl(5) > 0.0) ? hl(2)/hl(5) : hl(2))
                        << std::endl;
            }
          } else if (ckstsp_) {
            launch_ck_coef();                 // ck-store: the fluxes come below
          } else {
            launch_ck_full();
          }
          // ---- problem/ck_beam_par: the pseudo-spherical beam of the storing pass,
          // which rt_chain_ck (tm, FOP) skipped, as three kernels with a thread per
          // TARGET instead of a thread per column: (1) per beam chain the deepest target
          // still under the vertical depth 60 (the slant depth is never smaller, so every
          // target the in-kernel block deposits more than e^-60 into is kept); (2) per
          // (beam chain, face) the slant optical depth of the ray to that face, by the
          // in-kernel block's own chord walk (same shells, same order: the same number);
          // (3) per (chain block, target) the deposit, the block's arithmetic in its
          // order.  Scratch: the lin1 partials ck_lps (free on the storing pass).
          if (cbt_ && ckfst_ && ckbsph_ && cksph_ && ckform_ == 1) {
            const Real alb_ = albedo;
            const Real fst_ = Fstar;
            const int bst_ = (ck_nq_ == 2) ? 2 : 1;
            const int nbc_ = nch_/bst_;
            if (ck_bstop_ptr == nullptr) {
              ck_bstop_ptr = new DvceArray4D<int>("ck_bstop", nmb1+1, nbc_, n3, n2);
            }
            auto bsp_ = *ck_bstop_ptr;
            auto tsl_ = *ck_lps_ptr;         // (m, beam chain, face, k, j)
            par_for("ck_beam_stop", DevExeSpace(), 0, nmb1, 0, nbc_-1, ks, ke, js, je,
            KOKKOS_LAMBDA(const int m, const int p, const int k, const int j) {
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              const int icut = icc_g(m,k,j);
              const int c = p*bst_;
              Real vt = 0.0;
              int fs = icut;
              for (int i=ie; i>=icut; --i) {
                if (vt > 60.0) {
                  fs = i + 1;
                  break;
                }
                vt += static_cast<Real>(static_cast<RtF>(ckkro_g(m,c,i,k,j)))
                      *(X1F(m,i+1) - X1F(m,i));
              }
              bsp_(m,p,k,j) = fs;
            });
            // ck-next: CK_BTG beam chains per thread.  The chord lengths of the ray to
            // face f (the sqrt per shell) depend on the column geometry alone, so one
            // thread forms them once for its group and adds each chain's kappa rho
            // along them, every chain in its own accumulator in the old order: the
            // same slant depth, bitwise.  A chain whose stop face lies above f skips
            // its loads (its slot is not read, as before).
            constexpr int CK_BTG = 8;
            const int nbg_ = (nbc_ + CK_BTG - 1)/CK_BTG;
            // ck-lin2: A TEAM PER (block, beam group, column).  Every face's thread of a
            // column walked the same kappa rho through L2 (the kernel was bound by its
            // loads: L2 hit 98 %, the memory pipes 76 % busy).  The team stages the
            // group's kappa rho (cells icut .. ie, with the RtF cast the walk applies)
            // and the face radii in team scratch once; its threads then walk the faces
            // from there, each face by the old body in the old order: bitwise.  Faces go
            // to threads in snake order (f and the mirrored face of the next round), so
            // the long walks of the low faces are spread over the team.
            const int nbk_ = ke - ks + 1, nbj_ = je - js + 1;
            const int nbf_ = ie - is + 2;                    // faces is .. ie+1
            const size_t btsz_ = ScrArray2D<Real>::shmem_size(CK_BTG, nbf_)
                               + ScrArray1D<Real>::shmem_size(nbf_);
            const int btlv_ = TeamScratchLevel(btsz_, 0);
            const int btlg_ = (nmb1 + 1)*nbg_*nbk_*nbj_;
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
            Kokkos::TeamPolicy<> btpol(DevExeSpace(), btlg_, 128);
#else
            Kokkos::TeamPolicy<> btpol(DevExeSpace(), btlg_, Kokkos::AUTO);
#endif
            Kokkos::parallel_for("ck_beam_tau",
                                 btpol.set_scratch_size(btlv_, Kokkos::PerTeam(btsz_)),
            KOKKOS_LAMBDA(TeamMember_t tm) {
              const int nkj = nbk_*nbj_;
              const int lr = tm.league_rank();
              const int m = lr/(nbg_*nkj);
              const int lq = lr - m*nbg_*nkj;
              const int pg = lq/nkj;
              const int k = (lq - pg*nkj)/nbj_ + ks;
              const int j = (lq - pg*nkj) % nbj_ + js;
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              const int icut = icc_g(m,k,j);
              if (icut > ie) return;
              ScrArray2D<Real> kt(tm.team_scratch(btlv_), CK_BTG, nbf_);
              ScrArray1D<Real> xf(tm.team_scratch(btlv_), nbf_);
              const int ncl = ie + 1 - icut;
              Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ncl + 1),
              [&](const int t) {
                xf(icut - is + t) = X1F(m,icut + t);
              });
              Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, CK_BTG*ncl),
              [&](const int t) {
                const int q = t/ncl;
                const int i = icut + (t - q*ncl);
                const int p = pg*CK_BTG + q;
                if (p < nbc_) {
                  const RtF kr = static_cast<RtF>(ckkro_g(m,p*bst_,i,k,j));
                  kt(q,i - is) = static_cast<Real>(kr);
                }
              });
              int sb[CK_BTG];
              for (int q=0; q<CK_BTG; ++q) {
                const int p = pg*CK_BTG + q;
                sb[q] = (p < nbc_) ? bsp_(m,p,k,j) : (ie + 2);
              }
              const Real mu0 = cf_g(m,k,j,3);
              const Real sinz = sqrt((mu0*mu0 < 1.0) ? (1.0 - mu0*mu0) : 0.0);
              tm.team_barrier();
              const Real rcut = xf(icut - is);
              const int nf = ie + 2 - icut;
              const int ts = tm.team_size(), tr = tm.team_rank();
              for (int r=0; r*ts<nf; ++r) {
                const int o = (r & 1) ? (r*ts + ts - 1 - tr) : (r*ts + tr);
                if (o >= nf) continue;
                const int f = icut + o;
                bool on[CK_BTG];
                bool any = false;
                for (int q=0; q<CK_BTG; ++q) {
                  on[q] = (f >= sb[q]);
                  any = any || on[q];
                }
                if (!any) continue;
                const Real bb = xf(f - is)*sinz;
                const Real b2 = bb*bb;
                int jlo = f;
                bool dark = false;
                if (mu0 < 0.0) {
                  if (bb <= rcut) {
                    dark = true;
                  } else {
                    jlo = icut;
                    for (int jj=f-1; jj>=icut; --jj) {
                      if (xf(jj - is) <= bb) {
                        jlo = jj;
                        break;
                      }
                    }
                  }
                }
                Real tl[CK_BTG];
                for (int q=0; q<CK_BTG; ++q) tl[q] = 0.0;
                if (!dark) {
                  const Real rl = xf(jlo - is);
                  Real prev = (rl*rl > b2) ? sqrt(rl*rl - b2) : 0.0;
                  for (int jj=jlo; jj<ie+1; ++jj) {
                    const Real ru = xf(jj + 1 - is);
                    const Real cur = sqrt(ru*ru - b2);
                    const Real ds = (jj < f) ? 2.0*(cur - prev) : (cur - prev);
                    prev = cur;
                    for (int q=0; q<CK_BTG; ++q) {
                      if (on[q]) tl[q] += ds*kt(q,jj - is);
                    }
                  }
                }
                for (int q=0; q<CK_BTG; ++q) {
                  if (on[q]) tsl_(m,pg*CK_BTG + q,f,k,j) = dark ? 1.0e30 : tl[q];
                }
              }
            });
            par_for("ck_beam_dep", DevExeSpace(), 0, nmb1, 0, nblk-1, is, ie, ks, ke,
                    js, je,
            KOKKOS_LAMBDA(const int m, const int blk, const int i, const int k,
                          const int j) {
              constexpr int NC = RT_NB;
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              const int icut = icc_g(m,k,j);
              if (icut > ie || i < icut) return;
              Real q = 0.0;
              bool any = false;
              for (int cc=0; cc<NC; cc+=bst_) {
                const int c = blk*NC + cc;
                const int p = c/bst_;
                if (i < bsp_(m,p,k,j)) continue;
                int g, b;
                if (ck_nq_ == 1) {
                  g = c % CK_NG;
                  b = c/CK_NG;
                } else {
                  g = (c/2) % CK_NG;
                  b = c/(2*CK_NG);
                }
                const Real wgc = (ck_nq_ == 1) ? ckgw(g)/static_cast<Real>(ck_nq_)
                                               : ckgw(g);
                const Real th = tsl_(m,p,i+1,k,j);
                const Real tl = tsl_(m,p,i,k,j);
                const bool dark = (tl >= 1.0e30);
                const RtF thi = RT_EXP(-static_cast<RtF>(th));
                const Real dtl = dark ? 1.0e30 : (tl - th);
                const RtF tlo = dark ? static_cast<RtF>(0.0)
                                     : RT_EXP(-static_cast<RtF>(tl));
                const Real dif = static_cast<Real>(thi) - static_cast<Real>(tlo);
                const Real fac = (fabs(dtl) > 1.0e-3)
                    ? (dif/dtl) : (static_cast<Real>(thi)*(1.0 - 0.5*dtl));
                q += (1.0-alb_)*fst_*ckswf(b)*wgc*fac
                     *static_cast<Real>(static_cast<RtF>(ckkro_g(m,c,i,k,j)));
                any = true;
              }
              if (any) Qb_g(m,blk,i,k,j) += q;
            });
          }
          // the factorisation, on the storing pass, from what that pass stored
          // ck-jlin: on a pass that also builds the tridiagonal (ck_impl_jac_lin), pass 1
          // of rt_chain_ck_jlin runs here, on the values this kernel has in hand, and
          // parks its window in lpf/lps/lpj (the top one at ie+1) for pass 2
          const bool ckjfus_ = cklbuild_ && ckjlp_;
          if (cklbuild_) {
            ck_jlc_ok = false;       // ck-jlin: a new factorisation, new coefficients
            auto lpjb_g = ckjfus_ ? *ck_lpj_ptr : CkDum<DvceArray5D<Real>>("ck_lpj_d");
            // ck-next: the column geometry, by the expressions the tm body spells, in a
            // kernel of its own before ck_lin_build, which then reads beta and dz from
            // lG (j fastest, coalesced) instead of every chain re-forming beta from the
            // (k, j, i) volume/area Views (strided by n1 across a warp) -- the same
            // numbers, so bitwise.  lG slot 4 is dz itself.
            par_for("ck_lin_geom", DevExeSpace(), 0, nmb1, ks, ke, js, je,
            KOKKOS_LAMBDA(const int m, const int k, const int j) {
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              const int icut = icc_g(m,k,j);
              if (icut > ie) return;
              for (int i=icut; i<ie+2; ++i) {
                lG_g(m,0,i,k,j) = BTF(m,k,j,i,icut);
                if (i <= ie) {
                  lG_g(m,1,i,k,j) = (i == icut)
                      ? AFC(m,k,j,icut)/ACC(m,k,j,icut)
                      : ACC(m,k,j,i-1)/ACC(m,k,j,i);
                  lG_g(m,2,i,k,j) = (i == icut)
                      ? 1.0 : ACC(m,k,j,i-1)/AFC(m,k,j,i);
                  lG_g(m,3,i,k,j) = 1.0/dx1(m,k,j,i);
                  lG_g(m,4,i,k,j) = dx1(m,k,j,i);
                } else {
                  lG_g(m,2,i,k,j) = ACC(m,k,j,ie)/AFC(m,k,j,ie+1);
                }
              }
            });
            if (vef6_) {
              // ck_sph_face = 6: the ladder of the passive network (see ck_sph_face =
              // 6), one thread per (band, g) pair pp = (c0, c0 + 1), Norton from the
              // top: at the upper node of cell i Y_hi = Y_up + g_p; through the series
              // conductance b_s = 1/(g_s + Y_hi), a_s = g_s b_s, Y_lo = a_s Y_hi + g_p;
              // through the contact resistance of face i a_c = 1/(1 + R_c Y_lo), Y_up =
              // a_c Y_lo.  Stored in lP per cell index i: (slot, chain) (0..5, c0) =
              // g_p, K_eq/B, a_s, b_s, a_c, R_c (face i), (0..2, c0+1) = Y_lo, Y_hi,
              // g_s; the top's conductance y_t at (3, c0+1) of ie+1, the wall's
              // source factor cw and conductance y_w at (3, c0+1), (4, c0+1) of icut.
              const int npr6 = nch_/2;
              par_for("ck_f6_build", DevExeSpace(), 0, nmb1, 0, npr6-1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
                if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                const int icut = icc_g(m,k,j);
                const int c0 = 2*pp, c1 = c0 + 1;
                const int b = pp/CK_NG;
                if (m == 0 && k == ks && j == js) {
                  for (int q=0; q<2; ++q) {
                    const Real wfc = 2.0*M_PI*wg[q]*mug[q]*ckgw(pp % CK_NG);
                    lC_g(0,c0+q) = wfc;
                    lC_g(2,c0+q) = 2.0*(wfc/mug[q]);
                    lC_g(1,c0+q) = (int_at_cut ? boltz_sigma/M_PI*Tint4
                                    * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
                  }
                }
                if (icut > ie) return;
                auto P = [&](const int sl, const int c, const int i) -> Real & {
                  return lP_g(m,sl*nch_+c,i,k,j);
                };
                // the top: Auer's H = h J, h from zeta = (sqrt f_t - h)/(sqrt f_t + h),
                // G = y_t K with Phi = 1 at the top face
                const Real ft = SQR(vfs_g(m,pp,ie+2,k,j))/3.0;
                const Real zt = vfg_g(m,pp,ie+2,k,j);
                Real ht = sqrt(ft)*(1.0 - zt)/(1.0 + zt);
                ht = (ht > 1.0e-30) ? ht : 1.0e-30;
                const Real rt = X1F(m,ie+1);
                const Real yt = ht*rt*rt/ft;
                P(3,c1,ie+1) = yt;
                CkF6Cell cc = CkF6Make(X1F(m,ie), rt, ckkro_g(m,c0,ie,k,j),
                                       SQR(vfs_g(m,pp,ie+1,k,j))/3.0, 0.0);
                Real yup = yt;
                for (int i=ie; i>=icut; --i) {
                  const Real yhi = yup + cc.gp;
                  const Real bs = 1.0/(cc.gs + yhi);
                  const Real as = cc.gs*bs;
                  const Real ylo = as*yhi + cc.gp;
                  P(0,c0,i) = cc.gp;
                  P(1,c0,i) = cc.keq;
                  P(2,c0,i) = as;
                  P(3,c0,i) = bs;
                  P(0,c1,i) = ylo;
                  P(1,c1,i) = yhi;
                  P(2,c1,i) = cc.gs;
                  if (i > icut) {
                    const CkF6Cell cd = CkF6Make(X1F(m,i-1), X1F(m,i),
                                                 ckkro_g(m,c0,i-1,k,j),
                                                 SQR(vfs_g(m,pp,i,k,j))/3.0, cc.lpl);
                    const Real rc = CkF6Rc(cd, cc);
                    const Real ac = 1.0/(1.0 + rc*ylo);
                    P(4,c0,i) = ac;
                    P(5,c0,i) = rc;
                    yup = ac*ylo;
                    cc = cd;
                  } else {
                    // the wall: J + beta H = Ib (1/2 + beta/4), as a conductance y_w to
                    // the source potential Phi_w f_w (1/2 + beta/4) Ib
                    const Real fw = SQR(vfs_g(m,pp,icut,k,j))/3.0;
                    const Real bw = vfg_g(m,pp,icut-1,k,j);
                    const Real rw = X1F(m,icut);
                    P(4,c0,i) = 1.0;
                    P(5,c0,i) = 0.0;
                    P(3,c1,icut) = rw*rw*(0.5 + 0.25*bw)/bw;
                    P(4,c1,icut) = rw*rw/(bw*exp(cc.lpl)*fw);
                  }
                }
              });
              if (ck_f6_report > 0) {
                // problem/ck_f6_report: the signs of the elements and the dominance of
                // the node rows (lo_i: g_p + g_s + 1/R_c of face i; hi_i: g_p + g_s +
                // 1/R_c of face i+1, the top's y_t, the wall's y_w)
                DvceArray1D<Real> rr("ck_f6_rep", 4);
                Kokkos::deep_copy(rr, 0.0);
                Kokkos::parallel_for("ck_f6_rep_init",
                                     Kokkos::RangePolicy<DevExeSpace>(0, 1),
                                     KOKKOS_LAMBDA(const int) { rr(3) = 1.0; });
                par_for("ck_f6_rep", DevExeSpace(), 0, nmb1, 0, npr6-1, ks, ke, js, je,
                KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
                  if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                  const int icut = icc_g(m,k,j);
                  if (icut > ie) return;
                  const int c0 = 2*pp, c1 = c0 + 1;
                  auto P = [&](const int sl, const int c, const int i) -> Real {
                    return lP_g(m,sl*nch_+c,i,k,j);
                  };
                  Real nel = 0.0, nneg = 0.0, nnf = 0.0, rmin = 1.0;
                  for (int i=icut; i<=ie; ++i) {
                    const Real gp = P(0,c0,i), kq = P(1,c0,i), gs = P(2,c1,i);
                    const Real rc = (i > icut) ? P(5,c0,i) : 0.0;
                    const Real el[4] = {gp, kq, gs, rc};
                    for (int q=0; q<4; ++q) {
                      nel += 1.0;
                      if (!(el[q] == el[q]) || fabs(el[q]) > 1.0e300) nnf += 1.0;
                      else if (el[q] < 0.0) nneg += 1.0;
                    }
                    const Real cl = (i > icut) ? ((rc > 0.0) ? 1.0/rc : 1.0e300)
                                               : P(4,c1,icut);
                    const Real rcu = (i < ie) ? P(5,c0,i+1) : 0.0;
                    const Real cu = (i < ie) ? ((rcu > 0.0) ? 1.0/rcu : 1.0e300)
                                             : P(3,c1,ie+1);
                    const Real xl = (i > icut) ? gp : gp + cl;
                    const Real xh = (i < ie) ? gp : gp + cu;
                    const Real rl = xl/(gp + gs + cl), rh = xh/(gp + gs + cu);
                    rmin = fmin(rmin, fmin(rl, rh));
                  }
                  Kokkos::atomic_add(&rr(0), nel);
                  Kokkos::atomic_add(&rr(1), nneg);
                  Kokkos::atomic_add(&rr(2), nnf);
                  Kokkos::atomic_min(&rr(3), rmin);
                });
                auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rr);
                if (global_variable::my_rank == 0) {
                  std::printf("### f6rep ncycle=%d elements=%.0f negative=%.0f "
                              "nonfinite=%.0f min(row excess/diag)=%.3e\n",
                              static_cast<int>(pm->ncycle), hr(0), hr(1), hr(2), hr(3));
                }
              }
            } else if (vefl_) {   // ck_sph_face = 5: the VEF factorisation
              par_for("ck_lin_build", DevExeSpace(), 0, nmb1, 0, nch_-1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
                if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                const int icut = icc_g(m,k,j);
                // the per-chain constants: the flux weight and the internal-flux datum of
                // the chain's band at the cut (T-independent), written once per chain
                int b;
                Real wfc, muc;
                if (ck_nq_ == 1) {
                  b = c/CK_NG;
                  muc = 1.0/CK_DIFFUSIVITY;
                  wfc = M_PI*ckgw(c % CK_NG);
                } else {
                  const int nq = c % 2;
                  b = c/(2*CK_NG);
                  muc = mug[nq];
                  wfc = 2.0*M_PI*wg[nq]*mug[nq]*ckgw((c/2) % CK_NG);
                }
                if (m == 0 && k == ks && j == js) {
                  lC_g(0,c) = wfc;
                  lC_g(2,c) = 2.0*(wfc/muc);   // ck-next: the emission weight (CkEmW)
                  lC_g(1,c) = (int_at_cut ? boltz_sigma/M_PI*Tint4
                              * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
                }
                if (icut > ie) return;
                // the tm Moebius map of pass 1, on R alone (FRM = 1 of rt_chain_ck);
                // ck_dif_dtau: R = 1 at a handover column's cut (the flux datum)
                const Real gdf = ckdif_ ? difg_g(m,c,k,j) : -1.0;
                RtF rr = static_cast<RtF>(0.0);
                if (ckdif_ && gdf >= 0.0) rr = static_cast<RtF>(1.0);
                CkJlP1 w = CkJlP1Init(gdf);
                // ck_sph_face = 5 (ck_nquad = 2: the pair c/2 is the (band, g) of the
                // formal solution): the Marshak-type wall, R0 = -(1 - beta mu_eff)/(1 +
                // beta mu_eff), dSc_cut/dB_cut = (1 + beta/2)/(1 + beta mu_eff)
                const int cv = c/2;
                if (vefl_) {
                  const Real bwb = vfg_g(m,cv,icut-1,k,j);
                  const Real mw = vfs_g(m,cv,icut,k,j)/SQRT3_VEF;
                  rr = static_cast<RtF>(-(1.0 - bwb*mw)/(1.0 + bwb*mw));
                  w.s2 = (1.0 + 0.5*bwb)/(1.0 + bwb*mw);
                }
                // ck-next: the loads of cell i+1 are issued before cell i is worked
                // (the stores below keep the compiler from hoisting them itself)
                Real btr_n = lG_g(m,0,icut,k,j), dz_n = lG_g(m,4,icut,k,j);
                // ck-lin2: cin and cout are formed (CkCinCout), not read
                // ck-lin2: and so is e0 (CkE0)
                Real kro_n = ckkro_g(m,c,icut,k,j);
                for (int i=icut; i<ie+1; ++i) {
                  const Real btr = btr_n, dz = dz_n;
                  const Real kro = kro_n;
                  // ck_sph_face = 5: cell i at mu_eff = s/sqrt 3; rp and gamma of face i
                  Real mue = muc, vrp = 1.0, vgm = 1.0;
                  if (vefl_) {
                    const Real vsa = vfs_g(m,cv,i+1,k,j);
                    mue = vsa/SQRT3_VEF;
                    vrp = (1.0 - btr)/(1.0 + btr)*vfs_g(m,cv,i,k,j)/vsa;
                    vgm = vfg_g(m,cv,i,k,j);
                  }
                  const Real e0 = CkE0(kro, dz, mue);
                  if (i < ie) {
                    btr_n = lG_g(m,0,i+1,k,j);
                    dz_n = lG_g(m,4,i+1,k,j);
                    kro_n = ckkro_g(m,c,i+1,k,j);
                  }
                  const RtF bb = static_cast<RtF>(btr);
                  const RtF dn = static_cast<RtF>(1.0) + rr*bb;
                  const Real rj = static_cast<Real>(rr);
                  const Real idn = CkIdn(rj, btr);   // ck-next: as the readers form it
                  lP_g(m,0*nch_+c,i,k,j) = rj;
                  RtF rn = (rr + bb)/dn;
                  Real vr1 = 0.0, vc1 = 0.0;
                  if (vefl_) {
                    CkVefMap(rj, vrp, vgm, vr1, vc1);
                    rn = static_cast<RtF>(vr1);
                  }
                  // ck-jlin: the half-layer triple is NOT copied (ck-store: lP has no
                  // slots for it any more): the linear kernels read ckc0/ckci/ckco,
                  // which this pass stored and which stay frozen with the rest of the
                  // factorisation
                  const RtF tr = static_cast<RtF>(1.0) - static_cast<RtF>(e0);
                  rr = tr*tr*rn;
                  rr = tr*tr*rr;
                  // ck-nq2: slots 2-4 depend on kappa rho alone, which is bitwise the
                  // same for the two angles of a pair; the q = 1 chain reads the q = 0
                  // one
                  // (ck-jlin: and forms the same numbers itself for a fused pass 1)
                  // ck-lin2: the paired linear kernel and rt_chain_ck_jlin form them
                  // (CkLayW); stored only for the unpaired tiers (!cklp_)
                  const bool own = !(ck_nq_ == 2 && (c & 1)) && !cklp_;
                  Real wl = 0.0, wu = 0.0, ffj = 0.0;
                  if (i < ie && (own || ckjfus_)) {
                    CkLayW(kro, kro_n, dz, dz_n, bface_on, wl, wu, ffj);
                    if (own) {
                      lP_g(m,2*nch_+c,i,k,j) = wl;
                      lP_g(m,3*nch_+c,i,k,j) = wu;
                      lP_g(m,4*nch_+c,i,k,j) = ffj;
                    }
                  }
                  if (ckjfus_) {
                    lpf_g(m,c,i,k,j) = w.s0;
                    lps_g(m,c,i,k,j) = w.s1;
                    lpjb_g(m,c,i,k,j) = w.s2;
                    Real ci, co;
                    CkCinCout(kro, dz, mue, e0, ci, co);
                    if (vefl_) {
                      CkJlP1StepRC(w, i < ie, wl, wu, ffj, ci, co, e0, vr1, vc1);
                    } else {
                      CkJlP1Step(w, i < ie, wl, wu, ffj, ci, co, e0, rj, idn, btr);
                    }
                  }
                }
                if (ckjfus_) {
                  lpf_g(m,c,ie+1,k,j) = w.s0;
                  lps_g(m,c,ie+1,k,j) = w.s1;
                  lpjb_g(m,c,ie+1,k,j) = w.s2;
                }
                lP_g(m,0*nch_+c,ie+1,k,j) = static_cast<Real>(rr);
              });
            } else {
              par_for("ck_lin_build", DevExeSpace(), 0, nmb1, 0, nch_-1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
                if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                const int icut = icc_g(m,k,j);
                // the per-chain constants: the flux weight and the internal-flux datum of
                // the chain's band at the cut (T-independent), written once per chain
                int b;
                Real wfc, muc;
                if (ck_nq_ == 1) {
                  b = c/CK_NG;
                  muc = 1.0/CK_DIFFUSIVITY;
                  wfc = M_PI*ckgw(c % CK_NG);
                } else {
                  const int nq = c % 2;
                  b = c/(2*CK_NG);
                  muc = mug[nq];
                  wfc = 2.0*M_PI*wg[nq]*mug[nq]*ckgw((c/2) % CK_NG);
                }
                if (m == 0 && k == ks && j == js) {
                  lC_g(0,c) = wfc;
                  lC_g(2,c) = 2.0*(wfc/muc);   // ck-next: the emission weight (CkEmW)
                  lC_g(1,c) = (int_at_cut ? boltz_sigma/M_PI*Tint4
                              * ck_planck_frac(ckpf, pfl0, pfid, Tint, b) : 0.0);
                }
                if (icut > ie) return;
                // the tm Moebius map of pass 1, on R alone (FRM = 1 of rt_chain_ck);
                // ck_dif_dtau: R = 1 at a handover column's cut (the flux datum)
                const Real gdf = ckdif_ ? difg_g(m,c,k,j) : -1.0;
                RtF rr = static_cast<RtF>(0.0);
                if (ckdif_ && gdf >= 0.0) rr = static_cast<RtF>(1.0);
                CkJlP1 w = CkJlP1Init(gdf);
                // ck-next: the loads of cell i+1 are issued before cell i is worked
                // (the stores below keep the compiler from hoisting them itself)
                Real btr_n = lG_g(m,0,icut,k,j), dz_n = lG_g(m,4,icut,k,j);
                // ck-lin2: cin and cout are formed (CkCinCout), not read
                // ck-lin2: and so is e0 (CkE0)
                Real kro_n = ckkro_g(m,c,icut,k,j);
                for (int i=icut; i<ie+1; ++i) {
                  const Real btr = btr_n, dz = dz_n;
                  const Real kro = kro_n;
                  const Real e0 = CkE0(kro, dz, muc);
                  if (i < ie) {
                    btr_n = lG_g(m,0,i+1,k,j);
                    dz_n = lG_g(m,4,i+1,k,j);
                    kro_n = ckkro_g(m,c,i+1,k,j);
                  }
                  const RtF bb = static_cast<RtF>(btr);
                  const RtF dn = static_cast<RtF>(1.0) + rr*bb;
                  const Real rj = static_cast<Real>(rr);
                  const Real idn = CkIdn(rj, btr);   // ck-next: as the readers form it
                  lP_g(m,0*nch_+c,i,k,j) = rj;
                  const RtF rn = (rr + bb)/dn;
                  // ck-jlin: the half-layer triple is NOT copied (ck-store: lP has no
                  // slots for it any more): the linear kernels read ckc0/ckci/ckco,
                  // which this pass stored and which stay frozen with the rest of the
                  // factorisation
                  const RtF tr = static_cast<RtF>(1.0) - static_cast<RtF>(e0);
                  rr = tr*tr*rn;
                  rr = tr*tr*rr;
                  // ck-nq2: slots 2-4 depend on kappa rho alone, which is bitwise the
                  // same for the two angles of a pair; the q = 1 chain reads the q = 0
                  // one
                  // (ck-jlin: and forms the same numbers itself for a fused pass 1)
                  // ck-lin2: the paired linear kernel and rt_chain_ck_jlin form them
                  // (CkLayW); stored only for the unpaired tiers (!cklp_)
                  const bool own = !(ck_nq_ == 2 && (c & 1)) && !cklp_;
                  Real wl = 0.0, wu = 0.0, ffj = 0.0;
                  if (i < ie && (own || ckjfus_)) {
                    CkLayW(kro, kro_n, dz, dz_n, bface_on, wl, wu, ffj);
                    if (own) {
                      lP_g(m,2*nch_+c,i,k,j) = wl;
                      lP_g(m,3*nch_+c,i,k,j) = wu;
                      lP_g(m,4*nch_+c,i,k,j) = ffj;
                    }
                  }
                  if (ckjfus_) {
                    lpf_g(m,c,i,k,j) = w.s0;
                    lps_g(m,c,i,k,j) = w.s1;
                    lpjb_g(m,c,i,k,j) = w.s2;
                    Real ci, co;
                    CkCinCout(kro, dz, muc, e0, ci, co);
                    CkJlP1Step(w, i < ie, wl, wu, ffj, ci, co, e0, rj, idn, btr);
                  }
                }
                if (ckjfus_) {
                  lpf_g(m,c,ie+1,k,j) = w.s0;
                  lps_g(m,c,ie+1,k,j) = w.s1;
                  lpjb_g(m,c,ie+1,k,j) = w.s2;
                }
                lP_g(m,0*nch_+c,ie+1,k,j) = static_cast<Real>(rr);
              });
            }
          }
          // ---- problem/ck_impl_jac_lin: THE TRIDIAGONAL FROM THE FACTORISATION ------
          // The JAC recurrences of the tm body (see the JAC blocks in rt_chain_ck), term
          // for term, on what ck_lin_build parked: pass 1 carries dSc_f/dB_{f-2,f-1,f}
          // up and parks it at every face, pass 2 carries dd^+/dB_{i-1,i,i+1} and the
          // reflection scalar W down and forms row i.  One thread per chain; each chain
          // writes its three entries (times dB_b/dT, negative off-diagonal parts dropped
          // as in the JAC assembly) to per-chain partials, and ck_jlin_sum adds them in
          // chain order.  The dSc window parked at each face lives in the same three
          // partials (pass 2 reads face i's window just before it writes row i there),
          // so the kernel carries no per-thread column and no scratch.
          // ck-jlin: none of this depends on B or T -- only the final factor dB_b/dT
          // does -- so the kernel stores the three coefficients BEFORE that factor
          // (ljm, lj0, lpj) and runs once per stored factorisation (ck_jlc_ok, cleared
          // by ck_lin_build); ck_jlin_sum multiplies by dB/dT (a lone rounded product,
          // CkMulRn: bitwise the product the kernel used to store) and adds in chain
          // order.  Every column is swept, done or not: a column that is done at this
          // pass is live again at the next call's, which may re-apply the operator.
          if (ckjlp_) {
            auto jac_g = ckjac_g;
            auto db_g = ckdb_g;
            auto lpj_g = *ck_lpj_ptr;
            auto ljm_g = *ck_ljm_ptr;
            auto lj0_g = *ck_lj0_ptr;
            const bool jneg_ = ck_impl_jneg;
            const bool jnet_ = ck_impl_jnet;   // keep a negative NET off-diagonal
            if (ckdif_) ck_jlc_ok = false;   // the flux datum is not cached
            if (!ck_jlc_ok) {
              const bool jfus_ = ckjfus_;
              if (vef6_) {
                // ck_sph_face = 6: the exact (i-1, i, i+1) part of dSrc_i/dB of the
                // ladder (ck_f6_build).  A unit current into node m raises node m by
                // 1/D_m (D = Y looking up + Y looking down, the node's shunt once) and
                // every other node by the product of the voltage ratios of the elements
                // in between (up: a_s, a_c of the build; down: dn_s = g_s/(g_s + Yd_lo
                // + g_p), dn_c = 1/(1 + R_c (Yd_hi + g_p)), Yd from a sweep up from the
                // wall).  Cell j's B drives its two face nodes with g_p K_eq/B each (the
                // cut cell's also the wall, cw); dep_i = g_p (K_lo + K_hi) - 2 g_p K_eq.
                // Every factor is >= 0: the off-diagonals are >= 0 and the diagonal
                // <= 0 (passive network, cell-local emission).  Pair rows in the first
                // chain of the pair, zeros in the second.
                const int npr6 = nch_/2;
                CkParFor4("rt_chain_ck_jlin6", cklw_, 0, nmb1, 0, npr6-1, ks, ke, js, je,
                KOKKOS_LAMBDA(const int m, const int pp, const int k, const int j) {
                  // ck_f6_build skipped the columns done on its pass
                  if (jfus_ && ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                  const int icut = icc_g(m,k,j);
                  if (icut > ie) return;
                  const int c0 = 2*pp, c1 = c0 + 1;
                  auto P = [&](const int sl, const int c, const int i) -> Real {
                    return lP_g(m,sl*nch_+c,i,k,j);
                  };
                  const Real cw = P(3,c1,icut);
                  const Real cnv = (lC_g(2,c0) + lC_g(2,c1))
                                 *AFC(m,k,j,ie+1)/SQR(X1F(m,ie+1));
                  // cell i (current) and its lower neighbour's driving points
                  Real gp = P(0,c0,icut), sr = gp*P(1,c0,icut), as = P(2,c0,icut);
                  Real gs = P(2,c1,icut);
                  Real ydl = P(4,c1,icut);                       // the wall
                  Real ydh = gs*(ydl + gp)/(gs + ydl + gp);
                  Real dlo = P(0,c1,icut) + ydl, dhi = P(1,c1,icut) + ydh;
                  Real dns = gs/(gs + ydl + gp);
                  Real dlop = 1.0, dhip = 1.0, asp = 0.0, srp = 0.0;
                  for (int i=icut; i<=ie; ++i) {
                    const Real ac = P(4,c0,i);
                    // look-ahead: cell i+1
                    Real gpn = 0.0, srn = 0.0, asn = 0.0, ydln = 0.0, ydhn = 0.0;
                    Real dlon = 1.0, dhin = 1.0, dncn = 0.0, dnsn = 0.0;
                    if (i < ie) {
                      gpn = P(0,c0,i+1);
                      srn = gpn*P(1,c0,i+1);
                      asn = P(2,c0,i+1);
                      const Real gsn = P(2,c1,i+1);
                      const Real rcn = P(5,c0,i+1);
                      const Real yb = ydh + gp;
                      dncn = 1.0/(1.0 + rcn*yb);
                      ydln = yb*dncn;
                      dnsn = gsn/(gsn + ydln + gpn);
                      ydhn = (ydln + gpn)*dnsn;
                      dlon = P(0,c1,i+1) + ydln;
                      dhin = P(1,c1,i+1) + ydhn;
                    }
                    // own: sources at lo_i and hi_i (and the wall for the cut cell)
                    const Real wl = (i == icut) ? cw : 0.0;
                    const Real klo0 = (sr + wl)/dlo + sr*dns/dhi;
                    const Real khi0 = (sr + wl)*as/dlo + sr/dhi;
                    // below: cell i-1 (and the wall through lo_{i-1} if it is the cut)
                    Real klom = 0.0, khim = 0.0;
                    if (i > icut) {
                      const Real wlm = (i - 1 == icut) ? cw : 0.0;
                      klom = ac*(srp/dhip + (srp + wlm)*asp/dlop);
                      khim = as*klom;
                    }
                    // above: cell i+1
                    Real klop = 0.0, khip = 0.0;
                    if (i < ie) {
                      khip = dncn*(srn/dlon + srn*dnsn/dhin);
                      klop = dns*khip;
                    }
                    const Real cv = cnv/VLS(m,k,j,i);
                    lj0_g(m,c0,i,k,j) = cv*(gp*(klo0 + khi0) - 2.0*sr);
                    lpj_g(m,c0,i,k,j) = cv*gp*(klop + khip);
                    ljm_g(m,c0,i,k,j) = cv*gp*(klom + khim);
                    lj0_g(m,c1,i,k,j) = 0.0;
                    lpj_g(m,c1,i,k,j) = 0.0;
                    ljm_g(m,c1,i,k,j) = 0.0;
                    // shift
                    dlop = dlo;
                    dhip = dhi;
                    asp = as;
                    srp = sr;
                    gp = gpn;
                    sr = srn;
                    as = asn;
                    ydh = ydhn;
                    dlo = dlon;
                    dhi = dhin;
                    dns = dnsn;
                  }
                });
              } else if (vefl_) {   // ck_sph_face = 5: the VEF rows
                CkParFor4("rt_chain_ck_jlin", cklw_, 0, nmb1, 0, nch_-1, ks, ke, js, je,
                KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
                  // fused: ck_lin_build left the columns it skipped without a window
                  if (jfus_ && ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                  const int icut = icc_g(m,k,j);
                  if (icut > ie) return;
                  const Real wfc = lC_g(0,c);
                  // PASS 1: dSc/dB in the window (B_{i-2}, B_{i-1}, B_i) below face i
                  // (ck-jlin: parked by ck_lin_build when it ran on this pass, with the
                  // top window at ie+1)
                  // ck_dif_dtau: the flux datum g_c (B_{cut-1} - B_cut) of a handover
                  // column, whose row at the cut then couples to the deep cell below
                  const Real gdf = ckdif_ ? difg_g(m,c,k,j) : -1.0;
                  // ck_sph_face = 5: the pair c/2 is the formal solution's (band, g)
                  const int cv = c/2;
                  Real jS0, jS1, jS2;
                  if (jfus_) {
                    jS0 = lpf_g(m,c,ie+1,k,j);
                    jS1 = lps_g(m,c,ie+1,k,j);
                    jS2 = lpj_g(m,c,ie+1,k,j);
                  } else {
                    CkJlP1 w = CkJlP1Init(gdf);
                    if (vefl_) {             // ck_sph_face = 5: the Marshak-type wall
                      const Real bwb = vfg_g(m,cv,icut-1,k,j);
                      const Real mw = vfs_g(m,cv,icut,k,j)/SQRT3_VEF;
                      w.s2 = (1.0 + 0.5*bwb)/(1.0 + bwb*mw);
                    }
                    for (int i=icut; i<ie+1; ++i) {
                      lpf_g(m,c,i,k,j) = w.s0;
                      lps_g(m,c,i,k,j) = w.s1;
                      lpj_g(m,c,i,k,j) = w.s2;
                      const bool up = (i < ie);
                      // ck-lin2: the layer weights formed (CkLayW), not read
                      Real wl = 0.0, wu = 0.0, ffj = 0.0;
                      if (up) CkLayW(ckkro_g(m,c,i,k,j), ckkro_g(m,c,i+1,k,j),
                                     lG_g(m,4,i,k,j), lG_g(m,4,i+1,k,j), bface_on,
                                     wl, wu, ffj);
                      // ck-lin2: cin and cout formed (CkCinCout), not read
                      Real mu1 = (ck_nq_ == 1) ? 1.0/CK_DIFFUSIVITY : mug[c % 2];
                      if (vefl_) mu1 = vfs_g(m,cv,i+1,k,j)/SQRT3_VEF;
                      const Real e01 = CkE0(ckkro_g(m,c,i,k,j), lG_g(m,4,i,k,j), mu1);
                      Real ci1, co1;
                      CkCinCout(ckkro_g(m,c,i,k,j), lG_g(m,4,i,k,j), mu1, e01, ci1, co1);
                      if (vefl_) {
                        const Real btv = lG_g(m,0,i,k,j);
                        Real vr1, vc1;
                        CkVefMap(lP_g(m,0*nch_+c,i,k,j),
                                 (1.0 - btv)/(1.0 + btv)*vfs_g(m,cv,i,k,j)
                                 /vfs_g(m,cv,i+1,k,j), vfg_g(m,cv,i,k,j), vr1, vc1);
                        CkJlP1StepRC(w, up, wl, wu, ffj, ci1, co1, e01, vr1, vc1);
                      } else {
                        CkJlP1Step(w, up, wl, wu, ffj, ci1, co1, e01,
                                   lP_g(m,0*nch_+c,i,k,j),
                                   CkIdn(lP_g(m,0*nch_+c,i,k,j), lG_g(m,0,i,k,j)),
                                   lG_g(m,0,i,k,j));
                      }
                    }
                    jS0 = w.s0;
                    jS1 = w.s1;
                    jS2 = w.s2;
                  }
                  // PASS 2: the top datum is fixed, so d^+ carries nothing and W = 0
                  Real jQ, jDu0, jDu1, jDu2, jdm0, jdm1, jdm2;
                  Real jfu0 = 1.0, jfu1 = 0.0, jlc0 = 1.0, jlc1 = 0.0;
                  {
                    const Real bt = lG_g(m,0,ie+1,k,j);
                    const Real idn = CkIdn(lP_g(m,0*nch_+c,ie+1,k,j), lG_g(m,0,ie+1,k,j));
                    Real al = (1.0 + bt)*idn;
                    Real be = bt*idn;
                    if (vefl_) {
                      // ck_sph_face = 5: the reflected datum d_a = (zeta Sc' + datum)/
                      // (1 - zeta R') depends on Sc_top (the chain kernel's JAC top
                      // block):
                      // dd_b/dSc = al kp - be, dD/dSc = aD (1 - gamma (1 - R) kp)
                      const Real rj = lP_g(m,0*nch_+c,ie+1,k,j);
                      const Real vsb = vfs_g(m,cv,ie+1,k,j);
                      const Real vrp = (1.0 - bt)/(1.0 + bt)*vsb/vfs_g(m,cv,ie+2,k,j);
                      const Real vgm = vfg_g(m,cv,ie+1,k,j);
                      const Real zt = vfg_g(m,cv,ie+2,k,j);
                      const Real ev = (1.0 + rj) + vgm*vrp*(1.0 - rj);
                      const Real rup = 1.0 - 2.0*vgm*vrp*(1.0 - rj)/ev;
                      const Real vkp = zt*(2.0*vrp/ev)/(1.0 - zt*rup);
                      al = 2.0*vgm/ev;
                      be = (1.0 - vgm*vrp)/ev - al*vkp;
                      const Real aD = 2.0/ev*(1.0 - vgm*(1.0 - rj)*vkp);
                      jDu0 = aD*jS0;
                      jDu1 = aD*jS1;
                      jDu2 = aD*jS2;
                    } else {
                      jDu0 = al*jS0;
                      jDu1 = al*jS1;
                      jDu2 = al*jS2;
                    }
                    jdm0 = -be*jS0;
                    jdm1 = -be*jS1;
                    jdm2 = -be*jS2;
                    jQ = be;
                  }
                  // ck-next: cell i-1's loads are issued before cell i is worked
                  // ck-next: cin, cout and 1/(1 + R beta) formed (CkCinCout, CkIdn)
                  const Real muj = (ck_nq_ == 1) ? 1.0/CK_DIFFUSIVITY : mug[c % 2];
                  Real kr_p = ckkro_g(m,c,ie,k,j), dz_p = lG_g(m,4,ie,k,j);
                  Real rj_p = lP_g(m,0*nch_+c,ie,k,j);   // ck-lin2: e0 formed (CkE0)
                  Real bt_p = lG_g(m,0,ie,k,j);
                  Real ra_p = lG_g(m,1,ie,k,j), g3_p = lG_g(m,3,ie,k,j);
                  Real s0_p = lpf_g(m,c,ie,k,j), s1_p = lps_g(m,c,ie,k,j);
                  Real s2_p = lpj_g(m,c,ie,k,j);
                  // ck-lin2: the layer (i-1, i) weights formed (CkLayW): kappa rho and
                  // dz of cell i-1 loaded two cells ahead
                  Real kr_l = 0.0, dz_l = 0.0;
                  if (ie > icut) {
                    kr_l = ckkro_g(m,c,ie-1,k,j);
                    dz_l = lG_g(m,4,ie-1,k,j);
                  }
                  for (int i=ie; i>icut-1; --i) {
                    // ck_sph_face = 5: cell i at mu_eff, s below face i, gamma of face i
                    Real mui = muj, vsb = 1.0, vsa = 1.0, vgm = 1.0;
                    if (vefl_) {
                      vsb = vfs_g(m,cv,i,k,j);
                      vsa = vfs_g(m,cv,i+1,k,j);
                      vgm = vfg_g(m,cv,i,k,j);
                      mui = vsa/SQRT3_VEF;
                    }
                    const Real c0_p = CkE0(kr_p, dz_p, mui);
                    Real ci, co;
                    CkCinCout(kr_p, dz_p, mui, c0_p, ci, co);
                    const Real tj = 1.0 - c0_p;
                    const Real rj = rj_p, bt = bt_p, rat = ra_p, g3 = g3_p;
                    const Real idn = CkIdn(rj, bt);
                    const Real s0 = s0_p, s1 = s1_p, s2 = s2_p;
                    Real wc = 0.0, wa = 0.0, ffj = 0.0;
                    if (i > icut) CkLayW(kr_l, kr_p, dz_l, dz_p, bface_on, wc, wa, ffj);
                    if (i > icut) {
                      kr_p = kr_l;
                      dz_p = dz_l;
                      rj_p = lP_g(m,0*nch_+c,i-1,k,j);
                      bt_p = lG_g(m,0,i-1,k,j);
                      ra_p = lG_g(m,1,i-1,k,j);
                      g3_p = lG_g(m,3,i-1,k,j);
                      s0_p = lpf_g(m,c,i-1,k,j);
                      s1_p = lps_g(m,c,i-1,k,j);
                      s2_p = lpj_g(m,c,i-1,k,j);
                      if (i-1 > icut) {
                        kr_l = ckkro_g(m,c,i-2,k,j);
                        dz_l = lG_g(m,4,i-2,k,j);
                      }
                    }
                    // the layer joining cells i-1 and i: slots (B_{i-1}, B_i); the third
                    // (B_{i+1}) is 0 for all three
                    Real duv0 = 0.0, duv1 = 1.0, dnl0 = 0.0, dnl1 = 1.0;
                    Real dfw0 = 0.0, dfw1 = 1.0;
                    if (i > icut) {
                      const Real sw = (1.0 - ffj)*wc + ffj*(1.0 - wa);
                      duv0 = 1.0 - wa;
                      duv1 = wa;
                      dnl0 = wc;
                      dnl1 = 1.0 - wc;
                      dfw0 = sw;
                      dfw1 = 1.0 - sw;
                    }
                    // dd_i^+: the down ray's two half steps through cell i
                    Real v0 = tj*jdm0;
                    v0 = tj*v0 + ci*duv0 + co*dfw0;
                    Real v1 = tj*jdm1 + ci*jfu0 + co*jlc0;
                    v1 = tj*v1 + ci*duv1 + co*dfw1;
                    Real v2 = tj*jdm2 + ci*jfu1 + co*jlc1;
                    v2 = tj*v2;
                    Real al = (1.0 + bt)*idn;
                    Real be = bt*idn;
                    Real wr = tj*tj*tj*tj*(1.0 - bt)*idn*jQ;
                    Real omr = 1.0 - rj;
                    Real aD = al, ratv = rat, wfi = wfc;
                    if (vefl_) {
                      // ck_sph_face = 5: D = (2/E) (Sc - gamma (1 - R) d^+), d crosses
                      // the face with al = 2 gamma/E, be = (1 - gamma rho')/E
                      ratv = rat*vsb/vsa;
                      const Real ev = (1.0 + rj) + vgm*ratv*(1.0 - rj);
                      al = 2.0*vgm/ev;
                      be = (1.0 - vgm*ratv)/ev;
                      aD = 2.0/ev;
                      wr = tj*tj*tj*tj*(2.0*ratv/ev)*jQ;
                      omr = vgm*(1.0 - rj);
                      wfi = wfc*vsa/(SQRT3_VEF*muj);
                    }
                    const Real dd0 = aD*(1.0 + omr*wr)*s0;
                    const Real dd1 = aD*(s1 - omr*v0);
                    const Real dd2 = aD*(s2 - omr*v1);
                    const Real dd3 = -aD*omr*v2;
                    const Real wj = wfi*g3;
                    const Real jm = wj*(ratv*dd1 - jDu0);
                    const Real j0 = wj*(ratv*dd2 - jDu1);
                    const Real jp = wj*(ratv*dd3 - jDu2);
                    lj0_g(m,c,i,k,j) = j0;
                    ljm_g(m,c,i,k,j) = ((i > icut || gdf >= 0.0) && (jneg_ || jm > 0.0))
                                     ? jm : 0.0;
                    lpj_g(m,c,i,k,j) = (i < ie && (jneg_ || jp > 0.0)) ? jp : 0.0;
                    jDu0 = dd0;
                    jDu1 = dd1;
                    jDu2 = dd2;
                    const Real qn = al*wr + be;
                    jdm0 = -qn*s0;
                    jdm1 = al*v0 - be*s1;
                    jdm2 = al*v1 - be*s2;
                    jQ = qn;
                    jfu0 = dfw0;
                    jfu1 = dfw1;
                    jlc0 = dnl0;
                    jlc1 = dnl1;
                  }
                });
              } else {
                CkParFor4("rt_chain_ck_jlin", cklw_, 0, nmb1, 0, nch_-1, ks, ke, js, je,
                KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
                  // fused: ck_lin_build left the columns it skipped without a window
                  if (jfus_ && ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                  const int icut = icc_g(m,k,j);
                  if (icut > ie) return;
                  const Real wfc = lC_g(0,c);
                  // PASS 1: dSc/dB in the window (B_{i-2}, B_{i-1}, B_i) below face i
                  // (ck-jlin: parked by ck_lin_build when it ran on this pass, with the
                  // top window at ie+1)
                  // ck_dif_dtau: the flux datum g_c (B_{cut-1} - B_cut) of a handover
                  // column, whose row at the cut then couples to the deep cell below
                  const Real gdf = ckdif_ ? difg_g(m,c,k,j) : -1.0;
                  Real jS0, jS1, jS2;
                  if (jfus_) {
                    jS0 = lpf_g(m,c,ie+1,k,j);
                    jS1 = lps_g(m,c,ie+1,k,j);
                    jS2 = lpj_g(m,c,ie+1,k,j);
                  } else {
                    CkJlP1 w = CkJlP1Init(gdf);
                    for (int i=icut; i<ie+1; ++i) {
                      lpf_g(m,c,i,k,j) = w.s0;
                      lps_g(m,c,i,k,j) = w.s1;
                      lpj_g(m,c,i,k,j) = w.s2;
                      const bool up = (i < ie);
                      // ck-lin2: the layer weights formed (CkLayW), not read
                      Real wl = 0.0, wu = 0.0, ffj = 0.0;
                      if (up) CkLayW(ckkro_g(m,c,i,k,j), ckkro_g(m,c,i+1,k,j),
                                     lG_g(m,4,i,k,j), lG_g(m,4,i+1,k,j), bface_on,
                                     wl, wu, ffj);
                      // ck-lin2: cin and cout formed (CkCinCout), not read
                      const Real mu1 = (ck_nq_ == 1) ? 1.0/CK_DIFFUSIVITY : mug[c % 2];
                      const Real e01 = CkE0(ckkro_g(m,c,i,k,j), lG_g(m,4,i,k,j), mu1);
                      Real ci1, co1;
                      CkCinCout(ckkro_g(m,c,i,k,j), lG_g(m,4,i,k,j), mu1, e01, ci1, co1);
                      CkJlP1Step(w, up, wl, wu, ffj, ci1, co1, e01,
                                 lP_g(m,0*nch_+c,i,k,j),
                                 CkIdn(lP_g(m,0*nch_+c,i,k,j), lG_g(m,0,i,k,j)),
                                 lG_g(m,0,i,k,j));
                    }
                    jS0 = w.s0;
                    jS1 = w.s1;
                    jS2 = w.s2;
                  }
                  // PASS 2: the top datum is fixed, so d^+ carries nothing and W = 0
                  Real jQ, jDu0, jDu1, jDu2, jdm0, jdm1, jdm2;
                  Real jfu0 = 1.0, jfu1 = 0.0, jlc0 = 1.0, jlc1 = 0.0;
                  {
                    const Real bt = lG_g(m,0,ie+1,k,j);
                    const Real idn = CkIdn(lP_g(m,0*nch_+c,ie+1,k,j), lG_g(m,0,ie+1,k,j));
                    const Real al = (1.0 + bt)*idn;
                    const Real be = bt*idn;
                    jDu0 = al*jS0;
                    jDu1 = al*jS1;
                    jDu2 = al*jS2;
                    jdm0 = -be*jS0;
                    jdm1 = -be*jS1;
                    jdm2 = -be*jS2;
                    jQ = be;
                  }
                  // ck-next: cell i-1's loads are issued before cell i is worked
                  // ck-next: cin, cout and 1/(1 + R beta) formed (CkCinCout, CkIdn)
                  const Real muj = (ck_nq_ == 1) ? 1.0/CK_DIFFUSIVITY : mug[c % 2];
                  Real kr_p = ckkro_g(m,c,ie,k,j), dz_p = lG_g(m,4,ie,k,j);
                  Real rj_p = lP_g(m,0*nch_+c,ie,k,j);   // ck-lin2: e0 formed (CkE0)
                  Real bt_p = lG_g(m,0,ie,k,j);
                  Real ra_p = lG_g(m,1,ie,k,j), g3_p = lG_g(m,3,ie,k,j);
                  Real s0_p = lpf_g(m,c,ie,k,j), s1_p = lps_g(m,c,ie,k,j);
                  Real s2_p = lpj_g(m,c,ie,k,j);
                  // ck-lin2: the layer (i-1, i) weights formed (CkLayW): kappa rho and
                  // dz of cell i-1 loaded two cells ahead
                  Real kr_l = 0.0, dz_l = 0.0;
                  if (ie > icut) {
                    kr_l = ckkro_g(m,c,ie-1,k,j);
                    dz_l = lG_g(m,4,ie-1,k,j);
                  }
                  for (int i=ie; i>icut-1; --i) {
                    const Real c0_p = CkE0(kr_p, dz_p, muj);
                    Real ci, co;
                    CkCinCout(kr_p, dz_p, muj, c0_p, ci, co);
                    const Real tj = 1.0 - c0_p;
                    const Real rj = rj_p, bt = bt_p, rat = ra_p, g3 = g3_p;
                    const Real idn = CkIdn(rj, bt);
                    const Real s0 = s0_p, s1 = s1_p, s2 = s2_p;
                    Real wc = 0.0, wa = 0.0, ffj = 0.0;
                    if (i > icut) CkLayW(kr_l, kr_p, dz_l, dz_p, bface_on, wc, wa, ffj);
                    if (i > icut) {
                      kr_p = kr_l;
                      dz_p = dz_l;
                      rj_p = lP_g(m,0*nch_+c,i-1,k,j);
                      bt_p = lG_g(m,0,i-1,k,j);
                      ra_p = lG_g(m,1,i-1,k,j);
                      g3_p = lG_g(m,3,i-1,k,j);
                      s0_p = lpf_g(m,c,i-1,k,j);
                      s1_p = lps_g(m,c,i-1,k,j);
                      s2_p = lpj_g(m,c,i-1,k,j);
                      if (i-1 > icut) {
                        kr_l = ckkro_g(m,c,i-2,k,j);
                        dz_l = lG_g(m,4,i-2,k,j);
                      }
                    }
                    // the layer joining cells i-1 and i: slots (B_{i-1}, B_i); the third
                    // (B_{i+1}) is 0 for all three
                    Real duv0 = 0.0, duv1 = 1.0, dnl0 = 0.0, dnl1 = 1.0;
                    Real dfw0 = 0.0, dfw1 = 1.0;
                    if (i > icut) {
                      const Real sw = (1.0 - ffj)*wc + ffj*(1.0 - wa);
                      duv0 = 1.0 - wa;
                      duv1 = wa;
                      dnl0 = wc;
                      dnl1 = 1.0 - wc;
                      dfw0 = sw;
                      dfw1 = 1.0 - sw;
                    }
                    // dd_i^+: the down ray's two half steps through cell i
                    Real v0 = tj*jdm0;
                    v0 = tj*v0 + ci*duv0 + co*dfw0;
                    Real v1 = tj*jdm1 + ci*jfu0 + co*jlc0;
                    v1 = tj*v1 + ci*duv1 + co*dfw1;
                    Real v2 = tj*jdm2 + ci*jfu1 + co*jlc1;
                    v2 = tj*v2;
                    const Real al = (1.0 + bt)*idn;
                    const Real be = bt*idn;
                    const Real wr = tj*tj*tj*tj*(1.0 - bt)*idn*jQ;
                    const Real omr = 1.0 - rj;
                    const Real dd0 = al*(1.0 + omr*wr)*s0;
                    const Real dd1 = al*(s1 - omr*v0);
                    const Real dd2 = al*(s2 - omr*v1);
                    const Real dd3 = -al*omr*v2;
                    const Real wj = wfc*g3;
                    const Real jm = wj*(rat*dd1 - jDu0);
                    const Real j0 = wj*(rat*dd2 - jDu1);
                    const Real jp = wj*(rat*dd3 - jDu2);
                    lj0_g(m,c,i,k,j) = j0;
                    ljm_g(m,c,i,k,j) = ((i > icut || gdf >= 0.0) && (jneg_ || jm > 0.0))
                                     ? jm : 0.0;
                    lpj_g(m,c,i,k,j) = (i < ie && (jneg_ || jp > 0.0)) ? jp : 0.0;
                    jDu0 = dd0;
                    jDu1 = dd1;
                    jDu2 = dd2;
                    const Real qn = al*wr + be;
                    jdm0 = -qn*s0;
                    jdm1 = al*v0 - be*s1;
                    jdm2 = al*v1 - be*s2;
                    jQ = qn;
                    jfu0 = dfw0;
                    jfu1 = dfw1;
                    jlc0 = dnl0;
                    jlc1 = dnl1;
                  }
                });
              }
              ck_jlc_ok = true;
            }
            CkParFor4("ck_jlin_sum", cklw_, 0, nmb1, is, ie, ks, ke, js, je,
            KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
              if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
              const int icut = icc_g(m,k,j);
              if (i < icut || icut > ie) return;
              // a masked (zero) coefficient never reads dB/dT outside icut .. ie
              Real s0 = 0.0, s1 = 0.0, s2 = 0.0;
              for (int c=0; c<nch_; ++c) {
                const int b = (ck_nq_ == 1) ? (c/CK_NG) : (c/(2*CK_NG));
                const Real jm = ljm_g(m,c,i,k,j);
                const Real jp = lpj_g(m,c,i,k,j);
                s0 += (jm != 0.0) ? CkMulRn(jm, db_g(m,b,i-1,k,j)) : 0.0;
                s1 += CkMulRn(lj0_g(m,c,i,k,j), db_g(m,b,i,k,j));
                s2 += (jp != 0.0) ? CkMulRn(jp, db_g(m,b,i+1,k,j)) : 0.0;
              }
              // ck_impl_jneg: a negative NET off-diagonal is dropped (the M-matrix);
              // without it every part is already >= 0
              jac_g(m,0,k,j,i) = (jnet_ || s0 > 0.0) ? s0 : 0.0;
              jac_g(m,1,k,j,i) = s1;
              jac_g(m,2,k,j,i) = (jnet_ || s2 > 0.0) ? s2 : 0.0;
            });
          }
          // ck-store: the storing pass's fluxes, by the linear re-apply on the
          // factorisation just built (after rt_chain_ck_jlin, whose windows share the
          // lps/lpf partials the linear kernels write)
          // ck_sph_face = 6: a storing pass that ran the chain kernel (on the face-5
          // map) takes its Src, Fb and Em from the face-6 network too
          if (ckstsp_ || (vef6_ && cklbuild_)) launch_ck_lin_tier();
          // ---- problem/ck_impl_jfd (diagnostic): the tridiagonal of this pass against
          // a finite difference of the linear re-apply.  For every cell index i0 in turn
          // B_b(i0) of EVERY column moves by dB_b/dT dT (dT = 1e-3 T(i0)); the linear
          // kernel is affine in B, so (Src(i0) - Src0(i0))/dT is dSrc_i0/dT_i0 to
          // round-off.  Prints the worst relative mismatch of the diagonal and the
          // counts of positive diagonals (Jacobian / difference).  The pass's own Src,
          // Fb, Em and B are restored (the linear kernel rerun on the saved B).
          if (ck_impl_jfd > 0 && ckjlp_ && cklon_) {
            DvceArray5D<Real> bsv("ck_jfd_b", Bb_g.extent(0), Bb_g.extent(1),
                                  Bb_g.extent(2), Bb_g.extent(3), Bb_g.extent(4));
            Kokkos::deep_copy(bsv, Bb_g);
            DvceArray4D<Real> s0v("ck_jfd_s0", nmb1+1, n1+1, ke+1, je+1);
            DvceArray1D<Real> rv("ck_jfd_r", 10);
            Kokkos::deep_copy(rv, 0.0);
            auto jac_g = ckjac_g;
            auto db_g = ckdb_g;
            const int nbk = nblk;
            par_for("ck_jfd_s0", DevExeSpace(), 0, nmb1, is, ie, ks, ke, js, je,
            KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
              Real sm = 0.0;
              for (int bk=0; bk<nbk; ++bk) sm += Src_g(m,bk,i,k,j);
              s0v(m,i,k,j) = sm;
            });
            for (int i0=is; i0<=ie; ++i0) {
              par_for("ck_jfd_pert", DevExeSpace(), 0, nmb1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int k, const int j) {
                if (i0 < icc_g(m,k,j)) return;
                const Real dT = 1.0e-3*T_g(m,k,j,i0);
                for (int b=0; b<CK_NB; ++b) Bb_g(m,b,i0,k,j) += db_g(m,b,i0,k,j)*dT;
              });
              launch_ck_lin_tier();
              par_for("ck_jfd_cmp", DevExeSpace(), 0, nmb1, ks, ke, js, je,
              KOKKOS_LAMBDA(const int m, const int k, const int j) {
                if (i0 < icc_g(m,k,j)) return;
                if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
                const Real dT = 1.0e-3*T_g(m,k,j,i0);
                Real sm = 0.0;
                for (int bk=0; bk<nbk; ++bk) sm += Src_g(m,bk,i0,k,j);
                const Real fd = (sm - s0v(m,i0,k,j))/dT;
                const Real jj = jac_g(m,1,k,j,i0);
                const Real rel = fabs(fd - jj)/(fabs(fd) + fabs(jj) + 1.0e-300);
                Kokkos::atomic_max(&rv(0), rel);
                Kokkos::atomic_add(&rv(1), 1.0);
                if (jj > 0.0) Kokkos::atomic_add(&rv(2), 1.0);
                if (fd > 0.0) Kokkos::atomic_add(&rv(3), 1.0);
                if (jj > 0.0 && fd > 0.0) Kokkos::atomic_add(&rv(4), 1.0);
                if (fd > 0.0) Kokkos::atomic_max(&rv(5), fd/(fabs(jj) + 1.0e-300));
                if (rel > 1.0e-6) Kokkos::atomic_add(&rv(6), 1.0);
                // the off-diagonals: dSrc_{i0+1}/dT_i0 (row i0+1, slot 0) and
                // dSrc_{i0-1}/dT_i0 (row i0-1, slot 2), against the (clamped) entries
                for (int sd=-1; sd<=1; sd+=2) {
                  const int ir = i0 + sd;
                  if (ir < icc_g(m,k,j) || ir > ie) continue;
                  Real so = 0.0;
                  for (int bk=0; bk<nbk; ++bk) so += Src_g(m,bk,ir,k,j);
                  const Real fo = (so - s0v(m,ir,k,j))/dT;
                  const Real jo = jac_g(m,(sd > 0) ? 0 : 2,k,j,ir);
                  const Real jd = fabs(jac_g(m,1,k,j,ir)) + 1.0e-300;
                  if (fo < 0.0) {
                    Kokkos::atomic_add(&rv(7), 1.0);
                    Kokkos::atomic_max(&rv(8), -fo/jd);
                  } else {
                    Kokkos::atomic_max(&rv(9), fabs(fo - jo)/(fabs(fo) + fabs(jo)
                                                              + 1.0e-300));
                  }
                }
              });
              Kokkos::deep_copy(Bb_g, bsv);
            }
            launch_ck_lin_tier();
            auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), rv);
            std::cout << "### ck_jfd rank=" << global_variable::my_rank << " ncycle="
                      << pm->ncycle << " pass=" << ck_impl_pass << " cells=" << hr(1)
                      << " maxrel_diag=" << hr(0) << " n(rel>1e-6)=" << hr(6)
                      << " n(J>0)=" << hr(2) << " n(fd>0)=" << hr(3) << " both="
                      << hr(4) << " max fd/|J| where fd>0=" << hr(5)
                      << " offdiag: n(fd<0)=" << hr(7) << " max -fd/|J_ii| there="
                      << hr(8) << " maxrel where fd>=0=" << hr(9) << std::endl;
          }
        }
        // ---- ck-fast2 lever 1 (problem/ck_dif_dtau): THE DEEP DIFFUSION -----------
        // Every pass, after the chains: the cells ic .. ich-1 of a handover column get
        // the flux-form divergence of F(f) = sum_b G_b(f) (B_b(f-1) - B_b(f)) (the
        // wall face: route B's datum sum_c wfc I_int,b; face ich: the same F the
        // chains' datum carries), the per-area face fluxes, the emission and, on a
        // pass that builds the tridiagonal, their exact rows (after ck_jlin_sum, which
        // leaves the cells below the chain cut alone).  Written to chain block 0, whose
        // Src/Fb/Em below ich the chain kernels left 0.
        if (ckdif_) {
          auto difG_ = *ck_difG_ptr;
          auto difE_ = *ck_difE_ptr;
          const bool djac_ = ckjlp_;
          const int nchd_ = nblk*RT_NB;
          par_for("ck_dif_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je,
          KOKKOS_LAMBDA(const int m, const int k, const int j) {
            if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
            const int ic = icut_g(m,k,j);
            const int ich = icc_g(m,k,j);
            if (ich <= ic || ich > ie) return;
            Real flo = 0.0;          // per-area flux through the lower face of cell i
            for (int c=0; c<nchd_; ++c) flo += lC_g(0,c)*lC_g(1,c);
            for (int i=ic; i<ich; ++i) {
              Real fup = 0.0, em = 0.0;
              for (int b=0; b<CK_NB; ++b) {
                fup += difG_(m,b,i+1,k,j)*(Bb_g(m,b,i,k,j) - Bb_g(m,b,i+1,k,j));
                em += difE_(m,b,i,k,j)*Bb_g(m,b,i,k,j);
              }
              const Real vol = ACC(m,k,j,i)*dx1(m,k,j,i);
              const Real alo = AFC(m,k,j,i);
              const Real aup = AFC(m,k,j,i+1);
              Src_g(m,0,i,k,j) = (alo*flo - aup*fup)/vol;
              Fb_g(m,0,i,k,j) = flo;
              Em_g(m,0,i,k,j) = em;
              if (djac_) {
                Real glo = 0.0, gdl = 0.0, gdu = 0.0, gup = 0.0;
                for (int b=0; b<CK_NB; ++b) {
                  const Real gu = difG_(m,b,i+1,k,j);
                  gdu += gu*ckdb_g(m,b,i,k,j);
                  gup += gu*ckdb_g(m,b,i+1,k,j);
                  if (i > ic) {
                    const Real gl = difG_(m,b,i,k,j);
                    glo += gl*ckdb_g(m,b,i-1,k,j);
                    gdl += gl*ckdb_g(m,b,i,k,j);
                  }
                }
                ckjac_g(m,0,k,j,i) = alo*glo/vol;
                ckjac_g(m,1,k,j,i) = -(alo*gdl + aup*gdu)/vol;
                ckjac_g(m,2,k,j,i) = aup*gup/vol;
              }
              flo = fup;
            }
          });
        }
      } else {
      par_for("rt_chain", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int blk, const int k, const int j) {
        constexpr int NN = RT_NNC;
        auto tau_down_r_f = Kokkos::subview(tau_g, m, k, j, Kokkos::ALL);
        auto B            = Kokkos::subview(B_g,   m, k, j, Kokkos::ALL);
        auto F_ir_f       = Kokkos::subview(Fb_g,  m, blk, k, j, Kokkos::ALL);
        const Real gamir1 = cf_g(m,k,j,0);
        const Real gamir2 = cf_g(m,k,j,1);
        const Real beta   = cf_g(m,k,j,2);
        // this block owns its own flux slot, so it starts from zero exactly as the
        // serial accumulator did
        for (int i=is; i<ie+2; ++i) {
          Fb_g(m,blk,i,k,j) = 0.0;
          Em_g(m,blk,i,k,j) = 0.0;
          Src_g(m,blk,i,k,j) = 0.0;
        }
          Real gamirc[NC], fbc[NC], muggc[NC], wggc[NC];
          for (int cc=0; cc<NC; ++cc) {
            const int c = blk*NC + cc;
            const int n = (c/2) % 2;
            const int vir = c % 2;
            muggc[cc] = mug[n];
            wggc[cc] = wg[n];
            gamirc[cc] = (vir == 0) ? gamir1 : gamir2;
            fbc[cc] = (vir == 0) ? beta : (1.0-beta);
          }
          Real I_ir_down_c[NC][NN];
#if RT_CACHE
          Real e0c[NC][NN], alpc[NC][NN], betc[NC][NN];
#endif

          // top.  Everything carried from here on is the AREA-WEIGHTED intensity
          // J = A I of the dilution note above, so the incoming beam enters scaled by
          // the top face's area.
          const Real aft_ = AFC(m,k,j,ie+1);
          for (int cc=0; cc<NC; ++cc) {
            Real dtauir = gamirc[cc]*tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/muggc[cc]);
            I_ir_down_c[cc][ie+1] = (1.0-trans)*(fbc[cc]*B[ie+1]*aft_);
          }
          if (layer_legacy) {
            // down-sweep
            for (int i=ie; i>is-1; --i) {
              Real dtau_i = tau_down_r_f[i]-tau_down_r_f[i+1];
              // emissivity-weighted far endpoint of the layer's source function: dtau/dr
              // is kappa rho, so the two layers' optical thicknesses per unit length are
              // the weights (see BFace).  The ghost above ie has no dtau in the array, so
              // the topmost layer keeps the old endpoint exactly.
              const Real kro_g = dtau_i/dx1(m,k,j,i);
              const Real krf_g = (i < ie) ? (tau_down_r_f[i+1]-tau_down_r_f[i+2])
                                          / dx1(m,k,j,i+1) : kro_g;
              const Real bfr_g = BFace(kro_g, krf_g, B[i], B[i+1], bface_on);
              // the source function of the J equation is A B at the endpoint's own
              // position; this layer runs between the two cell centres' mean areas
              const Real bfa_g = bfr_g*ACC(m,k,j,i+1);
              const Real bwn_g = B[i]*ACC(m,k,j,i);
              for (int cc=0; cc<NC; ++cc) {
                Real dtauir = gamirc[cc]*dtau_i;
                Real x = dtauir/muggc[cc];
                Real e0 = -expm1(-x);
                Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
                Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
                // direct source: what this stream leaves in cell i, absorbed minus
                // emitted
                Src_g(m,blk,i,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i)
                                    *(e0*I_ir_down_c[cc][i+1]
                                      - fbc[cc]*(alp*bfa_g + bet*bwn_g));
                I_ir_down_c[cc][i] = (1.0-e0)*I_ir_down_c[cc][i+1]
                                   + alp*fbc[cc]*bfa_g + bet*fbc[cc]*bwn_g;
#if RT_CACHE
                e0c[cc][i] = e0;
                alpc[cc][i] = alp;
                betc[cc][i] = bet;
#endif
              }
            }

            // bottom
            Real I_ir_up_c[NC];
            const Real afb_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              I_ir_up_c[cc] = Iint*afb_ + I_ir_down_c[cc][is];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is]/afb_;
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afb_;
              Fb_g(m,blk,is,k,j) += (F_ir_up_f - F_ir_down_f);
            }
            // up-sweep, accumulating the band flux as it goes
            for (int i=is+1; i<ie+2; ++i) {
#if !RT_CACHE
              Real dtau_i = tau_down_r_f[i-1]-tau_down_r_f[i];
#endif
              // the same emissivity-weighted far endpoint for the upward stream and for
              // Em; the topmost layer keeps the old endpoint (see the down-sweep)
              const Real krou = (tau_down_r_f[i-1]-tau_down_r_f[i])/dx1(m,k,j,i-1);
              const Real krfu = (i < ie+1) ? (tau_down_r_f[i]-tau_down_r_f[i+1])
                                           / dx1(m,k,j,i) : krou;
              const Real bfru = BFace(krou, krfu, B[i-1], B[i], bface_on);
              const Real bfa_u = bfru*ACC(m,k,j,i);
              const Real bwn_u = B[i-1]*ACC(m,k,j,i-1);
              const Real afi_ = AFC(m,k,j,i);
              for (int cc=0; cc<NC; ++cc) {
#if RT_CACHE
                // layer i-1, already solved on the way down
                const Real e0 = e0c[cc][i-1];
                const Real bet = betc[cc][i-1];
                const Real gm = alpc[cc][i-1];
#else
                Real dtauir = gamirc[cc]*dtau_i;
                Real x = dtauir/muggc[cc];
                Real e0 = -expm1(-x);
                Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
                Real gm = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
#endif
                const Real Iup_in = I_ir_up_c[cc];
                Src_g(m,blk,i-1,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i-1)
                                      *(e0*Iup_in - fbc[cc]*(bet*bfa_u + gm*bwn_u));
                I_ir_up_c[cc] = (1.0-e0)*Iup_in
                              + bet*fbc[cc]*bfa_u + gm*fbc[cc]*bwn_u;
                Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][i]/afi_;
                Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afi_;
                Fb_g(m,blk,i,k,j) += (F_ir_up_f - F_ir_down_f);
                // the cell's own emission per unit volume, both hemispheres.  dtau_i is
                // kappa rho dr for this layer, so dtau_i/dr = kappa rho and the layer
                // thickness cancels out of the volumetric rate.
                const Real dtau_ly = tau_down_r_f[i-1]-tau_down_r_f[i];
                Em_g(m,blk,i-1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_ly
                                     / dx1(m,k,j,i-1)*fbc[cc]*0.5*(bfru+B[i-1]);
              }
            }
          } else {
            // ---- the centre-to-centre layers (see rt_layer_legacy) ----
            // Same construction as the grey and correlated-k sweeps: a layer runs from
            // one cell CENTRE to the next, split at the face into two halves, each half
            // lying entirely inside one cell.  Here the optical depths come ready-made
            // from the face-cumulative tau array, so a half layer is simply half a
            // cell's dtau and nothing extra is looked up.  The coefficient cache is not
            // used on this path: both halves of a cell have the SAME thickness, but the
            // cache is indexed by the whole-cell layer the legacy branch walks.
            Real Idn[NC];
            // one half-layer step; dsrc is absorbed minus emitted, which is exactly the
            // change in the intensity, so Src_g and the divergence of Fb_g agree
            auto step = [&](const Real dtau, const Real mu, const Real s_in,
                            const Real s_out, Real &I, Real &dsrc) {
              const Real x = dtau/mu;
              const Real e0 = -expm1(-x);
              const Real cin = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              const Real cout = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              const Real em = cin*s_in + cout*s_out;
              dsrc = e0*I - em;
              I = (1.0 - e0)*I + em;
            };
            for (int cc=0; cc<NC; ++cc) {
              Idn[cc] = I_ir_down_c[cc][ie+1];
            }
            // down-sweep, recording the FACE intensity in between the two halves
            for (int i=ie; i>is-1; --i) {
              const Real dtau_l = tau_down_r_f[i]-tau_down_r_f[i+1];
              const Real dzl = dx1(m,k,j,i);
              if (i == ie) {
                // the top half layer: the upper half of cell ie, its source held at the
                // cell's own Planck function -- there is nothing above it
                // the half layer runs from the TOP FACE to the cell CENTRE, so its
                // two source endpoints take those two areas
                const Real act_ = ACC(m,k,j,ie);
                for (int cc=0; cc<NC; ++cc) {
                  const Real s = fbc[cc]*B[ie];
                  Real dsrc;
                  step(0.5*gamirc[cc]*dtau_l, muggc[cc], s*aft_, s*act_, Idn[cc], dsrc);
                  Src_g(m,blk,ie,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,ie)*dsrc;
                }
              } else {
                // the layer between the centres of cells i+1 and i.  dtau/dr is kappa
                // rho, so the two cells' optical thicknesses per unit length are the
                // emissivity weights BFace wants, applied symmetrically to both
                // endpoints.
                const Real dtau_u = tau_down_r_f[i+1]-tau_down_r_f[i+2];
                const Real dzu = dx1(m,k,j,i+1);
                const Real krl = dtau_l/dzl;
                const Real kru = dtau_u/dzu;
                const Real s_l = BFace(kru, krl, B[i+1], B[i], bface_on);
                const Real s_u = BFace(krl, kru, B[i], B[i+1], bface_on);
                const Real dtc = dtau_l + dtau_u;
                const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dtau_l/dtc))
                                             : (0.5*(s_l + s_u));
                const Real acu_ = ACC(m,k,j,i+1), afm_ = AFC(m,k,j,i+1);
                const Real acl_ = ACC(m,k,j,i);
                for (int cc=0; cc<NC; ++cc) {
                  Real dsrc;
                  step(0.5*gamirc[cc]*dtau_u, muggc[cc], fbc[cc]*s_u*acu_,
                       fbc[cc]*s_f*afm_, Idn[cc], dsrc);
                  Src_g(m,blk,i+1,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i+1)
                                        *dsrc;
                  I_ir_down_c[cc][i+1] = Idn[cc];
                  step(0.5*gamirc[cc]*dtau_l, muggc[cc], fbc[cc]*s_f*afm_,
                       fbc[cc]*s_l*acl_, Idn[cc], dsrc);
                  Src_g(m,blk,i,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,i)*dsrc;
                }
              }
            }
            // the lower half of cell is, down to the bottom face
            const Real dtau_is = tau_down_r_f[is]-tau_down_r_f[is+1];
            const Real acs_ = ACC(m,k,j,is), afs_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              const Real s = fbc[cc]*B[is];
              Real dsrc;
              step(0.5*gamirc[cc]*dtau_is, muggc[cc], s*acs_, s*afs_, Idn[cc], dsrc);
              Src_g(m,blk,is,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,is)*dsrc;
              I_ir_down_c[cc][is] = Idn[cc];
            }

            // bottom, now AT the bottom face
            Real I_ir_up_c[NC];
            const Real afb_ = AFC(m,k,j,is);
            for (int cc=0; cc<NC; ++cc) {
              I_ir_up_c[cc] = Iint*afb_ + I_ir_down_c[cc][is];
              Real F_ir_down_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_down_c[cc][is]/afb_;
              Real F_ir_up_f = 2.0*M_PI*wggc[cc]*muggc[cc]*I_ir_up_c[cc]/afb_;
              Fb_g(m,blk,is,k,j) += (F_ir_up_f - F_ir_down_f);
            }
            // up-sweep: the bottom half layer first, then centre to centre, accumulating
            // the band flux as it goes.  Em_g is the cell's OWN emission now, since each
            // cell owns both of its halves.
            for (int cc=0; cc<NC; ++cc) {
              const Real s = fbc[cc]*B[is];
              Real dsrc;
              step(0.5*gamirc[cc]*dtau_is, muggc[cc], s*afs_, s*acs_, I_ir_up_c[cc],
                   dsrc);
              Src_g(m,blk,is,k,j) += 2.0*M_PI*wggc[cc]*muggc[cc]/VLS(m,k,j,is)*dsrc;
              Em_g(m,blk,is,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_is
                                  / dx1(m,k,j,is)*fbc[cc]*B[is];
            }
            for (int i=is; i<ie; ++i) {
              const Real dtau_l = tau_down_r_f[i]-tau_down_r_f[i+1];
              const Real dtau_u = tau_down_r_f[i+1]-tau_down_r_f[i+2];
              const Real dzl = dx1(m,k,j,i);
              const Real dzu = dx1(m,k,j,i+1);
              const Real krl = dtau_l/dzl;
              const Real kru = dtau_u/dzu;
              const Real s_l = BFace(kru, krl, B[i+1], B[i], bface_on);
              const Real s_u = BFace(krl, kru, B[i], B[i+1], bface_on);
              const Real dtc = dtau_l + dtau_u;
              const Real s_f = (dtc > 0.0) ? (s_l + (s_u - s_l)*(dtau_l/dtc))
                                           : (0.5*(s_l + s_u));
              const Real acl2_ = ACC(m,k,j,i), afm2_ = AFC(m,k,j,i+1);
              const Real acu2_ = ACC(m,k,j,i+1);
              for (int cc=0; cc<NC; ++cc) {
                const Real pref = 2.0*M_PI*wggc[cc]*muggc[cc];
                Real dsrc;
                step(0.5*gamirc[cc]*dtau_l, muggc[cc], fbc[cc]*s_l*acl2_,
                     fbc[cc]*s_f*afm2_, I_ir_up_c[cc], dsrc);
                Src_g(m,blk,i,k,j) += pref/VLS(m,k,j,i)*dsrc;
                Fb_g(m,blk,i+1,k,j) += pref*(I_ir_up_c[cc]
                                             - I_ir_down_c[cc][i+1])/afm2_;
                step(0.5*gamirc[cc]*dtau_u, muggc[cc], fbc[cc]*s_f*afm2_,
                     fbc[cc]*s_u*acu2_, I_ir_up_c[cc], dsrc);
                Src_g(m,blk,i+1,k,j) += pref/VLS(m,k,j,i+1)*dsrc;
                Em_g(m,blk,i+1,k,j) += 2.0*(2.0*M_PI*wggc[cc])*gamirc[cc]*dtau_u
                                     / dzu*fbc[cc]*B[i+1];
              }
            }
            // the top half of cell ie, back out through the top face
            {
              const Real dtau_ie = tau_down_r_f[ie]-tau_down_r_f[ie+1];
              const Real act2_ = ACC(m,k,j,ie);
              for (int cc=0; cc<NC; ++cc) {
                const Real pref = 2.0*M_PI*wggc[cc]*muggc[cc];
                const Real s = fbc[cc]*B[ie];
                Real dsrc;
                step(0.5*gamirc[cc]*dtau_ie, muggc[cc], s*act2_, s*aft_,
                     I_ir_up_c[cc], dsrc);
                Src_g(m,blk,ie,k,j) += pref/VLS(m,k,j,ie)*dsrc;
                Fb_g(m,blk,ie+1,k,j) += pref*(I_ir_up_c[cc]
                                              - I_ir_down_c[cc][ie+1])/aft_;
              }
            }
          }
      });
      }

      // ---- C: reduce over blocks in order, then apply ------------------------------
      int nclip = 0;
      const Real demax = rt_de_max;
      const bool semilin = rt_semi_lin;
      const bool newton_on = rt_newton;
      // see rt_ali_diag.  Only the band paths carry a per-cell continuum opacity in
      // kc_g, so the diagonal is available there and nowhere else.
      const bool ali_on = rt_ali_diag && band_on;
      const int ali_nq = (ck_nq_ > 1) ? 2 : 1;
      if (rt_efix_ptr == nullptr) {
        rt_efix_ptr = new DvceArray1D<int>("rt_efix", 3);
        Kokkos::deep_copy(*rt_efix_ptr, 0);
        rt_efix_rec = new DvceArray1D<Real>("rt_efix_rec", 16);
        Kokkos::deep_copy(*rt_efix_rec, 0.0);
      }
      auto efix_g = *rt_efix_ptr;
      auto efrec_g = *rt_efix_rec;
      if (rt_desum_ptr == nullptr) {
        rt_desum_ptr = new DvceArray1D<Real>("rt_desum", 2);
      }
      auto dsum_g = *rt_desum_ptr;
      if (ck_impl_nosync) {
        Kokkos::deep_copy(DevExeSpace(), dsum_g, 0.0);   // stream-ordered, no fence
      } else {
        Kokkos::deep_copy(dsum_g, 0.0);
      }
      const int efix_cyc = pm->ncycle;
      const bool resc_eq = rt_rescue_eq;
      // the sub-cycled local relaxation; see rt_relax_sub
      const int nsub_in = rt_relax_sub;
      const int nsub_max = (rt_relax_submax > 1) ? rt_relax_submax : 1;
      const Real xcrit = rt_relax_xcrit;
      // rt_cell_report: claimed once per RT call, so only the FIRST rescued cell prints
      // (problem/ck_impl_nosync: one cached View, zeroed in stream order)
      static DvceArray1D<int> *repc_c = nullptr;
      DvceArray1D<int> repc_g;
      if (ck_impl_nosync) {
        if (repc_c == nullptr) repc_c = new DvceArray1D<int>(std::string("rt_repc"), 1);
        repc_g = *repc_c;
        Kokkos::deep_copy(DevExeSpace(), repc_g, 0);
      } else {
        repc_g = DvceArray1D<int>(std::string("rt_repc"), 1);
        Kokkos::deep_copy(repc_g, 0);
      }
      const bool fixed_on = report_on && (rt_report_every > 0) &&
                            (pm->ncycle % rt_report_every == 0);
      const Real rep_r = rt_report_r;
      const int rep_k = ks, rep_j = js;
      const int rep_cyc = pm->ncycle;
      const Real rep_time = pm->time;
      const bool direct_on = rt_src_direct && (ck_on || grey_on);
      const bool semi_imp = rt_semi_implicit;
      // problem/rt_floor_consistent: bound the cooling step by the EOS floor state
      const bool efc_ = rt_floor_consistent;
      // see rt_apply_debug: which column, and how many calls are left to print
      const bool dbg_on = (rt_apply_debug > 0);
      const int dbg_m = rt_dump_m;
      const int dbg_j = (rt_dump_j >= 0) ? rt_dump_j : (js + je)/2;
      const int dbg_k = (rt_dump_k >= 0) ? rt_dump_k : (ks + ke)/2;
      const int dbg_n = rt_apply_debug_n;
      if (dbg_on) --rt_apply_debug;
      // per-cycle diagnostic: an empty View captures fine, so the lambda needs no
      // branch on the pointer itself
      const bool diag = rt_diag;
      DvceArray5D<Real> dg;
      if (diag) dg = *rt_diag_ptr;
      // ---- the radiative momentum source (problem/rt_rad_force) --------------------
      const bool radforce = rt_rad_force && grey_on;
      const Real inv_c = 1.0/2.99792458e10;      // cgs: this path runs in cgs code units
      const Real arad_f = eos.tbl.arad;
      const Real xlo_f = eos.tbl.rad_lrho_lo, xhi_f = eos.tbl.rad_lrho_hi;
      // the taper's TEMPERATURE GATE (see rad_taper::WeightGated).  The force must use
      // the SAME w the EOS uses, or the two stop being complementary and the total force
      // is neither -grad Prad nor kappa F/c; that is the whole point of sharing the
      // weight function between the two modules.
      const Real ylo_f = eos.tbl.rad_lt_lo, yhi_f = eos.tbl.rad_lt_hi;
      const bool tg_f = eos.tbl.rad_tgate;
      const bool md_f = pm->multi_d, td_f = pm->three_d;
      // problem/rt_force_tau_gate: the optical-depth gate s(tau) on the force above.
      // OFF (or s = 0) executes the old expressions verbatim -- see rt_force_tau_gate.
      const bool ftg_ = radforce && rt_force_tau_gate;
      const Real ftlg_lo = log10(rt_force_tau_lo);
      const Real ftlg_hi = log10(rt_force_tau_hi);
      const Real ft_hi = rt_force_tau_hi;
      const bool fverb = radforce && (rt_force_verbose > 0);
      if (fverb) --rt_force_verbose;
      // problem/rt_force_verbose normalises the cell balance by g.  rt_force_grav is a
      // CONSTANT (the pgen's GM/r_in^2), which on a global radial domain is not the
      // gravity any cell feels -- at the He star's photosphere it is 4x the local value,
      // so both the printed accelerations and the residual were mis-scaled.  Use the
      // LOCAL g(r) that the force itself works against, and keep rt_force_grav only as
      // the fallback for a mesh with no radius (the plane-parallel box).
      const Real gver = rt_force_grav;
      auto cfv_ = (rt_cf_ptr != nullptr) ? *rt_cf_ptr
                : CkDum<DvceArray4D<Real>>("rt_cf_dummy");
      const bool cfv_on = (rt_cf_ptr != nullptr);
      const Real apo_ = ap;             // the orbit's a; the print shadows `ap` below
      const int vcyc = pm->ncycle;
      auto eos_f = eos;
      // problem/rt_budget_verbose: the v.f work accumulator (see rt_bud_ptr)
      const bool budg_ = (rt_bud_ptr != nullptr) && radforce;
      auto bud_ = budg_ ? *rt_bud_ptr : CkDum<DvceArray1D<Real>>("rtbuddummy");
      const bool cklean_ = ck_impl_nosync && ckimp_ && !report_on && !diag && !radforce
                           && !dbg_on;
      par_reduce_clip4("rt_apply", 0, nmb1, ks, ke, js, je, is, ie, nclip,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, int &nc) {
        // ck_impl_colskip: a converged column is left exactly as its last pass left it
        if (ckskip_ && ckdone_g(m,k,j) > 0.0) return;
        Real Ft = 0.0, Fb = 0.0;
        for (int b=0; b<nblk; ++b) {
          Ft += Fb_g(m,b,i+1,k,j);
          Fb += Fb_g(m,b,i,k,j);
        }
        // the face fluxes are PER UNIT AREA; every divergence below is therefore
        // (A_top F_top - A_bot F_bot)/V (the dilution note above), which under
        // rt_plane_parallel is A = 1, V = DX1 and bitwise the old expression
        // ---- THE CORRELATED-K APPLY IS PLANE-PARALLEL UNLESS ck_spherical IS ON ------
        // The ck sweep builds Fb_g with the PLANE-PARALLEL divergence unless
        // ck_spherical is set, so the apply must difference it the same way: areas 1,
        // volume dx1, which is what the production binary (27ca5b13) did.  1159a8f3 made
        // this kernel area-weighted for the GREY centre-to-centre sweep -- which does
        // carry the r^-2 dilution in its face conditions -- and that form then reached
        // the ck path too, where it is inconsistent: an undiluted F differenced with
        // areas manufactures luminosity at the area ratio (column budget +0.63 on the
        // production grid, a 36-40 % sweep-to-gas gap, the top collapsing onto pfloor).
        // The grey and picket-fence paths are UNTOUCHED here.  See
        // tests_ck_implicit/README.md, section 0.
        const bool ckpp_ = ck_on && !cksph_;
        const Real aft_a = ckpp_ ? 1.0 : AFC(m,k,j,i+1);
        const Real afb_a = ckpp_ ? 1.0 : AFC(m,k,j,i);
        const Real vla_a = ckpp_ ? DX1(m,k,j,i) : VLA(m,k,j,i);
        // the two-stream's share of each face in the tau blend
        Real src;
        // THE BLEND HANDOVER IS NOT AN EMISSION, AND MUST NOT BE RELAXED.  src below is
        // the divergence of the two-stream's BLENDED share, and inside the tau ramp it
        // is dominated by the handover term d[w F]/dr: the flux the two-stream gives up
        // to the conduction operator (and the conduction operator, face by face, takes
        // up in full).  The semi-implicit step further down damps whatever it is handed
        // by (1 - e^-x)/x per cell, which is the right thing to do to the LOCAL
        // absorption-emission balance -- it is stiff -- and the wrong thing to do to a
        // handover: one half of the exchange is then damped and the other is not, the
        // two no longer cancel, and the residual is deposited in the blend layer for as
        // long as the run lasts.  MEASURED on the 1-D grey atmosphere with tau_bottom
        // = 1000 and rad_tau_lo/hi = 10/100, started FROM the analytic steady state:
        // +0.098 F of spurious heating, 68 % of it inside tau 10-100, against +0.006 F
        // with problem/rt_semi_implicit = false.
        //
        // So the source is split.  src_relax, the (1-w)-weighted local balance, goes
        // through the relaxation and is paired with the (1-w)-weighted Em it must
        // balance; src_ex = src - src_relax, which is the handover and nothing else, is
        // applied EXPLICITLY and therefore exactly.  The explicit limit of the whole
        // step is src*bdt either way.  Outside the blend (w = 0 on both faces, i.e.
        // every cell the two-stream owns alone, and every run with no tau blend at all)
        // src_relax == src and src_ex == 0 identically, so this is inert there.
        Real src_relax;
        const Real wbar = taublend ? 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1)) : 0.0;
        if (direct_on) {
          // The direct source is the divergence of the FULL two-stream flux.  What the
          // two-stream actually deposits is the divergence of its blended share,
          //     -d[(1-w)F]/dr = -dF/dr + d[wF]/dr,
          // and the second term is the flux it hands to the conduction operator across
          // the blend layer -- an O(F/dr) exchange, not a correction.  Keep it, formed
          // from the raw face fluxes.  Where w = 0 (the thin layers this rewrite is for)
          // it vanishes identically, so no round-off of F - F re-enters there; where it
          // is non-zero the cell is at tau > 10 and its energy dwarfs that round-off.
          src = 0.0;
          for (int b=0; b<nblk; ++b) src += Src_g(m,b,i,k,j);
          src_relax = src;
          if (taublend) {
            src += (w_g(m,k,j,i+1)*(Ft*aft_a)
                    - w_g(m,k,j,i)*(Fb*afb_a))/vla_a;
            src_relax *= (1.0 - wbar);
          }
        } else {
          const Real dvg = -(Ft*aft_a-Fb*afb_a)/vla_a;
          if (taublend) {
            Ft *= (1.0 - w_g(m,k,j,i+1));
            Fb *= (1.0 - w_g(m,k,j,i));
          }
          src = -(Ft*aft_a-Fb*afb_a)/vla_a;
          src_relax = taublend ? (1.0 - wbar)*dvg : src;
        }
        Real Qs_d = 0.0;   // the stellar heating that entered src, for the diagnostic
        if (band_on) {
          // deeper than the cut nothing radiative is applied: that region is optically
          // thick and convective, and the stellar beam died decades of optical depth
          // above
          if (i < icut_g(m,k,j)) {
            src = 0.0;
            src_relax = 0.0;
          } else {
            Real Qs = 0.0;
            for (int b=0; b<nblk; ++b) Qs += Qb_g(m,b,i,k,j);
            src += Qs;
            src_relax += Qs;
            Qs_d = Qs;
          }
        } else {
          Qs_d = Qv_g(m,k,j,i);
          src += Qs_d;
          src_relax += Qs_d;
        }
        // ---- problem/ck_implicit: hand the column solve the residual it needs -------
        // S_i is the FULL source -- the (1-w)-weighted flux divergence, the tau-blend
        // handover and the direct beam -- so the fixed point of the Newton is the exact
        // backward-Euler balance.  e^n is latched on the FIRST pass only; every later
        // pass re-reads the running e out of u0 through eiN, which is why rt_use_cons is
        // required.  NOTHING is applied to the gas here; see the skip below.
        if (ckimp_) {
          cksrc_g(m,k,j,i) = src;
          const Real eic = eiN(m,k,j,i);
          ckei_g(m,k,j,i) = eic;
          // e^n, recorded NET of any warm-start seed already applied to the gas (see
          // ck_impl_warm; the seed array is 0 with the switch off, so this is bitwise the
          // old expression there)
          if (ckpass0_) ckest_g(m,k,j,i) = eic - cksd_g(m,k,j,i);
          // the cell's OWN emission, with EXACTLY the weights the semi-implicit block
          // below gives it, so that S = A - E holds with the same E the two schemes
          // mean.  The thin-cell solve and the fallback both need the split.
          Real Emc = 0.0;
          for (int b=0; b<nblk; ++b) Emc += Em_g(m,b,i,k,j);
          if (taublend) Emc *= 1.0 - wbar;
          if (band_on && i < icut_g(m,k,j)) Emc = 0.0;
          ckem_g(m,k,j,i) = (Emc > 0.0) ? Emc : 0.0;
          // ---- THE TWO-LEVEL SPLIT'S VERDICT, see ck_impl_arat.  A = S + E; a cell is
          // safe to linearise while the field it absorbs is within ckarat_ of its own
          // emission, and a COOLING cell (A < E) always is.  Read by rt_pre_opac of the
          // NEXT pass, which is what freezes it for that pass's matrix and solve.
          ckthk_g(m,k,j,i) =
              ((Emc > 0.0) && (src + Emc <= ckarat_*Emc)) ? 1.0 : 0.0;
          // problem/ck_impl_nosync: nothing below reaches the gas or the solve under
          // ck_implicit (skip_de); with no diagnostic asked for, stop here
          if (cklean_) return;
        }
        // SEMI-IMPLICIT APPLICATION.  The source splits as src = A - E(T), A being the
        // absorption of the field from elsewhere, fixed on this step, and E the cell's
        // own emission, which is 4 sigma kappa_P rho T^4 and so scales as T^4.  The step
        // has to be right in two limits at once: it must reduce to the explicit rate
        // src*bdt as bdt -> 0, and it must land on the state where the cell stops
        // exchanging energy, E(T) = A, as bdt -> infinity.
        //
        // WHY NOT LINEARIZE.  Treating E implicitly by linearizing about the CURRENT
        // state gives a relaxation rate lambda = dE/de = (4E/T)/(rho c_v) ~ 4E/e and the
        // step (src/lambda)(1 - exp(-lambda bdt)).  That is only valid while |src| <~ E.
        // Its asymptote is e(A-E)/(4E), which for a cold optically thin cell in a hotter
        // field (A >> E) DIVERGES as E -> 0: the emission is damped but the absorption
        // stays explicit and unbounded.  A 400 K ambient cell sitting in a field whose
        // equilibrium temperature is 1500 K has A/E = (1500/400)^4 = 197 and is handed
        // de/e ~ 49 in one step, where the true answer is T_eq/T - 1 = 2.75 -- an
        // overshoot of ~18x, and the reason rt_de_max used to clip essentially every cell
        // above the photosphere.  Cooling was safe (A = 0 gives -e/4); heating was not.
        //
        // WHAT THIS DOES INSTEAD.  Relax toward the TRUE fixed point.  E ~ T^4 and, under
        // the same e ~ rho c_v T the old rate assumed, e ~ T, so E(T_eq) = A puts the
        // equilibrium at e_eq = e (A/E)^(1/4) and deq = e_eq - e.  Then take the rate
        // from the initial slope rather than from lambda -- x = src*bdt/deq -- so that
        //     de = deq (1 - exp(-x))
        // is exactly src*bdt for small bdt and exactly deq for large bdt.  src and deq
        // always share a sign (A > E <=> e_eq > e), so x >= 0 and the step is monotone
        // and can never cross the equilibrium from either side.  Where the departure is
        // small it agrees with the old linearized form, since deq -> src/lambda there.
        //
        // The e ~ T behind e_eq is the same approximation the old rate made, and it
        // UNDER-estimates e_eq wherever H2 or H is partly dissociated, i.e. it errs
        // toward a smaller step.  problem/rt_semi_lin recovers the old linearization.
        // ck_implicit deposits nothing here: the column solve owns the energy update,
        // and it takes the CONVERGED flux divergence with no relaxation factor and no
        // rt_de_max clip.  de below is kept only for the diagnostics.
        const bool skip_de = mode3_ || ckimp_;
        // ---- rt_ali_diag: 1 - Lambda*_ii, the local escape factor -----------------
        Real ome_ = 1.0;
        if (ali_on) {
          const Real dtau_a = kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*DX1(m,k,j,i);
          if (dtau_a > 0.0) {
            Real lst = 0.0;
            if (ali_nq > 1) {
              const Real mq[2] = {0.21132487, 0.78867513};
              for (int q=0; q<2; ++q) {
                const Real xq = dtau_a/mq[q];
                lst += 0.5*((xq > 1.0e-3) ? (1.0 - (-expm1(-xq))/xq)
                                          : (0.5*xq - xq*xq/6.0));
              }
            } else {
              const Real xq = dtau_a*1.66;
              lst = (xq > 1.0e-3) ? (1.0 - (-expm1(-xq))/xq) : (0.5*xq - xq*xq/6.0);
            }
            ome_ = 1.0 - lst;
            if (!(ome_ > 1.0e-8)) ome_ = 1.0e-8;   // never divide by zero below
            if (ome_ > 1.0) ome_ = 1.0;
          }
        }
        Real de = src*bdt;
        // rt_cell_report bookkeeping: the pieces of the step, kept for the report below
        Real dg_A = 0.0, dg_Em = 0.0, dg_deq = 0.0;
        Real dg_it[8];
        int dg_nit = 0;
        bool dg_resc = false;
        int dg_nsub = 1;   // the sub-cycle count actually used
        for (int q=0; q<8; ++q) dg_it[q] = 0.0;
        // problem/rt_semi_implicit = false: nothing else touches de.
        if (semi_imp && !skip_de) {
          Real Em = 0.0;
          for (int b=0; b<nblk; ++b) Em += Em_g(m,b,i,k,j);
          if (taublend) {
            Em *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (band_on && i < icut_g(m,k,j)) Em = 0.0;
          // e^n, the energy this stage STARTED from.  Without the outer iteration that
          // is just the state the solver reads; with it, u0 already carries de_{k-1}
          // and has to be walked back, because de below is the TOTAL stage increment
          // and not an addition to what the previous pass left.  T_g and Em, by
          // contrast, are deliberately the UPDATED state's: Em(T^{k-1}) together with
          // t0 = T^{k-1} is what makes the Newton's Em*(T(e^n+de)/t0)^4 the emission at
          // the NEW energy, i.e. F(de) = de - bdt(A - Em(e^n+de)) exactly.
          const Real ei = eiN(m,k,j,i);
          if (Em > 0.0 && ei > 0.0) {
            const Real src_ex = src - src_relax;     // the handover, applied exactly
            const Real sdt = src_relax*bdt;
            if (semilin) {
              const Real lam = ome_*4.0*Em/ei;     // see rt_ali_diag
              const Real x = lam*bdt;
              de = (x > 1.0e-4) ? (src_relax/lam)*(-expm1(-x)) : sdt;
              de += src_ex*bdt;
            } else {
              // SUB-CYCLED RELAXATION (problem/rt_relax_sub).  The step above relaxes
              // toward the equilibrium of the UN-relaxed column: A is one Jacobi pass,
              // formed from the intensities of the state the sweep saw, and the cell is
              // then damped toward it by a per-cell factor (1 - exp(-x))/x.  That factor
              // is what breaks the problem: the non-local exchange between two cells
              // cancels exactly only while both are damped by the SAME amount, and two
              // neighbours at different x are not.  What survives the near-cancellation
              // is a pressure perturbation proportional to dt times the local compression
              // rate, i.e. a term that pumps rather than damps, and in a closed box it
              // feeds the organ-pipe modes: the He-star box's saturated v_rms came out
              // LINEAR in dt, which no converged solution can be.
              //
              // The fix is to stop taking one large damped step.  A, the absorbed field,
              // stays frozen for the stage -- it is what the sweep computed and re-doing
              // the sweep is what this is trying to avoid -- but the LOCAL balance is
              // integrated in nsub sub-steps of bdt/nsub, re-forming the cell's own
              // emission from the running energy after each one (Em(T) ~ T^4, with T(e)
              // from the EOS where there is one).  Each sub-step is then taken at small
              // x, where (1 - exp(-x))/x -> 1 and the damping that broke the cancellation
              // is gone, while the sum still lands on the same fixed point E(T) = A.
              //
              // nsub is chosen per cell from the stiffness of the first sub-step,
              // ceil(x/rt_relax_xcrit), floored at rt_relax_sub and capped at
              // rt_relax_submax.  rt_relax_sub <= 1 disables the whole thing and is
              // BITWISE the single-step form above.
              // A_eff: the equilibrium the cell is relaxed toward is the one where the
              // NET exchange vanishes, (1 - Lambda*)(E(T) - E) = src, i.e. E_eq = E +
              // src/(1 - Lambda*).  With rt_ali_diag off ome_ is 1 and this is the old
              // A = src + Em, bit for bit.  See rt_ali_diag.
              const Real absn = src_relax/ome_ + Em;  // A, held fixed over the step
              // e_eq - e.  With nothing arriving the equilibrium is T = 0, i.e. -e.
              const Real deq0 = (absn > 0.0) ? ei*(sqrt(sqrt(absn/Em)) - 1.0) : -ei;
              int nsub = 1;
              if (nsub_in > 1) {
                nsub = nsub_in;
                if (deq0 != 0.0 && xcrit > 0.0) {
                  const Real xn = ceil((sdt/deq0)/xcrit);
                  if (xn > static_cast<Real>(nsub)) {
                    nsub = (xn >= static_cast<Real>(nsub_max)) ? nsub_max
                                                              : static_cast<int>(xn);
                  }
                }
                if (nsub > nsub_max) nsub = nsub_max;
                if (nsub < 1) nsub = 1;
              }
              dg_nsub = nsub;
              const Real sub_bdt = bdt/static_cast<Real>(nsub);
              const Real sub_sdt = src_relax*sub_bdt;
              const Real t0 = T_g(m,k,j,i);
              const Real d0 = rhoN(m,k,j,i);
              Real ec = ei;                            // the running internal energy
              Real Emc = Em;                           // its emission, re-formed below
              de = 0.0;
              for (int s=0; s<nsub; ++s) {
                const Real deq = (absn > 0.0) ? ec*(sqrt(sqrt(absn/Emc)) - 1.0) : -ec;
                dg_A = absn; dg_Em = Emc; dg_deq = deq;
                Real des;
                if (deq != 0.0) {
                  const Real x = sub_sdt/deq;          // >= 0: src and deq share a sign
                  des = (x > 1.0e-4) ? deq*(-expm1(-x)) : sub_sdt;
                } else {
                  des = sub_sdt;
                }
                // NEWTON REFINEMENT, seeded by the closed form above.  That estimate is
                // already the right asymptote, so this only has to correct the e ~ T it
                // assumed -- with H2 dissociating or H ionizing, e(T) is far steeper than
                // linear and the equilibrium moves.  Solve the exact backward-Euler
                // balance with the emission at the NEW temperature,
                //     F(de) = de - A dt + E(T_old) (T_new/T_old)^4 dt = 0,
                // which is the same E ~ T^4 the band solver emits with, but with T(e) and
                // c_v(e) taken from the EOS.  F is monotone in de (both terms increase),
                // so Newton from a bracketing-quality guess converges in a step or two;
                // F(0) = -src dt recovers the explicit answer if it stops immediately.
                // Under sub-cycling T_old and E(T_old) stay the values of the cell at the
                // START of the stage, so r4 is measured against the same fixed A.
                if (newton_on && eos.IsGeneral()) {
                  // both terms carry the same (1 - Lambda*), so F(de) = de -
                  // bdt*src_relax + ome*bdt*Em*(r4 - 1): explicit at r4 = 1, zero at the
                  // ALI equilibrium E(T) = A_eff.  ome_ = 1 is the old expression.
                  const Real abdt = ome_*absn*sub_bdt, embdt = ome_*Em*sub_bdt;
                  // SAFEGUARDED, because a bare Newton here does not converge.  F is
                  // MONOTONE INCREASING in de -- both de and E(T(e+de)) rise with de --
                  // so its root is unique and bracketing is available for free: every
                  // iterate with F < 0 is a lower bound on the root and every one with F
                  // > 0 an upper bound.  Without that, in an optically thin
                  // radiation-dominated cell the step is very nearly the EXPLICIT one
                  // (the damping denominator dfx -> 1 as the cell's own emission stops
                  // controlling its temperature) and it overshoots the root by orders of
                  // magnitude.  MEASURED in RG_v4/out.txt, where rt_cell_report printed
                  // the iterates: on a cell with e = 6.02e-3 the sequence was de =
                  // -1.63e+1, +1.06e+1, -3.33e+1 -- oscillating, each iterate 3-4 decades
                  // past e -- and on another, e = 3.06e-1 was handed de = -3.37e-1 on the
                  // FIRST step.  Both left e + de <= 0 and fell through to the rescue,
                  // which can only leave the cell at ~1e-3 e, from which the table EOS
                  // returns a floor temperature and the sound speed collapses the
                  // timestep.  That is the mechanism behind the He-star top-cell "dt
                  // COLLAPSE" events and the red giant's 152 eos_tclamp family.
                  //
                  // The lower bracket nlo is the 99.9 % floor the rescue used, and F(nlo)
                  // < 0 whenever the absorption A is non-negative, so it is a true lower
                  // bound on the root.  ONLY a step that would cross zero energy is
                  // replaced, by a bisection of [nlo, des]; every other step is taken
                  // exactly as before, so this is inert on every cell the bare iteration
                  // handled.  A wider safeguard was tried first -- bracketing from both
                  // sides and clamping the converged answer into the bracket -- and is
                  // NOT in the tree: it perturbed healthy cells (the 1-D He-star column's
                  // timestep fell 0.107 -> 0.034 s with no rescue anywhere in the run),
                  // because the Newton legitimately walks past the e ~ T equilibrium
                  // wherever the EOS is far from that, which is the reason the refinement
                  // exists at all.
                  const Real nlo = -(1.0 - 1.0e-3)*ec;
                  if (t0 > 0.0 && d0 > 0.0) {
                    for (int it=0; it<8; ++it) {
                      const Real e1 = ec + des;
                      if (!(e1 > 0.0)) break;
                      const Real tc = eos.Temperature(d0, e1);
                      const Real t1 = tc*eos.temp_cgs;
                      if (!(t1 > 0.0)) break;
                      const Real cv = d0*eos.SpecificHeatCv(d0, e1, tc);
                      if (!(cv > 0.0)) break;
                      const Real r4 = SQR(SQR(t1/t0));
                      const Real fx = des - abdt + embdt*r4;
                      // d(r4)/de = 4 r4/T dT/de, with dT/de = temp_cgs/(d c_v)
                      const Real dfx = 1.0 + embdt*4.0*r4/t1*(eos.temp_cgs/cv);
                      if (!(dfx > 0.0)) break;
                      const Real step = fx/dfx;
                      const Real dn = des - step;
                      // THE ONE INTERVENTION.  A step that would leave a non-positive
                      // internal energy is replaced by a bisection of [nlo, des], which
                      // brackets the root: F is monotone and F(des) > 0 is what makes the
                      // step negative in the first place, while F(nlo) < 0 for any
                      // non-negative absorption.  Every other step is taken EXACTLY as
                      // before -- dn is `des - step`, the same expression `des -= step`
                      // evaluated -- so a cell the bare iteration handled is untouched,
                      // bit for bit, and `ec + des > 1e-3 ec > 0` now holds
                      // unconditionally after the loop.
                      des = (ec + dn > 0.0) ? dn : 0.5*(nlo + des);
                      if (dg_nit < 8) dg_it[dg_nit++] = des;
                      if (fabs(step) <= 1.0e-8*(fabs(des) + fabs(ec))) break;
                    }
                  }
                  // THE LOOP CHECKS e1 > 0 AT THE TOP, NOT AT THE BOTTOM.  The last `des
                  // -= step` is never validated, so Newton can exit having pushed the
                  // cell to ec + des <= 0 -- a one-step NaN with no counterpart in the
                  // closed form, which guarantees e1 = ec*exp(-x) > 0.  R9's all-column
                  // NaN at t = 1.99e5 is under bisection; this closes the only path in
                  // this branch to a non-positive energy.  Falling back to a 99.9 % drop
                  // keeps the cell cooling hard without ever crossing zero.
                  if (!(ec + des > 0.0)) {
                    const Real de_nt = des;        // what Newton left, for the report
                    // problem/rt_rescue_eq: land on the equilibrium the cell actually
                    // sees rather than on a fixed 99.9 % drop.  deq is e_eq - e from the
                    // closed form above, i.e. Em(T_eq) = A = src + Em; the old floor is
                    // kept underneath it, so this can only make the rescue gentler.
                    const Real defl = -(1.0 - 1.0e-3)*ec;
                    if (resc_eq && deq < 0.0 && deq > defl) {
                      des = deq;
                      Kokkos::atomic_fetch_add(&efix_g(1), 1);
                    } else {
                      des = defl;
                      Kokkos::atomic_fetch_add(&efix_g(2), 1);
                    }
                    if (Kokkos::atomic_fetch_add(&efix_g(0), 1) == 0) {
                      efrec_g(0) = static_cast<Real>(m);
                      efrec_g(1) = static_cast<Real>(k);
                      efrec_g(2) = static_cast<Real>(j);
                      efrec_g(3) = static_cast<Real>(i);
                      efrec_g(4) = rhoN(m,k,j,i);
                      efrec_g(5) = ei;
                      efrec_g(6) = t0;
                      efrec_g(7) = src;
                      efrec_g(8) = src_relax;
                      efrec_g(9) = Em;
                      efrec_g(10) = absn;
                      efrec_g(11) = deq;
                      efrec_g(12) = de_nt;
                      efrec_g(13) = bdt;
                      efrec_g(14) = des;
                      efrec_g(15) = static_cast<Real>(efix_cyc);
                    }
                    dg_resc = true;
                  }
                }
                de += des;
                ec = ei + de;
                // Re-form the emission of the UPDATED state for the next sub-step.  Em
                // scales as T^4 and T comes from the EOS wherever the Newton block has
                // one; with e ~ T (ideal gas, and the closed form's own assumption) the
                // ratio is just ec/ei.  This is the whole point of the sub-cycle: without
                // it every sub-step would relax toward the same stale equilibrium.
                if (s + 1 < nsub) {
                  if (!(ec > 0.0)) break;
                  Real trat;
                  if (newton_on && eos.IsGeneral() && t0 > 0.0 && d0 > 0.0) {
                    const Real t1 = eos.Temperature(d0, ec)*eos.temp_cgs;
                    trat = (t1 > 0.0) ? t1/t0 : 0.0;
                  } else {
                    trat = ec/ei;
                  }
                  Emc = Em*SQR(SQR(trat));
                  if (!(Emc > 0.0)) break;
                }
              }
              // ...and the handover the relaxation was deliberately not shown.  Added
              // after the Newton refinement and the rescue, both of which are statements
              // about the LOCAL balance alone; src_ex is zero unless this cell sits
              // inside the tau ramp.  It is OUTSIDE the sub-cycle: it is an exact
              // explicit term, not part of the local balance being relaxed.
              de += src_ex*bdt;
            }
          }
        }
        const Real de_pre = de;
        if (demax > 0.0 && !skip_de) {
          const Real dl = LimitRTSource(de,
              eiN(m,k,j,i), demax);
          if (dl != de) { ++nc; de = dl; }
        }
        // ---- problem/rt_floor_consistent: STOP AT THE FLOOR, DO NOT CROSS IT --------
        // See the flag's note.  A cooling step is bounded by the state the EOS floors
        // would restore this cell to anyway, so the sweep and ConsToPrim stop trading the
        // same energy back and forth once per stage.  Only COOLING cells pay the EOS
        // evaluation, and only a step that would actually cross the floor is changed --
        // everything else is the identity, so with the switch off (and on, away from the
        // floors) the arithmetic is untouched.
        if (efc_ && !skip_de && de < 0.0) {
          const Real dfl_ = rhoN(m,k,j,i);
          const Real ei0_ = eiN(m,k,j,i);
          Real emin_ = eos_f.EnergyFromTemperature(dfl_, eos_f.tfloor);
          // the pressure floor, taken only where the cheap closed-form bound says the
          // cell could be near it (that bound is exact for an ideal gas)
          if (ei0_ + de < eos_f.EnergyFloorBound(dfl_)) {
            const Real ep_ = eos_f.EnergyFromPressure(dfl_, eos_f.pfloor);
            if (ep_ > emin_) emin_ = ep_;
          }
          if (ei0_ + de < emin_) {
            de = (ei0_ > emin_) ? (emin_ - ei0_) : 0.0;
          }
        }
        if (report_on) {
          const Real dxb = DX1(m,k,j,i);
          Kokkos::atomic_add(&dsum_g(0), src*bdt*dxb);
          Kokkos::atomic_add(&dsum_g(1), de*dxb);
        }
        if (diag) {
          Real Em_d = 0.0;
          for (int b=0; b<nblk; ++b) Em_d += Em_g(m,b,i,k,j);
          if (taublend) {
            Em_d *= 1.0 - 0.5*(w_g(m,k,j,i) + w_g(m,k,j,i+1));
          }
          if (band_on && i < icut_g(m,k,j)) Em_d = 0.0;
          dg(m,0,k,j,i) = src;
          dg(m,1,k,j,i) = de;
          dg(m,2,k,j,i) = (de != de_pre) ? 1.0 : 0.0;
          dg(m,3,k,j,i) = Ft;
          dg(m,4,k,j,i) = Fb;
          dg(m,5,k,j,i) = Qs_d;
          dg(m,6,k,j,i) = Em_d;
        }
        if (dbg_on && m == dbg_m && k == dbg_k && j == dbg_j &&
            i > ie - dbg_n) {
          Real Em = 0.0, Qs = 0.0;
          for (int b=0; b<nblk; ++b) { Em += Em_g(m,b,i,k,j); Qs += Qb_g(m,b,i,k,j); }
          // the raw direct source, before the blend handover and the beam are added
          Real srcraw = 0.0;
          for (int b=0; b<nblk; ++b) srcraw += Src_g(m,b,i,k,j);
          const Real ei = eiN(m,k,j,i);
          Kokkos::printf("rt_apply i=%d T=%.4e d=%.4e e=%.4e Fb=%.6e Ft=%.6e "
                         "divF=%.14e srcraw=%.14e src=%.4e Qs=%.4e Em=%.4e "
                         "lamdt=%.4e de=%.4e "
                         "de/e=%.4e tau=%.4e w=%.4e dx=%.4e "
                         "Ab=%.6e At=%.6e Ac=%.6e V=%.6e\n",
                         i, T_g(m,k,j,i), rhoN(m,k,j,i), ei, Fb, Ft,
                         -(Ft*aft_a-Fb*afb_a)/vla_a, srcraw, src, Qs, Em,
                         (Em > 0.0 && ei > 0.0) ? 4.0*Em/ei*bdt : 0.0,
                         de, de/ei, taublend ? tauf_g(m,k,j,i) : 0.0,
                         taublend ? w_g(m,k,j,i) : 0.0, DX1(m,k,j,i),
                         afb_a, aft_a, ACC(m,k,j,i), vla_a);
        }
        // ---- rt_cell_report ---------------------------------------------------
        if (report_on) {
          const bool fixedcell = fixed_on && (k == rep_k) && (j == rep_j) &&
              (fabs(X1V(m,i) - rep_r) < 0.5*DX1(m,k,j,i));
          bool doprint = fixedcell;
          if (dg_resc) {
            if (Kokkos::atomic_fetch_add(&repc_g(0), 1) == 0) doprint = true;
          }
          if (doprint && i >= icut_g(m,k,j)) {
            const Real kap = kc_g(m,0,i,k,j);
            const Real rho = rhoN(m,k,j,i);
            const Real dtc = kap*rho*DX1(m,k,j,i);
            Real tautop = 0.0;
            for (int i2=i; i2<=ie; ++i2) {
              tautop += kc_g(m,0,i2,k,j)*rhoN(m,k,j,i2)*DX1(m,k,j,i2);
            }
            // the hemispheric-mean weights the sweeps used (ck_nquad = 1); with
            // ck_nquad = 2 the split into absorption and emission below is indicative
            const Real mu1 = 1.0/CK_DIFFUSIVITY, wf = M_PI;
            const Real xq = dtc/mu1;
            const Real e0 = -expm1(-xq);
            const Real alp = (xq > 1.0e-3) ? (e0 - 1.0 + e0/xq) : (xq/2.0 - SQR(xq)/3.0);
            const Real bet = (xq > 1.0e-3) ? (1.0 - e0/xq) : (xq/2.0 - SQR(xq)/6.0);
            const Real gmq = alp;
            const Real bi = Bb_g(m,0,i,k,j);
            const Real bip = Bb_g(m,0,i+1,k,j);
            const Real bim = (i > is) ? Bb_g(m,0,i-1,k,j) : bi;
            const Real idn = idn_g(m,k,j,i+1), iup = iup_g(m,k,j,i);
            const Real absdn = wf/DX1(m,k,j,i)*e0*idn;
            const Real absup = wf/DX1(m,k,j,i)*e0*iup;
            const Real emidn = wf/DX1(m,k,j,i)*(alp*bip + bet*bi);
            const Real emiup = wf/DX1(m,k,j,i)*(bet*bip + gmq*bi);
            const Real ei = eiN(m,k,j,i);
            // the staleness itself, so a report says whether w0 and u0 had parted
            const Real ei_w0_r = w0_uc_(m,IEN,k,j,i);
            const Real ei_u0_r = EintFromCons(u0_uc_, m, k, j, i,
                                              cs_uc_ ? cosc_uc_(m,k,j) : 0.0, cs_uc_,
                                              etg_uc_,
                                              etg_uc_ ? phicc_uc_(m,k,j,i) : 0.0);
            const Real rho_w0_r = w0_uc_(m,IDN,k,j,i);
            const Real rho_u0_r = u0_uc_(m,IDN,k,j,i);
            Kokkos::printf(
              "### rt_cell_report %s ncycle=%d t=%.6e m=%d k=%d j=%d i=%d r=%.6e\n"
              "    rho=%.4e T=%.6e kap=%.4e kaprho=%.4e dtau=%.4e tau_to_top=%.4e "
              "dt=%.4e\n"
              "    I_dn(i+1)=%.6e I_dn(i)=%.6e I_up(i)=%.6e I_up(i+1)=%.6e "
              "F(i)=%.6e F(i+1)=%.6e\n"
              "    B(i-1)=%.6e B(i)=%.6e B(i+1)=%.6e T(i-1)=%.4e T(i+1)=%.4e "
              "kap(i-1)=%.3e kap(i+1)=%.3e\n"
              "    abs_dn=%.6e abs_up=%.6e emi_dn=%.6e emi_up=%.6e src=%.6e Em=%.6e "
              "A=%.6e\n"
              "    e=%.6e deq=%.6e e_eq=%.6e de=%.6e de/e=%.6e rescued=%d nit=%d\n"
              "    ei_w0=%.6e ei_u0=%.6e stale=%.4e rho_w0=%.6e rho_u0=%.6e\n"
              "    newton de: %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
              dg_resc ? "RESCUE" : "FIXED", rep_cyc, rep_time, m, k, j, i, X1V(m,i),
              rho, T_g(m,k,j,i), kap, kap*rho, dtc, tautop, bdt,
              idn, idn_g(m,k,j,i), iup, iup_g(m,k,j,i+1),
              Fb_g(m,0,i,k,j), Fb_g(m,0,i+1,k,j),
              bim, bi, bip, (i > is) ? T_g(m,k,j,i-1) : 0.0, T_g(m,k,j,i+1),
              (i > is) ? kc_g(m,0,i-1,k,j) : 0.0, kc_g(m,0,i+1,k,j),
              absdn, absup, emidn, emiup, src, dg_Em, dg_A,
              ei, dg_deq, ei + dg_deq, de, de/ei, dg_resc ? 1 : 0, dg_nit,
              ei_w0_r, ei_u0_r,
              (ei_w0_r != 0.0) ? (ei_u0_r - ei_w0_r)/ei_w0_r : 0.0,
              rho_w0_r, rho_u0_r,
              dg_it[0], dg_it[1], dg_it[2], dg_it[3],
              dg_it[4], dg_it[5], dg_it[6], dg_it[7]);
          }
        }
        const Real de_app = de;
        if (skip_de) {
          // NOTHING is applied here: the tridiagonal owns the energy update, and it is
          // called by the wrapper as soon as this kernel is done.  de is kept only for
          // the diagnostics above, where it is the EXPLICIT rate this state would give.
        } else {
          u0(m,IEN,k,j,i) += de_app;
        }
        // ---- radiative momentum source, see rt_rad_force --------------------------
        if (radforce) {
          const Real rho = rhoN(m,k,j,i);
          const Real tkc = T_g(m,k,j,i);
          Real wr, dwdx, dwdy;
          rad_taper::WeightGated(log10(rho), xlo_f, xhi_f, log10(tkc), ylo_f, yhi_f,
                                 tg_f, wr, dwdx, dwdy);
          // the cell's net two-stream flux, positive upward
          const Real fnet = 0.5*(Ft + Fb);
          // problem/rt_force_tau_gate: s(tau) with tau the grey optical depth to the top
          // of this column.  The ingredients are the ones the `tautop` diagnostic uses,
          // kappa_grey rho dr summed over the cells above, plus half of this cell's own
          // dtau; the sum is formed here because no per-cell tau-to-top array exists
          // outside the tau blend (pcond_rt->rad_tauf, which is a FACE tau and only
          // allocated when taublend is on).  It runs only with the gate on.
          Real sgt_ = 0.0;
          if (ftg_) {
            Real tau_ = 0.5*kc_g(m,0,i,k,j)*rhoN(m,k,j,i)*DX1(m,k,j,i);
            // s = 0 once tau passes tau_hi: deep cells leave after a cell or two
            for (int i2=i+1; i2<=ie && tau_ < ft_hi; ++i2) {
              tau_ += kc_g(m,0,i2,k,j)*rhoN(m,k,j,i2)*DX1(m,k,j,i2);
            }
            sgt_ = RTForceTauGate(tau_, ftlg_lo, ftlg_hi);
          }
          Real f1 = 0.0, f2 = 0.0, f3 = 0.0, prgw = 0.0;
          if (sgt_ > 0.0) {
            // ---- the GATED form (problem/rt_force_tau_gate) ------------------------
            //   f = s [ rho kappa F/c + grad(w Prad) ] + (1-s) [ old ],
            // with grad(w Prad) a centred difference of w(rho,T) a T^4/3 over the same
            // neighbours and the same spacings the old Prad grad w term differences:
            // X1V(i+1)-X1V(i-1) radially, the arc lengths 2 DX2 / 2 DX3 transversely.
            const Real kfc = rho*kc_g(m,0,i,k,j)*fnet*inv_c;
            const Real wpm = RTWPrad(rhoN(m,k,j,i-1), T_g(m,k,j,i-1), xlo_f, xhi_f,
                                     ylo_f, yhi_f, tg_f, arad_f);
            const Real wpp = RTWPrad(rhoN(m,k,j,i+1), T_g(m,k,j,i+1), xlo_f, xhi_f,
                                     ylo_f, yhi_f, tg_f, arad_f);
            Real g1 = kfc + (wpp - wpm)/(X1V(m,i+1) - X1V(m,i-1));
            Real g2 = 0.0, g3 = 0.0;
            if (md_f) {
              g2 = (RTWPrad(rhoN(m,k,j+1,i), T_g(m,k,j+1,i), xlo_f, xhi_f,
                            ylo_f, yhi_f, tg_f, arad_f)
                  - RTWPrad(rhoN(m,k,j-1,i), T_g(m,k,j-1,i), xlo_f, xhi_f,
                            ylo_f, yhi_f, tg_f, arad_f))/(2.0*DX2(m,k,j,i));
            }
            if (td_f) {
              g3 = (RTWPrad(rhoN(m,k+1,j,i), T_g(m,k+1,j,i), xlo_f, xhi_f,
                            ylo_f, yhi_f, tg_f, arad_f)
                  - RTWPrad(rhoN(m,k-1,j,i), T_g(m,k-1,j,i), xlo_f, xhi_f,
                            ylo_f, yhi_f, tg_f, arad_f))/(2.0*DX3(m,k,j,i));
            }
            // the ungated expressions, for the (1-s) share; identical algebra to the
            // branch below, kept separate so that branch stays untouched
            Real o1 = (1.0 - wr)*kfc;
            Real o2 = 0.0, o3 = 0.0;
            if (dwdx != 0.0 || dwdy != 0.0) {
              const Real tk = tkc;
              const Real cg = (arad_f*tk*tk*tk*tk/3.0)*dwdx*M_LOG10E/rho;
              const Real ct = (dwdy != 0.0)
                            ? (arad_f*tk*tk*tk*tk/3.0)*dwdy*M_LOG10E/tk : 0.0;
              prgw = cg*(rhoN(m,k,j,i+1) - rhoN(m,k,j,i-1))/(X1V(m,i+1) - X1V(m,i-1));
              if (ct != 0.0) {
                prgw += ct*(T_g(m,k,j,i+1) - T_g(m,k,j,i-1))/(X1V(m,i+1) - X1V(m,i-1));
              }
              o1 += prgw;
              if (md_f) {
                o2 = cg*(rhoN(m,k,j+1,i) - rhoN(m,k,j-1,i))/(2.0*DX2(m,k,j,i));
                if (ct != 0.0) {
                  o2 += ct*(T_g(m,k,j+1,i) - T_g(m,k,j-1,i))/(2.0*DX2(m,k,j,i));
                }
              }
              if (td_f) {
                o3 = cg*(rhoN(m,k+1,j,i) - rhoN(m,k-1,j,i))/(2.0*DX3(m,k,j,i));
                if (ct != 0.0) {
                  o3 += ct*(T_g(m,k+1,j,i) - T_g(m,k-1,j,i))/(2.0*DX3(m,k,j,i));
                }
              }
            }
            const Real oms = 1.0 - sgt_;
            f1 = sgt_*g1 + oms*o1;
            f2 = sgt_*g2 + oms*o2;
            f3 = sgt_*g3 + oms*o3;
          } else {
          // ---- the UNGATED form.  Left at its original indentation and untouched,
          // so that with rt_force_tau_gate off the arithmetic executed is bitwise the
          // arithmetic of every run before the gate existed.
          f1 = (1.0 - wr)*rho*kc_g(m,0,i,k,j)*fnet*inv_c;
          if (dwdx != 0.0 || dwdy != 0.0) {
            // Prad grad w, with grad w = w'(x) grad rho/(rho ln10) -- the SAME derivative
            // the EOS put into chi_rho, so the two cancel in the continuum limit
            const Real tk = tkc;
            const Real cg = (arad_f*tk*tk*tk*tk/3.0)*dwdx*M_LOG10E/rho;
            // ...and, under the temperature gate, w'(y) grad T/(T ln10), the partner of
            // the dw/dlnT the EOS put into chi_T.  Exactly one of cg and ct is non-zero,
            // because max() picks one branch of the weight.
            const Real ct = (dwdy != 0.0)
                          ? (arad_f*tk*tk*tk*tk/3.0)*dwdy*M_LOG10E/tk : 0.0;
            prgw = cg*(rhoN(m,k,j,i+1) - rhoN(m,k,j,i-1))/(X1V(m,i+1) - X1V(m,i-1));
            if (ct != 0.0) {
              prgw += ct*(T_g(m,k,j,i+1) - T_g(m,k,j,i-1))/(X1V(m,i+1) - X1V(m,i-1));
            }
            f1 += prgw;
            // the transverse components of Prad grad w.  The derivative is taken
            // against the physical ARC LENGTH (DX2/DX3), not the coordinate spacing:
            // on the cubed sphere x2/x3 are the gnomonic ANGLES and size.dx2 is
            // dimensionless, which made this term ~r = 2e11 times too large and drove
            // a grid-scale transverse instability from round-off in a few cycles.
            // What u0(IM2), u0(IM3) hold on the cubed sphere is the COVARIANT angular
            // momentum (coordinates.hpp:151-161), and the covariant component of a
            // gradient along the unit tangent e_xi IS dphi/dl_xi, so no metric
            // inverse enters here -- the arc-length derivative is already the right
            // component to add.
            if (md_f) {
              f2 = cg*(rhoN(m,k,j+1,i) - rhoN(m,k,j-1,i))/(2.0*DX2(m,k,j,i));
              if (ct != 0.0) {
                f2 += ct*(T_g(m,k,j+1,i) - T_g(m,k,j-1,i))/(2.0*DX2(m,k,j,i));
              }
            }
            if (td_f) {
              f3 = cg*(rhoN(m,k+1,j,i) - rhoN(m,k-1,j,i))/(2.0*DX3(m,k,j,i));
              if (ct != 0.0) {
                f3 += ct*(T_g(m,k+1,j,i) - T_g(m,k-1,j,i))/(2.0*DX3(m,k,j,i));
              }
            }
          }
          }
          const Real dinv = 1.0/u0(m,IDN,k,j,i);
          const Real v1 = u0(m,IM1,k,j,i)*dinv;
          const Real v2 = u0(m,IM2,k,j,i)*dinv;
          const Real v3 = u0(m,IM3,k,j,i)*dinv;
          u0(m,IM1,k,j,i) += f1*bdt;
          if (md_f) u0(m,IM2,k,j,i) += f2*bdt;
          if (td_f) u0(m,IM3,k,j,i) += f3*bdt;
          u0(m,IEN,k,j,i) += (v1*f1 + v2*f2 + v3*f3)*bdt;
          if (budg_) {
            const Real dvb = DX1(m,k,j,i)*size.d_view(m).dx2*size.d_view(m).dx3;
            Kokkos::atomic_add(&bud_(10), (v1*f1 + v2*f2 + v3*f3)*bdt*dvb);
          }
          if (fverb && wr > 0.0 && wr < 1.0) {
            // the hydrostatic balance of this cell: pressure gradient, gravity, source
            EOSThermoState sm, sp;
            eos_f.tbl.EvalNoMu(rhoN(m,k,j,i-1), T_g(m,k,j,i-1), sm);
            eos_f.tbl.EvalNoMu(rhoN(m,k,j,i+1), T_g(m,k,j,i+1), sp);
            const Real ap = -(sp.p - sm.p)/((X1V(m,i+1) - X1V(m,i-1))*rho);
            const Real gloc = !pp_
                ? EffGravAt(grav, apo_, X1V(m,i), grav_pmass, omega,
                            cfv_on ? cfv_(m,k,j,3) : 0.0, tide)
                  *(cfv_on ? cfv_(m,k,j,4) : 1.0)
                : gver;
            const Real gg = (gloc > 0.0) ? gloc : ((gver > 0.0) ? gver : 1.0);
            Kokkos::printf("### rt_force ncycle=%d i=%d rho=%.4e T=%.4e w=%.4f "
                           "a_p/g=%.6e a_f/g=%.6e a_prgw/g=%.6e resid/g=%.6e "
                           "F=%.4e kap=%.4e\n",
                           vcyc, i, rho, T_g(m,k,j,i), wr, ap/gg, f1/(rho*gg),
                           prgw/(rho*gg), (ap - gg + f1/rho)/gg,
                           0.5*(Ft + Fb), kc_g(m,0,i,k,j));
          }
        }
      }, ck_impl_nosync && ckimp_);   // ck_implicit never clips here (skip_de)
      rt_nclip_last = nclip;
      RTSourceLimiterWarn(nclip);
      // ---- the per-cycle clip/rescue CENSUS (problem/rt_outer_verbose) -------------
      // How many cells the per-STEP caps actually touched this stage: nclip is the
      // LimitRTSource count reduced over the apply kernel, efix(0..2) the running Newton
      // positivity rescue totals.  Diagnostic only; nothing here changes an answer.
      if (rt_outer_verbose && global_variable::my_rank == 0 && pm->ncycle % 200 == 0) {
        auto hcen = Kokkos::create_mirror_view(efix_g);
        Kokkos::deep_copy(hcen, efix_g);
        std::cout << "### rt_clip ncycle=" << pm->ncycle << " t=" << pm->time
                  << " nclip=" << nclip << " efix_tot=" << hcen(0)
                  << " resc_eq=" << hcen(1) << " resc_floor=" << hcen(2) << std::endl;
      }
      // The Newton positivity rescue should never fire.  Say so the first time it does,
      // with the running total, and stay quiet afterwards.  (problem/ck_impl_nosync:
      // checked on the first pass of a call only.)
      if (!(ck_impl_nosync && ck_impl_pass > 0)) {
        static bool efix_warned = false;
        static int efix_seen = 0;
        auto he = Kokkos::create_mirror_view(efix_g);
        Kokkos::deep_copy(he, efix_g);
        if (he(0) > efix_seen && !efix_warned) {
          std::cout << "### two_stream_rt: the Newton step left e <= 0 in " << he(0)
                    << " cell(s): " << he(1) << " rescued to the radiative equilibrium "
                    << "(problem/rt_rescue_eq), " << he(2) << " to the 99.9 % floor; "
                    << "this is reported once" << std::endl;
          auto hr = Kokkos::create_mirror_view(efrec_g);
          Kokkos::deep_copy(hr, efrec_g);
          std::cout << "    first such cell: cycle " << static_cast<int>(hr(15))
                    << " (m,k,j,i) = (" << static_cast<int>(hr(0)) << ","
                    << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
                    << static_cast<int>(hr(3)) << ")  rho = " << hr(4)
                    << "  e = " << hr(5) << "  T = " << hr(6) << " K" << std::endl
                    << "      src = " << hr(7) << "  src_relax = " << hr(8)
                    << "  Em = " << hr(9) << "  A = src_relax+Em = " << hr(10)
                    << "  bdt = " << hr(13) << std::endl
                    << "      src*bdt/e = " << (hr(5) != 0.0 ? hr(7)*hr(13)/hr(5) : 0.0)
                    << "  Em*bdt/e = " << (hr(5) != 0.0 ? hr(9)*hr(13)/hr(5) : 0.0)
                    << "  deq/e = " << (hr(5) != 0.0 ? hr(11)/hr(5) : 0.0)
                    << "  de_newton/e = " << (hr(5) != 0.0 ? hr(12)/hr(5) : 0.0)
                    << "  de_used/e = " << (hr(5) != 0.0 ? hr(14)/hr(5) : 0.0)
                    << std::endl;
          efix_warned = true;
        }
        efix_seen = he(0);
      }
      // the source's energy budget: what the flux divergence asked for against what the
      // (sub-cycled, relaxed, limited) application actually deposited.  See rt_desum_ptr.
      if (report_on && fixed_on) {
        auto hd = Kokkos::create_mirror_view(dsum_g);
        Kokkos::deep_copy(hd, dsum_g);
        std::cout << "### rt_desum ncycle=" << pm->ncycle << " sum(src*dt*dx) = " << hd(0)
                  << "  sum(de*dx) = " << hd(1) << "  rel = "
                  << ((hd(0) != 0.0) ? (hd(1) - hd(0))/hd(0) : 0.0) << std::endl;
      }

      // ---- one-shot column dump, for cross-code comparison -------------------------
      if (band_on && !rt_dump_file.empty() && !rt_dump_done) {
        rt_dump_done = true;
        const int jd = (rt_dump_j >= 0) ? rt_dump_j : (js + je)/2;
        const int kd = (rt_dump_k >= 0) ? rt_dump_k : (ks + ke)/2;
        const int md = (rt_dump_m <= nmb1) ? rt_dump_m : 0;
        // ck-nq2: the correlated-k dump also carries rho and, per BAND, the net
        // longwave flux F_b (face) and the cell-centre Planck intensity B_b, so the
        // thick limit F_b = (4 pi/3) dB_b/dtau can be checked band by band.  A chain
        // block never straddles a band (RT_NB divides CK_NG), so F_b is the sum of the
        // blocks of that band.  Appended columns: the first 12 are unchanged.
        // ... and five Rosseland diagnostics (ck_build_rosseland_table vs the opacity
        // the sweep actually sees): 35 kappa_R of the cell from the sweep's own
        // (xT, xP, kc_g) inputs, 36 the table's kappa_R at the cell (T, p), 37 the
        // exact NON-GREY discrete diffusion flux at the face, sum_bg (4 pi/3) gw dB_b /
        // dtau_bg over the same centre-to-centre layer the sweep uses, 38 the smallest
        // optical depth from the top over all (band, g) at the face, 39 the Rosseland
        // (grey) diffusion flux at the face from column 35.
        const int ncold0 = grey_on ? 12 : (18 + 2*CK_NB);
        const int ncold = (grey_on || !rt_dump_kap) ? ncold0
                        : (ncold0 + CK_NB*CK_NG + CK_NG);
        const bool dtab = (pcond_rt != nullptr) && pcond_rt->rad_kappa_tab &&
                          pcond_rt->rad_kr_nT > 0 && !pcond_rt->rad_kappa_rho;
        auto dkt = dtab ? pcond_rt->rad_kr_tab : CkDum<DvceArray2D<Real>>("d");
        auto dklT = dtab ? pcond_rt->rad_kr_lT : CkDum<DvceArray1D<Real>>("d");
        auto dklP = dtab ? pcond_rt->rad_kr_lP : CkDum<DvceArray1D<Real>>("d");
        const int dnT = dtab ? pcond_rt->rad_kr_nT : 0;
        const int dnP = dtab ? pcond_rt->rad_kr_nP : 0;
        const int bpb = (ck_nq_ == 1 ? CK_NG : 2*CK_NG)/RT_NB;   // blocks per band
        DvceArray2D<Real> col("rt_col", n1, ncold);
        par_for("rt_dumpcol", DevExeSpace(), is, ie+1, KOKKOS_LAMBDA(const int i) {
          Real Fs = 0.0;
          Real Qs = 0.0;
          Real Ss = 0.0;
          for (int b=0; b<nblk; ++b) {
            Fs += Fb_g(md,b,i,kd,jd);
            Qs += Qb_g(md,b,i,kd,jd);
            Ss += Src_g(md,b,i,kd,jd);
          }
          if (ncold > 12) {
            col(i,12) = rhoN(md,kd,jd,(i <= ie) ? i : ie);
            for (int bb=0; bb<CK_NB; ++bb) {
              Real fb = 0.0;
              for (int b=bb*bpb; b<(bb+1)*bpb && b<nblk; ++b) fb += Fb_g(md,b,i,kd,jd);
              col(i,13+bb) = fb;
              col(i,13+CK_NB+bb) = Bb_g(md,bb,i,kd,jd);
            }
            // the sweep's kappa of (band, g) at cell ii, cm^2/g
            auto kbg = [&](const int ii, const int b, const int g) {
              const Real xTv = xT_g(md,kd,jd,ii), xPv = xP_g(md,kd,jd,ii);
              const int iT = static_cast<int>(xTv), iP = static_cast<int>(xPv);
              return ck_kappa(cklk, iT, xTv - iT, iP, xPv - iP, b, g)
                     + kc_g(md,b,ii,kd,jd);
            };
            auto kros = [&](const int ii) {
              const Real TT = T_g(md,kd,jd,ii);
              const Real sp = boltz_sigma*SQR(SQR(1.01*TT));
              const Real sm = boltz_sigma*SQR(SQR(0.99*TT));
              Real num = 0.0, den = 0.0;
              for (int b=0; b<CK_NB; ++b) {
                const Real wb = sp*ck_planck_frac(ckpf, pfl0, pfid, 1.01*TT, b)
                              - sm*ck_planck_frac(ckpf, pfl0, pfid, 0.99*TT, b);
                if (!(wb > 0.0)) continue;
                Real inv = 0.0;
                for (int g=0; g<CK_NG; ++g) inv += ckgw(g)/kbg(ii, b, g);
                num += wb*inv;
                den += wb;
              }
              return den/num;
            };
            const int ic = (i <= ie) ? i : ie;
            col(i,13+2*CK_NB) = kros(ic);
            col(i,14+2*CK_NB) = dtab ? RosselandTable(dkt, dklT, dklP, dnT, dnP,
                                                      T_g(md,kd,jd,ic),
                                                      pb_g(md,kd,jd,ic)*1.0e6) : 0.0;
            Real fd = 0.0, fr = 0.0, tmn = 1.0e300;
            if (i > is && i <= ie) {
              const Real hl = 0.5*rhoN(md,kd,jd,i-1)*DX1(md,kd,jd,i-1);
              const Real hu = 0.5*rhoN(md,kd,jd,i)*DX1(md,kd,jd,i);
              Real bsl = 0.0, bsu = 0.0;
              for (int b=0; b<CK_NB; ++b) {
                const Real db = Bb_g(md,b,i,kd,jd) - Bb_g(md,b,i-1,kd,jd);
                bsl += Bb_g(md,b,i-1,kd,jd);
                bsu += Bb_g(md,b,i,kd,jd);
                for (int g=0; g<CK_NG; ++g) {
                  const Real dtc = kbg(i-1, b, g)*hl + kbg(i, b, g)*hu;
                  fd -= (4.0*M_PI/3.0)*ckgw(g)*db/dtc;
                  Real tt = 0.0;
                  for (int i2=i; i2<=ie; ++i2) {
                    tt += kbg(i2, b, g)*rhoN(md,kd,jd,i2)*DX1(md,kd,jd,i2);
                  }
                  tmn = (tt < tmn) ? tt : tmn;
                }
              }
              fr = -(4.0*M_PI/3.0)*(bsu - bsl)/(kros(i-1)*hl + kros(i)*hu);
            }
            if (ncold > ncold0) {
              for (int b=0; b<CK_NB; ++b) {
                for (int g=0; g<CK_NG; ++g) col(i,ncold0+b*CK_NG+g) = kbg(ic, b, g);
              }
              for (int g=0; g<CK_NG; ++g) col(i,ncold0+CK_NB*CK_NG+g) = ckgw(g);
            }
            col(i,15+2*CK_NB) = fd;
            col(i,16+2*CK_NB) = (tmn < 1.0e300) ? tmn : 0.0;
            col(i,17+2*CK_NB) = fr;
          }
          // THE ENERGY BUDGET OF THE COLUMN, in the solver's own geometry: the thermal
          // deposit per unit volume, the cell's FRAME area A_i = V_i/dx_i and its
          // volume.  With these three the telescoping identity of the spherical form,
          //     sum_i V_i Src_i = -(A_top F_top - A_bot F_bot),
          // is checkable from the dump alone, on any mesh, with no analytic
          // reconstruction of the geometry (see tests_ck_sph/README.md).
          col(i,9)  = Ss;
          col(i,10) = ACC(md,kd,jd,i);
          col(i,11) = VLA(md,kd,jd,i);
          col(i,0) = pb_g(md,kd,jd,i);
          col(i,1) = T_g(md,kd,jd,i);
          col(i,2) = Fs;
          col(i,3) = Qs;
          col(i,4) = X1F(md,i);
          // Gamma_1 and grad_ad of the CURRENT state. Both collapse to the ideal values
          // when the EOS is ideal; under the tabulated EOS they dip hard through the H2
          // dissociation and H ionization bands, and how steeply they vary across a cell
          // is what the reconstruction has to cope with.
          const Real dd = rhoN(md,kd,jd,i);
          const Real ee = eiN(md,kd,jd,i);
          col(i,5) = eos.IsGeneral() ? eos.Gamma1(dd, ee) : eos.gamma;
          col(i,6) = GradAd(eos, eos.gamma, Rgas, pb_g(md,kd,jd,i)*1.0e6,
                            T_g(md,kd,jd,i));
          col(i,7) = taublend ? tauf_g(md,kd,jd,i) : 0.0;
          col(i,8) = taublend ? w_g(md,kd,jd,i) : 0.0;
        });
        auto hc = Kokkos::create_mirror_view(col);
        Kokkos::deep_copy(hc, col);
        auto hcut = Kokkos::create_mirror_view(icut_g);
        Kokkos::deep_copy(hcut, icut_g);
        auto hcf = Kokkos::create_mirror_view(cf_g);
        Kokkos::deep_copy(hcf, cf_g);
        if (global_variable::my_rank == 0) {
          std::ofstream f(rt_dump_file);
          f.precision(10);
          f << std::scientific;
          f << (grey_on ? "# two_stream_rt GREY column dump\n"
                        : "# deep_hot_jupiter_rt correlated-k column dump\n")
            << "# meshblock " << md << ", k = " << kd << ", j = " << jd
            << ", mu0 = " << hcf(md,kd,jd,3) << ", icut = " << hcut(md,kd,jd)
            << " (is = " << is << ", ie = " << ie << ")\n"
            << "# T_int = " << Tint << " K, T_irr = " << Tirr
            << " K, grav = " << grav << " cm/s^2\n"
            << "# fluxes are cgs: 1 erg/s/cm^2 = 1e-3 W/m^2. F_lw is NET (up minus "
                << "down)\n"
            << "# i  r_face[cm]  p[bar]  T[K]  F_lw_net[erg/s/cm2]  Q_sw[erg/s/cm3]"
            << "  Gamma_1  grad_ad  tau_R(face)  w_diff(face)"
            << "  Src_lw[erg/s/cm3]  A_cell[cm2]  V_cell[cm3]"
            << ((ncold > 12) ? ("  rho[g/cm3]  F_b(b=0..10)  B_b(b=0..10)  kR_sweep  "
                                "kR_table  F_diff_nongrey  tau_min_bg  F_diff_R") : "")
            << "\n";
          if (ncold > 12) f << "# ck_nquad = " << ck_nq_ << "\n";
          for (int i=is; i<ie+2; ++i) {
            f << i << " " << hc(i,4) << " " << hc(i,0) << " " << hc(i,1)
              << " " << hc(i,2) << " " << hc(i,3)
              << " " << hc(i,5) << " " << hc(i,6)
              << " " << hc(i,7) << " " << hc(i,8)
              << " " << hc(i,9) << " " << hc(i,10) << " " << hc(i,11);
            for (int q=12; q<ncold; ++q) f << " " << hc(i,q);
            f << "\n";
          }
          f.close();
          std::cout << "two_stream_rt: wrote "
                    << (grey_on ? "grey" : "correlated-k") << " column dump to '"
                    << rt_dump_file << "' (k = " << kd << ", j = " << jd
                    << ", mu0 = " << hcf(md,kd,jd,3) << ")" << std::endl;
        }
      }
      return;
    }

//    size_t scr_size = 8 * ScrArray1D<Real>::shmem_size(n1);
//    int scr_level = 0;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, scr_level,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    // The private radial arrays here were a flat 270, with no check. Production runs at
    // nx1 = 256, i.e. n1 = 260, so they were ten cells from silently overrunning per-
    // thread
    // memory -- and nghost = 4 with the same nx1 would have gone over. Dispatched at run
    // time over compile-time tiers instead, as the correlated-k chain kernel is.
    int nclip_grey = 0;
    const Real demax_grey = rt_de_max;
    // --- problem/nan_report (see rt_nan_report): the in-kernel catcher.  One int
    // counter and one 16-slot record, zeroed per call; the FIRST offending cell on the
    // rank wins the record.  Everything it needs is recomputed inside the `bad` branch,
    // so the fast path pays nothing but a bool test.
    const bool nanrep_g = rt_nan_report;
    if (nanrep_g && rt_nanrep_cnt == nullptr) {
      rt_nanrep_cnt = new DvceArray1D<int>("rt_nanrep_cnt", 1);
      rt_nanrep_rec = new DvceArray1D<Real>("rt_nanrep_rec", 16);
    }
    auto nrcnt = nanrep_g ? *rt_nanrep_cnt : CkDum<DvceArray1D<int>>("d");
    auto nrrec = nanrep_g ? *rt_nanrep_rec : CkDum<DvceArray1D<Real>>("d");
    if (nanrep_g) {
      Kokkos::deep_copy(nrcnt, 0);
      Kokkos::deep_copy(nrrec, 0.0);
    }
    // nvcc: no extended (device) lambda may sit inside a generic lambda, so the
    // body of this launcher is the function template TsrtGreyRt above
    // picket_fence_two_stream_RT_pass; grey_rt_ctx hands it, by reference, the
    // locals it reads (unpacked there in the same order).
    auto grey_rt_ctx = std::forward_as_tuple(albedo, ap, area1, bdt, boltz_sigma,
                                             demax_grey, dx1, eiN, eos, Fstar, gm1, grav,
                                             grav_pmass, iap, ie, Iint, is, je, js, ke,
                                             ks, mbpanel_, met, mug, nanrep_g, nchain_rt,
                                             nclip_grey, nmb1, nrcnt, nrrec, omega, Rgas,
                                             rhoN, stcl_g, Teq, tide, Tint4, Tirr4, u0,
                                             use_cubed_sphere_, use_spherical_polar,
                                             volume, wg, wtemp_, x1f_, x1v_, x2v_, x3v_);
    // tusplit-1102: TsrtGreyRt is instantiated for this type (two_stream_rt.hpp)
    static_assert(std::is_same_v<decltype(grey_rt_ctx), TsrtGreyRtCtx>,
                  "grey_rt_ctx != TsrtGreyRtCtx: update two_stream_rt.hpp");
    auto launch_grey_rt = [&](auto nn_tag) {
      TsrtGreyRt(grey_rt_ctx, nn_tag);
    };
    if (n1 <= 72) {
      launch_grey_rt(std::integral_constant<int, 72>{});
    } else if (n1 <= 136) {
      launch_grey_rt(std::integral_constant<int, 136>{});
    } else if (n1 <= 264) {
      launch_grey_rt(std::integral_constant<int, 264>{});
    } else if (n1 <= 520) {
      launch_grey_rt(std::integral_constant<int, 520>{});
    } else if (n1 <= 1032) {
      launch_grey_rt(std::integral_constant<int, 1032>{});
    } else {
      std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: n1 = " << n1
                << " exceeds the largest radial tier (1032). Add a tier to the dispatch "
                << "in picket_fence_two_stream_RT." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    RTSourceLimiterWarn(nclip_grey);
    if (nanrep_g && rt_nanrep_lines < rt_nanrep_maxlines) {
      auto hgc = Kokkos::create_mirror_view(nrcnt);
      Kokkos::deep_copy(hgc, nrcnt);
      if (hgc(0) > 0) {
        auto hgr = Kokkos::create_mirror_view(nrrec);
        Kokkos::deep_copy(hgr, nrrec);
        ++rt_nanrep_lines;
        std::cout << "### rg nan_report [RT_grey_apply] rank "
                  << global_variable::my_rank << " cycle " << pm->ncycle
                  << " t = " << pm->time << ": " << hgc(0) << " cell(s); first (m,k,j,i)"
                  << " = (" << static_cast<int>(hgr(0)) << ","
                  << static_cast<int>(hgr(1)) << "," << static_cast<int>(hgr(2)) << ","
                  << static_cast<int>(hgr(3)) << ")"
                  << " d = " << hgr(4) << " e_pre = " << hgr(5) << " T = " << hgr(6)
                  << " rho_kap = " << hgr(7) << " du = " << hgr(8)
                  << " du_flux = " << hgr(9) << " u_new = " << hgr(10)
                  << " e[i-1] = " << hgr(11) << " T[i-1] = " << hgr(12)
                  << " e[i+1] = " << hgr(13) << " T[i+1] = " << hgr(14)
                  << " tau = " << hgr(15) << std::endl;
      }
    }

    return;
}

}  // namespace two_stream_rt
