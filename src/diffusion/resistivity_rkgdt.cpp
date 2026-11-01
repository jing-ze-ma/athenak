//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_rkgdt.cpp
//! \brief rkgdt-1101: the RKG super-step's diffusive dt from the forward-Euler limit of the
//! resistive operator (<mhd>/rkg_dt_safety), and a power-iteration check of that limit on
//! the actual operator (<mhd>/rkg_lmax_diag).
//!
//! WHY. RKG2 (resistivity_rkg_tasks.cpp TotStage: s^2 + s - 2 = 4 dt/dt_diff) is stable
//! when dt_diff is the forward-Euler limit 2/lambda_max of the discrete operator.  The old
//! dt_diff = cfl_number x (1/6) min_d dx_d^2/eta is the isotropic cubic-cell value times
//! the hyperbolic CFL number.  In thin radial cells (dr << r dtheta, WASP-121b C32:
//! dh/dr ~ 34) the limit is ~dr^2/(2 eta), so the old value is ~10x small and the stage
//! count ~3.2x large.  For a fixed dt the RKG2 result is second order whatever s is; s only
//! sets the stability margin.
//!
//! THE BOUND (per cell, maximised over the mesh):
//!   lambda_c = 4 eta_r (1/dx1)(1/dc_- + 1/dc_+)/2 + 4 eta_h (1/dx2^2 + 1/dx3^2)/s^2
//! dc_-/+ = centroid spacings to the radial neighbours (uniform grid: the first term is
//! 4 eta/dr^2, the 1-D Nyquist eigenvalue).  eta_r = max of eta_b over the cells whose
//! x2/x3 edges carry the radial derivative (3 radial x 5 transverse; the edge eta is the
//! mean of the four cells around the edge), eta_h = max over the 3x3x3 block.  The
//! diagonal (corner) ghosts -- at a cube vertex no neighbour fills them -- enter only
//! eta_h: in eta_r they pinned the bound at ~max_eta/dr_min^2 on every state (a t = 0
//! state with true lambda 0.108 gave 14.4).  s = sin of the cs cell angle (horizontal
//! metric, conservative); on sp dx3 is taken at the theta face nearer the pole (pole
//! cells: sin(dtheta/2)), since the modes of the pole rows exceed the centroid value.
//! The cs 4-point radial stencil (RadialFaceWeights)
//! reduces to the 2-point one on a uniform grid; on the stretched grid the power
//! iteration measures the difference (WASP-121b r1m: bound/lambda = 1.07).

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "bvals/bvals.hpp"
#include "pgen/pgen.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn Real Resistivity::LambdaBoundMax
//! \brief max over the mesh (all ranks) of the per-cell bound lambda_c above; optionally
//! the radius and eta of the cell where it is attained.

Real Resistivity::LambdaBoundMax(Real *rloc, Real *etaloc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  auto dx1_ = pmy_pack->pcoord->dx1;
  auto dx2_ = pmy_pack->pcoord->dx2;
  auto dx3_ = pmy_pack->pcoord->dx3;
  auto x1v_ = pmy_pack->pcoord->x1v;
  auto scell = pmy_pack->pcoord->sin_cell;
  auto x2f_ = pmy_pack->pcoord->xx2f;
  auto x3f_ = pmy_pack->pcoord->xx3f;
  const bool sp = pmy_pack->pmesh->use_spherical_polar;
  auto eta_ = eta_b;
  // ambipolar-1101 merge: eta_O + eta_A per cell (eta_A = min(ad_coef |B|^2, max_eta_ad),
  // the coefficient AddEMFAmbipolar uses).  The AD operator eta_A (J - (J.B)B/B^2) is the
  // Ohmic one with J projected perpendicular to B, so its frozen-coefficient spectrum is
  // bounded by the isotropic one with eta_A added.
  const bool ad_on = (ad_type != 0);
  auto coef_ = ad_coef;
  auto bcc_ = pmy_pack->pmhd->bcc0;
  const Real cap_ad = max_eta_ad;

  using ML = Kokkos::MaxLoc<Real, int>;
  ML::value_type res;
  Kokkos::parallel_reduce("res_rkgdt_bound", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, ML::value_type &v) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    const int k0 = three_d ? k-1 : k, k1 = three_d ? k+1 : k;
    const int j0 = multi_d ? j-1 : j, j1 = multi_d ? j+1 : j;
    // eta_r: the cells of the x2/x3 edges that carry the radial derivative (one transverse
    // offset at a time: the edge eta is the mean of the 4 cells around the edge); eta_h:
    // the whole 3x3x3 block, for the horizontal term (the radial x1 edges at a cube vertex
    // also average the vertex ghost, which no neighbour fills)
    Real er = 0.0, eh = 0.0;
    for (int kk = k0; kk <= k1; ++kk) {
      for (int jj = j0; jj <= j1; ++jj) {
        const bool diag = (kk != k) && (jj != j);
        for (int ii = i-1; ii <= i+1; ++ii) {
          Real e = eta_(m,kk,jj,ii);
          if (ad_on) {
            const Real b2 = SQR(bcc_(m,IBX,kk,jj,ii)) + SQR(bcc_(m,IBY,kk,jj,ii))
                          + SQR(bcc_(m,IBZ,kk,jj,ii));
            e += fmin(coef_(m,kk,jj,ii)*b2, cap_ad);
          }
          eh = fmax(eh, e);
          if (!diag) {er = fmax(er, e);}
        }
      }
    }
    const Real dcm = x1v_(m,i) - x1v_(m,i-1);
    const Real dcp = x1v_(m,i+1) - x1v_(m,i);
    Real lam = 2.0*(1.0/dcm + 1.0/dcp)/dx1_(m,k,j,i);
    Real h = 0.0;
    if (multi_d) {h += 4.0/SQR(dx2_(m,k,j,i));}
    if (three_d) {
      Real d3 = dx3_(m,k,j,i);
      if (sp) {
        // the phi length at the theta face nearer the pole, not at the centroid; a cell
        // on the pole takes sin(dtheta/2) (sp_test iprob 11: the pole cells' mode is 1.43x
        // the centroid value 4 eta/dx3^2, this bound 1.36x above it)
        const Real s0 = fabs(sin(x2f_(m,j))), s1 = fabs(sin(x2f_(m,j+1)));
        const Real sh = fabs(sin(0.5*(x2f_(m,j+1) - x2f_(m,j))));
        const Real se = fmax(fmin(s0, s1), sh);
        d3 = x1v_(m,i)*(x3f_(m,k+1) - x3f_(m,k))*se;
      }
      h += 4.0/SQR(d3);
    }
    if (cs) {h /= SQR(scell(m,k,j));}
    lam = er*lam + eh*h;
    if (lam > v.val) {v.val = lam; v.loc = idx;}
  }, ML(res));

  Real lmax = res.val;
  Real r = 0.0, e = 0.0;
  if (rloc != nullptr || etaloc != nullptr) {
    if (nmkji > 0) {
      const int idx = res.loc;
      int m = idx/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(x1v_, m, i));
      auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(eta_, m, k, j, i));
      r = hr();
      e = he();
    }
  }
#if MPI_PARALLEL_ENABLED
  struct {double v; int r;} in, out;
  in.v = static_cast<double>(lmax);
  in.r = global_variable::my_rank;
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  lmax = static_cast<Real>(out.v);
  if (rloc != nullptr || etaloc != nullptr) {
    double re[2] = {static_cast<double>(r), static_cast<double>(e)};
    MPI_Bcast(re, 2, MPI_DOUBLE, out.r, MPI_COMM_WORLD);
    r = re[0];
    e = re[1];
  }
#endif
  if (rloc != nullptr) {*rloc = r;}
  if (etaloc != nullptr) {*etaloc = e;}
  return lmax;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::NewTimeStepRKGFE
//! \brief dtnew = rkg_dt_safety x 2/lambda_bound (local max; mesh.cpp takes it as dt_diff
//! without the CFL factor, and the MPI min over ranks as before).

void Resistivity::NewTimeStepRKGFE() {
  const Real lam = LambdaBoundMax();
  dtnew = (lam > 0.0) ? rkg_dt_safety*2.0/lam
                      : static_cast<Real>(std::numeric_limits<float>::max());
}

namespace {
KOKKOS_INLINE_FUNCTION
Real RkgdtHash(const int m, const int k, const int j, const int i, const int c) {
  uint32_t h = 2166136261u;
  const int q[5] = {m, k, j, i, c};
  for (int n = 0; n < 5; ++n) {
    h ^= static_cast<uint32_t>(q[n]) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h *= 16777619u;
  }
  h ^= h >> 15;
  h *= 0x2c1b3c6du;
  h ^= h >> 12;
  return static_cast<Real>(h & 0xffffffu)/static_cast<Real>(0x800000u) - 1.0;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::LmaxDiag
//! \brief power iteration v <- L v/|L v| with L v = curl_CT(eta J(v)): the face field v is
//! held in pmhd->b0 so that the ordinary face-field halo (cs seam resample included), the
//! physical and user B boundary conditions, the resistive EMF (AddEMFGeneralResist ->
//! gnomonic two-pass curl on cs), the conducting/stress-free wall zeroing and the seam/edge
//! EMF synchronisation (PackAndSendFluxFC) act on it exactly as in a super-stage.  The curl
//! is resistivity_ct.cpp's curvilinear branch with dt = 1 and no recurrence.  b0, bcc0,
//! u0, w0 and efld are saved and restored.

void Resistivity::LmaxDiag() {
  auto *pmhd = pmy_pack->pmhd;
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pm->multi_d;
  const bool three_d = pm->three_d;
  const int niter = rkg_lmax_diag;

  auto &b0 = pmhd->b0;
  auto &ef = pmhd->efld;
  // saved state
  auto sb1 = Kokkos::create_mirror(DevMemSpace(), b0.x1f);
  auto sb2 = Kokkos::create_mirror(DevMemSpace(), b0.x2f);
  auto sb3 = Kokkos::create_mirror(DevMemSpace(), b0.x3f);
  auto se1 = Kokkos::create_mirror(DevMemSpace(), ef.x1e);
  auto se2 = Kokkos::create_mirror(DevMemSpace(), ef.x2e);
  auto se3 = Kokkos::create_mirror(DevMemSpace(), ef.x3e);
  auto su0 = Kokkos::create_mirror(DevMemSpace(), pmhd->u0);
  auto sw0 = Kokkos::create_mirror(DevMemSpace(), pmhd->w0);
  auto sbc = Kokkos::create_mirror(DevMemSpace(), pmhd->bcc0);
  Kokkos::deep_copy(sb1, b0.x1f);
  Kokkos::deep_copy(sb2, b0.x2f);
  Kokkos::deep_copy(sb3, b0.x3f);
  Kokkos::deep_copy(se1, ef.x1e);
  Kokkos::deep_copy(se2, ef.x2e);
  Kokkos::deep_copy(se3, ef.x3e);
  Kokkos::deep_copy(su0, pmhd->u0);
  Kokkos::deep_copy(sw0, pmhd->w0);
  Kokkos::deep_copy(sbc, pmhd->bcc0);
  // L v
  auto a1 = Kokkos::create_mirror(DevMemSpace(), b0.x1f);
  auto a2 = Kokkos::create_mirror(DevMemSpace(), b0.x2f);
  auto a3 = Kokkos::create_mirror(DevMemSpace(), b0.x3f);

  // the old dt_diff formula, for the printout
  const bool fe_save = rkg_dt_fe;
  const Real dtn_save = dtnew;
  rkg_dt_fe = false;
  NewTimeStepGeneralResist(pmhd->w0, pmhd->peos->eos_data);
  Real dt_old = pm->cfl_no*dtnew;
  rkg_dt_fe = fe_save;
  dtnew = dtn_save;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &dt_old, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif
  Real rb = 0.0, eb = 0.0;
  const Real lbound = LambdaBoundMax(&rb, &eb);

  // random start on the active faces
  {
    auto v1 = b0.x1f;
    auto v2 = b0.x2f;
    auto v3 = b0.x3f;
    Kokkos::deep_copy(v1, 0.0);
    Kokkos::deep_copy(v2, 0.0);
    Kokkos::deep_copy(v3, 0.0);
    par_for("rkgdt_v1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v1(m,k,j,i) = RkgdtHash(m,k,j,i,1);});
    par_for("rkgdt_v2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v2(m,k,j,i) = RkgdtHash(m,k,j,i,2);});
    par_for("rkgdt_v3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v3(m,k,j,i) = RkgdtHash(m,k,j,i,3);});
  }

  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  auto &area1 = pmy_pack->pcoord->area.x1f;
  auto &area2 = pmy_pack->pcoord->area.x2f;
  auto &area3 = pmy_pack->pcoord->area.x3f;
  auto &dxe1 = pmy_pack->pcoord->dxedge.x1e;
  auto &dxe2 = pmy_pack->pcoord->dxedge.x2e;
  auto &dxe3 = pmy_pack->pcoord->dxedge.x3e;
  auto *pbb = pmhd->pbval_b;

  Real lam_rq = 0.0, lam_n = 0.0, lam_rq_prev = 0.0;
  const int nprint = std::max(1, niter/20);
  for (int it = 0; it < niter; ++it) {
    // ghosts of v: halo (with the cs seam resample) + physical/user BCs
    while (pbb->InitRecv(3) != TaskStatus::complete) {}
    while (pbb->InitFluxRecv(3) != TaskStatus::complete) {}
    while (pbb->PackAndSendFC(b0, pmhd->coarse_b0) != TaskStatus::complete) {}
    while (pbb->RecvAndUnpackFC(b0, pmhd->coarse_b0) != TaskStatus::complete) {}
    if (!pm->strictly_periodic) {
      MeshBoundaryValuesFC::BFieldBCs(pmy_pack, pbb->b_in, b0);
      if (pm->pgen->user_bcs) {(pm->pgen->user_bcs_func)(pm);}
    }
    // E = eta J(v), walls, then the edge sync, as in a super-stage
    ClearResistiveEMFs(efld_resist);
    AddEMFGeneralResist(b0, efld_resist);
    // ambipolar: J(v) with eta_A and the B direction of the saved state (sbc), i.e. the
    // frozen-coefficient linearisation of the AD EMF
    if (ad_type != 0) {AddEMFAmbipolar(b0, sbc, efld_resist);}
    if (wall_cond_reflect || wall_cond_ix1) {
      pmhd->WallEmfZero(efld_resist, wall_cond_reflect, wall_cond_ix1);
    }
    AddEMFDirect(efld_resist, ef);
    while (pbb->PackAndSendFluxFC(ef) != TaskStatus::complete) {}
    while (pbb->RecvAndUnpackFluxFC(ef) != TaskStatus::complete) {}
    while (pbb->ClearSend() != TaskStatus::complete) {}
    while (pbb->ClearRecv() != TaskStatus::complete) {}
    while (pbb->ClearFluxSend() != TaskStatus::complete) {}
    while (pbb->ClearFluxRecv() != TaskStatus::complete) {}

    // a = curl_CT E (b_new = b - dt a in resistivity_ct.cpp)
    auto e1 = ef.x1e;
    auto e2 = ef.x2e;
    auto e3 = ef.x3e;
    Kokkos::deep_copy(a1, 0.0);
    Kokkos::deep_copy(a2, 0.0);
    Kokkos::deep_copy(a3, 0.0);
    if (multi_d) {
      par_for("rkgdt_a1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        Real a = (dxe3(m,k,j+1,i)*e3(m,k,j+1,i) - dxe3(m,k,j,i)*e3(m,k,j,i))/area1(m,k,j,i);
        if (three_d) {
          a -= (dxe2(m,k+1,j,i)*e2(m,k+1,j,i) - dxe2(m,k,j,i)*e2(m,k,j,i))/area1(m,k,j,i);
        }
        a1(m,k,j,i) = a;
      });
    }
    par_for("rkgdt_a2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      const bool do_pole =
          (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) ||
          (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1);
      Real ar = (do_pole) ? 1.0 : area2(m,k,j,i);
      Real dxe31 = (do_pole) ? 0.0 : dxe3(m,k,j,i+1);
      Real dxe30 = (do_pole) ? 0.0 : dxe3(m,k,j,i);
      Real a = -(dxe31*e3(m,k,j,i+1) - dxe30*e3(m,k,j,i))/ar;
      if (three_d) {
        Real dxe11 = (do_pole) ? 0.0 : dxe1(m,k+1,j,i);
        Real dxe10 = (do_pole) ? 0.0 : dxe1(m,k,j,i);
        a += (dxe11*e1(m,k+1,j,i) - dxe10*e1(m,k,j,i))/ar;
      }
      a2(m,k,j,i) = a;
    });
    par_for("rkgdt_a3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      Real a = (dxe2(m,k,j,i+1)*e2(m,k,j,i+1) - dxe2(m,k,j,i)*e2(m,k,j,i))/area3(m,k,j,i);
      if (multi_d) {
        a -= (dxe1(m,k,j+1,i)*e1(m,k,j+1,i) - dxe1(m,k,j,i)*e1(m,k,j,i))/area3(m,k,j,i);
      }
      a3(m,k,j,i) = a;
    });

    // dots over the active faces (block-boundary faces counted on both blocks)
    auto v1 = b0.x1f;
    auto v2 = b0.x2f;
    auto v3 = b0.x3f;
    Real s[3] = {0.0, 0.0, 0.0};
    {
      Real vv = 0.0, va = 0.0, aa = 0.0;
      const int n1 = ie - is + 2, n2 = je - js + 2, n3 = ke - ks + 2;
      const int ntot = (nmb1 + 1)*n3*n2*n1;
      Kokkos::parallel_reduce("rkgdt_dot", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
      KOKKOS_LAMBDA(const int idx, Real &svv, Real &sva, Real &saa) {
        int m = idx/(n3*n2*n1);
        int r = idx - m*(n3*n2*n1);
        int k = r/(n2*n1) + ks;
        r -= (k - ks)*(n2*n1);
        int j = r/n1 + js;
        int i = r - (j - js)*n1 + is;
        if (k <= ke && j <= je) {
          svv += SQR(v1(m,k,j,i)); sva += v1(m,k,j,i)*a1(m,k,j,i); saa += SQR(a1(m,k,j,i));
        }
        if (k <= ke && i <= ie) {
          svv += SQR(v2(m,k,j,i)); sva += v2(m,k,j,i)*a2(m,k,j,i); saa += SQR(a2(m,k,j,i));
        }
        if (j <= je && i <= ie) {
          svv += SQR(v3(m,k,j,i)); sva += v3(m,k,j,i)*a3(m,k,j,i); saa += SQR(a3(m,k,j,i));
        }
      }, vv, va, aa);
      s[0] = vv; s[1] = va; s[2] = aa;
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, s, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    lam_rq_prev = lam_rq;
    lam_rq = s[1]/s[0];
    lam_n = std::sqrt(s[2]/s[0]);
    if (global_variable::my_rank == 0 && ((it+1) % nprint == 0 || it == 0)) {
      std::cout << "rkg_lmax_diag it=" << it+1 << " lambda_rq=" << lam_rq
                << " lambda_norm=" << lam_n << " lambda_bound=" << lbound
                << " bound/rq=" << lbound/lam_rq << std::endl;
    }
    // v <- a/|a|
    const Real inv = 1.0/std::sqrt(s[2]);
    Kokkos::deep_copy(v1, 0.0);
    Kokkos::deep_copy(v2, 0.0);
    Kokkos::deep_copy(v3, 0.0);
    par_for("rkgdt_n1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v1(m,k,j,i) = inv*a1(m,k,j,i);});
    par_for("rkgdt_n2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v2(m,k,j,i) = inv*a2(m,k,j,i);});
    par_for("rkgdt_n3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {v3(m,k,j,i) = inv*a3(m,k,j,i);});
  }

  // where the converged vector lives: the tangential x2-face entry of largest |v|
  Real rv = 0.0, ev = 0.0;
  {
    auto v2 = b0.x2f;
    auto eta_ = eta_b;
    auto x1v_ = pmy_pack->pcoord->x1v;
    const int n1 = ie - is + 1, n2 = je - js + 2, n3 = ke - ks + 1;
    const int ntot = (nmb1 + 1)*n3*n2*n1;
    using ML = Kokkos::MaxLoc<Real, int>;
    ML::value_type res;
    Kokkos::parallel_reduce("rkgdt_vloc", Kokkos::RangePolicy<>(DevExeSpace(), 0, ntot),
    KOKKOS_LAMBDA(const int idx, ML::value_type &v) {
      int m = idx/(n3*n2*n1);
      int r = idx - m*(n3*n2*n1);
      int k = r/(n2*n1) + ks;
      r -= (k - ks)*(n2*n1);
      int j = r/n1 + js;
      int i = r - (j - js)*n1 + is;
      const Real a = fabs(v2(m,k,j,i));
      if (a > v.val) {v.val = a; v.loc = idx;}
    }, ML(res));
    if (ntot > 0) {
      int idx = res.loc;
      int m = idx/(n3*n2*n1);
      int r = idx - m*(n3*n2*n1);
      int k = r/(n2*n1) + ks;
      r -= (k - ks)*(n2*n1);
      int j = r/n1 + js;
      int i = r - (j - js)*n1 + is;
      j = std::min(j, je);
      auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(x1v_, m, i));
      auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(eta_, m, k, j, i));
      rv = hr();
      ev = he();
    }
#if MPI_PARALLEL_ENABLED
    struct {double v; int r;} in, out;
    in.v = static_cast<double>(res.val);
    in.r = global_variable::my_rank;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    double re[2] = {static_cast<double>(rv), static_cast<double>(ev)};
    MPI_Bcast(re, 2, MPI_DOUBLE, out.r, MPI_COMM_WORLD);
    rv = re[0];
    ev = re[1];
#endif
  }

  if (global_variable::my_rank == 0) {
    const Real dt = pm->dt;
    auto stages = [](Real tau) {
      Real sr = -0.5 + std::sqrt(2.25 + 4.0*tau);
      int s = static_cast<int>(std::ceil(sr));
      return (s < 3) ? 3 : s;
    };
    std::cout << "rkg_lmax_diag RESULT iters=" << niter
              << " lambda_rq=" << lam_rq << " (prev " << lam_rq_prev << ")"
              << " lambda_norm=" << lam_n << std::endl
              << "rkg_lmax_diag RESULT lambda_bound=" << lbound
              << " at r=" << rb << " eta=" << eb
              << "  bound/lambda_rq=" << lbound/lam_rq
              << "  vector max |v2| at r=" << rv << " eta=" << ev << std::endl
              << "rkg_lmax_diag RESULT dt=" << dt
              << " dt_diff_old=" << dt_old << " (s=" << stages(dt/dt_old) << ")"
              << "  2/lambda_rq=" << 2.0/lam_rq << " (s=" << stages(dt*lam_rq/2.0) << ")"
              << "  2/lambda_bound=" << 2.0/lbound
              << " (s at safety 1/0.8/0.5: " << stages(dt*lbound/2.0) << "/"
              << stages(dt*lbound/1.6) << "/" << stages(dt*lbound) << ")"
              << "  this run dt_diff=" << pm->dt_diff << std::endl;
  }

  // restore
  Kokkos::deep_copy(b0.x1f, sb1);
  Kokkos::deep_copy(b0.x2f, sb2);
  Kokkos::deep_copy(b0.x3f, sb3);
  Kokkos::deep_copy(ef.x1e, se1);
  Kokkos::deep_copy(ef.x2e, se2);
  Kokkos::deep_copy(ef.x3e, se3);
  Kokkos::deep_copy(pmhd->u0, su0);
  Kokkos::deep_copy(pmhd->w0, sw0);
  Kokkos::deep_copy(pmhd->bcc0, sbc);
}
