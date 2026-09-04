//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sp_test.cpp
//! \brief SPHERICAL-POLAR counterpart of cs_test iprob=8: a uniform gas AT REST threaded
//! by a UNIFORM Cartesian field B = b0 zhat.  curl B = 0, so the exact evolution is
//! nothing; after one step the momentum per cell is the spurious force of the geometric
//! source term against the flux divergence (the cancellation the cubed sphere's
//! well-balanced source was built for).  Measure it from two consecutive dumps.
//!
//! The face-normal components are FACE AVERAGES, not point values: on an r-face the
//! area-weighted mean of b0 cos(theta) is b0 (cos th_l + cos th_r)/2, and on a theta-face
//! -b0 sin(theta_f) is constant.  With those the discrete divergence is zero EXACTLY
//! (both terms reduce to b0 (r_r^2 - r_l^2) dphi (sin^2 th_r - sin^2 th_l)/2), so no
//! monopole force contaminates the measurement.

#include <cmath>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen.hpp"
#include "hydro/hydro.hpp"

// FLUX PROBE (problem/flux_probe = 1, needs problem/user_srcs = true): on the first stage
// of the first cycle, print the code's momentum flux at one polar cell's outer theta-face
// and two radial faces next to the EXACT stress projection T.n of the analytic field, with
// each face's area/volume weight.  Tells which face carries the polar-row residual.
namespace {
Real sp_b0 = 1.0, sp_p0 = 1.0;
int sp_bdir = 2, sp_probe_done = 0;
}  // namespace
void SPTestFluxProbe(Mesh *pm, const Real bdt);

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (!pmy_mesh_->use_spherical_polar || pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "sp_test requires mesh/use_spherical_polar = true and an <mhd> block"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  const Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  const Real p0 = pin->GetOrAddReal("problem", "p0", 1.0);
  const Real b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  // field direction: bdir = 2 (default) is b0 zhat, bdir = 0 is b0 xhat.  For xhat the
  // face averages are B_r = b0 <sin th> cos-average, B_th = b0 cos(th_f) <cos ph>,
  // B_ph = -b0 sin(ph_f); each is the exact area-weighted mean over its face.
  const int bdir = pin->GetOrAddInteger("problem", "bdir", 2);
  sp_b0 = b0; sp_p0 = p0; sp_bdir = bdir;
  if (pin->GetOrAddInteger("problem", "flux_probe", 0) != 0) user_srcs_func = SPTestFluxProbe;
  // econsistent = true: the magnetic energy is built from the SAME cell-centred bcc that
  // ConsToPrim subtracts, so the recovered pressure is exactly uniform and only the
  // geometric-source/flux-divergence mismatch survives.  false (cs_test's choice): the
  // analytic b0^2/2, so the |bcc|^2 - b0^2 defect appears as an O(dtheta^2) pressure
  // ripple, largest in the polar rows -- which then dominates the one-step force.
  const bool econs = pin->GetOrAddBoolean("problem", "econsistent", true);
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0f = pmbp->pmhd->b0;
  auto &bcc = pmbp->pmhd->bcc0;
  auto &x2f_ = pmbp->pcoord->xx2f;   // theta at the theta-faces, (m, j)
  auto &x3f_ = pmbp->pcoord->xx3f;   // phi at the phi-faces, (m, k)
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x3v_ = pmbp->pcoord->x3v;

  par_for("sp_test_b", DevExeSpace(), 0,nmb1, ks,ke+1, js,je+1, is,ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    if (bdir == 2) {
      if (j <= je && k <= ke) {
        b0f.x1f(m,k,j,i) = b0*0.5*(std::cos(x2f_(m,j)) + std::cos(x2f_(m,j+1)));
      }
      if (i <= ie && k <= ke) {
        b0f.x2f(m,k,j,i) = -b0*std::sin(x2f_(m,j));
      }
      if (i <= ie && j <= je) {
        b0f.x3f(m,k,j,i) = 0.0;
      }
    } else {
      // B = b0 xhat: B_r = b0 sin(th) cos(ph), B_th = b0 cos(th) cos(ph),
      // B_ph = -b0 sin(ph)
      const Real tl = x2f_(m,j), tr = x2f_(m,(j <= je) ? j+1 : j);
      const Real pl = x3f_(m,k), pr = x3f_(m,(k <= ke) ? k+1 : k);
      if (j <= je && k <= ke) {
        // <sin th>_area = Int sin^2 / Int sin ; <cos ph> = (sin pr - sin pl)/(pr - pl)
        const Real isin2 = 0.5*(tr - tl) - 0.25*(std::sin(2.0*tr) - std::sin(2.0*tl));
        const Real isin  = std::cos(tl) - std::cos(tr);
        b0f.x1f(m,k,j,i) = b0*(isin2/isin)*(std::sin(pr) - std::sin(pl))/(pr - pl);
      }
      if (i <= ie && k <= ke) {
        b0f.x2f(m,k,j,i) = b0*std::cos(tl)*(std::sin(pr) - std::sin(pl))/(pr - pl);
      }
      if (i <= ie && j <= je) {
        b0f.x3f(m,k,j,i) = -b0*std::sin(pl);
      }
    }
  });

  par_for("sp_test_u", DevExeSpace(), 0,nmb1, ks,ke, js,je, is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    u0(m,IDN,k,j,i) = d0;
    u0(m,IM1,k,j,i) = 0.0;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    // the SAME position-weighted interpolation to the volume centroid that ConsToPrim
    // uses on this grid (ideal_mhd.cpp), so econsistent really is consistent
    Real lw, rw;
    lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
    rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
    bcc(m,IBX,k,j,i) = lw*b0f.x1f(m,k,j,i) + rw*b0f.x1f(m,k,j,i+1);
    lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
    rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
    bcc(m,IBY,k,j,i) = lw*b0f.x2f(m,k,j,i) + rw*b0f.x2f(m,k,j+1,i);
    lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
    rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
    bcc(m,IBZ,k,j,i) = lw*b0f.x3f(m,k,j,i) + rw*b0f.x3f(m,k+1,j,i);
    const Real bsq = econs ? (SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i))
                             + SQR(bcc(m,IBZ,k,j,i))) : b0*b0;
    u0(m,IEN,k,j,i) = p0/gm1 + 0.5*bsq;
  });
  return;
}

void SPTestFluxProbe(Mesh *pm, const Real bdt) {
  if (sp_probe_done++ != 0) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  auto x1f = Kokkos::create_mirror_view(pmbp->pcoord->xx1f); Kokkos::deep_copy(x1f, pmbp->pcoord->xx1f);
  auto x2f = Kokkos::create_mirror_view(pmbp->pcoord->xx2f); Kokkos::deep_copy(x2f, pmbp->pcoord->xx2f);
  auto x2v = Kokkos::create_mirror_view(pmbp->pcoord->x2v);  Kokkos::deep_copy(x2v, pmbp->pcoord->x2v);
  auto a1 = Kokkos::create_mirror_view(pmbp->pcoord->area.x1f); Kokkos::deep_copy(a1, pmbp->pcoord->area.x1f);
  auto a2 = Kokkos::create_mirror_view(pmbp->pcoord->area.x2f); Kokkos::deep_copy(a2, pmbp->pcoord->area.x2f);
  auto vol = Kokkos::create_mirror_view(pmbp->pcoord->volume); Kokkos::deep_copy(vol, pmbp->pcoord->volume);
  auto f1 = Kokkos::create_mirror_view(pmbp->pmhd->uflx.x1f); Kokkos::deep_copy(f1, pmbp->pmhd->uflx.x1f);
  auto f2 = Kokkos::create_mirror_view(pmbp->pmhd->uflx.x2f); Kokkos::deep_copy(f2, pmbp->pmhd->uflx.x2f);
  auto bcc = Kokkos::create_mirror_view(pmbp->pmhd->bcc0); Kokkos::deep_copy(bcc, pmbp->pmhd->bcc0);
  auto b2f = Kokkos::create_mirror_view(pmbp->pmhd->b0.x2f); Kokkos::deep_copy(b2f, pmbp->pmhd->b0.x2f);
  auto w = Kokkos::create_mirror_view(pmbp->pmhd->w0); Kokkos::deep_copy(w, pmbp->pmhd->w0);
  const int m = 0, k = ks, i = is + 3;
  const Real b0 = sp_b0, p0 = sp_p0;
  std::printf("### SP FLUX PROBE (bdir=%d, b0=%g, p0=%g) block %d k=%d i=%d, cell j=js\n",
              sp_bdir, b0, p0, m, k, i);
  std::printf("###   V=%.6e  theta faces: %.6e %.6e  x2v=%.6e\n", vol(m,k,js,i),
              x2f(m,js), x2f(m,js+1), x2v(m,js));
  std::printf("###   cell bcc: Br=%.8e Bth=%.8e Bph=%.8e  p(w)=%.8e\n", bcc(m,IBX,k,js,i),
              bcc(m,IBY,k,js,i), bcc(m,IBZ,k,js,i), w(m,IEN,k,js,i)*(5.0/3.0-1.0));
  // exact stress projections for B = b0 zhat, v = 0:
  //   theta-face (n = thhat):  F_r = -B_th B_r = b0^2 sin cos,  F_th = p + b0^2 cos(2th)/2, F_ph = 0
  //   r-face     (n = rhat):   F_r = p + b0^2 (sin^2 - cos^2)/2 = p - b0^2 cos(2th)/2,  F_th = -B_r B_th = b0^2 sin cos
  for (int jf = js; jf <= js + 2; ++jf) {
    const Real th = x2f(m,jf);
    const Real exr = b0*b0*std::sin(th)*std::cos(th), exth = p0 + 0.5*b0*b0*std::cos(2.0*th);
    std::printf("###   THETA face j=%d (th=%.6f) A/V=%.6e  B_n(face)=%.8e\n"
                "###       code F_r=%.8e F_th=%.8e F_ph=%.3e | exact F_r=%.8e F_th=%.8e\n"
                "###       (code-exact)*A/V: r %.3e  th %.3e   in units (b0^2/2)/r: r %.3e th %.3e\n",
                jf, th, a2(m,k,jf,i)/vol(m,k,js,i), b2f(m,k,jf,i),
                f2(m,IM1,k,jf,i), f2(m,IM2,k,jf,i), f2(m,IM3,k,jf,i), exr, exth,
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,js,i),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,js,i),
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)));
  }
  // r-faces of the polar cell: the exact flux is the AREA MEAN over the face's theta range
  const Real tl = x2f(m,js), tr = x2f(m,js+1);
  const Real isin = std::cos(tl) - std::cos(tr);
  const Real mcos2 = (0.5*(std::cos(tl)*std::cos(tl)*std::cos(tl) - std::cos(tr)*std::cos(tr)*std::cos(tr))*2.0/3.0
                      - 0.5*isin*0.0);  // placeholder replaced below
  // <cos 2th> over the r-face with weight sin th: Int (2cos^2-1) sin dth / Int sin dth
  const Real icos2sin = (std::cos(tl)*std::cos(tl)*std::cos(tl) - std::cos(tr)*std::cos(tr)*std::cos(tr))/3.0;
  const Real avg_cos2 = 2.0*icos2sin/isin - 1.0;
  // <sin cos> = Int sin^2 cos / Int sin = (sin^3 tr - sin^3 tl)/3 / isin
  const Real avg_sc = (std::pow(std::sin(tr),3) - std::pow(std::sin(tl),3))/3.0/isin;
  (void)mcos2;
  for (int ii = i; ii <= i+1; ++ii) {
    const Real exr = p0 - 0.5*b0*b0*avg_cos2, exth = b0*b0*avg_sc;
    std::printf("###   R face i=%d (r=%.4f) A/V=%.6e\n"
                "###       code F_r=%.8e F_th=%.8e | exact(area-mean) F_r=%.8e F_th=%.8e\n"
                "###       (code-exact)*A/V in (b0^2/2)/r: r %.3e th %.3e\n",
                ii, x1f(m,ii), a1(m,k,js,ii)/vol(m,k,js,i),
                f1(m,IM1,k,js,ii), f1(m,IM2,k,js,ii), exr, exth,
                (f1(m,IM1,k,js,ii)-exr)*a1(m,k,js,ii)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)),
                (f1(m,IM2,k,js,ii)-exth)*a1(m,k,js,ii)/vol(m,k,js,i)/(0.5*b0*b0/x1f(m,i)));
  }
  // and the interior cell j = js+4 for reference (its outer theta-face)
  {
    const int jc = js + 4, jf = jc + 1; const Real th = x2f(m,jf);
    const Real exr = b0*b0*std::sin(th)*std::cos(th), exth = p0 + 0.5*b0*b0*std::cos(2.0*th);
    std::printf("###   INTERIOR ref: theta face j=%d (th=%.6f) A/V=%.6e  (code-exact)*A/V in (b0^2/2)/r: r %.3e th %.3e\n",
                jf, th, a2(m,k,jf,i)/vol(m,k,jc,i),
                (f2(m,IM1,k,jf,i)-exr)*a2(m,k,jf,i)/vol(m,k,jc,i)/(0.5*b0*b0/x1f(m,i)),
                (f2(m,IM2,k,jf,i)-exth)*a2(m,k,jf,i)/vol(m,k,jc,i)/(0.5*b0*b0/x1f(m,i)));
  }
}
