//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_closure_rays.cpp
//! \brief the closure of a cell (D, h+-, k+, C1) from the moments of the ray sweeps.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_ring.cpp, rays/sc_source.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_init.cpp

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <type_traits>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_moments.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RayCut
//! \brief clean_ray_cut: first swept shell vray_icut (global min over columns); below
//! it the closure is the diffusion limit (ThickI, error O(R^2),
//! R = |grad S|/(chi S)).
//! A column's cut: 10 optical depths below its lowest cell with R >= ray_cut_r
//! or thinner than one optical depth.  Moves down at once, up after ray_cut_wait
//! builds (both in the restart file).  sc rays only; no cut with several MeshBlocks
//! along x1

void RadiationVet::RayCut() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  vray_icut = 0;
  if (!cln_rcut_on || !cln_init || pmy_pack->pmesh->mesh_indcs.nx1 != indcs.nx1) return;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, n1 = indcs.nx1;
  const int nk = ke - ks + 1, nj = je - js + 1, nf = pmy_pack->pmesh->three_d ? 3 :
                                                  (pmy_pack->pmesh->multi_d ? 2 : 1);
  SolveArgs x = MakeArgs();
  auto cs_ = vray_cs;
  const Real rr = cln_rcut_r;
  int lmin = n1;
  Kokkos::parallel_reduce("cln_raycut", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          pmy_pack->nmb_thispack*nk*nj),
  KOKKOS_LAMBDA(const int idx, int &lm) {
    const int m = idx/(nk*nj), k = ks + (idx/nj) % nk, j = js + idx % nj;
    auto tc = [&](const int i) {return exp(cs_(m,0,k,j,i))*(x.x1f(m,i+1) - x.x1f(m,i));};
    int i = is;
    for (; i <= ie && !(tc(i) < 1.0); ++i) {
      Real r2 = 0.0;
      for (int dir = 0; dir < nf; ++dir) {
        const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
        Real g = 0.0;
        for (int up = 0; up < 2; ++up) {
          const int s = up ? 1 : -1;
          if (IsBoundaryFace(x, dir, m, k + up*dk, j + up*dj, i + up*di)) continue;
          g = fmax(g, fabs(cs_(m,1,k+s*dk,j+s*dj,i+s*di) - cs_(m,1,k,j,i))/
                      x.dxf[dir](m,k+up*dk,j+up*dj,i+up*di));
        }
        r2 += g*g;
      }
      if (!(sqrt(r2) < rr*exp(cs_(m,0,k,j,i)))) break;
    }
    int lc = i - is;
    for (Real t = 0.0; lc > 0 && t < 10.0;) t += tc(is + --lc);
    lm = (lc < lm) ? lc : lm;
  }, Kokkos::Min<int>(lmin));
  lmin = -static_cast<int>(GlobalReduce(-lmin, true));
  lmin = std::max(0, std::min(lmin, n1 - 2));
  const int prev = cln_rcut;
  cln_rcut_n = (lmin > cln_rcut && cln_rcut >= 0) ? cln_rcut_n + 1 : 0;
  if (lmin < cln_rcut || cln_rcut < 0 || cln_rcut_n >= cln_rcut_wait) {
    cln_rcut = lmin;
    cln_rcut_n = 0;
  }
  if (global_variable::my_rank == 0 && cln_verbose > 0 && prev != cln_rcut) {
    std::cout << "CLEAN ray cut cycle " << pmy_pack->pmesh->ncycle << " shell "
              << prev << " -> " << cln_rcut << " (want " << lmin << ")" << std::endl;
  }
  pin_rcut_->SetInteger("rad_vet", "ray_cut_now", cln_rcut);
  pin_rcut_->SetInteger("rad_vet", "ray_cut_count", cln_rcut_n);
  vray_icut = cln_rcut;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureFromSc
//! \brief closure from the vet_gd intensities, mu_b = n . e_b in the cell's basis:
//! J = <I>, K_bc = <I mu_b mu_c>, H+-_b = <I mu_b> over mu_b >< 0, K+^d_b = <I mu_b mu_d>
//! over mu_d > 0, C1 = <|I - S|>/J.  Cubed sphere: Cartesian D, half ranges about the
//! family normals

void RadiationVet::ClosureFromSc() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n, og = vgd_w - indcs.ng, ok = vgd_w - ks;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cw_ = vgd_cs;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;
  const bool hf = cln_rfh;   // (CLN_HO)
  const bool kp = cs.extent_int(1) > CLN_KP + 8;   // (k+, not CLN_HO)
  auto gn = cln_gn;
  // Cartesian slab: the raw moments J, K, H also into vet_cell (the vet_sc outputs)
  const bool slb = (vgd_slab == 1);
  auto vc_ = vet_cell;
  SolveArgs x = MakeArgs();
  ScSlabArgs(x);
  // (in device memory, RayFaceSums)
  cln_gdx = DevArgs(x, 3);
  const SolveArgs *px = cln_gdx;
  // clean_ray_cut: below the cut vgd_i holds the diffusion limit
  if (vgd_on && vray_icut > 0) {
    const Real ic = 1.0/c_light;
    // (one team per cell, lanes along the directions: contiguous stores)
    const int nc = vray_icut, nj = je - js + 1, nk = ke - ks + 1;
    const int tsz =
        std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace> ? 1 : 32;
    Kokkos::parallel_for("cln_gd_thk",
                         Kokkos::TeamPolicy<>(DevExeSpace(), (nmb1 + 1)*nk*nj*nc, tsz),
    KOKKOS_LAMBDA(TeamMember_t tm) {
      const int t = tm.league_rank(), i = is + t % nc, j = js + (t/nc) % nj;
      const int k = ks + (t/(nc*nj)) % nk, m = t/(nc*nj*nk);
      const SolveArgs &x = *px;
      Real e[3][3];   // (the cell's frame, F and E once, not per direction)
      ThickFrame(x, x3v, m, k, j, i, e);
      const Real f[3] = {iw_(m,VET_IW_F1,k,j,i), iw_(m,VET_IW_F2,k,j,i),
                         iw_(m,VET_IW_F3,k,j,i)}, en = iw_(m,VET_IW_EN,k,j,i);
      const DirRow r = ScDirRow(vi_, m, k+ok, j+og, i);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, n), [&](const int d) {
        const Real nv[3] = {dir_(d,0), dir_(d,1), dir_(d,2)};
        r.p[d*r.s] = ThickI(e, f, en, ic, nv);
      });
    });
  }
  SolveForK("cln_gd", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const SolveArgs &x = *px;
    // (spherical polar only: other meshes hold no angle arrays)
    const Real st = x.sph ? sin(x2v(m,j)) : 0.0, ct = x.sph ? cos(x2v(m,j)) : 0.0;
    const Real sp = x.sph ? sin(x3v(m,k)) : 0.0, cp = x.sph ? cos(x3v(m,k)) : 0.0;
    const Real s0 = exp(cw_(m,1,k+ok,j+og,i));
    // cubed sphere: Cartesian D, half ranges about the families n_d
    Real nd[3][3];
    if (x.cub) CsFamily(x, m, k, j, i, nd);
    Real q[CLN_N] = {0.0}, hp[3] = {0.0}, hm[3] = {0.0};   // (q in registers)
    Real fh[3][3] = {{0.0}}, fl[3][3] = {{0.0}}, oh[3] = {0.0}, ol[3] = {0.0};
    for (int a = 0; a < 9 && hf; ++a) {   // (the normals of the high and low faces)
      fh[a/3][a%3] = gn.At<CLG_N>(m,a,k+(a/3 == 2),j+(a/3 == 1),i+(a/3 == 0));
      fl[a/3][a%3] = gn.At<CLG_N>(m,a,k,j,i);
    }
    const DirRow r = ScDirRow(vi_, m, k+ok, j+og, i);
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      Real mu[3] = {nx*st*cp + ny*st*sp + nz*ct, nx*ct*cp + ny*ct*sp - nz*st,
                    -nx*sp + ny*cp}, nu[3];
      for (int b = 0; b < 3; ++b) {
        if (!x.sph) mu[b] = dir_(d,b);
        nu[b] = x.cub ? nx*nd[b][0] + ny*nd[b][1] + nz*nd[b][2] : mu[b];
      }
      const Real nv[3] = {nx, ny, nz};
      const Real iv = r[d];
      const Real a = dir_(d,3)*iv;
      q[CLN_J] += a;
      for (int b = 0; b < 3; ++b) {
        q[CLN_D11+b] += a*mu[b]*mu[b];
        if (!(nu[b] > 0.0)) {
          hm[b] += a*nu[b];
          continue;
        }
        hp[b] += a*nu[b];
        for (int c = 0; c < 3 && kp; ++c) q[CLN_KP+3*b+c] += a*mu[c]*nu[b];
      }
      for (int b = 0; b < 3 && hf; ++b) {
        Real mh = 0.0, ml = 0.0;
        for (int c = 0; c < 3; ++c) {
          mh += dir_(d,c)*fh[b][c];
          ml += dir_(d,c)*fl[b][c];
        }
        if (mh > 0.0) oh[b] += a*mh;
        if (!(ml > 0.0)) ol[b] += a*ml;
      }
      q[CLN_D11+3] += a*mu[0]*mu[1];
      q[CLN_D11+4] += a*mu[0]*mu[2];
      q[CLN_D11+5] += a*mu[1]*mu[2];
      q[CLN_C1] += dir_(d,3)*fabs(iv - s0);
    }
    for (int b = 0; b < 3; ++b) {
      q[CLN_HP+b] = hp[b];
      q[CLN_HM+b] = hm[b];
    }
    q[CLN_JR] = q[CLN_J];
    // (Cartesian slab: the face sums the sweep made, ScSlabFaceSum)
    for (int c = 0; c < 4 && slb; ++c) {
      q[CLN_F2+c] = cs(m,CLN_F2+c,k,j,i);
      q[CLN_FX+c] = cs(m,CLN_FX+c,k,j,i);
      q[CLN_F3+c] = cs(m,CLN_F3+c,k,j,i);
      if (CLN_F1 + c < cs.extent_int(1)) {q[CLN_F1+c] = cs(m,CLN_F1+c,k,j,i);}
    }
    for (int c = 0; c < 10 && slb; ++c) {
      vc_(m,VETRAY_J+c,k,j,i) = (c == 0) ? q[CLN_J] : ((c < 7) ? q[CLN_D11+c-1]
                                                              : hp[c-7] + hm[c-7]);
    }
    if (q[CLN_J] < ray_min*iw_(m,VET_IW_EN,k,j,i)) q[CLN_J] = 0.0;
    StoreMoments(cs, valid, q, m, k, j, i);
    // vet_rotate_avg_ema_wb on the cubed sphere: FaceHalf's half ranges per cell,
    // each about the normal of the face it leaves through: CLN_HO + b = h+ out of the
    // high face of family b, CLN_HO + 3 + b = h- out of the low one (over the cell's J).
    // The ring carries them as cell moments, so a set change keeps a static state's; a
    // face takes its two cells' (a ghost's in its own axes, CsSlot)
    const bool ohk = cs(m,CLN_J,k,j,i) > 0.0 && q[CLN_JR] > 0.0;
    for (int b = 0; b < 3 && hf; ++b) {
      cs(m,CLN_HO+b,k,j,i) = ohk ? oh[b]/q[CLN_JR] : 0.25;
      cs(m,CLN_HO+3+b,k,j,i) = ohk ? ol[b]/q[CLN_JR] : -0.25;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureFromRays
//! \brief closure from the vet_sc moments where J >= closure_ray_min E (else isotropic,
//! marked for the fill); the rays' face fluxes copied everywhere

void RadiationVet::ClosureFromRays() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto vc = vet_cell;
  auto iw_ = iw;
  auto cs = cln_s;
  auto valid = cln_v;
  Real ray_min = cln_raymin;

  const bool f2 = cln_consist;
  DvceArray1D<Real> nb("cln_raynan", 1);   // (a ray J not finite: counted, no silent D)
  SolveFor("cln_rays", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = vc(m,VETRAY_J,k,j,i);
    if (!Kokkos::isfinite(jr)) Kokkos::atomic_add(&nb(0), 1.0);
    for (int n = 0; n < VETRAY_NFS; ++n) {   // (F2, F1, FX, F3)
      const int c = (n < 4) ? CLN_F2 + n : ((n < 8) ? CLN_F1 + n - 4 : CLN_FX + n - 8);
      if (c < cs.extent_int(1)) cs(m,c,k,j,i) = f2 ? vc(m,VETRAY_F2+n,k,j,i) : 0.0;
    }
    cs(m,CLN_JR,k,j,i) = jr;
    if (jr > 0.0 && jr >= ray_min*iw_(m,VET_IW_EN,k,j,i)) {
      for (int a = 0; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = vc(m,VETRAY_D11+a,k,j,i);
      for (int a = 0; a < 3; ++a) {
        const Real hp = vc(m,VETRAY_HP1+a,k,j,i)/jr;
        const Real h = vc(m,VETRAY_H1+a,k,j,i)/jr;
        cs(m,CLN_HP+a,k,j,i) = fmax(hp, 0.0);
        cs(m,CLN_HM+a,k,j,i) = fmin(h - hp, 0.0);
      }
      for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
        cs(m,CLN_KP+a,k,j,i) = vc(m,VETRAY_KP+a,k,j,i)/jr;
      }
      cs(m,CLN_C1,k,j,i) = fmax(vc(m,VETRAY_NL,k,j,i)/jr, 0.0);
      cs(m,CLN_J,k,j,i) = jr;
      valid(m,0,k,j,i) = 1.0;
    } else {
      // isotropic until the fill gives it the closure of its neighbours
      for (int a = 0; a < 3; ++a) {
        cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
        cs(m,CLN_D11+3+a,k,j,i) = 0.0;
        cs(m,CLN_HP+a,k,j,i) = 0.25;
        cs(m,CLN_HM+a,k,j,i) = -0.25;
      }
      for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
        cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
      }
      cs(m,CLN_C1,k,j,i) = 1.0;
      cs(m,CLN_J,k,j,i) = 0.0;
      valid(m,0,k,j,i) = 0.0;
    }
  });
  const Real nn = GlobalReduce(Kokkos::create_mirror_view_and_copy(HostMemSpace(), nb)(0),
                               false);
  if (nn > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: %.0f cells with a ray J not finite "
                "(isotropic)\n", pmy_pack->pmesh->ncycle, nn);
  }
  // Cartesian sweep: ray intensities < 0 set to 0 (none by construction: the PLM foot is
  // bounded, the segment weights >= 0)
  const Real nc = GlobalReduce(Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                                   vray_cnt)(0), false);
  if (nc > 0.0 && global_variable::my_rank == 0) {
    std::printf("CLEAN WARNING cycle %d: %.0f ray intensities < 0 set to 0\n",
                pmy_pack->pmesh->ncycle, nc);
  }
  Kokkos::deep_copy(vray_cnt, 0.0);
}

}  // namespace radvet
