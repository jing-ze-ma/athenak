//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_clean_cell.cpp
//! \brief Cell-centred layout of the clean implicit radiation solve
//! (<rad_m1>/clean_layout = cell): E and F are both cell unknowns of the stage solve.
//! The closure, the gas coupling (exchange Newton, implicit velocity, work) and the time
//! stepping are those of the staggered layout; see rad_m1_clean_cell.hpp for the
//! discrete equations.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1_clean.hpp"
#include "rad_m1/rad_m1_clean_cell.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellOldVector
//! \brief old vector (E, F/c) of the stage: E from CleanOldVector, F of the cells plus
//! the stage increment.  Flux components beyond the mesh dimension are zero.

void RadiationM1::CleanCellOldVector() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool stage = (t2_solve == M1_T2S_STAGE1 || t2_solve == M1_T2S_STAGE2);
  const int nf = pmy_pack->pmesh->three_d ? 3 : (pmy_pack->pmesh->multi_d ? 2 : 1);
  const Real ic = 1.0/c_light;
  auto u0_ = u0;
  auto inc = t2inc;
  auto w = cln_w;
  auto u = cln_x;

  par_for("clc_old", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    u(m,CLC_O,k,j,i) = w(m,CLW_EO,k,j,i);
    for (int b = 0; b < 3; ++b) {
      Real f = 0.0;
      if (b < nf) f = u0_(m,M1_F1+b,k,j,i) + (stage ? inc(m,M1_T2_F1+b,k,j,i) : 0.0);
      u(m,CLC_O+1+b,k,j,i) = ic*f;
    }
  });
  CleanGhostFill(cln_x, CLC_O, 4);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellFaceCoefficients
//! \brief per face: the dissipation weight al, the RC weight and the dissipation speed
//! (CleanCellFaceSpeeds).  The smoothness sensor psi (Jameson type, on E and F/c of the
//! old vector, the larger of the two cells) enters only clean_alpha = tau_sensor.

void RadiationM1::CleanCellFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int nf = three_d ? 3 : (multi_d ? 2 : 1);
  const Real eps = cln_psi_eps;
  CleanArgs x = CleanMakeArgs();
  // on the cubed sphere the opacities of seam ghosts are the neighbour's cells (not the
  // resampled gas ghosts), so that a seam face is single-valued
  if (pmy_pack->pmesh->use_cubed_sphere) {
    CleanGhostFill(opac, 0, opac.extent_int(1), false);
  }
  auto op = opac;
  auto u = cln_x;

  // smoothness sensor |u+ - 2u + u-|/(|u+ - u| + |u - u-|) of every cell and direction
  par_for("clc_psi", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int dir = 0; dir < 3; ++dir) {
      Real p = 0.0;
      if (dir == 0 || (dir == 1 && multi_d) || (dir == 2 && three_d)) {
        const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
        const Real scale = eps*(fabs(u(m,CLC_O,k-dk,j-dj,i-di)) + fabs(u(m,CLC_O,k,j,i))
                                + fabs(u(m,CLC_O,k+dk,j+dj,i+di))) + 1.0e-300;
        for (int q = 0; q <= nf; ++q) {
          const Real um = u(m,CLC_O+q,k-dk,j-dj,i-di);
          const Real u0 = u(m,CLC_O+q,k,j,i);
          const Real up = u(m,CLC_O+q,k+dk,j+dj,i+di);
          p = fmax(p, fabs(up - 2.0*u0 + um)/(fabs(up - u0) + fabs(u0 - um) + scale));
        }
      }
      u(m,CLC_PSI+dir,k,j,i) = p;
    }
  });
  CleanGhostFill(cln_x, CLC_PSI, 3);
  // clean_face_model = milne: B = kP a T^4/kE of the old vector, fixed during the stage
  if (cln_milne) {
    const Real arad_ = arad;
    par_for("clc_milneb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real kap = op(m,M1_OP_E,k,j,i);
      u(m,CLC_MB,k,j,i) = (kap > 0.0) ?
          op(m,M1_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,i)))/kap : 0.0;
    });
    CleanGhostFill(cln_x, CLC_MB, 1);
  }

  CleanCellFaceSpeeds();
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fc(m,CLF_LO,k,j,i) = 0.0;
      fc(m,CLF_UG,k,j,i) = 0.0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellFaceSpeeds
//! \brief dissipation weight al, RC weight (1 - al) theta and dissipation speed of every
//! face (after a closure update: the low-order flags are kept).  force_reference = clean
//! takes the steady RC coefficient (CleanDissipationWeight), so that the scheme's steady
//! states do not depend on h.

void RadiationM1::CleanCellFaceSpeeds() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real h = dt_sub, c = c_light;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const int mode = cln_alpha, pw = cln_dpow;
  const bool wb = (force_ref == M1_FREF_CLEAN);
  CleanArgs x = CleanMakeArgs();
  auto cs = cln_s;
  auto op = opac;
  auto u = cln_x;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_facespeed", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real c1 = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      const Real sray = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      const Real psi = fmax(u(m,CLC_PSI+dir,k,j,i), u(m,CLC_PSI+dir,kl,jl,il));
      const Real dd = fmax(CleanDnn(x, dir, m, k, j, i, k, j, i),
                           CleanDnn(x, dir, m, k, j, i, kl, jl, il));
      const Real dxd = CleanCellFaceDistance(x, dir, m, k, j, i);
      Real rc;
      const Real al = CleanDissipationWeight(mode, chi*dxd, psi, sray, rc, pw);
      fc(m,CLF_AL,k,j,i) = al;
      fc(m,CLF_TH,k,j,i) = wb ? rc*dxd/(h*c) : (1.0 - al)/(1.0 + h*c*chi);
      fc(m,CLF_AF,k,j,i) = c*std::sqrt(fmax(dd, 0.0));
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellMilneSetup
//! \brief clean_face_model = milne: the weight mu of every face and the terms
//! CleanCellMilne adds to the derivatives of the E column (CLF_ML + 1 + row for L,
//! + 5 + row for R), fixed while the closure, the face flags and B are

void RadiationM1::CleanCellMilneSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_milne", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real jl[4][4], jr[4][4], bl[4], br[4];
      CleanCellM1Jacobian(x, dir, m, k, j, i, jl, jr);
      for (int r = 0; r < 4; ++r) {
        bl[r] = jl[r][0];
        br[r] = jr[r][0];
      }
      const Real mu = CleanCellMilneWeight(x, dir, m, k, j, i);
      CleanCellMilne(x, dir, m, k, j, i, jl, jr, nullptr);
      x.fc[dir](m,CLF_ML,k,j,i) = mu;
      for (int r = 0; r < 4; ++r) {
        x.fc[dir](m,CLF_ML+1+r,k,j,i) = jl[r][0] - (1.0 - mu)*bl[r];
        x.fc[dir](m,CLF_ML+5+r,k,j,i) = jr[r][0] - (1.0 - mu)*br[r];
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellConsistency
//! \brief consistency bracket of every interior and open x1 face (clean_consistency).
//! An interior face takes the share w (the rays' streaming weight) of its flux from the
//! rays: F = (1 - w) f(U) + w [F_ray + F_kin(U - Ur)], F_ray the rays' face flux of the
//! sweep (CLN_FX, CLN_F2), whose cell balance is the exchange of the half-segments inside
//! the cell, and F_kin the rays' kinetic upwind flux (their half-range factors on E and
//! y), so that the deviation from the rays' moments Ur = (J, J h) moves only as the rays
//! do; the bracket is w (F_ray - F_kin(Ur)), w in CLF_WK.  An open x1 boundary face
//! carries the rays' flux CLN_F1 less the boundary flux of Ur.  Reflecting, imposed-flux
//! and x2, x3 boundary faces carry none; x3 faces take CLN_F3 as x2 faces CLN_F2.  A
//! physical boundary face takes w of its inner cell (C1 copied to the ghost) for the
//! momentum deviation it carries out.  On the cubed sphere the high x2, x3 faces of a
//! MeshBlock take CLN_FH of its last cell and F_kin the half ranges about the face normal
//! from the two cells' FV intensities (CLF_HK; with clean_face_model = milne also the
//! stress half ranges, CLF_KK): every term of a seam face is computed alike on both
//! sides.  Without clean_consistency (cubed sphere) only these half ranges are filled.

void RadiationM1::CleanCellConsistency() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool cwt = cln_conswt;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  CleanArgs x = CleanMakeArgs();
  const bool ml = x.milne;
  x.milne = false;   // the bracket is that of the moment scheme
  auto fi = cln_fvi;
  auto dir_ = vgd_dir;
  const int nd = x.cub ? vgd_n : 0;

  // the rays' cell moments Ur (y in its frame on the cubed sphere)
  par_for("clc_raymom", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = x.cs(m,CLN_J,k,j,i);
    Real s[3], h[3], q[3][3];
    for (int b = 0; b < 3; ++b) s[b] = x.cs(m,CLN_HP+b,k,j,i) + x.cs(m,CLN_HM+b,k,j,i);
    if (x.cub) {
      CleanRayFlux(x, m, k, j, i, s, h);
      CleanCsFrame(x, m, k, j, i, q);
      for (int b = 0; b < 3; ++b) s[b] = q[b][0]*h[0] + q[b][1]*h[1] + q[b][2]*h[2];
    }
    x.u(m,CLC_U,k,j,i) = jr;
    for (int b = 0; b < 3; ++b) x.u(m,CLC_U+1+b,k,j,i) = jr*s[b];
  });
  CleanGhostFill(cln_x, CLC_U, 4);

  // the bracket of every interior face
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_bracket", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // an open x1 boundary face: the rays' flux through it (CLN_F1 of the inner cell),
      // against the boundary flux of Ur with its source
      for (int r = 0; r < 4; ++r) x.fc[dir](m,CLF_CB+r,k,j,i) = 0.0;
      const bool bdry = CleanBoundaryFace(x, dir, m, k, j, i), blo = (i == x.is);
      const int kl = k - dk, jl = j - dj, il = i - di;
      // clean_consistency_weight: the rays' streaming weight of the face
      Real wb = 1.0;
      if (cwt) {
        const Real c1 = 0.5*(x.cs(m,CLN_C1,kl,jl,il) + x.cs(m,CLN_C1,k,j,i));
        wb *= CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      }
      x.fc[dir](m,CLF_WK,k,j,i) = wb;
      Real nv[3], hj[4] = {0.0, 0.0, 0.0, 0.0};   // J_L, J_R, H+_L, H-_R
      Real kk[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // K+_L, K-_R
      for (int b = 0; b < 3; ++b) nv[b] = CleanNf(x, b, dir, m, k, j, i);
      for (int d = 0; d < nd; ++d) {
        const Real mu = dir_(d,0)*nv[0] + dir_(d,1)*nv[1] + dir_(d,2)*nv[2];
        const Real wl = dir_(d,3)*fi(m,d,kl,jl,il), wr = dir_(d,3)*fi(m,d,k,j,i);
        hj[0] += wl;
        hj[1] += wr;
        hj[(mu > 0.0) ? 2 : 3] += ((mu > 0.0) ? wl : wr)*mu;
        for (int b = 0; b < 3 && ml; ++b) {
          kk[(mu > 0.0) ? b : 3 + b] += ((mu > 0.0) ? wl : wr)*mu*dir_(d,b);
        }
      }
      if (nd > 0) {
        x.fc[dir](m,CLF_HK,k,j,i) = (hj[0] > 0.0) ? hj[2]/hj[0] : 0.25;
        x.fc[dir](m,CLF_HK+1,k,j,i) = (hj[1] > 0.0) ? hj[3]/hj[1] : -0.25;
        for (int b = 0; b < 3 && ml; ++b) {
          x.fc[dir](m,CLF_KK+b,k,j,i) = (hj[0] > 0.0) ? kk[b]/hj[0] : nv[b]/6.0;
          x.fc[dir](m,CLF_KK+3+b,k,j,i) = (hj[1] > 0.0) ? kk[3+b]/hj[1] : nv[b]/6.0;
        }
      }
      if (!x.consist) return;
      if (bdry && (dir != 0 || !CleanOpen(x, m, blo ? 0 : 1))) return;
      // the rays' flux: of the face's own cell, an x1 boundary's or the high face of a
      // MeshBlock on the cubed sphere from the inner cell
      const bool up = bdry ? !blo : (x.cub && dir > 0 && (j > x.je || k > x.ke));
      const int n = bdry ? CLN_F1 : (up ? CLN_FH + 4*(dir - 1) : ((dir == 0) ? CLN_FX :
                                     ((dir == 1) ? CLN_F2 : CLN_F3)));
      Real f[4] = {0.0, 0.0, 0.0, 0.0}, fr[4] = {0.0, 0.0, 0.0, 0.0};
      if (bdry) {
        CleanCellFaceFlux(x, CLC_U, dir, m, k, j, i, true, f);
      } else {
        // the rays' kinetic flux of their own moments Ur: the part w takes in the solve
        Real hp, hm, ul[4], ur[4];
        CleanCellKinetic(x, dir, m, k, j, i, hp, hm);
        CleanCellState(x, CLC_U, m, kl, jl, il, ul);
        CleanCellState(x, CLC_U, m, k, j, i, ur);
        for (int r = 0; r <= x.nf; ++r) f[r] = x.c*(hp*ul[r] + hm*ur[r]);
      }
      for (int r = 0; r <= x.nf; ++r) {
        fr[r] = x.c*x.cs(m,n+r,up ? kl : k,up ? jl : j,up ? il : i);
      }
      // causality: the rays, a steady solution, run ahead of the moments' field; with the
      // FV rays (true face fluxes) the bracket fades out where a cell holds less than
      // half of the rays' J (gone below 0.1 J; the short characteristics' bracket is
      // defined only up to a divergence-free part and stays)
      auto share = [&](const int kk, const int jj, const int ii) {
        return x.fvrays ? CleanSmoothStep(2.5*fmax(x.u(m,CLC_O,kk,jj,ii), 0.0)/
                                          fmax(x.cs(m,CLN_J,kk,jj,ii), 1.0e-300) - 0.25)
                        : 1.0;
      };
      const Real q = fmin(bdry && blo ? 1.0 : share(kl, jl, il),
                          bdry && !blo ? 1.0 : share(k, j, i));
      for (int r = 0; r < 4; ++r) x.fc[dir](m,CLF_CB+r,k,j,i) = q*wb*(fr[r] - f[r]);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellSlopeWeights
//! \brief lagged weights of the van Leer slopes of the cell vector at slot v (the Newton
//! iterate, or Ur for the consistency bracket):
//! pl = dr/(dl + dr), pr = dl/(dl + dr) where dl dr > 0, else 0, so that the half-slope
//! dl dr/(dl + dr) = (pl dl + pr dr)/2 is linear in the state for fixed weights

void RadiationM1::CleanCellSlopeWeights(const int v) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  CleanGhostFill(cln_x, v, 4);
  auto u = cln_x;
  CleanArgs x = CleanMakeArgs();

  // (y in Cartesian components on the cubed sphere, as in the face fluxes)
  par_for("clc_slopes", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int dir = 0; dir < 3; ++dir) {
      const bool on = (dir == 0 || (dir == 1 && multi_d) || (dir == 2 && three_d));
      const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
      Real sm[4], s0[4], sp[4];
      if (on) {
        CleanCellState(x, v, m, k-dk, j-dj, i-di, sm);
        CleanCellState(x, v, m, k, j, i, s0);
        CleanCellState(x, v, m, k+dk, j+dj, i+di, sp);
      }
      for (int q = 0; q < 4; ++q) {
        Real pl = 0.0, pr = 0.0;
        if (on) {
          const Real dl = s0[q] - sm[q];
          const Real dr = sp[q] - s0[q];
          if (dl*dr > 0.0) {
            pl = dr/(dl + dr);
            pr = dl/(dl + dr);
          }
        }
        u(m,CLC_PL+4*dir+q,k,j,i) = pl;
        u(m,CLC_PR+4*dir+q,k,j,i) = pr;
      }
    }
  });
  CleanGhostFill(cln_x, CLC_PL, 24);
  if (!x.cub) return;
  // a seam ghost's weights along this block's axis from the ghost's own index line (the
  // exchanged ones belong to the neighbour's axes), so that the seam faces are second
  // order (both layers of a seam ghost are the neighbour's cells)
  const int no = std::max(indcs.nx2, indcs.nx3);
  par_for("clc_slopes_seam", DevExeSpace(), 0, nmb1, 0, 3, 0, no - 1, is, ie,
  KOKKOS_LAMBDA(const int m, const int sd, const int o, const int i) {
    const int dir = 1 + sd/2, hi = sd%2;
    if (x.mbc(m,2*dir+hi) != BoundaryFlag::panel) return;
    const int g = (dir == 1) ? (hi ? je + 1 : js - 1) : (hi ? ke + 1 : ks - 1);
    const int j = (dir == 1) ? g : js + o, k = (dir == 2) ? g : ks + o;
    if (j > je + 1 || k > ke + 1) return;
    const int dk = (dir == 2), dj = (dir == 1);
    Real sm[4], s0[4], sp[4];
    CleanCellState(x, v, m, k-dk, j-dj, i, sm);
    CleanCellState(x, v, m, k, j, i, s0);
    CleanCellState(x, v, m, k+dk, j+dj, i, sp);
    for (int q = 0; q < 4; ++q) {
      const Real dl = s0[q] - sm[q], dr = sp[q] - s0[q];
      const bool on = dl*dr > 0.0;
      u(m,CLC_PL+4*dir+q,k,j,i) = on ? dr/(dl + dr) : 0.0;
      u(m,CLC_PR+4*dir+q,k,j,i) = on ? dl/(dl + dr) : 0.0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellFirstIterate
//! \brief first iterate: E of CleanFirstIterate, F of the old vector

void RadiationM1::CleanCellFirstIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  auto u = cln_x;
  par_for("clc_iter0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    u(m,CLC_K,k,j,i) = w(m,CLW_EK,k,j,i);
    for (int b = 1; b < 4; ++b) u(m,CLC_K+b,k,j,i) = u(m,CLC_O+b,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellRightHandSide
//! \brief B_E = E_old + h c kP a T^4 - r h c kE E_k - W + h esrc, B_y = y_old, each less
//! h div of the E-independent face fluxes (exchange linearised about the iterate as in
//! CleanRightHandSide; esrc the explicit energy source of <rad_m1>, an external one),
//! and the initial guess X = the iterate; with gas coupling B_y adds the iterate's part
//! of the implicit gas velocity (CleanCellDrag)

void RadiationM1::CleanCellRightHandSide() {
  if (cln_feedback) CleanGhostFill(cln_x, CLC_K, 4);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  const Real arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  const bool eso = esrc_on;
  auto es = esrc;

  const int vfreeze = cln_vfreeze;
  if (cln_milne) CleanCellMilneSetup();
  par_for("clc_rhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real t = x.w(m,CLW_T,k,j,i);
    const Real hckp = x.h*x.c*op(m,M1_OP_P,k,j,i);
    const Real hcke = x.h*x.c*op(m,M1_OP_E,k,j,i);
    Real rr = 0.0, work = 0.0;
    if (x.feedback) {
      const Real d = uh(m,IDN,k,j,i);
      Real e, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e, p, chi_rho, chi_t, cv);
      const Real s = 4.0*hckp*arad_*t*t*t;
      rr = s/(d*cv + s);
      // work of the radiation force on the gas: the kinetic energy the gas gains
      Real pg[3];
      CleanGasMomentum(x, uh, m, k, j, i, pg);
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (pg[n] + (vfreeze == 2 ? 0.0 : 0.5*dm))*dm/d;
      }
    }
    Real sdiv[4];
    CleanCellSourceDivergence(x, m, k, j, i, sdiv);
    x.w(m,CLW_RR,k,j,i) = rr;
    x.w(m,CLW_WK,k,j,i) = work;
    x.u(m,CLC_B,k,j,i) = x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
                         - rr*hcke*fmax(x.w(m,CLW_EK,k,j,i), 0.0) - work - sdiv[0]
                         + (eso ? x.h*es(m,k,j,i) : 0.0);
    for (int b = 1; b < 4; ++b) x.u(m,CLC_B+b,k,j,i) = x.u(m,CLC_O+b,k,j,i) - sdiv[b];
    if (x.feedback) {
      Real g[3][3], dk[4];
      CleanCellDrag(x, op(m,M1_OP_T,k,j,i), uh(m,IDN,k,j,i), m, k, j, i, g);
      CleanCellDivergence(x, CLC_K, m, k, j, i, false, dk);
      for (int b = 0; b < x.nf; ++b) {
        for (int c = 0; c < x.nf; ++c) {
          x.u(m,CLC_B+1+b,k,j,i) += g[b][c]*(x.u(m,CLC_K+1+c,k,j,i) + dk[1+c]);
        }
      }
    }
    for (int b = 0; b < 4; ++b) x.u(m,CLC_X+b,k,j,i) = x.u(m,CLC_K+b,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellApply
//! \brief u(dst) = A u(src) over the active cells: the local exchange and flux sources
//! plus h div of the face fluxes; flux components beyond the mesh dimension are identity.
//! The cells at least ng from the block faces are done while the ghosts are exchanged.

void RadiationM1::CleanCellApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  auto op = opac;
  auto uh = FluidRef::Get(pmy_pack).u0;
  const int ng = indcs.ng;
  const bool md = x.multi_d, td = x.three_d;

  for (int part = 1; part <= 2; ++part) {
  CleanGhostFill(cln_x, src, 4, true, part);
  par_for("clc_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (CleanInner(k, j, i, ng, is, ie, js, je, ks, ke, md, td) != (part == 1)) return;
    Real div[4];
    CleanCellDivergence(x, src, m, k, j, i, false, div);
    const Real rr = x.feedback ? x.w(m,CLW_RR,k,j,i) : 0.0;
    const Real e = x.u(m,src,k,j,i);
    const Real chi = op(m,M1_OP_T,k,j,i);
    x.u(m,dst,k,j,i) = e*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr)) + div[0];
    Real g[3][3];
    if (x.feedback) CleanCellDrag(x, chi, uh(m,IDN,k,j,i), m, k, j, i, g);
    for (int b = 0; b < 3; ++b) {
      Real y = x.u(m,src+1+b,k,j,i);
      if (b < x.nf) {
        y = y*(1.0 + x.h*x.c*chi) + div[1+b];
        if (x.feedback) {
          y -= x.h*chi*x.w(m,CLW_A1+b,k,j,i)*e;
          for (int c = 0; c < x.nf; ++c) y += g[b][c]*(x.u(m,src+1+c,k,j,i) + div[1+c]);
        }
      }
      x.u(m,dst+1+b,k,j,i) = y;
    }
  });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellPreconditionerSetup
//! \brief 4x4 blocks of the preconditioner, from A with first-order jumps: the diagonal
//! block of every cell (all directions) and the blocks of its two neighbours along each
//! direction (states of cells i-1, i, i+1).  The RC terms in cells i+-2 and the
//! reconstruction are left to the Krylov solve.  Then the blocks of the coarse levels.

void RadiationM1::CleanCellPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.recon = false;
  auto op = opac;
  auto uh = FluidRef::Get(pmy_pack).u0;

  par_for("clc_pcsetup", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real a0[4][4];
    Real jlp[4][4], jrp[4][4], jlm[4][4], jrm[4][4];
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) a0[r][c] = 0.0;
    }
    const Real rr = x.feedback ? x.w(m,CLW_RR,k,j,i) : 0.0;
    const Real chi = op(m,M1_OP_T,k,j,i);
    a0[0][0] = 1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr);
    for (int b = 0; b < 3; ++b) {
      a0[1+b][1+b] = (b < x.nf) ? 1.0 + x.h*x.c*chi : 1.0;
      if (b < x.nf && x.feedback) a0[1+b][0] = -x.h*chi*x.w(m,CLW_A1+b,k,j,i);
    }
    if (x.feedback) {
      Real g[3][3];
      CleanCellDrag(x, chi, uh(m,IDN,k,j,i), m, k, j, i, g);
      for (int b = 0; b < x.nf; ++b) {
        for (int c = 0; c < x.nf; ++c) a0[1+b][1+c] += g[b][c];
      }
    }
    for (int dir = 0; dir < 3; ++dir) {
      if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
      const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
      const int sm = CleanLineSlotM(dir), sp = CleanLineSlotP(dir);
      CleanCellFaceJacobian(x, dir, m, k+dk, j+dj, i+di, jlp, jrp);
      CleanCellFaceJacobian(x, dir, m, k, j, i, jlm, jrm);
      if (x.sph || x.cub) {
        Real bp[4][4], bm[4][4], bpn[4][4], bmn[4][4];
        CleanCellBlock(x, dir, m, k+dk, j+dj, i+di, k, j, i, k, j, i, jlp, bp);
        CleanCellBlock(x, dir, m, k, j, i, k, j, i, k, j, i, jrm, bm);
        CleanCellBlock(x, dir, m, k+dk, j+dj, i+di, k, j, i, k+dk, j+dj, i+di, jrp, bpn);
        CleanCellBlock(x, dir, m, k, j, i, k, j, i, k-dk, j-dj, i-di, jlm, bmn);
        for (int r = 0; r < 4; ++r) {
          for (int c = 0; c < 4; ++c) {
            a0[r][c] += bp[r][c] - bm[r][c];
            x.u(m,sp+4*r+c,k,j,i) = bpn[r][c];
            x.u(m,sm+4*r+c,k,j,i) = -bmn[r][c];
          }
        }
        continue;
      }
      const Real s = x.h*x.idx[dir];
      for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
          a0[r][c] += s*(jlp[r][c] - jrm[r][c]);
          x.u(m,sp+4*r+c,k,j,i) = s*jrp[r][c];
          x.u(m,sm+4*r+c,k,j,i) = -s*jlm[r][c];
        }
      }
    }
    if (x.sph) {
      Real g[3];
      CleanCellSphSource(x, m, k, j, i, g);
      for (int b = 0; b < x.nf; ++b) a0[1+b][0] += g[b];
    }
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) x.u(m,CLC_A0+4*r+c,k,j,i) = a0[r][c];
    }
  });
  CleanCellCoarseSetup();
}

//----------------------------------------------------------------------------------------
//! \fn CleanLevel RadiationM1::CleanCellLevel
//! \brief level l of the preconditioner: l = 0 is the mesh (blocks in cln_x), l > 0 the
//! coarse levels of the multigrid (cln_mg[l-1])

CleanLevel RadiationM1::CleanCellLevel(const int l) {
  CleanLevel lv;
  if (l == 0) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    lv.u = cln_x;
    lv.o[0] = indcs.is;
    lv.o[1] = indcs.js;
    lv.o[2] = indcs.ks;
    lv.n[0] = indcs.nx1;
    lv.n[1] = indcs.nx2;
    lv.n[2] = indcs.nx3;
    lv.a0 = CLC_A0;
    for (int d = 0; d < 3; ++d) {
      lv.sm[d] = CleanLineSlotM(d);
      lv.sp[d] = CleanLineSlotP(d);
    }
  } else {
    lv.u = cln_mg[l-1];
    for (int d = 0; d < 3; ++d) {
      lv.o[d] = 0;
      lv.n[d] = cln_mgn[l-1][d];
      lv.sm[d] = CLM_M + 32*d;
      lv.sp[d] = CLM_P + 32*d;
    }
    lv.a0 = CLM_A0;
  }
  return lv;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellLineSolve
//! \brief u(dst) = L^{-1} u(src) on a level, L the block-tridiagonal part of its operator
//! along direction dir (diagonal blocks and the neighbour blocks along dir), one team
//! per line in team scratch: rounds of block parallel cyclic reduction split the line
//! into P = 2^k interleaved systems of at most cln_chunk rows, each solved by block
//! Thomas (cln_chunk >= the line: plain Thomas; 1: full PCR)

void RadiationM1::CleanCellLineSolve(const CleanLevel &lv, const int dir, const int src,
                                     const int dst) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int is = lv.o[0], js = lv.o[1], ks = lv.o[2];
  const int a0 = lv.a0, sm = lv.sm[dir], sp = lv.sp[dir];
  // lines run along dir; (p, q) label them
  const int np = (dir == 2) ? lv.n[1] : lv.n[2];
  const int nq = (dir == 0) ? lv.n[1] : lv.n[0];
  const int nl = lv.n[dir];
  auto u = lv.u;
  int nsys = 1;
  while (nl > cln_chunk*nsys) nsys *= 2;
  // row n: A x_n-P + B x_n + C x_n+P = r (52 values) in buffer 0 or 1
  const size_t scr = ScrArray1D<Real>::shmem_size(104*nl);
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb1 + 1)*np*nq, Kokkos::AUTO);
  const int lev = TeamScratchLevel(scr, 0);
  Kokkos::parallel_for("clc_line", policy.set_scratch_size(lev, Kokkos::PerTeam(scr)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/(np*nq);
    const int p = (tm.league_rank() - m*np*nq)/nq, q = tm.league_rank() % nq;
    auto cell = [&](const int n, int &k, int &j, int &i) {
      k = (dir == 2) ? ks + n : ks + p;
      j = (dir == 1) ? js + n : ((dir == 2) ? js + p : js + q);
      i = (dir == 0) ? is + n : is + q;
    };
    ScrArray1D<Real> w(tm.team_scratch(lev), 104*nl);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nl), [&](const int n) {
      int k, j, i;
      cell(n, k, j, i);
      for (int e = 0; e < 16; ++e) {
        w(52*n + e) = (n > 0) ? u(m,sm+e,k,j,i) : 0.0;
        w(52*n + 16 + e) = u(m,a0+e,k,j,i);
        w(52*n + 32 + e) = (n < nl - 1) ? u(m,sp+e,k,j,i) : 0.0;
      }
      for (int r = 0; r < 4; ++r) w(52*n + 48 + r) = u(m,src+r,k,j,i);
    });
    tm.team_barrier();
    // PCR rounds: row n takes rows n -+ st out, f = -X Bnb^-1 (X = A or C of row n) from
    // Bnb^T f^T = -X^T
    int o = 0;
    for (int st = 1; st < nsys; st *= 2) {
      const int d = 52*nl - o;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nl), [&](const int n) {
        Real an[4][4], bn[4][4], cn[4][4], rn[4];
        for (int r = 0; r < 4; ++r) {
          for (int c = 0; c < 4; ++c) {
            an[r][c] = 0.0;
            bn[r][c] = w(o + 52*n + 16 + 4*r + c);
            cn[r][c] = 0.0;
          }
          rn[r] = w(o + 52*n + 48 + r);
        }
        for (int sd = -1; sd <= 1; sd += 2) {
          const int nb = n + sd*st;
          if (nb < 0 || nb >= nl) continue;
          const int ox = o + 52*n + ((sd < 0) ? 0 : 32), ob = o + 52*nb;
          Real bt[4][4], f[4][5];
          for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
              bt[r][c] = w(ob + 16 + 4*c + r);
              f[r][c] = -w(ox + 4*c + r);
            }
          }
          CleanBlockSolve(bt, f, 4);   // f[l][r] = (-X Bnb^-1)[r][l]
          for (int r = 0; r < 4; ++r) {
            for (int l = 0; l < 4; ++l) {
              const Real g = f[l][r];
              for (int c = 0; c < 4; ++c) {
                if (sd < 0) {
                  an[r][c] += g*w(ob + 4*l + c);
                } else {
                  cn[r][c] += g*w(ob + 32 + 4*l + c);
                }
                bn[r][c] += g*w(ob + ((sd < 0) ? 32 : 0) + 4*l + c);
              }
              rn[r] += g*w(ob + 48 + l);
            }
          }
        }
        for (int r = 0; r < 4; ++r) {
          for (int c = 0; c < 4; ++c) {
            w(d + 52*n + 4*r + c) = an[r][c];
            w(d + 52*n + 16 + 4*r + c) = bn[r][c];
            w(d + 52*n + 32 + 4*r + c) = cn[r][c];
          }
          w(d + 52*n + 48 + r) = rn[r];
        }
      });
      tm.team_barrier();
      o = d;
    }
    // block Thomas on each system s (rows s, s + P, ...): forward D = B - A C_n-P,
    // C_n = D^-1 C, d_n = D^-1 (r - A d_n-P) (C_n, d_n in the other buffer); backward
    // x_n = d_n - C_n x_n+P
    const int t = 52*nl - o;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nsys), [&](const int sy) {
      for (int n = sy; n < nl; n += nsys) {
        Real a[4][4], b[4][5];
        for (int r = 0; r < 4; ++r) {
          Real rhs = w(o + 52*n + 48 + r);
          for (int c = 0; c < 4; ++c) {
            a[r][c] = w(o + 52*n + 16 + 4*r + c);
            b[r][c] = w(o + 52*n + 32 + 4*r + c);
          }
          if (n >= nsys) {
            for (int l = 0; l < 4; ++l) {
              const Real aml = w(o + 52*n + 4*r + l);
              if (aml == 0.0) continue;
              for (int c = 0; c < 4; ++c) a[r][c] -= aml*w(t + 52*(n-nsys) + 4*l + c);
              rhs -= aml*w(t + 52*(n-nsys) + 48 + l);
            }
          }
          b[r][4] = rhs;
        }
        CleanBlockSolve(a, b, 5);
        for (int r = 0; r < 4; ++r) {
          for (int c = 0; c < 4; ++c) w(t + 52*n + 4*r + c) = b[r][c];
          w(t + 52*n + 48 + r) = b[r][4];
        }
      }
      for (int n = sy + nsys*((nl - 1 - sy)/nsys) - nsys; n >= 0; n -= nsys) {
        for (int r = 0; r < 4; ++r) {
          for (int c = 0; c < 4; ++c) {
            w(t + 52*n + 48 + r) -= w(t + 52*n + 4*r + c)*w(t + 52*(n+nsys) + 48 + c);
          }
        }
      }
    });
    tm.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nl), [&](const int n) {
      int k, j, i;
      cell(n, k, j, i);
      for (int r = 0; r < 4; ++r) u(m,dst+r,k,j,i) = w(t + 52*n + 48 + r);
    });
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellLevelApply
//! \brief u(dst) = A_l u(src) with the blocks of a level: diagonal block plus the
//! neighbour blocks of the neighbours inside the level

void RadiationM1::CleanCellLevelApply(const CleanLevel &lv, const int src,
                                      const int dst) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel l = lv;

  par_for("clc_lvapply", DevExeSpace(), 0, nmb1, l.o[2], l.o[2]+l.n[2]-1,
          l.o[1], l.o[1]+l.n[1]-1, l.o[0], l.o[0]+l.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real y[4];
    for (int r = 0; r < 4; ++r) {
      y[r] = 0.0;
      for (int c = 0; c < 4; ++c) y[r] += l.u(m,l.a0+4*r+c,k,j,i)*l.u(m,src+c,k,j,i);
    }
    const int pos[3] = {i - l.o[0], j - l.o[1], k - l.o[2]};
    for (int d = 0; d < 3; ++d) {
      const int dk = (d == 2), dj = (d == 1), di = (d == 0);
      for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
          if (pos[d] > 0) {
            y[r] += l.u(m,l.sm[d]+4*r+c,k,j,i)*l.u(m,src+c,k-dk,j-dj,i-di);
          }
          if (pos[d] < l.n[d] - 1) {
            y[r] += l.u(m,l.sp[d]+4*r+c,k,j,i)*l.u(m,src+c,k+dk,j+dj,i+di);
          }
        }
      }
    }
    for (int r = 0; r < 4; ++r) l.u(m,dst+r,k,j,i) = y[r];
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellCoarseSetup
//! \brief blocks of the coarse levels by aggregation (Galerkin product with
//! piecewise-constant prolongation and summing restriction): a coarse block is the sum of
//! the blocks of its fine cells, the couplings between fine cells of the same coarse cell
//! folded into the diagonal block.  On curvilinear meshes the mesh rows are summed with
//! their cell volumes, so that the fluxes between fine cells cancel as in the
//! conservative operator.

void RadiationM1::CleanCellCoarseSetup() {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nlev = static_cast<int>(cln_mg.size()) + 1;
  const bool curv = pmy_pack->pmesh->use_spherical_polar ||
                    pmy_pack->pmesh->use_cubed_sphere;
  auto vol = pmy_pack->pcoord->volume;
  for (int lc = 1; lc < nlev; ++lc) {
    const CleanLevel f = CleanCellLevel(lc-1), c = CleanCellLevel(lc);
    const int r[3] = {f.n[0]/c.n[0], f.n[1]/c.n[1], f.n[2]/c.n[2]};
    const bool vw = (lc == 1) && curv;

    par_for("clc_galerkin", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1,
            0, c.n[0]-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
      for (int s = 0; s < CLM_B; ++s) c.u(m,s,kc,jc,ic) = 0.0;
      for (int a2 = 0; a2 < r[2]; ++a2) {
        for (int a1 = 0; a1 < r[1]; ++a1) {
          for (int a0 = 0; a0 < r[0]; ++a0) {
            const int k = f.o[2] + r[2]*kc + a2;
            const int j = f.o[1] + r[1]*jc + a1;
            const int i = f.o[0] + r[0]*ic + a0;
            const int a[3] = {a0, a1, a2};
            const Real v = vw ? vol(m,k,j,i) : 1.0;
            for (int q = 0; q < 16; ++q) c.u(m,c.a0+q,kc,jc,ic) += v*f.u(m,f.a0+q,k,j,i);
            for (int d = 0; d < 3; ++d) {
              if (f.n[d] == 1) continue;
              const int dm = (a[d] > 0) ? c.a0 : c.sm[d];
              const int dp = (a[d] < r[d] - 1) ? c.a0 : c.sp[d];
              for (int q = 0; q < 16; ++q) {
                c.u(m,dm+q,kc,jc,ic) += v*f.u(m,f.sm[d]+q,k,j,i);
                c.u(m,dp+q,kc,jc,ic) += v*f.u(m,f.sp[d]+q,k,j,i);
              }
            }
          }
        }
      }
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellSmooth
//! \brief nsweep sweeps of line relaxation on level l for A_l u(dst) = u(src): along x1,
//! x2, x3 in turn, each a block-tridiagonal solve of the residual on all lines of that
//! direction (multiplicative).  first: u(dst) starts from the x1-line solve of u(src).
//! The residual on the mesh uses the full operator, on coarse levels their blocks.

void RadiationM1::CleanCellSmooth(const int l, const int src, const int dst,
                                  const bool first, const int nsweep) {
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel lv = CleanCellLevel(l);
  const int t1 = (l == 0) ? CLC_T1 : CLM_T1, t2 = (l == 0) ? CLC_T2 : CLM_T2;
  const int is = lv.o[0], ie = lv.o[0] + lv.n[0] - 1;
  const int js = lv.o[1], je = lv.o[1] + lv.n[1] - 1;
  const int ks = lv.o[2], ke = lv.o[2] + lv.n[2] - 1;
  auto u = lv.u;
  if (first) CleanCellLineSolve(lv, 0, src, dst);
  for (int n = first ? 1 : 0; n < 3*nsweep; ++n) {
    const int dir = n % 3;
    if (lv.n[dir] == 1) continue;
    if (l == 0) {
      CleanCellApply(dst, t1);
    } else {
      CleanCellLevelApply(lv, dst, t1);
    }
    par_for("clc_pcres", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < 4; ++c) u(m,t1+c,k,j,i) = u(m,src+c,k,j,i) - u(m,t1+c,k,j,i);
    });
    CleanCellLineSolve(lv, dir, t1, t2);
    par_for("clc_pccor", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < 4; ++c) u(m,dst+c,k,j,i) += u(m,t2+c,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellVCycle
//! \brief u(dst) = M^{-1} u(src) on level l by a multigrid V-cycle: line-relaxation
//! smoothing, the residual summed over each coarse cell (with the cell volumes from a
//! curvilinear mesh, as CleanCellCoarseSetup), the V-cycle on level l+1, its solution
//! added to every fine cell of its coarse cell, smoothing again.  The coarsest level (a
//! few cells per direction) is relaxed with four times the sweeps.

void RadiationM1::CleanCellVCycle(const int l, const int src, const int dst) {
  const int nlev = static_cast<int>(cln_mg.size()) + 1;
  if (l == nlev - 1) {
    CleanCellSmooth(l, src, dst, true, (l == 0) ? cln_pcsweeps : 4*cln_pcsweeps);
    return;
  }
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const CleanLevel f = CleanCellLevel(l), c = CleanCellLevel(l+1);
  const int t1 = (l == 0) ? CLC_T1 : CLM_T1;
  const int r[3] = {f.n[0]/c.n[0], f.n[1]/c.n[1], f.n[2]/c.n[2]};

  // pre-smoothing and the residual
  CleanCellSmooth(l, src, dst, true, cln_pcsweeps);
  if (l == 0) {
    CleanCellApply(dst, t1);
  } else {
    CleanCellLevelApply(f, dst, t1);
  }
  auto uf = f.u;
  auto uc = c.u;
  const bool vw = (l == 0) && (pmy_pack->pmesh->use_spherical_polar ||
                                pmy_pack->pmesh->use_cubed_sphere);
  auto vol = pmy_pack->pcoord->volume;
  par_for("clc_restrict", DevExeSpace(), 0, nmb1, 0, c.n[2]-1, 0, c.n[1]-1, 0, c.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int ic) {
    for (int q = 0; q < 4; ++q) {
      Real s = 0.0;
      for (int a2 = 0; a2 < r[2]; ++a2) {
        for (int a1 = 0; a1 < r[1]; ++a1) {
          for (int a0 = 0; a0 < r[0]; ++a0) {
            const int k = f.o[2] + r[2]*kc + a2;
            const int j = f.o[1] + r[1]*jc + a1;
            const int i = f.o[0] + r[0]*ic + a0;
            s += (vw ? vol(m,k,j,i) : 1.0)*(uf(m,src+q,k,j,i) - uf(m,t1+q,k,j,i));
          }
        }
      }
      uc(m,CLM_B+q,kc,jc,ic) = s;
    }
  });

  // coarse-level correction
  CleanCellVCycle(l+1, CLM_B, CLM_X);
  par_for("clc_prolong", DevExeSpace(), 0, nmb1, f.o[2], f.o[2]+f.n[2]-1,
          f.o[1], f.o[1]+f.n[1]-1, f.o[0], f.o[0]+f.n[0]-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int kc = (k - f.o[2])/r[2], jc = (j - f.o[1])/r[1], ic = (i - f.o[0])/r[0];
    for (int q = 0; q < 4; ++q) uf(m,dst+q,k,j,i) += uc(m,CLM_X+q,kc,jc,ic);
  });

  // post-smoothing
  CleanCellSmooth(l, src, dst, false, cln_pcsweeps);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellPrecondition
//! \brief u(dst) = M^{-1} u(src): one multigrid V-cycle from the mesh level

void RadiationM1::CleanCellPrecondition(const int src, const int dst) {
  CleanCellVCycle(0, src, dst);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanCellNewtonUpdate
//! \brief new iterate from the solution X: with gas coupling E and T by
//! CleanNewtonUpdate, F/c from X.  Returns the largest relative change: that of
//! CleanNewtonUpdate and max |dE|, |dy| over max E.

Real RadiationM1::CleanCellNewtonUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  int nmb1 = nmb - 1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const bool feedback = cln_feedback;
  auto w = cln_w;
  auto u = cln_x;

  Real de_newton = 0.0;
  if (feedback) {
    par_for("clc_xe", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_X,k,j,i) = u(m,CLC_X,k,j,i);
    });
    de_newton = CleanNewtonUpdate();
  }
  Real dmax = 0.0, emax = 0.0;
  Kokkos::parallel_reduce("clc_newton", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &d_max, Real &e_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real e = feedback ? w(m,CLW_EK,k,j,i) : u(m,CLC_X,k,j,i);
    Real d = fabs(e - u(m,CLC_K,k,j,i));
    u(m,CLC_K,k,j,i) = e;
    for (int b = 1; b < 4; ++b) {
      d = fmax(d, fabs(u(m,CLC_X+b,k,j,i) - u(m,CLC_K+b,k,j,i)));
      u(m,CLC_K+b,k,j,i) = u(m,CLC_X+b,k,j,i);
    }
    d_max = fmax(d_max, d);
    e_max = fmax(e_max, fabs(e));
  }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(emax));
  return fmax(de_newton, CleanReduce(dmax, true)/fmax(CleanReduce(emax, true), 1.0e-300));
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanCellVelocityUpdate
//! \brief momentum change of the iterate, -[(y - y_old) + h div fY]/c per cell (minus the
//! radiation momentum change besides its flux divergence) less h times the reference
//! force (CleanRefForce; force_reference = clean takes it from the first iterate of the
//! first coupled solve), and the advective coefficients of the new velocity.  Returns
//! max |dv| over max |v|.

Real RadiationM1::CleanCellVelocityUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CleanGhostFill(cln_x, CLC_K, 4);
  CleanArgs x = CleanMakeArgs();
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto op = opac;
  const Real arad_ = arad;
  const int fref = force_ref;
  const bool vzero = cln_vzero;
  const bool cap = (fref == M1_FREF_CLEAN) && !cln_fref_set && !vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;

  Real dvmax = 0.0, vmax = 0.0, femax = 0.0, fvmax = 0.0;
  Kokkos::parallel_reduce("clc_velocity", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fe_max, Real &fv_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real d = uh(m,IDN,k,j,i);
    const Real e = x.u(m,CLC_K,k,j,i), hc = x.h*x.c;
    const Real chi = op(m,M1_OP_T,k,j,i);
    Real div[4], adiv[4] = {0.0, 0.0, 0.0, 0.0};
    CleanCellDivergence(x, CLC_K, m, k, j, i, true, div, adiv);
    // round-off levels of E (its terms over E) and of dm (the terms of the flux rows,
    // with that of E in the stress)
    const Real ek = fabs(e)*(1.0 + hc*op(m,M1_OP_E,k,j,i)) + fabs(x.w(m,CLW_EO,k,j,i));
    const Real re = (ek + adiv[0] + hc*op(m,M1_OP_P,k,j,i)*arad_*
                     SQR(SQR(x.w(m,CLW_T,k,j,i))))/fmax(fabs(e), 1.0e-300);
    fe_max = fmax(fe_max, re);
    Real v2 = 0.0, p[3];
    CleanGasMomentum(x, uh, m, k, j, i, p);
    for (int b = 0; b < 3; ++b) {
      Real dm = 0.0;
      if (b < x.nf) {
        const Real y = x.u(m,CLC_K+1+b,k,j,i), yo = x.u(m,CLC_O+1+b,k,j,i);
        dm = -(y - yo + div[1+b])/x.c;
        if (cap || vzero) {
          if (cap) fr(m,b,k,j,i) = dm/x.h;
          dm = 0.0;
        } else {
          dm -= x.h*CleanRefForce(fref, fr, aref, d, b, m, k, j, i);
        }
        fv_max = fmax(fv_max, (fabs(y)*(1.0 + hc*chi) + fabs(yo) + adiv[1+b]*(1.0 + re)
                               + x.h*chi*fabs(x.w(m,CLW_A1+b,k,j,i)*e))/(x.c*d));
      }
      dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+b,k,j,i))/d);
      x.w(m,CLW_DM+b,k,j,i) = dm;
      v2 += SQR((p[b] + dm)/d);
    }
    v_max = fmax(v_max, std::sqrt(v2));
  }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(femax),
  Kokkos::Max<Real>(fvmax));
  if (cap) cln_fref_set = true;
  cln_fe = std::numeric_limits<Real>::epsilon()*femax;
  cln_fv = std::numeric_limits<Real>::epsilon()*fvmax/fmax(vmax, 1.0e-300);
  if (cln_vfreeze != 1) CleanAdvectiveCoefficients();
  return CleanReduce(dvmax, true)/fmax(CleanReduce(vmax, true), 1.0e-300);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellCopyIterate
//! \brief the solution is the last iterate

void RadiationM1::CleanCellCopyIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u = cln_x;
  par_for("clc_copy", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int b = 0; b < 4; ++b) u(m,CLC_X+b,k,j,i) = u(m,CLC_K+b,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanCellCountNegative
//! \brief CleanCountNegative on E of the solution; also counts the cells with
//! |F| > c E beyond the solver tolerance (cln_nnreal, diagnostic only).  Returns the
//! number of negative cells.

int RadiationM1::CleanCellCountNegative() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  int nmb1 = nmb - 1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w = cln_w;
  auto u = cln_x;

  par_for("clc_neg_in", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_X,k,j,i) = u(m,CLC_X,k,j,i);
    w(m,CLW_B,k,j,i) = u(m,CLC_B,k,j,i);
  });
  const int nneg = CleanCountNegative();
  Real bmax = 0.0;
  Kokkos::parallel_reduce("clc_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(w(m,CLW_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));
  const Real zero = 10.0*cln_tol*CleanReduce(bmax, true);
  Real nreal = 0.0;
  Kokkos::parallel_reduce("clc_neg_out", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real e = w(m,CLW_X,k,j,i);
    u(m,CLC_X,k,j,i) = e;
    const Real y = std::sqrt(SQR(u(m,CLC_X+1,k,j,i)) + SQR(u(m,CLC_X+2,k,j,i)) +
                             SQR(u(m,CLC_X+3,k,j,i)));
    if (y > (1.0 + 1.0e-12)*e + zero) n += 1.0;
  }, Kokkos::Sum<Real>(nreal));
  cln_nnreal = static_cast<int>(CleanReduce(nreal, false));
  return nneg;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanCellMarkFaces
//! \brief marks every face of the cells whose solution has E < 0 (beyond the solver
//! tolerance, 10 clean_tol max|B_E|) as low-order.  |F| > c E is left alone: a beam
//! overshoots |F| = c E slightly, and low-order faces there smear it into E < 0.  The
//! flags of these cells (cln_v slot 1) take the ghost exchange, so that a face shared by
//! two MeshBlocks is marked on both.  Returns the number of such cells.

int RadiationM1::CleanCellMarkFaces() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto u = cln_x;
  auto bad = cln_v;

  Real bmax = 0.0;
  Kokkos::parallel_reduce("clc_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(u(m,CLC_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));
  const Real zero = 10.0*cln_tol*CleanReduce(bmax, true);
  Real nbad = 0.0;
  Kokkos::parallel_reduce("clc_bad", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const bool b = !(u(m,CLC_X,k,j,i) >= -zero);
    bad(m,1,k,j,i) = b ? 1.0 : 0.0;
    if (b) n += 1.0;
  }, Kokkos::Sum<Real>(nbad));
  nbad = CleanReduce(nbad, false);
  if (nbad == 0.0) return 0;
  CleanGhostFill(cln_v, 1, 1);

  // faces on a physical boundary are marked only with clean_pos_boundary; with
  // clean_pos_ghost they take the upwind stress against the ghost instead
  const bool posbdry = cln_posbdry, ug = cln_posghost;
  CleanArgs x = CleanMakeArgs();
  for (int dir = 0; dir < x.nf; ++dir) {
    auto fc = x.fc[dir];
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_mark", DevExeSpace(), 0, nmb-1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real bl = bad(m,1,k-dk,j-dj,i-di), br = bad(m,1,k,j,i);
      if (!CleanBoundaryFace(x, dir, m, k, j, i)) {
        if (bl > 0.0 || br > 0.0) fc(m,CLF_LO,k,j,i) = 1.0;
        return;
      }
      const int pos = (dir == 0) ? i : ((dir == 1) ? j : k);
      const int lo = (dir == 0) ? is : ((dir == 1) ? js : ks);
      if (((pos == lo) ? br : bl) > 0.0 && (posbdry || ug)) {
        fc(m,posbdry ? CLF_LO : CLF_UG,k,j,i) = 1.0;
      }
    });
  }
  return static_cast<int>(nbad);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellAllLowOrder
//! \brief marks every face low-order (faces on a physical boundary only with
//! clean_pos_boundary; otherwise their stress takes the upwind form)

void RadiationM1::CleanCellAllLowOrder() {
  cln_all_low = true;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool posbdry = cln_posbdry;
  CleanArgs x = CleanMakeArgs();
  for (int dir = 0; dir < x.nf; ++dir) {
    auto fc = x.fc[dir];
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_alllow", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (posbdry || !CleanBoundaryFace(x, dir, m, k, j, i)) fc(m,CLF_LO,k,j,i) = 1.0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellDiagnose
//! \brief diagnostic (clean_diag_cells): for the cells with |F| > c E in the solution,
//! the smallest eigenvalue of D - f f^T (in-plane) for f = F/(cE) of the old vector and
//! of the solution, the rays' D, C1 and flux factor h = H/J against the moment f

void RadiationM1::CleanCellDiagnose() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks;
  auto u = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_x);
  auto cs = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_s);
  auto lmin = [](Real a, Real b, Real d) {
    return 0.5*(a + d) - std::sqrt(0.25*(a - d)*(a - d) + b*b);
  };
  int nbad = 0, nold = 0, nnew = 0, nray = 0, nprint = 0;
  for (int j = js; j <= je; ++j) {
    for (int i = is; i <= ie; ++i) {
      const Real e = u(0,CLC_X,ks,j,i), eo = u(0,CLC_O,ks,j,i);
      const Real y1 = u(0,CLC_X+1,ks,j,i), y2 = u(0,CLC_X+2,ks,j,i);
      const Real yo1 = u(0,CLC_O+1,ks,j,i), yo2 = u(0,CLC_O+2,ks,j,i);
      if (std::sqrt(y1*y1 + y2*y2) <= (1.0 + 1.0e-6)*e) continue;
      ++nbad;
      const Real d11 = cs(0,CLN_D11,ks,j,i), d22 = cs(0,CLN_D11+1,ks,j,i);
      const Real d12 = cs(0,CLN_D11+3,ks,j,i);
      const Real f1 = y1/e, f2 = y2/e;
      const Real g1 = (eo > 0.0) ? yo1/eo : 0.0, g2 = (eo > 0.0) ? yo2/eo : 0.0;
      const Real h1 = cs(0,CLN_HP,ks,j,i) + cs(0,CLN_HM,ks,j,i);
      const Real h2 = cs(0,CLN_HP+1,ks,j,i) + cs(0,CLN_HM+1,ks,j,i);
      const Real lo = lmin(d11 - g1*g1, d12 - g1*g2, d22 - g2*g2);
      const Real ln = lmin(d11 - f1*f1, d12 - f1*f2, d22 - f2*f2);
      const Real lr = lmin(d11 - h1*h1, d12 - h1*h2, d22 - h2*h2);
      if (lo < -1.0e-12) ++nold;
      if (ln < -1.0e-12) ++nnew;
      if (lr < -1.0e-12) ++nray;
      if (nprint < 25) {
        ++nprint;
        std::printf("CLDIAG i %d j %d E_old %.3e |f_old| %.3f E %.3e |f| %.3f "
                    "D %.3f %.3f %.3f C1 %.3f J/E %.3e f_ray %.3f %.3f "
                    "f_old %.3f %.3f f %.3f %.3f "
                    "lmin_old %.3e lmin_new %.3e lmin_ray %.3e\n", i - is, j - js, eo,
                    std::sqrt(g1*g1 + g2*g2), e, std::sqrt(f1*f1 + f2*f2), d11, d22, d12,
                    cs(0,CLN_C1,ks,j,i), cs(0,CLN_J,ks,j,i)/fmax(e, 1.0e-300), h1, h2,
                    g1, g2, f1, f2, lo, ln, lr);
      }
    }
  }
  std::printf("CLDIAG cycle %d t2s %d bad %d D-ff^T<0: old %d new %d ray %d\n",
              pmy_pack->pmesh->ncycle, t2_solve, nbad, nold, nnew, nray);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellWriteBack
//! \brief stores the solution: E and F of the cells, the face fluxes of E (f0x*), the gas
//! momentum and energy, and the stage slopes K = (U - U_old)/h.  The gas takes exactly
//! the momentum and energy the radiation gives up besides its face fluxes, so both
//! totals change only by the fluxes through the boundary and the external source esrc
//! (the radiation's; the slope of E carries it); with force_reference the reference force
//! is withheld from the gas momentum (an external source, like gravity), and the
//! radiation pays the work of the residual, exactly the kinetic energy the gas gains,
//! and with wb_arad (force_reference_work = split) also that of the reference force,
//! which the problem generator gives the gas.  cln_neint counts the cells left with
//! eint < 0 (the fluid's floor would create energy there).  On the cubed sphere F is stored in the frame of y (F_r, F_theta,
//! F_phi) and the gas momentum change (and its stage slope) in the hydro's components
//! (CleanGasChange): E and L_z of gas plus radiation are conserved to round-off.

void RadiationM1::CleanCellWriteBack() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool slopes = (t2_solve != M1_T2S_NONE);
  const bool feedback = cln_feedback;
  const Real h = dt_sub;
  CleanGhostFill(cln_x, CLC_X, 4);
  CleanArgs x = CleanMakeArgs();
  auto kk = (t2_solve == M1_T2S_STAGE1) ? t2k2 : t2k1;
  const int fref = force_ref;
  const bool vzero = cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;
  const bool eso = esrc_on;
  auto es = esrc;

  // face fluxes of E
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto f0 = (dir == 0) ? f0x1 : ((dir == 1) ? f0x2 : f0x3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_wbface", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real f[4];
      CleanCellFaceFlux(x, CLC_X, dir, m, k, j, i, true, f);
      f0(m,k,j,i) = f[0];
    });
  }

  // E, F, gas, slopes
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto u0_ = u0;
  Real neint = 0.0;
  Kokkos::parallel_reduce("clc_wbcell", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, ks, js, is}, {nmb1 + 1, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &nb) {
    Real div[4];
    CleanCellDivergence(x, CLC_X, m, k, j, i, true, div);
    Real e = x.u(m,CLC_X,k,j,i);
    const Real e_old = x.w(m,CLW_EO,k,j,i);
    Real dm[3] = {0.0, 0.0, 0.0}, dh[3];
    Real den = 0.0;
    if (feedback) {
      for (int b = 0; b < x.nf && !vzero; ++b) {
        dm[b] = -(x.u(m,CLC_X+1+b,k,j,i) - x.u(m,CLC_O+1+b,k,j,i) + div[1+b])/x.c
                - h*CleanRefForce(fref, fr, aref, uh(m,IDN,k,j,i), b, m, k, j, i);
      }
      // the radiation pays exactly the kinetic energy (p + dm/2).dm/rho the gas gains,
      // not the work W of the iterate the solve took
      // and, with wb_arad, the work of the reference force at the velocity after the
      // kick, which the problem generator's well-balanced source gives the gas
      Real pg[3], wk = 0.0, wr = 0.0;
      CleanGasMomentum(x, uh, m, k, j, i, pg);
      for (int b = 0; b < 3; ++b) {
        const Real d = uh(m,IDN,k,j,i);
        wk += (pg[b] + 0.5*dm[b])*dm[b]/d;
        if (fref == M1_FREF_WB_ARAD) {
          wr += h*CleanRefForce(fref, fr, aref, d, b, m, k, j, i)*(pg[b] + dm[b])/d;
        }
      }
      e -= wk - x.w(m,CLW_WK,k,j,i) + wr;
      // (the hydro's components; dm.e_a covariant on the cubed sphere)
      CleanGasChange(x, m, k, j, i, dm, dh);
      for (int b = 0; b < x.nf; ++b) uh(m,IM1+b,k,j,i) += dh[b];
      den = -(e - e_old) - div[0] + (eso ? h*es(m,k,j,i) : 0.0) - wr;
      if (vzero) {
        // steady start: the internal energy of the gas temperature of the iterate, not
        // the residual of a step h -> infinity
        Real ei, p, chi_rho, chi_t, cv;
        eos.ThermoAt(uh(m,IDN,k,j,i), x.w(m,CLW_T,k,j,i), ei, p, chi_rho, chi_t, cv);
        den = ei - x.w(m,CLW_EI,k,j,i);
      } else {
        // the gas eint changes by the exchange of the solve only
        const Real ke = 0.5*(SQR(pg[0]) + SQR(pg[1]) + SQR(pg[2]))/uh(m,IDN,k,j,i);
        if (uh(m,IEN,k,j,i) + den - ke - wk < 0.0) nb += 1.0;
      }
      uh(m,IEN,k,j,i) += den;
    }
    if (slopes) {
      kk(m,M1_T2_E,k,j,i) = (e - e_old)/h;
      for (int b = 0; b < 3; ++b) {
        kk(m,M1_T2_F1+b,k,j,i) = x.c*(x.u(m,CLC_X+1+b,k,j,i) - x.u(m,CLC_O+1+b,k,j,i))/h;
        kk(m,M1_T2_M1+b,k,j,i) = feedback ? dh[b]/h : 0.0;
      }
      kk(m,M1_T2_EN,k,j,i) = den/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    for (int b = 0; b < 3; ++b) u0_(m,M1_F1+b,k,j,i) = x.c*x.u(m,CLC_X+1+b,k,j,i);
  }, Kokkos::Sum<Real>(neint));
  cln_neint = static_cast<int>(CleanReduce(neint, false));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellClosureSlope
//! \brief adds sign times the closure-dependent part of the implicit rates of the
//! solution to the step's last slope K1: the flux divergence of E and F and, with gas
//! coupling, the advective source chi a E of F with a = v + D v (the gas takes the
//! opposite momentum).
//! Called with -1 before and +1 after a closure update, it re-evaluates K1 = f(U^{n+1})
//! with the new closure.

void RadiationM1::CleanCellClosureSlope(const Real sign) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real ih = 1.0/dt_sub;
  CleanArgs x = CleanMakeArgs();
  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool feedback = cln_feedback;
  auto uh = fl.u0;
  auto op = opac;
  auto kk = t2k1;
  par_for("clc_clslope", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real div[4];
    CleanCellDivergence(x, CLC_X, m, k, j, i, true, div);
    const Real e = x.u(m,CLC_X,k,j,i);
    kk(m,M1_T2_E,k,j,i) -= sign*div[0]*ih;
    // velocity of the gas after the write-back
    Real v[3] = {0.0, 0.0, 0.0}, a[3], dm[3], dh[3];
    if (feedback) {
      CleanGasMomentum(x, uh, m, k, j, i, v);
      for (int b = 0; b < 3; ++b) v[b] /= uh(m,IDN,k,j,i);
    }
    CleanAdvect(x, m, k, j, i, v, a);
    const Real chi = op(m,M1_OP_T,k,j,i);
    for (int b = 0; b < 3; ++b) dm[b] = (b < x.nf) ? -sign*chi*a[b]*e/x.c : 0.0;
    CleanGasChange(x, m, k, j, i, dm, dh);
    for (int b = 0; b < x.nf; ++b) {
      kk(m,M1_T2_F1+b,k,j,i) += sign*x.c*(chi*a[b]*e - div[1+b]*ih);
      kk(m,M1_T2_M1+b,k,j,i) += dh[b];
    }
  });
}

}  // namespace radm1
