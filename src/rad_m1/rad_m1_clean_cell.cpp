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

#include <cmath>
#include <cstdio>

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
//! \brief per face: the dissipation weight al, theta and the dissipation speed.  The
//! smoothness sensor psi (Jameson type, on E and F/c of the old vector, the larger of the
//! two cells) enters only clean_alpha = tau_sensor.

void RadiationM1::CleanCellFaceCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const int nf = three_d ? 3 : (multi_d ? 2 : 1);
  const Real h = dt_sub, c = c_light, eps = cln_psi_eps;
  const Real c1lo = cln_c0, c1hi = cln_c1;
  const int mode = cln_alpha;
  auto &size = pmy_pack->pmb->mb_size;
  const Real dx[3] = {size.h_view(0).dx1, size.h_view(0).dx2, size.h_view(0).dx3};
  auto cs = cln_s;
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

  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !multi_d) continue;
    if (dir == 2 && !three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const Real dxd = dx[dir];
    par_for("clc_facecoef", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int kl = k - dk, jl = j - dj, il = i - di;
      const Real chi = 0.5*(op(m,M1_OP_T,k,j,i) + op(m,M1_OP_T,kl,jl,il));
      const Real c1 = fmax(cs(m,CLN_C1,k,j,i), cs(m,CLN_C1,kl,jl,il));
      const Real sray = CleanSmoothStep((c1 - c1lo)/(c1hi - c1lo));
      const Real psi = fmax(u(m,CLC_PSI+dir,k,j,i), u(m,CLC_PSI+dir,kl,jl,il));
      const Real dd = fmax(cs(m,CLN_D11+dir,k,j,i), cs(m,CLN_D11+dir,kl,jl,il));
      fc(m,CLF_AL,k,j,i) = CleanDissipationWeight(mode, chi*dxd, psi, sray);
      fc(m,CLF_TH,k,j,i) = 1.0/(1.0 + h*c*chi);
      fc(m,CLF_AF,k,j,i) = c*std::sqrt(fmax(dd, 0.0));
      fc(m,CLF_LO,k,j,i) = 0.0;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellSlopeWeights
//! \brief lagged weights of the van Leer slopes of the Newton iterate:
//! pl = dr/(dl + dr), pr = dl/(dl + dr) where dl dr > 0, else 0, so that the half-slope
//! dl dr/(dl + dr) = (pl dl + pr dr)/2 is linear in the state for fixed weights

void RadiationM1::CleanCellSlopeWeights() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  CleanGhostFill(cln_x, CLC_K, 4);
  auto u = cln_x;

  par_for("clc_slopes", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int dir = 0; dir < 3; ++dir) {
      const bool on = (dir == 0 || (dir == 1 && multi_d) || (dir == 2 && three_d));
      const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
      for (int q = 0; q < 4; ++q) {
        Real pl = 0.0, pr = 0.0;
        if (on) {
          const Real dl = u(m,CLC_K+q,k,j,i) - u(m,CLC_K+q,k-dk,j-dj,i-di);
          const Real dr = u(m,CLC_K+q,k+dk,j+dj,i+di) - u(m,CLC_K+q,k,j,i);
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
//! \brief B_E = E_old + h c kP a T^4 - r h c kE E_k - W, B_y = y_old, each less h div of
//! the incoming boundary fluxes (exchange linearised about the iterate as in
//! CleanRightHandSide), and the initial guess X = the iterate

void RadiationM1::CleanCellRightHandSide() {
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

  const int vfreeze = cln_vfreeze;
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
      for (int n = 0; n < 3; ++n) {
        const Real dm = x.w(m,CLW_DM+n,k,j,i);
        work += (uh(m,IM1+n,k,j,i) + (vfreeze == 2 ? 0.0 : 0.5*dm))*dm/d;
      }
    }
    Real sdiv[4];
    CleanCellSourceDivergence(x, m, k, j, i, sdiv);
    x.w(m,CLW_RR,k,j,i) = rr;
    x.u(m,CLC_B,k,j,i) = x.w(m,CLW_EO,k,j,i) + hckp*arad_*SQR(SQR(t))
                         - rr*hcke*x.w(m,CLW_EK,k,j,i) - work - sdiv[0];
    for (int b = 1; b < 4; ++b) x.u(m,CLC_B+b,k,j,i) = x.u(m,CLC_O+b,k,j,i) - sdiv[b];
    for (int b = 0; b < 4; ++b) x.u(m,CLC_X+b,k,j,i) = x.u(m,CLC_K+b,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellApply
//! \brief u(dst) = A u(src) over the active cells: the local exchange and flux sources
//! plus h div of the face fluxes; flux components beyond the mesh dimension are identity

void RadiationM1::CleanCellApply(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanGhostFill(cln_x, src, 4);
  CleanArgs x = CleanMakeArgs();
  auto op = opac;

  par_for("clc_apply", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real div[4];
    CleanCellDivergence(x, src, m, k, j, i, false, div);
    const Real rr = x.feedback ? x.w(m,CLW_RR,k,j,i) : 0.0;
    const Real e = x.u(m,src,k,j,i);
    const Real chi = op(m,M1_OP_T,k,j,i);
    x.u(m,dst,k,j,i) = e*(1.0 + x.h*x.c*op(m,M1_OP_E,k,j,i)*(1.0 - rr)) + div[0];
    for (int b = 0; b < 3; ++b) {
      Real y = x.u(m,src+1+b,k,j,i);
      if (b < x.nf) {
        y = y*(1.0 + x.h*x.c*chi) + div[1+b];
        if (x.feedback) y -= x.h*chi*x.w(m,CLW_A1+b,k,j,i)*e;
      }
      x.u(m,dst+1+b,k,j,i) = y;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellPreconditionerSetup
//! \brief 4x4 blocks of the preconditioner, from A with first-order jumps: the diagonal
//! block of every cell (all directions) and the blocks of its two neighbours along each
//! direction (states of cells i-1, i, i+1).  The RC terms in cells i+-2 and the
//! reconstruction are left to the Krylov solve.

void RadiationM1::CleanCellPreconditionerSetup() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  CleanArgs x = CleanMakeArgs();
  x.recon = false;
  auto op = opac;

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
    for (int dir = 0; dir < 3; ++dir) {
      if ((dir == 1 && !x.multi_d) || (dir == 2 && !x.three_d)) continue;
      const int dk = (dir == 2), dj = (dir == 1), di = (dir == 0);
      const Real s = x.h*x.idx[dir];
      const int sm = CleanLineSlotM(dir), sp = CleanLineSlotP(dir);
      CleanCellFaceJacobian(x, dir, m, k+dk, j+dj, i+di, jlp, jrp);
      CleanCellFaceJacobian(x, dir, m, k, j, i, jlm, jrm);
      for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
          a0[r][c] += s*(jlp[r][c] - jrm[r][c]);
          x.u(m,sp+4*r+c,k,j,i) = s*jrp[r][c];
          x.u(m,sm+4*r+c,k,j,i) = -s*jlm[r][c];
        }
      }
    }
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) x.u(m,CLC_A0+4*r+c,k,j,i) = a0[r][c];
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellLineSolve
//! \brief u(dst) = L^{-1} u(src), L the block-tridiagonal part of the preconditioner
//! along direction dir (diagonal blocks and the neighbour blocks along dir): block Thomas
//! solve on every line

void RadiationM1::CleanCellLineSolve(const int dir, const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int sm = CleanLineSlotM(dir), sp = CleanLineSlotP(dir);
  // lines run along dir; (p, q) label them
  const int np = (dir == 2) ? nx[1] : nx[2];
  const int nq = (dir == 0) ? nx[1] : nx[0];
  const int nl = nx[dir];
  auto u = cln_x;

  par_for("clc_thomas", DevExeSpace(), 0, nmb1, 0, np-1, 0, nq-1,
  KOKKOS_LAMBDA(const int m, const int p, const int q) {
    auto cell = [&](const int n, int &k, int &j, int &i) {
      k = (dir == 2) ? ks + n : ks + p;
      j = (dir == 1) ? js + n : ((dir == 2) ? js + p : js + q);
      i = (dir == 0) ? is + n : is + q;
    };
    // forward: D = A0 - Am C_n-1, C_n = D^-1 Ap, d_n = D^-1 (r_n - Am d_n-1)
    int kp = 0, jp = 0, ip = 0;
    for (int n = 0; n < nl; ++n) {
      int k, j, i;
      cell(n, k, j, i);
      Real a[4][4], b[4][5];
      for (int r = 0; r < 4; ++r) {
        Real rhs = u(m,src+r,k,j,i);
        for (int c = 0; c < 4; ++c) {
          a[r][c] = u(m,CLC_A0+4*r+c,k,j,i);
          b[r][c] = (n < nl - 1) ? u(m,sp+4*r+c,k,j,i) : 0.0;
        }
        if (n > 0) {
          for (int l = 0; l < 4; ++l) {
            const Real aml = u(m,sm+4*r+l,k,j,i);
            if (aml == 0.0) continue;
            for (int c = 0; c < 4; ++c) a[r][c] -= aml*u(m,CLC_CP+4*l+c,kp,jp,ip);
            rhs -= aml*u(m,dst+l,kp,jp,ip);
          }
        }
        b[r][4] = rhs;
      }
      CleanBlockSolve(a, b, 5);
      for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) u(m,CLC_CP+4*r+c,k,j,i) = b[r][c];
        u(m,dst+r,k,j,i) = b[r][4];
      }
      kp = k;
      jp = j;
      ip = i;
    }
    // backward: x_n = d_n - C_n x_n+1
    for (int n = nl - 2; n >= 0; --n) {
      int k, j, i, kn, jn, in;
      cell(n, k, j, i);
      cell(n+1, kn, jn, in);
      Real d[4];
      for (int r = 0; r < 4; ++r) {
        d[r] = u(m,dst+r,k,j,i);
        for (int c = 0; c < 4; ++c) {
          d[r] -= u(m,CLC_CP+4*r+c,k,j,i)*u(m,dst+c,kn,jn,in);
        }
      }
      for (int r = 0; r < 4; ++r) u(m,dst+r,k,j,i) = d[r];
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellPrecondition
//! \brief u(dst) = M^{-1} u(src): x1-line solve, then line-solve corrections of the
//! residual along x2, x3 (multi-D), and x1, x2, x3 again in each further sweep
//! (clean_pc_sweeps; multiplicative line relaxation)

void RadiationM1::CleanCellPrecondition(const int src, const int dst) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto u = cln_x;
  CleanCellLineSolve(0, src, dst);
  for (int n = 1; n < 3*cln_pcsweeps; ++n) {
    const int dir = n % 3;
    if ((dir == 1 && !pmy_pack->pmesh->multi_d) ||
        (dir == 2 && !pmy_pack->pmesh->three_d)) continue;
    CleanCellApply(dst, CLC_T1);
    par_for("clc_pcres", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < 4; ++c) {
        u(m,CLC_T1+c,k,j,i) = u(m,src+c,k,j,i) - u(m,CLC_T1+c,k,j,i);
      }
    });
    CleanCellLineSolve(dir, CLC_T1, CLC_T2);
    par_for("clc_pccor", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int c = 0; c < 4; ++c) u(m,dst+c,k,j,i) += u(m,CLC_T2+c,k,j,i);
    });
  }
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
  return fmax(de_newton, dmax/fmax(emax, 1.0e-300));
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::CleanCellVelocityUpdate
//! \brief momentum change of the iterate, -[(y - y_old) + h div fY]/c per cell (minus the
//! radiation momentum change besides its flux divergence), and the advective
//! coefficients of the new velocity.  Returns max |dv| over max |v|.

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

  Real dvmax = 0.0, vmax = 0.0;
  Kokkos::parallel_reduce("clc_velocity", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real d = uh(m,IDN,k,j,i);
    Real div[4];
    CleanCellDivergence(x, CLC_K, m, k, j, i, true, div);
    Real v2 = 0.0;
    for (int b = 0; b < 3; ++b) {
      Real dm = 0.0;
      if (b < x.nf) {
        dm = -(x.u(m,CLC_K+1+b,k,j,i) - x.u(m,CLC_O+1+b,k,j,i) + div[1+b])/x.c;
      }
      dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+b,k,j,i))/d);
      x.w(m,CLW_DM+b,k,j,i) = dm;
      v2 += SQR((uh(m,IM1+b,k,j,i) + dm)/d);
    }
    v_max = fmax(v_max, std::sqrt(v2));
  }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax));
  if (cln_vfreeze != 1) CleanAdvectiveCoefficients();
  return dvmax/fmax(vmax, 1.0e-300);
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
  const Real zero = 10.0*cln_tol*bmax;
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
  cln_nnreal = static_cast<int>(nreal);
  return nneg;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::CleanCellMarkFaces
//! \brief marks every face of the cells whose solution has E < 0 or |F| > c E (beyond
//! the solver tolerance, 10 clean_tol max|B_E|) as low-order.  Returns the number of
//! such cells.

int RadiationM1::CleanCellMarkFaces() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto u = cln_x;
  auto fc1 = cln_fc1, fc2 = cln_fc2, fc3 = cln_fc3;

  Real bmax = 0.0;
  Kokkos::parallel_reduce("clc_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(u(m,CLC_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));
  const Real zero = 10.0*cln_tol*bmax;
  Real nbad = 0.0;
  Kokkos::parallel_reduce("clc_mark", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real e = u(m,CLC_X,k,j,i);
    const Real y = std::sqrt(SQR(u(m,CLC_X+1,k,j,i)) + SQR(u(m,CLC_X+2,k,j,i)) +
                             SQR(u(m,CLC_X+3,k,j,i)));
    if (e >= -zero && y <= (1.0 + 1.0e-12)*e + zero) return;
    n += 1.0;
    fc1(m,CLF_LO,k,j,i) = 1.0;
    fc1(m,CLF_LO,k,j,i+1) = 1.0;
    if (multi_d) {
      fc2(m,CLF_LO,k,j,i) = 1.0;
      fc2(m,CLF_LO,k,j+1,i) = 1.0;
    }
    if (three_d) {
      fc3(m,CLF_LO,k,j,i) = 1.0;
      fc3(m,CLF_LO,k+1,j,i) = 1.0;
    }
  }, Kokkos::Sum<Real>(nbad));
  return static_cast<int>(nbad);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CleanCellAllLowOrder
//! \brief marks every face low-order

void RadiationM1::CleanCellAllLowOrder() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  for (int dir = 0; dir < 3; ++dir) {
    if (dir == 1 && !pmy_pack->pmesh->multi_d) continue;
    if (dir == 2 && !pmy_pack->pmesh->three_d) continue;
    auto fc = (dir == 0) ? cln_fc1 : ((dir == 1) ? cln_fc2 : cln_fc3);
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    par_for("clc_alllow", DevExeSpace(), 0, nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fc(m,CLF_LO,k,j,i) = 1.0;
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
//! totals change only by the fluxes through the boundary.

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
  auto u0_ = u0;
  par_for("clc_wbcell", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real div[4];
    CleanCellDivergence(x, CLC_X, m, k, j, i, true, div);
    const Real e = x.u(m,CLC_X,k,j,i);
    const Real e_old = x.w(m,CLW_EO,k,j,i);
    Real dm[3] = {0.0, 0.0, 0.0};
    Real den = 0.0;
    if (feedback) {
      for (int b = 0; b < x.nf; ++b) {
        dm[b] = -(x.u(m,CLC_X+1+b,k,j,i) - x.u(m,CLC_O+1+b,k,j,i) + div[1+b])/x.c;
        uh(m,IM1+b,k,j,i) += dm[b];
      }
      den = -(e - e_old) - div[0];
      uh(m,IEN,k,j,i) += den;
    }
    if (slopes) {
      kk(m,M1_T2_E,k,j,i) = (e - e_old)/h;
      for (int b = 0; b < 3; ++b) {
        kk(m,M1_T2_F1+b,k,j,i) = x.c*(x.u(m,CLC_X+1+b,k,j,i) - x.u(m,CLC_O+1+b,k,j,i))/h;
        kk(m,M1_T2_M1+b,k,j,i) = dm[b]/h;
      }
      kk(m,M1_T2_EN,k,j,i) = den/h;
    }
    u0_(m,M1_E,k,j,i) = e;
    for (int b = 0; b < 3; ++b) u0_(m,M1_F1+b,k,j,i) = x.c*x.u(m,CLC_X+1+b,k,j,i);
  });
}

}  // namespace radm1
