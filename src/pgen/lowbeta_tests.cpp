//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file lowbeta_tests.cpp
//! \brief lowbeta-tests-1102: TEST-ONLY initial conditions for the low-plasma-beta
//! Cartesian MHD robustness study (no numerical method is changed by this file).
//!
//! iprob = 1  low-beta rotor (Zhang, Lani & Poedts 2026, sect. 6.5): [0,1]^2, ambient
//!            rho = rho_amb (0.1), rotor rho = rho_rot (1) for r < r0 = 0.1, linear taper
//!            to r1 = 0.115, rigid rotation v = v0 (-(y-0.5), (x-0.5))/r0 inside, the
//!            same times f in the taper; p = 0.05; B = (b0x, b0y) = 2.5/sqrt(4 pi) (1,1);
//!            gamma = 5/3; t = 0.12.
//! iprob = 2  magnetised blast: rho = 1, p = p_in for r < r_blast, p_amb outside, uniform
//!            B = b0 (sin th cos ph, sin th sin ph, cos th) (Balsara & Spicer 1999 /
//!            Mignone & Del Zanna 2021 sect. 6.4; Tremblin et al. 2024 sect. 8.2.3-4).
//!            Centre (xc, yc, zc).  2-D or 3-D.
//! iprob = 3  circularly polarised Alfven wave (Toth 2000; Gardiner & Stone 2008 sect.
//!            6.3) in 1-D or 2-D.  The built-in pgen cpaw.cpp is 3-D only: its vector
//!            potential arrays have one cell in x3 (and x2) in 2-D (1-D) but are indexed
//!            at k+1 (j+1), so 1-D/2-D runs start from a garbage field.  Same formulas,
//!            same wave vector (along the grid diagonal of a 2:1 box in 2-D, along x1 in
//!            1-D), rho = 1, B_par = b_par, B_perp = v_perp = b_perp, uniform p = pres.
//!            The exact solution after a whole number of periods is the initial state.
//!
//! <problem>/user_hist = true adds the history columns
//!   pmin   minimum gas pressure over the mesh (after the EOS floors)
//!   emin   minimum of e_int = E - rho v^2/2 - B^2/2 recomputed from u0 and bcc0
//!   bmin   minimum plasma beta 2p/B^2
//!   dmin   minimum density
//! The mins are taken as (-max(-x)) so that the history's MPI_SUM over ranks of a single
//! rank's value is right only in serial; this test binary is serial.

#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "pgen/pgen.hpp"

namespace {
void LowBetaHistory(HistoryData *pdata, Mesh *pm);

struct CpawPar {
  Real pol, b_par, b_perp, k_par, sin_a3, cos_a3;
};
// vector potential of the 2-D CPAW (cpaw.cpp A1..A3 with sin_a2 = 0); A_x = -ay sin_a3,
// A_y = ay cos_a3, A_z = az, all independent of x3
KOKKOS_INLINE_FUNCTION
void CpawA(const Real x1, const Real x2, const CpawPar &c, Real &a1, Real &a2, Real &a3) {
  Real x = x1*c.cos_a3 + x2*c.sin_a3;
  Real y = -x1*c.sin_a3 + x2*c.cos_a3;
  Real ay = c.pol*(c.b_perp/c.k_par)*sin(c.k_par*x);
  Real az = (c.b_perp/c.k_par)*cos(c.k_par*x) + c.b_par*y;
  a1 = -ay*c.sin_a3;
  a2 = ay*c.cos_a3;
  a3 = az;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void CpawInit()
//! \brief CPAW for 1-D / 2-D meshes (iprob = 3), face fields as the discrete curl of A

void CpawInit(ParameterInput *pin, Mesh *pmesh) {
  MeshBlockPack *pmbp = pmesh->pmb_pack;
  const Real b_par = pin->GetReal("problem", "b_par");
  const Real b_perp = pin->GetReal("problem", "b_perp");
  const Real pres = pin->GetReal("problem", "pres");
  const Real pol = pin->GetOrAddBoolean("problem", "right_polar", true) ? 1.0 : -1.0;
  if (pmesh->three_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "lowbeta_tests iprob=3 is the 1-D/2-D CPAW; use cpaw.cpp in 3-D" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  Real x1size = pmesh->mesh_size.x1max - pmesh->mesh_size.x1min;
  Real x2size = pmesh->mesh_size.x2max - pmesh->mesh_size.x2min;
  CpawPar c;
  c.pol = pol; c.b_par = b_par; c.b_perp = b_perp;
  c.cos_a3 = 1.0; c.sin_a3 = 0.0;
  if (pmesh->multi_d) {
    Real ang3 = std::atan(x1size/x2size);
    c.sin_a3 = std::sin(ang3);
    c.cos_a3 = std::cos(ang3);
  }
  Real lambda = x1size*c.cos_a3;
  if (c.sin_a3 > 0.0) lambda = std::min(lambda, x2size*c.sin_a3);
  c.k_par = 2.0*M_PI/lambda;
  const Real v_perp = b_perp;   // rho = 1
  if (global_variable::my_rank == 0) {
    std::cout << "lowbeta_tests CPAW: lambda = " << lambda << " period = " << lambda/b_par
              << " beta = " << pres/(0.5*(b_par*b_par + b_perp*b_perp)) << std::endl;
  }
  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  const Real gm1 = eos.gamma - 1.0;
  auto &indcs = pmesh->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mb_size;
  par_for("pgen_cpaw2d", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real x1f = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real x2f = LeftEdgeX(j-js, indcs.nx2, x2min, x2max);
    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real t1, t2, t3;
    // A3 at the four x1x2 corners, A2 on x1 faces at x2v, A1 on x2 faces at x1v
    Real a3ll, a3rl, a3lr, a3rr, a2l, a2r, a1l, a1r;
    CpawA(x1f, x2f, c, t1, t2, a3ll);
    CpawA(x1f+dx1, x2f, c, t1, t2, a3rl);
    CpawA(x1f, x2f+dx2, c, t1, t2, a3lr);
    CpawA(x1f+dx1, x2f+dx2, c, t1, t2, a3rr);
    CpawA(x1f, x2v, c, t1, a2l, t3);
    CpawA(x1f+dx1, x2v, c, t1, a2r, t3);
    CpawA(x1v, x2f, c, a1l, t2, t3);
    CpawA(x1v, x2f+dx2, c, a1r, t2, t3);
    Real b1l = (a3lr - a3ll)/dx2, b1r = (a3rr - a3rl)/dx2;
    Real b2l = -(a3rl - a3ll)/dx1, b2r = -(a3rr - a3lr)/dx1;
    Real b3 = (a2r - a2l)/dx1 - (a1r - a1l)/dx2;
    b0.x1f(m,k,j,i) = b1l;
    b0.x2f(m,k,j,i) = b2l;
    b0.x3f(m,k,j,i) = b3;
    if (i == ie) b0.x1f(m,k,j,i+1) = b1r;
    if (j == je) b0.x2f(m,k,j+1,i) = b2r;
    if (k == ke) b0.x3f(m,k+1,j,i) = b3;

    Real x = x1v*c.cos_a3 + x2v*c.sin_a3;
    Real sn = sin(c.k_par*x);
    Real cs = c.pol*cos(c.k_par*x);
    Real my = -c.pol*v_perp*sn;
    Real mz = -c.pol*v_perp*cs;
    Real m1 = -my*c.sin_a3;
    Real m2 = my*c.cos_a3;
    u0(m,IDN,k,j,i) = 1.0;
    u0(m,IM1,k,j,i) = m1;
    u0(m,IM2,k,j,i) = m2;
    u0(m,IM3,k,j,i) = mz;
    Real bc1 = 0.5*(b1l + b1r), bc2 = 0.5*(b2l + b2r);
    u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(bc1*bc1 + bc2*bc2 + b3*b3)
                    + 0.5*(m1*m1 + m2*m2 + mz*mz);
  });
}

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (pin->GetOrAddBoolean("problem", "user_hist", false)) {
    user_hist_func = LowBetaHistory;
  }
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "lowbeta_tests needs an <mhd> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const int iprob = pin->GetInteger("problem", "iprob");
  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  const Real gm1 = eos.gamma - 1.0;

  // parameters
  Real rho_amb = 1.0, rho_rot = 1.0, r0 = 0.1, r1 = 0.115, v0 = 1.0, pres = 1.0;
  Real p_in = 1.0, p_amb = 1.0, rb = 0.1, xc = 0.0, yc = 0.0, zc = 0.0;
  Real bx = 0.0, by = 0.0, bz = 0.0;
  if (iprob == 3) {
    CpawInit(pin, pmy_mesh_);
    return;
  }
  if (iprob == 1) {
    rho_amb = pin->GetOrAddReal("problem", "rho_amb", 0.1);
    rho_rot = pin->GetOrAddReal("problem", "rho_rot", 1.0);
    r0 = pin->GetOrAddReal("problem", "r0", 0.1);
    r1 = pin->GetOrAddReal("problem", "r1", 0.115);
    v0 = pin->GetOrAddReal("problem", "v0", 1.0);
    pres = pin->GetOrAddReal("problem", "pres", 0.05);
    bx = pin->GetOrAddReal("problem", "b0x", 2.5/std::sqrt(4.0*M_PI));
    by = pin->GetOrAddReal("problem", "b0y", 2.5/std::sqrt(4.0*M_PI));
    bz = 0.0;
    xc = pin->GetOrAddReal("problem", "xc", 0.5);
    yc = pin->GetOrAddReal("problem", "yc", 0.5);
  } else if (iprob == 2) {
    p_in = pin->GetReal("problem", "p_in");
    p_amb = pin->GetReal("problem", "p_amb");
    rb = pin->GetOrAddReal("problem", "r_blast", 0.1);
    Real b0 = pin->GetReal("problem", "b0");
    Real th = pin->GetOrAddReal("problem", "theta", 0.5*M_PI);
    Real ph = pin->GetOrAddReal("problem", "phi", 0.25*M_PI);
    bx = b0*std::sin(th)*std::cos(ph);
    by = b0*std::sin(th)*std::sin(ph);
    bz = b0*std::cos(th);
    xc = pin->GetOrAddReal("problem", "xc", 0.0);
    yc = pin->GetOrAddReal("problem", "yc", 0.0);
    zc = pin->GetOrAddReal("problem", "zc", 0.0);
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "lowbeta_tests: iprob = 1 (rotor) or 2 (blast)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (global_variable::my_rank == 0) {
    Real pb = 0.5*(bx*bx + by*by + bz*bz);
    Real pa = (iprob == 1) ? pres : p_amb;
    std::cout << "lowbeta_tests: iprob=" << iprob << " ambient beta = "
              << ((pb > 0.0) ? pa/pb : 0.0) << std::endl;
  }

  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mb_size;
  const bool three_d = pmy_mesh_->three_d;

  par_for("pgen_lowbeta", DevExeSpace(), 0, (pmbp->nmb_thispack-1), ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);

    Real d = 1.0, vx = 0.0, vy = 0.0, p = 1.0;
    if (iprob == 1) {
      Real dx = x1v - xc, dy = x2v - yc;
      Real r = std::sqrt(dx*dx + dy*dy);
      p = pres;
      if (r < r0) {
        d = rho_rot;
        vx = -v0*dy/r0;
        vy = v0*dx/r0;
      } else if (r < r1) {
        Real f = (r1 - r)/(r1 - r0);
        d = rho_amb + (rho_rot - rho_amb)*f;
        vx = -f*v0*dy/r0;
        vy = f*v0*dx/r0;
      } else {
        d = rho_amb;
      }
    } else {
      Real dx = x1v - xc, dy = x2v - yc, dz = three_d ? (x3v - zc) : 0.0;
      Real r = std::sqrt(dx*dx + dy*dy + dz*dz);
      d = 1.0;
      p = (r < rb) ? p_in : p_amb;
    }
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*vx;
    u0(m,IM2,k,j,i) = d*vy;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = p/gm1 + 0.5*d*(vx*vx + vy*vy) + 0.5*(bx*bx + by*by + bz*bz);

    // uniform face field: divergence-free exactly
    b0.x1f(m,k,j,i) = bx;
    b0.x2f(m,k,j,i) = by;
    b0.x3f(m,k,j,i) = bz;
    if (i == ie) b0.x1f(m,k,j,i+1) = bx;
    if (j == je) b0.x2f(m,k,j+1,i) = by;
    if (k == ke) b0.x3f(m,k+1,j,i) = bz;
  });
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void LowBetaHistory()
//! \brief minimum p, e_int (from u0), beta and density over the mesh

void LowBetaHistory(HistoryData *pdata, Mesh *pm) {
  pdata->nhist = 4;
  pdata->label[0] = "pmin";
  pdata->label[1] = "emin";
  pdata->label[2] = "bmin";
  pdata->label[3] = "dmin";
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &w0 = pmbp->pmhd->w0;
  auto &u0 = pmbp->pmhd->u0;
  auto &bcc = pmbp->pmhd->bcc0;
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  Real mp = -FLT_MAX, me = -FLT_MAX, mb = -FLT_MAX, md = -FLT_MAX;
  Kokkos::parallel_reduce("lowbeta_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &a, Real &b, Real &c, Real &e) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real p = gm1*w0(m,IEN,k,j,i);
    Real b2 = bcc(m,IBX,k,j,i)*bcc(m,IBX,k,j,i) + bcc(m,IBY,k,j,i)*bcc(m,IBY,k,j,i)
            + bcc(m,IBZ,k,j,i)*bcc(m,IBZ,k,j,i);
    Real d = u0(m,IDN,k,j,i);
    Real ek = 0.5*(u0(m,IM1,k,j,i)*u0(m,IM1,k,j,i) + u0(m,IM2,k,j,i)*u0(m,IM2,k,j,i)
                 + u0(m,IM3,k,j,i)*u0(m,IM3,k,j,i))/d;
    Real ei = u0(m,IEN,k,j,i) - ek - 0.5*b2;
    a = fmax(a, -p);
    b = fmax(b, -ei);
    c = fmax(c, -((b2 > 0.0) ? 2.0*p/b2 : FLT_MAX));
    e = fmax(e, -w0(m,IDN,k,j,i));
  }, Kokkos::Max<Real>(mp), Kokkos::Max<Real>(me), Kokkos::Max<Real>(mb),
     Kokkos::Max<Real>(md));
  pdata->hdata[0] = -mp;
  pdata->hdata[1] = -me;
  pdata->hdata[2] = -mb;
  pdata->hdata[3] = -md;
}
}  // namespace
