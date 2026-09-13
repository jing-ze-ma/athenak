//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file balsara_vortex.cpp
//! \brief Advected magnetic vortex of Balsara (2004), in the form used by Leidi et al.
//! (2022, Sect. 5.2) as the low-Mach-number benchmark for the LHLLD Riemann solver.
//!
//! The vortex is an exact stationary solution of the ideal 2D homogeneous MHD equations:
//! the centrifugal acceleration, the magnetic tension and the gas + magnetic pressure
//! gradients balance exactly.  Superposing a uniform advection velocity leaves it an
//! exact solution, and after one advective crossing time the vortex has returned to its
//! starting point, so THE EXACT SOLUTION AT t_adv IS THE INITIAL CONDITION and the L1
//! errors measure numerical dissipation alone.
//!
//! Initial conditions, quoted from Leidi et al. (2022) Eq. (53), on the periodic domain
//! (x,y) in [-5,5] x [-5,5], with r^2 = x^2 + y^2:
//!
//!   (Vx, Vy) = Vtil * exp((1 - r^2)/2) * (-y, x)
//!   (Bx, By) = Btil * exp((1 - r^2)/2) * (-y, x)
//!   p        = 1 + [ (Btil^2/2)*(1 - r^2) - Vtil^2/2 ] * exp(1 - r^2)
//!   rho      = 1
//!
//! Vtil is the maximum rotational velocity of the vortex (reached at r = 1) and Btil
//! sets the maximum Alfven speed on the grid.  The two amplitudes are parametrised by
//! the vortex strength Vtil and
//!
//!   beta_K = Btil^2 / Vtil^2,
//!
//! the ratio of the MAGNETIC to the (rotational) KINETIC energy of the vortex, which is
//! constant across the domain.  With gamma = 5/3 the maximum sonic Mach number at t = 0
//! is 1.55*Vtil, since the peak speed |V_adv| + Vtil = 2*Vtil sits on a background sound
//! speed sqrt(gamma) = 1.291.
//!
//! To make the problem harder the vortex is advected along the grid diagonal with
//! |V_adv| = Vtil, i.e. V_adv = (Vtil/sqrt(2), Vtil/sqrt(2)), for one advective crossing
//!
//!   t_adv = 10*sqrt(2)/Vtil,
//!
//! during which it rotates 2.25 times.  The field is initialised from the vector
//! potential
//!
//!   A_z = Btil * exp((1 - r^2)/2),
//!
//! because dA_z/dy = -Btil*y*exp((1-r^2)/2) = Bx and -dA_z/dx = Btil*x*exp((1-r^2)/2)
//! = By.  Taking the discrete curl of A_z on cell edges makes the face-centred field
//! divergence-free to round-off.
//!
//! At the end of the run the L1 error of every primitive variable is computed following
//! Leidi et al. (2022) Eq. (55),
//!
//!   L1(w_k) = (1/N^2) * sum_ij | w_k,ij(t_adv) - w_k,ij(0) |,
//!
//! along with the volume integrals of the rotational kinetic energy, Eq. (56),
//!
//!   E_R = 0.5*rho*[ (Vx - Vtil/sqrt(2))^2 + (Vy - Vtil/sqrt(2))^2 ],
//!
//! and of the magnetic energy 0.5*B^2.  Both are normalised by their value in the
//! DISCRETE initial state (recomputed here with exactly the expressions used to
//! initialise the run), so the dissipated fractions are reported directly and the
//! discretisation error of the initial data cancels.
//!
//! REFERENCES:
//! - D. S. Balsara, "Second-order-accurate schemes for magnetohydrodynamics with
//!   divergence-free reconstruction", ApJS, 151, 149 (2004)
//! - G. Leidi et al., "A finite-volume scheme for modeling compressible
//!   magnetohydrodynamic flows at low Mach numbers in stellar interiors", A&A, 668, A143
//!   (2022), Sect. 5.2
//! - T. Minoshima & T. Miyoshi, "A low-dissipation HLLD approximate Riemann solver for a
//!   very wide range of Mach numbers", JCP, 446, 110639 (2021)

// C++ headers
#include <math.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

namespace {

// vortex parameters, shared between the generator and the error routine below
Real vtil_, btil_, vadv_;

//----------------------------------------------------------------------------------------
//! \fn Real VortexAz(x, y, bt)
//! \brief the vector potential A_z = Btil*exp((1-r^2)/2) of the Balsara vortex

KOKKOS_INLINE_FUNCTION
Real VortexAz(const Real x, const Real y, const Real bt) {
  return bt*exp(0.5*(1.0 - x*x - y*y));
}

void BalsaraVortexErrors(ParameterInput *pin, Mesh *pm);

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::BalsaraVortex(ParameterInput *pin, const bool restart)
//! \brief Problem generator for the advected magnetic vortex.  Assumes a periodic 2D
//! domain centred on the origin (the standard choice is [-5,5] x [-5,5]).

void ProblemGenerator::BalsaraVortex(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Balsara vortex test can only be run in MHD, but no <mhd> block "
              << "in input file" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(pmy_mesh_->two_d)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Balsara vortex test requires a 2D grid" << std::endl;
    exit(EXIT_FAILURE);
  }

  // vortex strength and the magnetic-to-kinetic energy ratio set the two amplitudes
  vtil_ = pin->GetOrAddReal("problem", "vtil", 1.0e-2);
  Real beta_k = pin->GetOrAddReal("problem", "beta_k", 1.0);
  btil_ = std::sqrt(beta_k)*vtil_;
  vadv_ = vtil_/std::sqrt(2.0);

  // one advective crossing of the diagonal returns the vortex to its starting point
  Real ncross = pin->GetOrAddReal("problem", "ncross", 1.0);
  Real lx = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  pin->SetReal("time", "tlim", ncross*lx*std::sqrt(2.0)/vtil_);

  // errors are computed against the initial condition at the end of the run
  pgen_final_func = BalsaraVortexErrors;
  if (restart) return;

  EOS_Data &eos = pmbp->pmhd->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  auto &u0 = pmbp->pmhd->u0;
  auto &b0 = pmbp->pmhd->b0;
  auto &size = pmbp->pmb->mb_size;
  Real vt = vtil_, bt = btil_, va = vadv_;

  par_for("pgen_balsara", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1f = LeftEdgeX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2f = LeftEdgeX(j-js, nx2, x2min, x2max);

    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;

    // face-centred field as the discrete curl of A_z, so that div(B) = 0 to round-off
    Real a_ll = VortexAz(x1f,     x2f,     bt);
    Real a_rl = VortexAz(x1f+dx1, x2f,     bt);
    Real a_lr = VortexAz(x1f,     x2f+dx2, bt);
    Real a_rr = VortexAz(x1f+dx1, x2f+dx2, bt);
    Real b1l = (a_lr - a_ll)/dx2;
    Real b1r = (a_rr - a_rl)/dx2;
    Real b2l = -(a_rl - a_ll)/dx1;
    Real b2r = -(a_rr - a_lr)/dx1;

    b0.x1f(m,k,j,i) = b1l;
    b0.x2f(m,k,j,i) = b2l;
    b0.x3f(m,k,j,i) = 0.0;
    if (i==ie) {
      b0.x1f(m,k,j,i+1) = b1r;
    }
    if (j==je) {
      b0.x2f(m,k,j+1,i) = b2r;
    }
    if (k==ke) {
      b0.x3f(m,k+1,j,i) = 0.0;
    }

    // cell-centred state: rho = 1, so the momenta equal the velocities
    Real rsq = x1v*x1v + x2v*x2v;
    Real ex = exp(0.5*(1.0 - rsq));
    Real vx = va - vt*ex*x2v;
    Real vy = va + vt*ex*x1v;
    Real pgas = 1.0 + (0.5*bt*bt*(1.0 - rsq) - 0.5*vt*vt)*exp(1.0 - rsq);
    Real bccx = 0.5*(b1l + b1r);
    Real bccy = 0.5*(b2l + b2r);

    u0(m,IDN,k,j,i) = 1.0;
    u0(m,IM1,k,j,i) = vx;
    u0(m,IM2,k,j,i) = vy;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = pgas/gm1 + 0.5*(vx*vx + vy*vy) + 0.5*(bccx*bccx + bccy*bccy);
  });

  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void BalsaraVortexErrors(ParameterInput *pin, Mesh *pm)
//! \brief Compare the final state to the initial condition (which is the exact solution
//! after an integer number of advective crossings) and write the L1 error of each
//! primitive variable, plus the surviving fractions of the rotational kinetic and the
//! magnetic energy of the vortex, to "<basename>-vortex.dat".

void BalsaraVortexErrors(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;

  // capture class variables for kernel
  auto &indcs = pm->mb_indcs;
  int &nx1 = indcs.nx1;
  int &nx2 = indcs.nx2;
  int &nx3 = indcs.nx3;
  int &is = indcs.is;
  int &js = indcs.js;
  int &ks = indcs.ks;
  auto &size = pmbp->pmb->mb_size;
  auto &w0_ = pmbp->pmhd->w0;
  auto &bcc0_ = pmbp->pmhd->bcc0;
  Real vt = vtil_, bt = btil_, va = vadv_;
  // w0(IEN) holds the internal energy density, so p = (gamma-1)*w0(IEN)
  Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;

  const int nmkji = (pmbp->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("balsara-err", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    // compute m,k,j,i indices of thread
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
    Real x1f = LeftEdgeX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
    Real x2f = LeftEdgeX(j-js, nx2, x2min, x2max);

    Real dx1 = size.d_view(m).dx1;
    Real dx2 = size.d_view(m).dx2;
    Real vol = dx1*dx2*size.d_view(m).dx3;

    // rebuild the discrete initial state with the expressions used to initialise it
    Real a_ll = VortexAz(x1f,     x2f,     bt);
    Real a_rl = VortexAz(x1f+dx1, x2f,     bt);
    Real a_lr = VortexAz(x1f,     x2f+dx2, bt);
    Real a_rr = VortexAz(x1f+dx1, x2f+dx2, bt);
    Real bx_ex = 0.5*((a_lr - a_ll) + (a_rr - a_rl))/dx2;
    Real by_ex = -0.5*((a_rl - a_ll) + (a_rr - a_lr))/dx1;

    Real rsq = x1v*x1v + x2v*x2v;
    Real ex = exp(0.5*(1.0 - rsq));
    Real vx_ex = va - vt*ex*x2v;
    Real vy_ex = va + vt*ex*x1v;
    Real p_ex = 1.0 + (0.5*bt*bt*(1.0 - rsq) - 0.5*vt*vt)*exp(1.0 - rsq);

    array_sum::GlobalSum evars;
    evars.the_array[0] = fabs(w0_(m,IDN,k,j,i) - 1.0);
    evars.the_array[1] = fabs(w0_(m,IVX,k,j,i) - vx_ex);
    evars.the_array[2] = fabs(w0_(m,IVY,k,j,i) - vy_ex);
    evars.the_array[3] = fabs(w0_(m,IVZ,k,j,i));
    evars.the_array[4] = fabs(gm1*w0_(m,IEN,k,j,i) - p_ex);
    evars.the_array[5] = fabs(bcc0_(m,IBX,k,j,i) - bx_ex);
    evars.the_array[6] = fabs(bcc0_(m,IBY,k,j,i) - by_ex);
    evars.the_array[7] = fabs(bcc0_(m,IBZ,k,j,i));

    // rotational kinetic energy (Leidi Eq. 56) and magnetic energy, final and initial
    Real dvx = w0_(m,IVX,k,j,i) - va;
    Real dvy = w0_(m,IVY,k,j,i) - va;
    evars.the_array[8] = vol*0.5*w0_(m,IDN,k,j,i)*(dvx*dvx + dvy*dvy
                         + SQR(w0_(m,IVZ,k,j,i)));
    evars.the_array[9] = vol*0.5*(SQR(vx_ex - va) + SQR(vy_ex - va));
    evars.the_array[10] = vol*0.5*(SQR(bcc0_(m,IBX,k,j,i)) + SQR(bcc0_(m,IBY,k,j,i))
                          + SQR(bcc0_(m,IBZ,k,j,i)));
    evars.the_array[11] = vol*0.5*(bx_ex*bx_ex + by_ex*by_ex);

    for (int n=12; n<NREDUCTION_VARIABLES; ++n) {
      evars.the_array[n] = 0.0;
    }
    mb_sum += evars;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));

  constexpr int kNsum = 12;
  Real gsum[kNsum];
  for (int n=0; n<kNsum; ++n) {
    gsum[n] = sum_this_mb.the_array[n];
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, gsum, kNsum, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  // L1 errors are averages over cells (Leidi Eq. 55)
  Real ncells = static_cast<Real>(pm->nmb_total)*static_cast<Real>(nx1*nx2*nx3);
  for (int n=0; n<8; ++n) {
    gsum[n] /= ncells;
  }
  Real l1_v = std::sqrt(SQR(gsum[1]) + SQR(gsum[2]) + SQR(gsum[3]));
  Real l1_b = std::sqrt(SQR(gsum[5]) + SQR(gsum[6]) + SQR(gsum[7]));
  Real ek_frac = (gsum[9] > 0.0) ? gsum[8]/gsum[9] : 0.0;
  Real em_frac = (gsum[11] > 0.0) ? gsum[10]/gsum[11] : 0.0;

  if (global_variable::my_rank == 0) {
    std::string fname;
    fname.assign(pin->GetString("job", "basename"));
    fname.append("-vortex.dat");
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(), "r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(), "a", pfile)) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else {
      if ((pfile = std::fopen(fname.c_str(), "w")) == nullptr) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "Error output file could not be opened" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::fprintf(pfile, "# Nx1  Nx2   vtil          beta_K        d_L1          ");
      std::fprintf(pfile, "v1_L1         v2_L1         v3_L1         p_L1          ");
      std::fprintf(pfile, "b1_L1         b2_L1         b3_L1         v_L1          ");
      std::fprintf(pfile, "b_L1          Ekin/Ekin0    Emag/Emag0\n");
    }
    Real beta_k = (vtil_ > 0.0) ? SQR(btil_/vtil_) : 0.0;
    std::fprintf(pfile, "%04d  %04d", pm->mesh_indcs.nx1, pm->mesh_indcs.nx2);
    std::fprintf(pfile, "  %e  %e", vtil_, beta_k);
    for (int n=0; n<8; ++n) {
      std::fprintf(pfile, "  %e", gsum[n]);
    }
    std::fprintf(pfile, "  %e  %e  %e  %e\n", l1_v, l1_b, ek_frac, em_frac);
    std::fclose(pfile);

    std::cout << std::endl << "Balsara vortex: vtil=" << vtil_ << " beta_K=" << beta_k
              << std::endl;
    std::cout << "  L1(d)=" << gsum[0] << " L1(v)=" << l1_v << " L1(p)=" << gsum[4]
              << " L1(b)=" << l1_b << std::endl;
    std::cout << "  Ekin/Ekin0=" << ek_frac << " Emag/Emag0=" << em_frac << std::endl;
  }

  return;
}

}  // namespace
