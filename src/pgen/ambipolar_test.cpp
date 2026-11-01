//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ambipolar_test.cpp
//! \brief Cartesian tests of ambipolar diffusion (ambipolar-1101).
//! Build with -D PROBLEM=ambipolar_test.
//!
//! <problem>/iprob:
//!   1  force-free field, B parallel to J.  1-D (nx2 = 1): the helical field
//!      B = b0 (0, sin kx, cos kx), for which the discrete J is parallel to the averaged
//!      edge B, so the AD EMF must vanish to round-off.  3-D: the ABC (Beltrami) field
//!      B = b0 (sin kz + cos ky, sin kx + cos kz, sin ky + cos kx), J = k B, for which
//!      the AD EMF vanishes to truncation order.  <problem>/pert adds a uniform-direction
//!      perturbation pert*b0 (1,1,1)/sqrt3 so that J x B != 0 (energy test).
//!      The final function prints max|E_AD| / max|eta_A J| over the active edges.
//!   2  nonlinear AD diffusion of B_y(x) (B_x = B_z = 0, so J.B = 0 and eta_A J is the
//!      whole AD EMF): dB/dt = eta_ad d/dx (B^2 dB/dx) = (eta_ad/3) d2(B^3)/dx2, whose
//!      Barenblatt solution (m = 3, d = 1) at s = eta_ad t/3 + s0 is
//!        B = s^(-1/4) [C - x^2 s^(-1/2)/12]_+^(1/2).
//!      Heavy gas (dens) in total-pressure balance keeps v ~ 0.  The final function
//!      prints
//!      the L1 error of B_y against the exact profile (whole line and |x| < x_front/2;
//!      the front, where B ~ sqrt(distance), limits the whole-line order to ~1).
//!   3  the cs_test iprob = 11 field B = b0c (-y, x, 0) + (bvx, bvy, bvz), J = 2 b0c zhat
//!      uniform, J.B = 2 b0c bvz != 0: the final function compares the code's AD edge EMF
//!      with the exact eta_ad (B^2 J - (J.B) B) on every interior edge (run with nlim = 1
//!      and a tiny cfl_number, as the cs/sp checks are); same normalisation as cs_test.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"
#include "diffusion/current_density.hpp"
#include "pgen.hpp"

namespace {
int iprob;
Real b0, kwav, pert, cbar, s0, dens, ptot, eta_ad_t, b0c, bvx, bvy, bvz;

KOKKOS_INLINE_FUNCTION
Real Barenblatt(Real x, Real s, Real c) {
  Real q = c - x*x/(12.0*sqrt(s));
  return (q > 0.0) ? pow(s, -0.25)*sqrt(q) : 0.0;
}

void AmbipolarTestFinal(ParameterInput *pin, Mesh *pm);
}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  pgen_final_func = AmbipolarTestFinal;
  iprob = pin->GetOrAddInteger("problem", "iprob", 1);
  b0 = pin->GetOrAddReal("problem", "b0", 1.0);
  kwav = pin->GetOrAddReal("problem", "kwav", 2.0*M_PI);
  pert = pin->GetOrAddReal("problem", "pert", 0.0);
  cbar = pin->GetOrAddReal("problem", "cbar", 1.0);
  s0 = pin->GetOrAddReal("problem", "s0", 1.0e-3);
  dens = pin->GetOrAddReal("problem", "dens", 1.0);
  ptot = pin->GetOrAddReal("problem", "ptot", 10.0);
  eta_ad_t = pin->GetOrAddReal("mhd", "eta_ad", 0.0);
  b0c = pin->GetOrAddReal("problem", "b0c", 0.3);
  bvx = pin->GetOrAddReal("problem", "bvx", 0.3);
  bvy = pin->GetOrAddReal("problem", "bvy", 0.2);
  bvz = pin->GetOrAddReal("problem", "bvz", 0.5);
  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &size = pmbp->pmb->mb_size;
  const Real gm1 = pmbp->pmhd->peos->eos_data.gamma - 1.0;
  auto &w0 = pmbp->pmhd->w0;
  auto &b = pmbp->pmhd->b0;
  auto &bcc = pmbp->pmhd->bcc0;
  const bool one_d = pmy_mesh_->one_d;
  const int ip = iprob;
  const Real b0_ = b0, k_ = kwav, pe_ = pert, c_ = cbar, s0_ = s0, d_ = dens, pt_ = ptot;
  const Real ph_ = 1.0/sqrt(3.0);
  const Real b0c_ = b0c, bvx_ = bvx, bvy_ = bvy, bvz_ = bvz;

  // face fields; the 1-D helical and Barenblatt fields depend on x only, the ABC field
  // component along each face normal does not depend on that coordinate (div B = 0)
  par_for("ad_pgen_b", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke+1, js, je+1,
          is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
    const Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real x3v = CellCenterX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    const Real x1f = LeftEdgeX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real x2f = LeftEdgeX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    const Real x3f = LeftEdgeX(k-ks, nx3, size.d_view(m).x3min, size.d_view(m).x3max);
    const Real dx = size.d_view(m).dx1;
    Real b1, b2, b3;
    if (ip == 2) {
      // cell average of B_y over the x-cell (Simpson) so the initial error is O(dx^4)
      const Real s = s0_;
      b1 = 0.0;
      b2 = (Barenblatt(x1v-0.5*dx, s, c_) + 4.0*Barenblatt(x1v, s, c_)
            + Barenblatt(x1v+0.5*dx, s, c_))/6.0;
      b3 = 0.0;
    } else if (ip == 3) {
      // linear in the transverse coordinates: the centre value is the face average
      b1 = -b0c_*x2v + bvx_;
      b2 = b0c_*x1v + bvy_;
      b3 = bvz_;
    } else if (one_d) {
      // face averages of sin kx, cos kx over the x-cell
      const Real xl = x1v - 0.5*dx, xr = x1v + 0.5*dx;
      b1 = 0.0;
      b2 = b0_*(cos(k_*xl) - cos(k_*xr))/(k_*dx);
      b3 = b0_*(sin(k_*xr) - sin(k_*xl))/(k_*dx);
    } else {
      // ABC with A = B = C = 1; point values at the face centres (div B = 0 exactly:
      // each component is independent of its own coordinate)
      b1 = b0_*(sin(k_*x3v) + cos(k_*x2v)) + pe_*b0_*ph_;
      b2 = b0_*(sin(k_*x1v) + cos(k_*x3v)) + pe_*b0_*ph_;
      b3 = b0_*(sin(k_*x2v) + cos(k_*x1v)) + pe_*b0_*ph_;
    }
    (void) x1f; (void) x2f; (void) x3f;
    // each face array has the extra layer only along its own normal
    if (j <= je && k <= ke) b.x1f(m,k,j,i) = b1;
    if (i <= ie && k <= ke) b.x2f(m,k,j,i) = b2;
    if (i <= ie && j <= je) b.x3f(m,k,j,i) = b3;
  });
  par_for("ad_pgen_w", DevExeSpace(), 0, pmbp->nmb_thispack-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real bx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
    const Real by = one_d ? b.x2f(m,k,j,i) : 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
    const Real bz = (indcs.nx3 == 1) ? b.x3f(m,k,j,i)
                                     : 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
    bcc(m,IBX,k,j,i) = bx;
    bcc(m,IBY,k,j,i) = by;
    bcc(m,IBZ,k,j,i) = bz;
    w0(m,IDN,k,j,i) = d_;
    w0(m,IVX,k,j,i) = 0.0;
    w0(m,IVY,k,j,i) = 0.0;
    w0(m,IVZ,k,j,i) = 0.0;
    const Real pg = (ip == 2) ? pt_ - 0.5*(bx*bx + by*by + bz*bz) : pt_;
    w0(m,IEN,k,j,i) = pg/gm1;
  });
  pmbp->pmhd->peos->PrimToCons(w0, bcc, pmbp->pmhd->u0, is, ie, js, je, ks, ke);
  return;
}

namespace {
void AmbipolarTestFinal(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  auto &size = pmbp->pmb->mb_size;
  auto *pr = pmbp->pmhd->presist;
  if (iprob == 2) {
    const Real s = eta_ad_t*pm->time/3.0 + s0;
    const Real c_ = cbar;
    auto &bcc = pmbp->pmhd->bcc0;
    Real l1 = 0.0, l1ref = 0.0, l1in = 0.0;
    const Real xin = 0.5*sqrt(12.0*c_)*pow(s, 0.25);   // half the front position
    Kokkos::parallel_reduce("ad_l1", Kokkos::RangePolicy<>(DevExeSpace(), is, ie+1),
    KOKKOS_LAMBDA(const int i, Real &sum, Real &ref, Real &sin) {
      const Real dx = size.d_view(0).dx1;
      const Real x1v = CellCenterX(i-is, indcs.nx1, size.d_view(0).x1min,
                                   size.d_view(0).x1max);
      const Real ex = (Barenblatt(x1v-0.5*dx, s, c_) + 4.0*Barenblatt(x1v, s, c_)
                       + Barenblatt(x1v+0.5*dx, s, c_))/6.0;
      sum += fabs(bcc(0,IBY,ks,js,i) - ex)*dx;
      ref += fabs(ex)*dx;
      if (fabs(x1v) < xin) sin += fabs(bcc(0,IBY,ks,js,i) - ex)*dx;
    }, l1, l1ref, l1in);
    std::cout << "AMBIPOLAR_TEST iprob=2 nx1=" << indcs.nx1 << " t=" << pm->time
              << " s=" << s << " L1(By)=" << l1 << " L1rel=" << l1/l1ref
              << " L1(|x|<xf/2)=" << l1in
              << std::endl;
    return;
  }
  // max |div B| dx/|B| over the active cells (CT: round-off)
  {
    auto &b = pmbp->pmhd->b0;
    const bool md = pm->multi_d, td = pm->three_d;
    const int nx1 = indcs.nx1;
    const int nkji = indcs.nx3*indcs.nx2*indcs.nx1, nji = indcs.nx2*indcs.nx1;
    Real dmax = 0.0;
    Kokkos::parallel_reduce("ad_divb", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            pmbp->nmb_thispack*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dm) {
      int m = idx/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks; j += js;
      const Real dx = size.d_view(m).dx1;
      Real d = (b.x1f(m,k,j,i+1) - b.x1f(m,k,j,i))/dx;
      if (md) d += (b.x2f(m,k,j+1,i) - b.x2f(m,k,j,i))/size.d_view(m).dx2;
      if (td) d += (b.x3f(m,k+1,j,i) - b.x3f(m,k,j,i))/size.d_view(m).dx3;
      dm = fmax(dm, fabs(d)*dx);
    }, Kokkos::Max<Real>(dmax));
    std::cout << "AMBIPOLAR_TEST divB max|divB| dx = " << dmax << std::endl;
  }
  // iprob 3: the code's AD edge EMF (efld_resist of the last stage) against the exact
  if (iprob == 3 && pr != nullptr && pr->ad_type == 1) {
    auto e1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr->efld_resist.x1e);
    auto e2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr->efld_resist.x2e);
    auto e3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pr->efld_resist.x3e);
    const Real ead = pr->eta_ad, jz = 2.0*b0c;
    const Real xm = pm->mesh_size.x1max;
    const Real scale = fabs(ead*(SQR(b0c*xm) + bvx*bvx + bvy*bvy + bvz*bvz)*jz);
    auto exact = [&](Real x, Real y, int c) {
      const Real bx = -b0c*y + bvx, by = b0c*x + bvy, bz = bvz;
      const Real b2 = bx*bx + by*by + bz*bz, jb = jz*bz;
      const Real ev[3] = {-ead*jb*bx, -ead*jb*by, ead*(b2*jz - jb*bz)};
      return ev[c];
    };
    Real l1[3] = {0.0, 0.0, 0.0}, mx[3] = {0.0, 0.0, 0.0};
    int64_t nc[3] = {0, 0, 0};
    const int nb = 2;   // skip edges within nb cells of a physical boundary
    for (int m=0; m<pmbp->nmb_thispack; ++m) {
      auto sz = size.h_view(m);
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            const Real xf = LeftEdgeX(i-is, indcs.nx1, sz.x1min, sz.x1max);
            const Real yf = LeftEdgeX(j-js, indcs.nx2, sz.x2min, sz.x2max);
            const Real xc = xf + 0.5*sz.dx1, yc = yf + 0.5*sz.dx2;
            const bool in = (xf > pm->mesh_size.x1min + nb*sz.dx1 - 1e-12) &&
                            (xf < pm->mesh_size.x1max - nb*sz.dx1 + 1e-12) &&
                            (yf > pm->mesh_size.x2min + nb*sz.dx2 - 1e-12) &&
                            (yf < pm->mesh_size.x2max - nb*sz.dx2 + 1e-12);
            if (!in) continue;
            if (i <= ie) {
              const Real d = fabs(e1(m,k,j,i) - exact(xc, yf, 0))/scale;
              l1[0] += d; mx[0] = fmax(mx[0], d); ++nc[0];
            }
            if (j <= je) {
              const Real d = fabs(e2(m,k,j,i) - exact(xf, yc, 1))/scale;
              l1[1] += d; mx[1] = fmax(mx[1], d); ++nc[1];
            }
            if (k <= ke) {
              const Real d = fabs(e3(m,k,j,i) - exact(xf, yf, 2))/scale;
              l1[2] += d; mx[2] = fmax(mx[2], d); ++nc[2];
            }
          }
        }
      }
    }
    std::printf("AMBIPOLAR_TEST iprob=3 nx1=%d AD EMF vs exact (scale %.4e): x1e L1=%.4e "
                "max=%.4e | x2e L1=%.4e max=%.4e | x3e L1=%.4e max=%.4e\n", indcs.nx1,
                scale, l1[0]/nc[0], mx[0], l1[1]/nc[1], mx[1], l1[2]/nc[2], mx[2]);
    return;
  }
  // iprob 1: max |E_AD| over max |eta_A J|, evaluated on the final state
  if (pr == nullptr || pr->ad_type == 0) return;
  pr->ClearResistiveEMFs(pr->efld_resist);
  pr->AddEMFAmbipolar(pmbp->pmhd->b0, pmbp->pmhd->bcc0, pr->efld_resist);
  auto e1 = pr->efld_resist.x1e, e2 = pr->efld_resist.x2e, e3 = pr->efld_resist.x3e;
  auto j1 = pr->jcur.x1e, j2 = pr->jcur.x2e, j3 = pr->jcur.x3e;
  auto &bcc = pmbp->pmhd->bcc0;
  auto coef = pr->ad_coef;
  const bool one_d = pm->one_d;
  Real emax = 0.0, jmax = 0.0;
  const int nmkji = pmbp->nmb_thispack*indcs.nx3*indcs.nx2*indcs.nx1;
  const int nkji = indcs.nx3*indcs.nx2*indcs.nx1, nji = indcs.nx2*indcs.nx1;
  const int nx1 = indcs.nx1;
  Kokkos::parallel_reduce("ad_ff", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &em, Real &jm) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    const Real b2 = SQR(bcc(m,IBX,k,j,i)) + SQR(bcc(m,IBY,k,j,i)) + SQR(bcc(m,IBZ,k,j,i));
    const Real eta = coef(m,k,j,i)*b2;
    Real e = fmax(fabs(e2(m,k,j,i)), fabs(e3(m,k,j,i)));
    Real jj = fmax(fabs(j2(m,k,j,i)), fabs(j3(m,k,j,i)));
    if (!one_d) {
      e = fmax(e, fabs(e1(m,k,j,i)));
      jj = fmax(jj, fabs(j1(m,k,j,i)));
    }
    em = fmax(em, e);
    jm = fmax(jm, eta*jj);
  }, Kokkos::Max<Real>(emax), Kokkos::Max<Real>(jmax));
  std::cout << "AMBIPOLAR_TEST iprob=1 nx1=" << indcs.nx1 << " t=" << pm->time
            << " max|E_AD|=" << emax << " max|eta_A J|=" << jmax
            << " ratio=" << emax/jmax << std::endl;
}
}  // namespace
