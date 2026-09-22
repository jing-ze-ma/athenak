//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vet.cpp
//! \brief <rad_m1>/closure = vet_sc: a VARIABLE EDDINGTON TENSOR for the multi-D implicit
//! M1 solve from a grey formal solution of the transfer equation by SHORT
//! CHARACTERISTICS on the M1 grid.
//!
//! WHAT IS SOLVED.  Once per hydro step, at the start-of-step state,
//!     mu . grad I = chi (S - I),   chi = rho (kappa_F + kappa_s)  (= opac M1_OP_T)
//! with the grey source (E units, eps = 4 pi I / c, so that E = <eps>)
//!     S = eps_th a T^4 + (1 - eps_th) E,   eps_th = min(kappa_P / chi, 1),
//! i.e. (kappa_a B + kappa_s J)/chi with J taken from the current M1 E (no lambda
//! iteration: E is what the moment solve delivers and it is lagged with the tensor).
//! The stellar tables carry no separate scattering opacity (kappa_s = 0, kappa_F =
//! Rosseland incl. electron scattering), so the absorption fraction is estimated by the
//! Planck/extinction ratio and CAPPED: S is a convex blend of a T^4 and E, never
//! negative.
//! (The form E + (kappa_P a T^4 - kappa_E E)/chi, the module's own net emission, is NOT
//! usable: kappa_P >> kappa_R in the iron bump amplifies T_gas - T_rad by kappa_P/chi and
//! drove S to 0 at tau ~ 2-9, K/J = 0.48 there, in the first 200 s test.)
//! The Eddington tensor D = K/J of the formal solution is projected onto the uniaxial
//! form the solver takes (vet_axis = flux, the default):
//!     D ~ (1-chi)/2 delta + (3 chi - 1)/2 n n,   n = H/|H|,   chi = n.D.n
//! with H the flux of the FORMAL SOLUTION (not the M1 flux of the cell): the trace and
//! the component along the flux are kept exactly, the two perpendicular ones averaged
//! (exact for a plane-parallel or axisymmetric field).  chi is clamped to [1/3, 1].
//! vet_axis = eigen takes chi = the largest eigenvalue of D and n its vector instead;
//! in the nearly isotropic interior (|D - 1/3| ~ 1e-4) that axis FLIPS between x1, x2
//! and x3 from cell to cell, D_22 jumps by 1.5 (chi - 1/3) between neighbours, and
//! the huge deep E turns the jump into a horizontal force (first 200 s test: bottom-cell
//! rms v2 1.8e4 cm/s, KE_2 3.5e4 x the T3_tau arm).  The flux axis is smooth.
//!
//! ANGLES.  Double-Gauss in mu_1 = Omega.x1 (vet_nmu nodes on each hemisphere) times
//! vet_nphi uniform azimuths about x1; no ray is horizontal.
//!
//! SWEEP.  Layer by layer in x1: the ray arriving at a cell centre of layer i is traced
//! back to the PLANE of cell centres of the upwind layer i -+ 1, where I, chi and S are
//! (bi)linearly interpolated (periodic wrap in x2/x3); along the segment the source is
//! linear in tau (the positive, first-order short-characteristics weights).  Upward and
//! downward rays are swept in the same launches (layer is+l up, ie-l down), each launch
//! a par_for over the (k, j) plane with the rays looped inside, so the moments need no
//! atomics.  3-D is the same code with the x3 shift switched on.
//!
//! BOUNDARIES.  Bottom (x1min face), upward rays: the diffusion intensity
//! eps = E + 3 (F . Omega) / c of the bottom cell's M1 state, carried half a cell with
//! S of that cell.  Top (x1max face), downward rays: vacuum.
//!
//! LIMITS (checked at start-up).  The whole horizontal plane must live in one MeshBlock
//! and the sweep in x1 must not cross a block: ONE MeshBlock in the mesh, periodic x2
//! (and x3).  Several blocks would need the upwind plane from the neighbours (a halo of
//! width dx1 |mu_perp/mu_1| per layer, i.e. a pipelined sweep over the x1 stack and a
//! horizontal exchange per layer); MPI is not supported for the same reason.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

namespace {

void VetFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in rad_m1_vet.cpp" << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

// Gauss-Legendre nodes and weights on [0,1] (weights sum to 1)
void GaussLegendre01(const int n, std::vector<double> &x, std::vector<double> &w) {
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(M_PI*(i + 0.75)/(n + 0.5));
    double pp = 1.0;
    for (int it = 0; it < 200; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int jj = 1; jj <= n; ++jj) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*jj - 1.0)*z*p2 - (jj - 1.0)*p3)/jj;
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;
      if (std::fabs(z - z1) < 1.0e-15) break;
    }
    x[i] = 0.5*(1.0 + z);
    w[i] = 1.0/((1.0 - z*z)*pp*pp);
  }
}

// periodic index in [0, n)
KOKKOS_INLINE_FUNCTION
int VetWrap(const int a, const int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

// largest eigenvalue and its unit eigenvector of a symmetric 3x3 matrix (cyclic Jacobi)
KOKKOS_INLINE_FUNCTION
void VetEigMax(Real a[3][3], Real &lmax, Real &e1, Real &e2, Real &e3) {
  Real v[3][3];
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 6; ++sw) {
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  int im = 0;
  if (a[1][1] > a[im][im]) im = 1;
  if (a[2][2] > a[im][im]) im = 2;
  lmax = a[im][im];
  e1 = v[0][im];
  e2 = v[1][im];
  e3 = v[2][im];
}

// all eigenvalues lam[p] and unit eigenvectors v[.][p] of a symmetric 3x3 matrix (the
// cyclic Jacobi of VetEigMax, run to round-off); a is destroyed
KOKKOS_INLINE_FUNCTION
void VetEigSym(Real a[3][3], Real lam[3], Real v[3][3]) {
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 12; ++sw) {
    Real off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off <= 1.0e-18*(fabs(a[0][0]) + fabs(a[1][1]) + fabs(a[2][2]))) break;
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  for (int p = 0; p < 3; ++p) {lam[p] = a[p][p];}
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetInit

void RadiationM1::VetInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  if (!trans_on) {
    VetFatal("<rad_m1>/closure = vet_sc needs transport = implicit on a MULTI-D mesh");
  }
  if (pm->nmb_total != 1 || global_variable::nranks != 1) {
    VetFatal("<rad_m1>/closure = vet_sc: the short-characteristics sweep needs ONE "
             "MeshBlock (the whole horizontal plane and the whole x1 column) on one "
             "rank");
  }
  if (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
      (trans_x3 && pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic)) {
    VetFatal("<rad_m1>/closure = vet_sc needs periodic x2 (and x3) boundaries");
  }
  vet_nmu = pin->GetOrAddInteger("rad_m1", "vet_nmu", 4);
  vet_nphi = pin->GetOrAddInteger("rad_m1", "vet_nphi", 8);
  vet_milne = pin->GetOrAddBoolean("rad_m1", "vet_milne", false);
  vet_dump = pin->GetOrAddString("rad_m1", "vet_dump", "none");
  if (vet_dump.compare("none") == 0) {vet_dump = "";}
  vet_dump_every = pin->GetOrAddInteger("rad_m1", "vet_dump_every", 0);
  {std::string ax = pin->GetOrAddString("rad_m1", "vet_axis", "flux");
  if (ax.compare("flux") == 0) {
    vet_axis_flux = true;
  } else if (ax.compare("eigen") == 0) {
    vet_axis_flux = false;
  } else {
    VetFatal("<rad_m1>/vet_axis = '" + ax + "' not implemented (flux | eigen)");
  }
  }
  // vet_tensor = uniaxial (default) | full.  Read only when given, so that an input
  // without it writes the restart file the uniaxial path always wrote.
  vet_full = false;
  vet_eig_min = 0.0;
  if (pin->DoesParameterExist("rad_m1", "vet_tensor")) {
    std::string vt = pin->GetString("rad_m1", "vet_tensor");
    if (vt.compare("full") == 0) {
      vet_full = true;
    } else if (vt.compare("uniaxial") != 0) {
      VetFatal("<rad_m1>/vet_tensor = '" + vt + "' not implemented (uniaxial | full)");
    }
  }
  if (pin->DoesParameterExist("rad_m1", "vet_eig_min")) {
    vet_eig_min = pin->GetReal("rad_m1", "vet_eig_min");
  }
  if (vet_eig_min < 0.0 || vet_eig_min > 0.33) {
    VetFatal("<rad_m1>/vet_eig_min must lie in [0, 0.33]");
  }
  if (vet_nmu < 1 || vet_nphi < 1) {VetFatal("<rad_m1>/vet_nmu, vet_nphi must be >= 1");}
  vet_nray = 2*vet_nmu*vet_nphi;

  // the angle set: (mu1, mu2, mu3, weight), upward hemisphere first
  std::vector<double> xg, wg;
  GaussLegendre01(vet_nmu, xg, wg);
  Kokkos::realloc(vet_ang, vet_nray, 4);
  auto ah = Kokkos::create_mirror_view(vet_ang);
  int r = 0;
  for (int h = 0; h < 2; ++h) {
    for (int a = 0; a < vet_nmu; ++a) {
      double m1 = (h == 0) ? xg[a] : -xg[a];
      double st = std::sqrt(1.0 - m1*m1);
      for (int l = 0; l < vet_nphi; ++l) {
        double ph = 2.0*M_PI*(l + 0.5)/vet_nphi;
        ah(r, 0) = m1;
        ah(r, 1) = st*std::cos(ph);
        ah(r, 2) = st*std::sin(ph);
        ah(r, 3) = 0.5*wg[a]/vet_nphi;
        ++r;
      }
    }
  }
  Kokkos::deep_copy(vet_ang, ah);

  auto &indcs = pm->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pm->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(vet_cell, nmb, M1_VET_NC, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(vet_cell, 0.0);
  Kokkos::realloc(vet_ipl, nmb, 2, vet_nray, ncells3, ncells2);
  Kokkos::deep_copy(vet_ipl, 0.0);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> closure = vet_sc: short characteristics, " << vet_nray
              << " rays (double-Gauss nmu=" << vet_nmu << " x nphi=" << vet_nphi
              << "), lagged per hydro step" << (vet_milne ? ", MILNE source (diag)" : "")
              << (vet_full ? ", FULL tensor D = K/J (eigenvalue floor " : "")
              << (vet_full ? std::to_string(vet_eig_min) + ")" : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetShortChar

void RadiationM1::VetShortChar() {
  Kokkos::fence();
  Kokkos::Timer timer;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const int nray = vet_nray;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto ip_ = vet_ipl;
  auto ang_ = vet_ang;
  auto opac_ = opac;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real ar = arad;
  const Real efl = e_floor;
  const bool thermal = (pmy_pack->phydro != nullptr) && coupling && !opac_zero;
  const bool milne = vet_milne;
  const Real fmil = iflux_x1min;

  // (1) extinction and source, column by column (the Milne diagnostic needs tau)
  par_for("m1_vet_src", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real dx1 = mbsize.d_view(m).dx1;
    Real tau = 0.0;
    for (int i = ie; i >= is; --i) {
      Real chx = fmax(opac_(m,M1_OP_T,k,j,i), 1.0e-300);
      Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
      Real s = e;
      if (milne) {
        // exact grey Milne problem: S = J = 3 H (tau + q(tau)), Hopf q fitted
        Real tc = tau + 0.5*chx*dx1;
        Real qh = 0.710446 - 0.133054*exp(-3.4488*tc);
        s = 3.0*fmil/cl*(tc + qh);
      } else if (thermal) {
        Real t = iw_(m,M1_IW_TP,k,j,i);
        Real t2 = t*t;
        // thermal fraction: the Planck (emission) mean over the extinction, capped at 1
        Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chx, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
      tau += chx*dx1;
      vc_(m,M1_VET_CHX,k,j,i) = chx;
      vc_(m,M1_VET_SRC,k,j,i) = s;
      for (int n = M1_VET_J; n < M1_VET_CHI; ++n) {vc_(m,n,k,j,i) = 0.0;}
    }
  });

  // (2) the sweep: launch l does layer is+l for the upward rays and ie-l for the
  // downward ones, reading the intensity of launch l-1 from the other plane buffer
  for (int l = 0; l < nx1; ++l) {
    const int pw = l & 1, pr = pw ^ 1;
    par_for("m1_vet_sweep", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const Real dx1 = mbsize.d_view(m).dx1;
      const Real dx2 = mbsize.d_view(m).dx2;
      const Real dx3 = mbsize.d_view(m).dx3;
      const int iu = is + l, id = ie - l;
      Real acu[10], acd[10];
      for (int n = 0; n < 10; ++n) {acu[n] = 0.0; acd[n] = 0.0;}
      for (int r = 0; r < nray; ++r) {
        const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
        const bool up = (m1 > 0.0);
        const int i = up ? iu : id;
        const Real am1 = fabs(m1);
        const Real c0 = vc_(m,M1_VET_CHX,k,j,i);
        const Real s0 = vc_(m,M1_VET_SRC,k,j,i);
        Real iv;
        if (l == 0) {
          // half a cell from the boundary face; S at the face linearly extrapolated
          // from the two boundary cells (exact for the deep linear S of a diffusion
          // regime), and linear in tau along the segment
          const int iin = up ? (i + 1) : (i - 1);
          const Real su = fmax(1.5*s0 - 0.5*vc_(m,M1_VET_SRC,k,j,iin), 0.0);
          Real ib = 0.0;
          if (up) {
            if (milne) {
              // the exact deep Milne intensity 3 H (tau + q + mu) = S + 3 H mu
              ib = su + 3.0*fmil/cl*m1;
            } else {
              // diffusion: eps = E + 3 (F . Omega)/c from the bottom cell's M1 state
              Real e = fmax(iw_(m,M1_IW_EN,k,j,i), efl);
              Real fo = iw_(m,M1_IW_F1,k,j,i)*m1 + iw_(m,M1_IW_F2,k,j,i)*m2
                        + iw_(m,M1_IW_F3,k,j,i)*m3;
              ib = fmax(e + 3.0*fo/cl, 0.0);
            }
          }
          const Real dtau = 0.5*c0*dx1/am1;
          const Real ex = exp(-dtau);
          Real w0, wu;
          if (dtau < 1.0e-3) {
            w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
            wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
          } else {
            const Real g = (1.0 - ex)/dtau;
            w0 = 1.0 - g;
            wu = g - ex;
          }
          iv = ib*ex + wu*su + w0*s0;
        } else {
          const int iup = up ? (i - 1) : (i + 1);
          // the foot of the characteristic on the upwind plane, in cell units
          const Real sh2 = -dx1*m2/(am1*dx2);
          const Real fl2 = floor(sh2);
          const Real a2 = sh2 - fl2;
          const int o2 = static_cast<int>(fl2);
          const int ja = js + VetWrap(j - js + o2, nx2);
          const int jb = js + VetWrap(j - js + o2 + 1, nx2);
          int ka = k, kb = k;
          Real a3 = 0.0;
          if (thrd) {
            const Real sh3 = -dx1*m3/(am1*dx3);
            const Real fl3 = floor(sh3);
            a3 = sh3 - fl3;
            const int o3 = static_cast<int>(fl3);
            ka = ks + VetWrap(k - ks + o3, nx3);
            kb = ks + VetWrap(k - ks + o3 + 1, nx3);
          }
          const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
          const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
          const Real iup_v = waa*ip_(m,pr,r,ka,ja) + wab*ip_(m,pr,r,ka,jb)
                             + wba*ip_(m,pr,r,kb,ja) + wbb*ip_(m,pr,r,kb,jb);
          const Real cu = waa*vc_(m,M1_VET_CHX,ka,ja,iup)
                          + wab*vc_(m,M1_VET_CHX,ka,jb,iup)
                          + wba*vc_(m,M1_VET_CHX,kb,ja,iup)
                          + wbb*vc_(m,M1_VET_CHX,kb,jb,iup);
          const Real su = waa*vc_(m,M1_VET_SRC,ka,ja,iup)
                          + wab*vc_(m,M1_VET_SRC,ka,jb,iup)
                          + wba*vc_(m,M1_VET_SRC,kb,ja,iup)
                          + wbb*vc_(m,M1_VET_SRC,kb,jb,iup);
          const Real dtau = 0.5*(cu + c0)*dx1/am1;
          const Real ex = exp(-dtau);
          Real w0, wu;
          if (dtau < 1.0e-3) {
            w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
            wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
          } else {
            const Real g = (1.0 - ex)/dtau;
            w0 = 1.0 - g;
            wu = g - ex;
          }
          iv = iup_v*ex + wu*su + w0*s0;
        }
        iv = fmax(iv, 0.0);
        ip_(m,pw,r,k,j) = iv;
        Real *ac = up ? acu : acd;
        const Real a = wr*iv;
        ac[0] += a;
        ac[1] += a*m1*m1;
        ac[2] += a*m2*m2;
        ac[3] += a*m3*m3;
        ac[4] += a*m1*m2;
        ac[5] += a*m1*m3;
        ac[6] += a*m2*m3;
        ac[7] += a*m1;
        ac[8] += a*m2;
        ac[9] += a*m3;
      }
      for (int n = 0; n < 10; ++n) {
        vc_(m,M1_VET_J+n,k,j,iu) += acu[n];
        vc_(m,M1_VET_J+n,k,j,id) += acd[n];
      }
    });
  }

  // (3) the uniaxial projection (chi, n) of D = K/J
  const bool axflux = vet_axis_flux;
  par_for("m1_vet_proj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real jj = vc_(m,M1_VET_J,k,j,i);
    Real chi = 1.0/3.0, n1 = 1.0, n2 = 0.0, n3 = 0.0;
    if (jj > 0.0) {
      Real ij = 1.0/jj;
      Real a[3][3];
      a[0][0] = vc_(m,M1_VET_K11,k,j,i)*ij;
      a[1][1] = vc_(m,M1_VET_K11+1,k,j,i)*ij;
      a[2][2] = vc_(m,M1_VET_K11+2,k,j,i)*ij;
      a[0][1] = a[1][0] = vc_(m,M1_VET_K11+3,k,j,i)*ij;
      a[0][2] = a[2][0] = vc_(m,M1_VET_K11+4,k,j,i)*ij;
      a[1][2] = a[2][1] = vc_(m,M1_VET_K11+5,k,j,i)*ij;
      // axis = the formal solution's OWN flux direction, chi = n.D.n (vet_axis = flux)
      // or the principal axis of D (vet_axis = eigen, see the file header)
      Real h1 = vc_(m,M1_VET_H1,k,j,i), h2 = vc_(m,M1_VET_H1+1,k,j,i);
      Real h3 = vc_(m,M1_VET_H1+2,k,j,i);
      Real hm = sqrt(h1*h1 + h2*h2 + h3*h3);
      if (axflux && hm > 1.0e-8*jj) {
        n1 = h1/hm;
        n2 = h2/hm;
        n3 = h3/hm;
        chi = n1*(a[0][0]*n1 + a[0][1]*n2 + a[0][2]*n3)
              + n2*(a[1][0]*n1 + a[1][1]*n2 + a[1][2]*n3)
              + n3*(a[2][0]*n1 + a[2][1]*n2 + a[2][2]*n3);
      } else {
        VetEigMax(a, chi, n1, n2, n3);
      }
      Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
      if (nn > 0.0) {
        n1 /= nn;
        n2 /= nn;
        n3 /= nn;
      } else {
        n1 = 1.0;
        n2 = 0.0;
        n3 = 0.0;
      }
      chi = fmin(fmax(chi, 1.0/3.0), 1.0);
    }
    vc_(m,M1_VET_CHI,k,j,i) = chi;
    vc_(m,M1_VET_N1,k,j,i) = n1;
    vc_(m,M1_VET_N1+1,k,j,i) = n2;
    vc_(m,M1_VET_N1+2,k,j,i) = n3;
  });

  // (4) vet_tensor = full: the GUARDED D = K/J the solve reads in place of (chi, n)
  if (vet_full) {VetFullTensor();}
  Kokkos::fence();
  vet_time += timer.seconds();

  // column dumps (diagnostic, outside the timer)
  if (!vet_dump.empty()) {
    int nc = static_cast<int>(vet_ncall);
    if (nc == 0 || (vet_dump_every > 0 && (nc % vet_dump_every) == 0)) {
      char buf[32];
      std::snprintf(buf, sizeof(buf), ".%06d.txt", nc);
      VetDump(vet_dump + std::string(buf));
    }
  }
  vet_ncall += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetFullTensor
//! \brief vet_tensor = full: D_ab = K_ab/J of the last formal solution, REALIZABILITY
//! GUARDED, into vet_cell(M1_VET_D11 .. +5), active cells and ghosts.
//!
//! By construction K = sum_r w_r I_r Omega_r Omega_r^T with w_r > 0 and I_r >= 0 (the
//! sweep clips I at 0) and sum_r w_r Omega_r Omega_r^T has trace sum_r w_r = J/J, so the
//! raw D is symmetric, has trace 1 to round-off and is positive semidefinite with
//! eigenvalues in [0, 1].  The guard makes that EXACT and enforces a floor:
//!   (1) J <= 0 or a non-positive trace: D = delta/3 (counted);
//!   (2) D /= tr D (a round-off correction, not counted);
//!   (3) eigen-decomposition (cyclic Jacobi); if an eigenvalue lies below vet_eig_min
//!       (default 0) the spectrum is mapped affinely, lam -> emin + (lam' - emin)
//!       (1 - 3 emin)/(sum lam' - 3 emin) with lam' = max(lam, emin), which keeps the
//!       eigenvectors, the trace 1 and every eigenvalue in [emin, 1] (counted).
//! Consequences for the solve: D_dd >= emin >= 0 on every diagonal, so each x1 line row
//! gains nu df (D_11,i + D_11,i) >= 0 on its diagonal and loses nu df D_11,i+-1 <= 0 off
//! it -- the column-diagonally-dominant M-matrix of the uniaxial form, and the diagonal
//! stays >= 1.  The transverse diagonal TDIA >= 0 likewise.
//! Ghosts: periodic wrap in x2/x3 (the SC sweep already requires ONE MeshBlock with
//! periodic x2/x3), copy of the edge cell in x1 (never read at a physical x1 face).
//! Per-cell |D_guarded - D_raw| goes to M1_VET_GD for the run statistics.

void RadiationM1::VetFullTensor() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*(indcs.ng);
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  const bool thrd = trans_x3;
  const Real emin = vet_eig_min;
  auto vc_ = vet_cell;

  par_for("m1_vet_full", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jj = vc_(m,M1_VET_J,k,j,i);
    Real d[6];
    Real raw[6];
    for (int n = 0; n < 6; ++n) {raw[n] = 0.0;}
    Real tr = 0.0;
    if (jj > 0.0) {
      const Real ij = 1.0/jj;
      for (int n = 0; n < 6; ++n) {raw[n] = vc_(m,M1_VET_K11+n,k,j,i)*ij;}
      tr = raw[0] + raw[1] + raw[2];
    }
    bool hit = false;
    if (!(jj > 0.0) || !(tr > 0.0)) {
      d[0] = d[1] = d[2] = 1.0/3.0;
      d[3] = d[4] = d[5] = 0.0;
      hit = true;
    } else {
      const Real itr = 1.0/tr;
      for (int n = 0; n < 6; ++n) {d[n] = raw[n]*itr;}
      Real a[3][3], v[3][3], lam[3];
      a[0][0] = d[0];
      a[1][1] = d[1];
      a[2][2] = d[2];
      a[0][1] = a[1][0] = d[3];
      a[0][2] = a[2][0] = d[4];
      a[1][2] = a[2][1] = d[5];
      VetEigSym(a, lam, v);
      if (lam[0] < emin || lam[1] < emin || lam[2] < emin) {
        Real sl = 0.0;
        for (int p = 0; p < 3; ++p) {
          lam[p] = fmax(lam[p], emin);
          sl += lam[p];
        }
        const Real sc = (1.0 - 3.0*emin)/fmax(sl - 3.0*emin, 1.0e-300);
        for (int p = 0; p < 3; ++p) {lam[p] = emin + (lam[p] - emin)*sc;}
        Real b[3][3];
        for (int p = 0; p < 3; ++p) {
          for (int q = 0; q < 3; ++q) {
            b[p][q] = lam[0]*v[p][0]*v[q][0] + lam[1]*v[p][1]*v[q][1]
                      + lam[2]*v[p][2]*v[q][2];
          }
        }
        d[0] = b[0][0];
        d[1] = b[1][1];
        d[2] = b[2][2];
        d[3] = 0.5*(b[0][1] + b[1][0]);
        d[4] = 0.5*(b[0][2] + b[2][0]);
        d[5] = 0.5*(b[1][2] + b[2][1]);
        hit = true;
      }
    }
    Real dev = 0.0;
    for (int n = 0; n < 6; ++n) {
      vc_(m,M1_VET_D11+n,k,j,i) = d[n];
      dev = fmax(dev, fabs(d[n] - raw[n]));
    }
    // <= -1 marks a cell the guard left alone (trace normalisation only): -1 - dev
    vc_(m,M1_VET_GD,k,j,i) = hit ? dev : (-1.0 - dev);
  });

  // ghosts: periodic in x2 (x3), edge copy in x1
  par_for("m1_vet_full_gh", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
    if (act) return;
    const int ii = (i < is) ? is : ((i > ie) ? ie : i);
    const int jj = js + VetWrap(j - js, nx2);
    const int kk = thrd ? (ks + VetWrap(k - ks, nx3)) : k;
    for (int n = 0; n < 6; ++n) {
      vc_(m,M1_VET_D11+n,k,j,i) = vc_(m,M1_VET_D11+n,kk,jj,ii);
    }
  });

  // run statistics of the guard (host reductions over the active cells)
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nc = (nmb1 + 1)*nk*nj*ni;
  Real nhit = 0.0, gmax = 0.0, dmin = 1.0;
  Kokkos::parallel_reduce("m1_vet_full_st", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
  KOKKOS_LAMBDA(const int idx, Real &sh, Real &sg, Real &sd) {
    const int m = idx/(nk*nj*ni);
    int r = idx - m*(nk*nj*ni);
    const int k = ks + r/(nj*ni);
    r -= (k - ks)*(nj*ni);
    const int j = js + r/ni;
    const int i = is + r - (j - js)*ni;
    const Real g = vc_(m,M1_VET_GD,k,j,i);
    if (g >= 0.0) {sh += 1.0;}
    sg = fmax(sg, (g >= 0.0) ? g : (-1.0 - g));
    sd = fmin(sd, vc_(m,M1_VET_D11,k,j,i));
  }, Kokkos::Sum<Real>(nhit), Kokkos::Max<Real>(gmax), Kokkos::Min<Real>(dmin));
  vet_nguard += nhit;
  vet_ncell += static_cast<Real>(nc);
  vet_guard_max = std::max(vet_guard_max, gmax);
  vet_d11_min = std::min(vet_d11_min, dmin);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetDump

void RadiationM1::VetDump(const std::string &fname) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, ks = indcs.ks;
  auto vh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_cell);
  auto ih = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
  Real dx1 = pmy_pack->pmb->mb_size.h_view(0).dx1;
  Real x1m = pmy_pack->pmb->mb_size.h_view(0).x1min;
  {
    // the whole (j, i) plane at k = ks: raw K/J, the uniaxial projection and, under
    // vet_tensor = full, the guarded D the solve reads
    const int je = indcs.je;
    std::ofstream g(fname + ".plane");
    g << "# closure = vet_sc plane dump, m = 0, k = ks; call "
      << static_cast<int>(vet_ncall) << (vet_full ? " (full)" : " (uniaxial)") << "\n"
      << "# 1 i  2 j  3 J  4-9 K/J 11 22 33 12 13 23  10-12 H/J  13 chi  14-16 n  "
      << "17-22 D_guarded 11 22 33 12 13 23 (0 unless full)  23 E_m1  24 chi_ext\n";
    g << std::setprecision(10);
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real jj = vh(0,M1_VET_J,ks,j,i);
        Real ij = (jj > 0.0) ? 1.0/jj : 0.0;
        g << i - is << " " << j - js << " " << jj;
        for (int n = 0; n < 6; ++n) {g << " " << vh(0,M1_VET_K11+n,ks,j,i)*ij;}
        for (int n = 0; n < 3; ++n) {g << " " << vh(0,M1_VET_H1+n,ks,j,i)*ij;}
        for (int n = 0; n < 4; ++n) {g << " " << vh(0,M1_VET_CHI+n,ks,j,i);}
        for (int n = 0; n < 6; ++n) {
          g << " " << (vet_full ? vh(0,M1_VET_D11+n,ks,j,i) : 0.0);
        }
        g << " " << ih(0,M1_IW_EN,ks,j,i) << " " << vh(0,M1_VET_CHX,ks,j,i) << "\n";
      }
    }
  }
  std::ofstream f(fname);
  f << "# closure = vet_sc column dump, m = 0, k = ks, j = js; call "
    << static_cast<int>(vet_ncall) << "\n"
    << "# 1 i  2 x1  3 tau_top(centre)  4 chi_ext  5 S  6 J_sc  7 K11/J  8 K22/J  "
    << "9 K33/J  10 K12/J  11 H1/J (=F1/cE)  12 chi_proj  13 n1  14 n2  15 E_m1  "
    << "16 F1/(c E)_m1  17 K/J Milne-exact (tau+q_inf)/(3(tau+q(tau)))\n";
  f << std::setprecision(10);
  Real tau = 0.0;
  std::vector<Real> tc(ie + 1);
  for (int i = ie; i >= is; --i) {
    Real chx = vh(0,M1_VET_CHX,ks,js,i);
    tc[i] = tau + 0.5*chx*dx1;
    tau += chx*dx1;
  }
  for (int i = is; i <= ie; ++i) {
    Real jj = vh(0,M1_VET_J,ks,js,i);
    Real ij = (jj > 0.0) ? 1.0/jj : 0.0;
    Real t = tc[i];
    Real qh = 0.710446 - 0.133054*std::exp(-3.4488*t);
    Real exact = (t + 0.710446)/(3.0*(t + qh));
    Real e = ih(0,M1_IW_EN,ks,js,i);
    f << i - is << " " << x1m + (i - is + 0.5)*dx1 << " " << t << " "
      << vh(0,M1_VET_CHX,ks,js,i) << " " << vh(0,M1_VET_SRC,ks,js,i) << " " << jj
      << " " << vh(0,M1_VET_K11,ks,js,i)*ij << " " << vh(0,M1_VET_K11+1,ks,js,i)*ij
      << " " << vh(0,M1_VET_K11+2,ks,js,i)*ij << " " << vh(0,M1_VET_K11+3,ks,js,i)*ij
      << " " << vh(0,M1_VET_H1,ks,js,i)*ij << " " << vh(0,M1_VET_CHI,ks,js,i)
      << " " << vh(0,M1_VET_N1,ks,js,i) << " " << vh(0,M1_VET_N1+1,ks,js,i)
      << " " << e << " " << ih(0,M1_IW_F1,ks,js,i)/(c_light*e) << " " << exact << "\n";
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetReport

void RadiationM1::VetReport() {
  if (global_variable::my_rank != 0 || vet_ncall <= 0.0) return;
  std::cout << "<rad_m1> vet_sc: formal solutions=" << vet_ncall
            << " SC seconds=" << vet_time << " (" << vet_time/vet_ncall << " per call)"
            << " implicit-solve seconds incl. SC=" << vet_itime
            << " SC/(implicit solve without SC)="
            << vet_time/std::fmax(vet_itime - vet_time, 1.0e-300) << std::endl;
  if (vet_full) {
    std::cout << "<rad_m1> vet_tensor=full: guard triggered in " << vet_nguard << " of "
              << vet_ncell << " cell-calls (eigenvalue floor " << vet_eig_min
              << "); max |D_guarded - K/J|=" << vet_guard_max
              << " (incl. trace normalisation); min D_11=" << vet_d11_min << std::endl;
  }
}

} // namespace radm1
