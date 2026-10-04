//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd.cpp
//! \brief <rad_m1>/vet_gd = true (m1-vet-gd; design rt_design_1003/SC_PROPER_SP.md,
//! candidate b): the Eddington tensor on the spherical-polar wedge from a
//! short-characteristics (SC) formal solution along GLOBALLY FIXED directions.
//!
//! WHAT IS SOLVED.  At every vet_col build (the same state, source S and extinction chi
//! as vet_col and vet_col_lat: VetLatSweep(0) fills ln chi, ln S and the first shell
//! vlat_icut), the grey transfer equation  n . grad I = chi (S - I)  is swept shell by
//! shell for the N = 12 vet_gd_nside^2 HEALPix ring centres n_d, fixed in the CARTESIAN
//! frame (straight rays: the curvature of the sp grid enters only through the positions,
//! there is no angular redistribution and no angular numerical diffusion).  The weights
//! are equal up to one global factor in n_z^2 that makes the set integrate 1, n and nn
//! exactly, so f = 1/3 holds for an isotropic or linear intensity in every cell frame.
//! Two passes, one launch per shell, one thread per (block, cell, direction):
//!   inward, top-down, for the (cell, direction) pairs with n . r_hat < 0: the upwind
//!     point is where the backward ray meets the sphere of shell i+1 (the top shell from
//!     the vacuum top face);
//!   outward, bottom-up, for n . r_hat >= 0: the sphere of shell i-1, or, for a ray that
//!     misses it, the same shell beyond the tangent point (an inward value of the first
//!     pass), or the inner face of the first shell (the diffusion intensity, as
//!     vet_col_lat).
//! For straight rays this order has no cycles (r is monotone on either side of the one
//! tangent point).  At the upwind point I (the SAME direction index), ln chi and ln S are
//! bilinear in (theta, phi) index space; the source is linear in tau along the segment
//! (vet_col's positive first-order weights).  This is rt_design_1003/sc_proper/gd_sc.py,
//! against which the C++ is gated.
//!
//! WALLS.  At a periodic theta or phi wall of the wedge the ghost cells hold the copied
//! intensities of the other side; they are re-indexed after every exchange by the LOCAL
//! FRAME convention of the run's own periodic boundaries (the (r, theta, phi)
//! components of a direction are kept, i.e. n -> M n with M the rotation between the
//! two frames, and the nearest set direction is read).  For a phi period of pi/2 the
//! HEALPix set is invariant and the map is an exact permutation.
//!
//! SEVERAL MESHBLOCKS.  An upwind point inside the block reads the current sweep; one in
//! the lateral ghost band reads the neighbour's intensity of the LAST sweep (exchanged
//! after every sweep): a lagged inflow across block boundaries (block-Jacobi; the first
//! build iterates vet_col_lat_init_iter times).  Reads beyond the ghost band are clamped
//! and counted (very oblique rays through laterally thin cells).
//!
//! HAND-OVER: the vet_col_lat interface (tau_ten slots M1_TT_LAT0.., the mesh basis):
//!   LAT0 = D_rr - f_K(vet_col), folded into slot 0 by VetLatBuild (so slot 0 = D_rr of
//!          the sweep above the first shell, vet_col's f_K below it),
//!   LAT1, LAT2 = D_r,theta and D_r,phi (the lagged Picard term M1SphLat when
//!          vet_col_lat_offdiag),
//!   LAT3, LAT4, LAT5 = D_tt - (1 - f_K)/2, D_tp, D_pp - (1 - f_K)/2.
//! All vet_col_lat keys apply (taucut, every, init_iter, offdiag, dump).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {

namespace {
void VgdFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/vet_gd: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

KOKKOS_INLINE_FUNCTION
Real VgdLerp(const Real a, const Real b, const Real u) {return a + u*(b - a);}

// the positive first-order SC weights (VcolW of rad_m1_vetcol.cpp, term for term)
KOKKOS_INLINE_FUNCTION
void VgdW(const Real dtau, Real &ex, Real &w0, Real &wu) {
  ex = exp(-dtau);
  if (dtau < 1.0e-3) {
    w0 = dtau*(0.5 - dtau*(1.0/6.0 - dtau/24.0));
    wu = dtau*(0.5 - dtau*(1.0/3.0 - dtau/8.0));
  } else {
    const Real g = (1.0 - ex)/dtau;
    w0 = 1.0 - g;
    wu = g - ex;
  }
}

// HEALPix ring-scheme centres (12 ns^2) with the moment-fixed weights
void VgdDirections(const int ns, std::vector<double> &d) {
  d.clear();
  for (int i = 1; i < 4*ns; ++i) {
    double z, ph0, dph;
    int cnt;
    if (i < ns) {
      z = 1.0 - static_cast<double>(i*i)/(3.0*ns*ns);
      cnt = 4*i; dph = M_PI/(2.0*i); ph0 = -0.5*dph;
    } else if (i <= 3*ns) {
      z = 4.0/3.0 - 2.0*i/(3.0*ns);
      const int s = (i - ns + 1) % 2;
      cnt = 4*ns; dph = M_PI/(2.0*ns); ph0 = -0.5*s*dph;
    } else {
      const int ii = 4*ns - i;
      z = -(1.0 - static_cast<double>(ii*ii)/(3.0*ns*ns));
      cnt = 4*ii; dph = M_PI/(2.0*ii); ph0 = -0.5*dph;
    }
    const double st = std::sqrt(std::max(1.0 - z*z, 0.0));
    for (int j = 1; j <= cnt; ++j) {
      const double p = ph0 + j*dph;
      d.push_back(st*std::cos(p));
      d.push_back(st*std::sin(p));
      d.push_back(z);
      d.push_back(0.0);
    }
  }
  const int n = static_cast<int>(d.size()/4);
  double m2 = 0.0, m4 = 0.0;
  for (int q = 0; q < n; ++q) {
    const double z2 = d[4*q+2]*d[4*q+2];
    m2 += z2/n;
    m4 += z2*z2/n;
  }
  const double c = (1.0/3.0 - m2)/(m4 - m2*m2);
  for (int q = 0; q < n; ++q) {
    d[4*q+3] = (1.0 + c*(d[4*q+2]*d[4*q+2] - m2))/n;
  }
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdInit
//! \brief direction table, intensity array, exchange object, wall maps (from VetLatInit)

void RadiationM1::VetGdInit() {
  Mesh *pm = pmy_pack->pmesh;
  if (vgd_nside < 1 || vgd_nside > 6) {VgdFatal("vet_gd_nside must be 1..6");}
  if (!pm->three_d) {VgdFatal("needs a 3-D wedge");}
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int c2 = indcs.nx2 + 2*indcs.ng;
  const int c3 = indcs.nx3 + 2*indcs.ng;
  std::vector<double> d;
  VgdDirections(vgd_nside, d);
  vgd_n = static_cast<int>(d.size()/4);
  const int n = vgd_n;
  Kokkos::realloc(vgd_dir, n, 4);
  auto d_h = Kokkos::create_mirror_view(vgd_dir);
  for (int q = 0; q < n; ++q) {
    for (int c = 0; c < 4; ++c) {d_h(q,c) = d[4*q+c];}
  }
  Kokkos::deep_copy(vgd_dir, d_h);
  Kokkos::realloc(vgd_i, nmb, n, c3, c2, c1);
  Kokkos::deep_copy(vgd_i, 0.0);
  Kokkos::realloc(vgd_i_c, nmb, n, 1, 1, 1);
  pbval_gd = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_gd->InitializeBuffers(n);
  pbval_gd->SetVectorPairs(n, {});
  // wall maps: lateral ghost cells outside the mesh in theta or phi hold the other
  // side's intensities; new(d) = old(map(d)), map = nearest set direction to M n_d
  Kokkos::realloc(vgd_wall, nmb, c3, c2);
  Kokkos::realloc(vgd_map, nmb, c3, c2, n);
  auto w_h = Kokkos::create_mirror_view(vgd_wall);
  auto mp_h = Kokkos::create_mirror_view(vgd_map);
  auto x2_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x2v);
  auto x3_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x3v);
  const double t0 = pm->mesh_size.x2min, t1 = pm->mesh_size.x2max;
  const double p0 = pm->mesh_size.x3min, p1 = pm->mesh_size.x3max;
  int nwall = 0, nexact = 0, nmap = 0;
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        const double th = x2_h(m,j), ph = x3_h(m,k);
        double thw = th, phw = ph;
        if (th < t0) {thw = th + (t1 - t0);}
        if (th > t1) {thw = th - (t1 - t0);}
        if (ph < p0) {phw = ph + (p1 - p0);}
        if (ph > p1) {phw = ph - (p1 - p0);}
        const bool wall = (thw != th) || (phw != ph);
        w_h(m,k,j) = wall ? 1 : 0;
        for (int q = 0; q < n; ++q) {mp_h(m,k,j,q) = q;}
        if (!wall) {continue;}
        ++nwall;
        const double st = std::sin(th), ct = std::cos(th);
        const double sp = std::sin(ph), cp = std::cos(ph);
        const double sw = std::sin(thw), cw = std::cos(thw);
        const double spw = std::sin(phw), cpw = std::cos(phw);
        for (int q = 0; q < n; ++q) {
          const double nx = d[4*q], ny = d[4*q+1], nz = d[4*q+2];
          const double a = nx*st*cp + ny*st*sp + nz*ct;
          const double b = nx*ct*cp + ny*ct*sp - nz*st;
          const double c = -nx*sp + ny*cp;
          const double mx = a*sw*cpw + b*cw*cpw - c*spw;
          const double my = a*sw*spw + b*cw*spw + c*cpw;
          const double mz = a*cw - b*sw;
          int best = 0;
          double bd = -2.0;
          for (int r = 0; r < n; ++r) {
            const double dd = mx*d[4*r] + my*d[4*r+1] + mz*d[4*r+2];
            if (dd > bd) {bd = dd; best = r;}
          }
          mp_h(m,k,j,q) = best;
          ++nmap;
          if (bd > 1.0 - 1.0e-12) {++nexact;}
        }
      }
    }
  }
  // the branch of every stored value: n . r_hat of the direction actually stored at a
  // lateral cell (the mapped one, at the source position, for a wall ghost)
  Kokkos::realloc(vgd_mr, nmb, c3, c2, n);
  auto mr_h = Kokkos::create_mirror_view(vgd_mr);
  for (int m = 0; m < nmb; ++m) {
    for (int k = 0; k < c3; ++k) {
      for (int j = 0; j < c2; ++j) {
        double th = x2_h(m,j), ph = x3_h(m,k);
        if (th < t0) {th += (t1 - t0);}
        if (th > t1) {th -= (t1 - t0);}
        if (ph < p0) {ph += (p1 - p0);}
        if (ph > p1) {ph -= (p1 - p0);}
        const double st = std::sin(th), ct = std::cos(th);
        const double sp = std::sin(ph), cp = std::cos(ph);
        for (int q = 0; q < n; ++q) {
          const int r = mp_h(m,k,j,q);
          mr_h(m,k,j,q) = d[4*r]*st*cp + d[4*r+1]*st*sp + d[4*r+2]*ct;
        }
      }
    }
  }
  Kokkos::deep_copy(vgd_mr, mr_h);
  Kokkos::deep_copy(vgd_wall, w_h);
  Kokkos::deep_copy(vgd_map, mp_h);
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> vet_gd: global-direction SC closure on the sp wedge, nside "
              << vgd_nside << " (" << n << " HEALPix directions, moment-fixed weights);"
              << " rank 0 wall ghost columns " << nwall << ", direction maps exact "
              << nexact << " of " << nmap << " (the rest: nearest direction)"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdWall
//! \brief re-index the wall ghost intensities after an exchange (local-frame periodicity)

void RadiationM1::VetGdWall() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int c2 = indcs.nx2 + 2*indcs.ng, c3 = indcs.nx3 + 2*indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  auto vi_ = vgd_i;
  auto wl_ = vgd_wall;
  auto mp_ = vgd_map;
  const int ilo = is + vlat_icut;
  par_for("m1_vgd_wall", DevExeSpace(), 0, nmb1, 0, c3 - 1, 0, c2 - 1, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (wl_(m,k,j) == 0) {return;}
    Real buf[M1_VGD_NMAX];
    for (int q = 0; q < n; ++q) {buf[q] = vi_(m,q,k,j,i);}
    for (int q = 0; q < n; ++q) {vi_(m,q,k,j,i) = buf[mp_(m,k,j,q)];}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdSweep
//! \brief one sweep (inward then outward) of every direction

void RadiationM1::VetGdSweep() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int n1 = indcs.nx1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int lcut = vlat_icut;
  auto cs_ = vlat_cs;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto cnt_ = vlat_cnt;
  auto mr_ = vgd_mr;
  auto cx1v = pmy_pack->pcoord->x1v;
  auto cx2v = pmy_pack->pcoord->x2v;
  auto cx3v = pmy_pack->pcoord->x3v;
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const int jlo = js - ng, jhi = je + ng - 1;
  const int klo = ks - ng, khi = ke + ng - 1;
  const Real twopi = 2.0*M_PI;
  for (int pass = 0; pass < 2; ++pass) {
    const bool inw = (pass == 0);
    for (int q = 0; q < n1 - lcut; ++q) {
      const int l = inw ? (n1 - 1 - q) : (lcut + q);
      const int i = is + l;
      par_for("m1_vgd_shell", DevExeSpace(), 0, nmb1, ks, ke, js, je, 0, n - 1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int d) {
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        const Real th = cx2v(m,j), ph = cx3v(m,k);
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        const Real mr = nx*st*cp + ny*st*sp + nz*ct;
        if (inw == (mr >= 0.0)) {return;}
        const Real r = cx1v(m,i);
        const Real zr = r*fabs(mr);
        const Real p2 = fmax(r*r - zr*zr, 0.0);
        const Real ch0 = exp(cs_(m,0,k,j,i));
        const Real s0 = exp(cs_(m,1,k,j,i));
        Real iv;
        if (inw && l == n1 - 1) {
          // from the vacuum top face (the end cell's S extrapolated, as vet_col_lat)
          const Real rtop = cx1f(m,ie+1);
          const Real ds = sqrt(fmax(rtop*rtop - p2, 0.0)) - zr;
          const Real stp = fmax(1.5*s0 - 0.5*exp(cs_(m,1,k,j,i-1)), 1.0e-300);
          Real ex, w0, wu;
          VgdW(ch0*ds, ex, w0, wu);
          iv = wu*stp + w0*s0;
        } else if (!inw && zr < 1.0e-12*r) {
          iv = s0;                     // tangent at the cell centre (measure zero)
        } else {
          int typ = 1;                 // 1 upwind shell, 2 inner face, 3 own shell
          int iu = i + 1;
          Real ds;
          if (inw) {
            const Real ru = cx1v(m,i+1);
            ds = sqrt(ru*ru - p2) - zr;
          } else {
            const Real rd = (l > lcut) ? cx1v(m,i-1) : cx1f(m,i);
            if (l > lcut && p2 <= rd*rd*(1.0 + 1.0e-13)) {
              iu = i - 1;
              ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
            } else if (l == lcut && p2 <= rd*rd) {
              typ = 2;
              iu = i;
              ds = zr - sqrt(fmax(rd*rd - p2, 0.0));
            } else {
              typ = 3;
              iu = i;
              ds = 2.0*zr;
            }
          }
          // the upwind point, its lateral cell and bilinear fractions (index space)
          const Real xu = r*st*cp - ds*nx;
          const Real yu = r*st*sp - ds*ny;
          const Real zu = r*ct - ds*nz;
          const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
          const Real thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
          Real dph = atan2(yu, xu) - ph;
          dph -= twopi*floor((dph + M_PI)/twopi);
          const Real fj = j + (thu - th)/mbsize.d_view(m).dx2;
          const Real fk = k + dph/mbsize.d_view(m).dx3;
          int j0 = static_cast<int>(floor(fj));
          int k0 = static_cast<int>(floor(fk));
          Real uj = fj - j0, uk = fk - k0;
          bool clp = false;
          if (j0 < jlo) {j0 = jlo; uj = 0.0; clp = true;}
          if (j0 > jhi) {j0 = jhi; uj = 1.0; clp = true;}
          if (k0 < klo) {k0 = klo; uk = 0.0; clp = true;}
          if (k0 > khi) {k0 = khi; uk = 1.0; clp = true;}
          if (clp) {Kokkos::atomic_add(&cnt_(0), 1.0);}
          const int j1 = j0 + 1, k1 = k0 + 1;
          if (typ == 2) {
            // the diffusion intensity at the inner face of the first shell
            const Real rf = cx1f(m,i);
            const Real muf = sqrt(fmax(1.0 - p2/(rf*rf), 0.0));
            const int i1 = i + 1;
            auto ib = [&](const int kk, const int jj) {
              const Real sa = exp(cs_(m,1,kk,jj,i)), sb = exp(cs_(m,1,kk,jj,i1));
              const Real g0 = -(sb - sa)/(cx1v(m,i1) - cx1v(m,i))/exp(cs_(m,0,kk,jj,i));
              return sa + g0*muf;
            };
            auto sbt = [&](const int kk, const int jj) {
              return fmax(1.5*exp(cs_(m,1,kk,jj,i)) - 0.5*exp(cs_(m,1,kk,jj,i1)),
                          1.0e-300);
            };
            Real ex, w0, wu;
            VgdW(ch0*ds, ex, w0, wu);
            const Real bv = fmax(VgdLerp(VgdLerp(ib(k0,j0), ib(k0,j1), uj),
                                         VgdLerp(ib(k1,j0), ib(k1,j1), uj), uk), 0.0);
            const Real qv = VgdLerp(VgdLerp(sbt(k0,j0), sbt(k0,j1), uj),
                                    VgdLerp(sbt(k1,j0), sbt(k1,j1), uj), uk);
            iv = bv*ex + wu*qv + w0*s0;
          } else {
            // each (cell, direction) holds ONE intensity, of the pass its n_d . r_hat
            // selects; the ray read here is inward at the upwind point for the inward
            // pass and for a turned ray, outward otherwise.  Stencil cells of the other
            // pass hold the wrong branch (near the tangent band): they are dropped and
            // the bilinear weights renormalised (all valid: the plain bilinear)
            const bool want_in = inw || (typ == 3);
            bool ok[4];
            bool all = true;
            for (int q = 0; q < 4; ++q) {
              const int kk = (q < 2) ? k0 : k1;
              const int jj = (q % 2 == 0) ? j0 : j1;
              const Real mq = mr_(m,kk,jj,d);
              ok[q] = want_in ? (mq < 0.0) : (mq >= 0.0);
              all = all && ok[q];
            }
            Real ivu;
            if (all) {
              ivu = VgdLerp(VgdLerp(vi_(m,d,k0,j0,iu), vi_(m,d,k0,j1,iu), uj),
                            VgdLerp(vi_(m,d,k1,j0,iu), vi_(m,d,k1,j1,iu), uj), uk);
            } else {
              Real ws = 0.0, vs = 0.0;
              for (int q = 0; q < 4; ++q) {
                const int kk = (q < 2) ? k0 : k1;
                const int jj = (q % 2 == 0) ? j0 : j1;
                const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                *((q % 2 == 0) ? (1.0 - uj) : uj);
                if (ok[q] && wq > 0.0) {ws += wq; vs += wq*vi_(m,d,kk,jj,iu);}
              }
              ivu = (ws > 0.0) ? vs/ws : -1.0;   // none: the upwind source (below)
            }
            auto lin = [&](const int c) {
              return VgdLerp(VgdLerp(cs_(m,c,k0,j0,iu), cs_(m,c,k0,j1,iu), uj),
                             VgdLerp(cs_(m,c,k1,j0,iu), cs_(m,c,k1,j1,iu), uj), uk);
            };
            const Real chu = exp(lin(0));
            const Real su = exp(lin(1));
            Real ex, w0, wu;
            VgdW(0.5*(chu + ch0)*ds, ex, w0, wu);
            iv = fmax(((ivu < 0.0) ? su : ivu)*ex + wu*su + w0*s0, 0.0);
          }
        }
        vi_(m,d,k,j,i) = iv;
      });
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdMoments
//! \brief D = K/J in the local (r, theta, phi) frame -> tau_ten LAT slots (relative to
//! vet_col's uniaxial tensor of this build; VetLatBuild folds LAT0 into slot 0)

void RadiationM1::VetGdMoments() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n = vgd_n;
  const int ilo = is + vlat_icut;
  auto vi_ = vgd_i;
  auto dir_ = vgd_dir;
  auto tt_ = tau_ten;
  auto cx2v = pmy_pack->pcoord->x2v;
  auto cx3v = pmy_pack->pcoord->x3v;
  par_for("m1_vgd_mom", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int c = 0; c < M1_TT_NLAT; ++c) {tt_(m,M1_TT_LAT0+c,k,j,i) = 0.0;}
    if (i < ilo) {return;}
    const Real th = cx2v(m,j), ph = cx3v(m,k);
    const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
    Real jm = 0.0, rr = 0.0, rt = 0.0, rp = 0.0, tq = 0.0, tp = 0.0, pq = 0.0;
    for (int d = 0; d < n; ++d) {
      const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
      const Real a = nx*st*cp + ny*st*sp + nz*ct;
      const Real b = nx*ct*cp + ny*ct*sp - nz*st;
      const Real c = -nx*sp + ny*cp;
      const Real wi = dir_(d,3)*vi_(m,d,k,j,i);
      jm += wi;
      rr += wi*a*a;
      rt += wi*a*b;
      rp += wi*a*c;
      tq += wi*b*b;
      tp += wi*b*c;
      pq += wi*c*c;
    }
    if (jm > 0.0) {
      const Real fk = tt_(m,0,k,j,i);
      const Real ht = 0.5*(1.0 - fk);
      tt_(m,M1_TT_LAT0,k,j,i) = rr/jm - fk;
      tt_(m,M1_TT_LAT0+1,k,j,i) = rt/jm;
      tt_(m,M1_TT_LAT0+2,k,j,i) = rp/jm;
      tt_(m,M1_TT_LAT0+3,k,j,i) = tq/jm - ht;
      tt_(m,M1_TT_LAT0+4,k,j,i) = tp/jm;
      tt_(m,M1_TT_LAT0+5,k,j,i) = pq/jm - ht;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::VetGdBuild
//! \brief the gd build inside VetLatBuild: source and first shell (VetLatSweep(0)), the
//! sweep(s) with the lagged lateral inflow, the moments; timed per part (fenced)

void RadiationM1::VetGdBuild() {
  Kokkos::Timer tm;
  VetLatSweep(0);
  Kokkos::fence();
  vgd_tsrc += tm.seconds();
  const int nit = (vlat_nbuild == 0) ? vlat_iinit : 1;
  for (int it = 0; it < nit; ++it) {
    tm.reset();
    VetGdSweep();
    Kokkos::fence();
    vgd_tswp += tm.seconds();
    tm.reset();
    VetLatExchange(vgd_i, vgd_i_c, pbval_gd);
    VetGdWall();
    Kokkos::fence();
    vgd_texc += tm.seconds();
    vlat_ncall += 1.0;
  }
  tm.reset();
  VetGdMoments();
  Kokkos::fence();
  vgd_tmom += tm.seconds();
  auto cnt_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vlat_cnt);
  vlat_nclamp = cnt_h(0);
}

} // namespace radm1
