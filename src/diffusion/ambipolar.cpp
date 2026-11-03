//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ambipolar.cpp
//! \brief Ambipolar diffusion (ambipolar-1101): the edge EMF
//!     E_AD = eta_A (J - (J.B) B / B^2),   eta_A = ad_coef |B|^2,
//! added to the resistive edge EMF (efld_resist) before the walls and the Poynting flux
//! are formed from it.  The explicit stage, the RKG super-step (both go through
//! AddFluxGeneralResist), the conducting / stress-free walls and the dhj x1 energy
//! ledger therefore see the total Ohmic + ambipolar EMF without further changes.
//!
//! Stencils follow Athena++ FieldDiffusion::AmbipolarEMF: J along the edge from the
//! curl, the two other J components averaged over the 2 (2-D) or 4 neighbouring edges;
//! B along the edge from the face field averaged over the 2 adjacent faces, the two other
//! components from the face field (2 faces) and the cell-centred field (4 cells); eta_A
//! averaged over the 2 or 4 cells around the edge.
//!
//! Spherical polar: the faces carry orthonormal (r, theta, phi) components and the curl
//! of current_density.hpp returns orthonormal components too, so the same algebra holds.
//! The polar ghost row is the far side of the axis: the J2 average of the edges ON the
//! axis (x1 and x3 edges of the polar face) uses the in-domain row only, and J3 on the
//! phi edges of the polar face (zero length, phi undefined there) is replaced by J3 of
//! the next face in wherever an average would read it (x1 edges on the axis, x2 edges of
//! the polar row).  The ghost B and
//! eta_A across the axis are the polar mirror (B_theta, B_phi sign-flipped), i.e. the
//! same vector in the local basis, and enter the averages as usual.

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "resistivity.hpp"
#include "current_density.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::SetAmbipolar()
//! \brief <mhd>/ambipolar = eos: ad_coef = eta_A/B^2 from the EOS ionization.
//!
//! In Gaussian units eta_A = B_G^2/(4 pi gamma_in rho_i rho_n) [cm^2/s].  The code field
//! has B^2/2 = magnetic pressure, so B_G^2 = 4 pi pres_cgs B^2 and
//!     ad_coef = pres_cgs/(gamma_in rho_i rho_n)       (rho in g/cm^3),
//! returned in cm^2/s per code B^2, the same convention as the Ohmic eos path (eta in
//! cm^2/s, i.e. code units for length_cgs = time_cgs = 1, as in the dhj).
//!   x_e    = n_e/n_tot from the same Saha table the Ohmic eos path reads
//!   rho_i  = x_e n_tot m_ion,  m_ion = <mhd>/ad_mion_amu amu (default 30: a Na+/K+/Fe+
//!            mix; Na 23, K 39, Fe 56, Mg 24, Ca 40, Al 27, H+ would be 1)
//!   rho_n  = rho - rho_i,  m_n = mu m_u (mu of the whole gas; differs by O(x_e))
//!   gamma_in = <sigma v>_in/(m_ion + m_n), <sigma v>_in = <mhd>/ad_sigv_in (Langevin,
//!            default 1.9e-9 cm^3/s, Draine, Roberge & Dalgarno 1983)

void Resistivity::SetAmbipolar(const DvceArray5D<Real> &w, const EOS_Data &eos,
                               const Real &Rgas, const int il, const int iu,
                               const int jl, const int ju, const int kl, const int ku) {
  int &nmb = pmy_pack->nmb_thispack;
  const Real gm1 = eos.gamma - 1.0;
  const bool gen = eos.IsGeneral();
  if (!eos.HasElectronFraction()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/ambipolar = eos needs the EOS electron fraction (eos = general "
              << "with a table, or an ideal gas with the eos_* composition table), as "
              << "ohmic_resistivity = eos does" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!gen && !(Rgas > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/ambipolar = eos on an ideal gas needs <problem>/Rgas"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &wtemp_ = pmy_pack->pmhd->wtemp;
  auto coef_ = ad_coef;
  const Real amu = 1.66053906660e-24;
  const Real mion = ad_mion_amu*amu;
  const Real sigv = ad_sigv_in;
  const int ixe = xe_idx;   // xetracer-1101: advected y = n_e/n_H when >= 0
  const Real nhfac = xe_nhfac;
  par_for("ad_coef", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real rho = w(m,IDN,k,j,i);
    const Real e = w(m,IEN,k,j,i);
    const Real kb = 1.380649e-16;
    const Real mh = 1.67262192369e-24;
    Real T, ntot, mu, rhoc, pcgs;
    if (gen) {
      const Real tcode = wtemp_(m,k,j,i);
      T = tcode*eos.temp_cgs;
      mu = eos.MeanMolecularWeight(rho,e,tcode);
      ntot = rho*eos.pres_cgs/(mu*eos.temp_cgs*kb);
      rhoc = rho*eos.dens_cgs;
      pcgs = eos.pres_cgs;
    } else {
      T = gm1*e/Rgas/rho;
      mu = kb/mh/Rgas;
      ntot = rho/(mu*mh);
      rhoc = rho;
      pcgs = 1.0;
    }
    // n_e = y n_H with the tracer, x_e n_tot with local Saha
    const Real ne = (ixe >= 0) ? w(m,ixe,k,j,i)*nhfac*rhoc
                               : eos.ElectronFractionKelvin(rho, T)*ntot;
    Real rhoi = ne*mion;
    rhoi = fmin(rhoi, 0.5*rhoc);
    const Real rhon = rhoc - rhoi;
    const Real gin = sigv/(mion + mu*amu);
    const Real den = gin*rhoi*rhon;
    // x_e = 0 (fully neutral): eta_A -> infinity, left to max_eta_ad; 1e300 keeps the
    // product with B^2 = 0 at 0 rather than NaN
    coef_(m,k,j,i) = (den > pcgs*1.0e-300) ? pcgs/den : 1.0e300;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::AddEMFAmbipolar()
//! \brief adds eta_A (J - (J.B)B/B^2) to the edge EMF efld (efld_resist).  bc is the
//! cell-centred field consistent with b0 (bcc0).  Cartesian and spherical polar; the
//! cubed sphere has its own routine (AddEMFAmbipolarGnomonic).

void Resistivity::AddEMFAmbipolar(const DvceFaceFld4D<Real> &b0,
                                  const DvceArray5D<Real> &bc,
                                  DvceEdgeFld4D<Real> &efld) {
  if (pmy_pack->pmesh->use_cubed_sphere) {
    AddEMFAmbipolarGnomonic(b0, bc, efld);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ncells1 = indcs.nx1 + 2*(indcs.ng);
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool one_d = pmy_pack->pmesh->one_d;
  const bool two_d = pmy_pack->pmesh->two_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool sph = pmy_pack->pmesh->use_spherical_polar;
  const auto geom = MakeCurrentDensityGeom(pmy_pack);
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mb_bcs_ = pmy_pack->pmb->mb_bcs;
  auto &x2v_ = pmy_pack->pcoord->x2v;
  auto &x2f_ = pmy_pack->pcoord->xx2f;
  auto jc1 = jcur.x1e, jc2 = jcur.x2e, jc3 = jcur.x3e;

  // (1) J on every edge the stencils below touch: i in [is-1, ie+1], j in [js-1, je+1]
  //     (multi-D), k in [ks-1, ke+1] (3-D)
  const int kl = three_d ? ks-1 : ks, ku = three_d ? ke+1 : ks;
  const int jl = one_d ? js : js-1, ju = one_d ? js : je+1;
  {
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;
    int scr_level = TeamScratchLevel(scr_size, 0);
    par_for_outer("ad_j", DevExeSpace(), scr_size, scr_level, 0, nmb1, kl, ku, jl, ju,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);
      int pole_edge = 0;
      Real pole_c = 0.0, pole_c1 = 0.0, pole_s = 1.0, pole_a = 1.0, pole_l = 1.0;
      if (three_d) {
        // the pole x1 edge, exactly as AddEMFGeneralResist forms it
        const bool pn = (mb_bcs_.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar);
        const bool ps = (mb_bcs_.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
        pole_edge = (j == js && pn) ? 1 : ((j == je+1 && ps) ? 2 : 0);
        if (pole_edge != 0) {
          const Real av = (pole_edge == 1) ? x2v_(m,js) - x2f_(m,js)
                                           : x2f_(m,je+1) - x2v_(m,je);
          const Real dth = (pole_edge == 1) ? x2f_(m,js+1) - x2f_(m,js)
                                            : x2f_(m,je+1) - x2f_(m,je);
          const Real a = 0.5*dth;
          pole_c = (a*a/3.0 - dth*dth)/(3.0*dth*dth);
          pole_c1 = a/(4.0*dth);
          pole_s = sin(a)/sin(av);
          pole_a = (1.0 - cos(a))/(1.0 - cos(av));
          pole_l = a/av;
        }
      }
      CurrentDensity(geom, member, m, k, j, is-1, ie+1, b0, mbsize.d_view(m),
                     j1, j2, j3, pole_edge, pole_c, pole_c1, pole_s, pole_a, pole_l);
      par_for_inner(member, is-1, ie+1, [&](const int i) {
        jc1(m,k,j,i) = j1(i);
        jc2(m,k,j,i) = j2(i);
        jc3(m,k,j,i) = j3(i);
      });
    });
  }

  // (2) eta_A (J - (J.B)B/B^2) on the edges
  auto coef_ = ad_coef;
  const Real cap = max_eta_ad;
  auto e1 = efld.x1e, e2 = efld.x2e, e3 = efld.x3e;
  auto b1f = b0.x1f, b2f = b0.x2f, b3f = b0.x3f;

  if (one_d) {
    par_for("ad_emf1", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int i) {
      auto eta = [&](int ii) {
        const Real b2 = SQR(bc(m,IBX,ks,js,ii)) + SQR(bc(m,IBY,ks,js,ii))
                      + SQR(bc(m,IBZ,ks,js,ii));
        return fmin(coef_(m,ks,js,ii)*b2, cap);
      };
      const Real eta_a = 0.5*(eta(i-1) + eta(i));
      const Real bx = b1f(m,ks,js,i);
      const Real by = 0.5*(b2f(m,ks,js,i) + b2f(m,ks,js,i-1));
      const Real bz = 0.5*(b3f(m,ks,js,i) + b3f(m,ks,js,i-1));
      const Real bsq = SQR(bx) + SQR(by) + SQR(bz) + 1.0e-300;
      const Real jdb = jc2(m,ks,js,i)*by + jc3(m,ks,js,i)*bz;
      const Real de2 = eta_a*(jc2(m,ks,js,i) - jdb*by/bsq);
      const Real de3 = eta_a*(jc3(m,ks,js,i) - jdb*bz/bsq);
      e2(m,ks,  js,i) += de2;
      e2(m,ke+1,js,i) += de2;
      e3(m,ks,js,  i) += de3;
      e3(m,ks,je+1,i) += de3;
    });
    return;
  }

  // k-neighbour layers: in 2-D the single x3 layer stands for both
  const int dk = three_d ? 1 : 0;
  par_for("ad_emf", DevExeSpace(), 0, nmb1, ks, ks+(ke+1-ks)*dk, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    auto eta = [&](int kk, int jj, int ii) {
      const Real b2 = SQR(bc(m,IBX,kk,jj,ii)) + SQR(bc(m,IBY,kk,jj,ii))
                    + SQR(bc(m,IBZ,kk,jj,ii));
      return fmin(coef_(m,kk,jj,ii)*b2, cap);
    };
    const int km = k - dk;   // lower x3 cell / edge layer (= k in 2-D)
    const int kp = k + dk;   // upper x3 edge layer for the x3-edge averages
    // spherical polar: is this j-face the polar axis?  The J2 row across it is not used.
    const bool pin2 = sph &&
                      (mb_bcs_.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar);
    const bool pout2 = sph &&
                       (mb_bcs_.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar);
    const bool pin = pin2 && (j == js);
    const bool pout = pout2 && (j == je+1);
    // J2 at x2-edge row pair (j-1, j) at x3 layer kk and x1 face ii
    auto j2pair = [&](int kk, int ii) {
      if (pin) return jc2(m,kk,j,ii);
      if (pout) return jc2(m,kk,j-1,ii);
      return 0.5*(jc2(m,kk,j,ii) + jc2(m,kk,j-1,ii));
    };
    // J3 on the x3 edges of j-face jj; the phi edges ON the axis (polar face) carry no
    // usable J_phi (zero length, phi undefined), so the next face in is used instead
    // a PHYSICAL radial boundary face (not block, not periodic): the J1 of the x1 edges
    // in the radial ghost column depends on how faithfully the boundary condition fills
    // the corner ghosts (at a polar corner it did not: E2 off by 0.1 |eta_A J| on the
    // sp_test force-free field), so the x1-edge average on that face uses the active
    // column only
    auto phys = [&](int f) {
      const BoundaryFlag bf = mb_bcs_.d_view(m,f);
      return (bf != BoundaryFlag::block) && (bf != BoundaryFlag::periodic) &&
             (bf != BoundaryFlag::shear_periodic);
    };
    const int ia = (i == is && phys(BoundaryFace::inner_x1)) ? i : i-1;
    const int ib = (i == ie+1 && phys(BoundaryFace::outer_x1)) ? i-1 : i;
    auto jf3 = [&](int jj) {
      if (pin2 && jj == js) return js+1;
      if (pout2 && jj == je+1) return je;
      return jj;
    };

    // ---- x1 edge (k-face, j-face, i-centre): i <= ie
    if (i <= ie && (three_d ? true : k == ks)) {
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j-1,i) + eta(km,j,i) + eta(km,j-1,i));
      const Real jx = jc1(m,k,j,i);
      const Real jy = 0.5*(j2pair(k,i) + j2pair(k,i+1));
      const int j3 = jf3(j);
      const Real jz = 0.25*(jc3(m,k,j3,i) + jc3(m,k,j3,i+1)
                          + jc3(m,km,j3,i) + jc3(m,km,j3,i+1));
      const Real bx = 0.25*(bc(m,IBX,k,j,i) + bc(m,IBX,k,j-1,i)
                          + bc(m,IBX,km,j,i) + bc(m,IBX,km,j-1,i));
      const Real by = 0.5*(b2f(m,k,j,i) + b2f(m,km,j,i));
      const Real bz = 0.5*(b3f(m,k,j,i) + b3f(m,k,j-1,i));
      const Real bsq = SQR(bx) + SQR(by) + SQR(bz) + 1.0e-300;
      const Real jdb = jx*bx + jy*by + jz*bz;
      const Real de = eta_a*(jx - jdb*bx/bsq);
      e1(m,k,j,i) += de;
      if (two_d) e1(m,ke+1,j,i) += de;
    }
    // ---- x2 edge (k-face, j-centre, i-face): j <= je
    if (j <= je) {
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j,i-1) + eta(km,j,i) + eta(km,j,i-1));
      const Real jx = 0.25*(jc1(m,k,j,ib) + jc1(m,k,j,ia)
                          + jc1(m,k,j+1,ib) + jc1(m,k,j+1,ia));
      const Real jy = jc2(m,k,j,i);
      const int ja = jf3(j), jb = jf3(j+1);
      const Real jz = 0.25*(jc3(m,k,ja,i) + jc3(m,k,jb,i)
                          + jc3(m,km,ja,i) + jc3(m,km,jb,i));
      const Real bx = 0.5*(b1f(m,k,j,i) + b1f(m,km,j,i));
      const Real by = 0.25*(bc(m,IBY,k,j,i) + bc(m,IBY,k,j,i-1)
                          + bc(m,IBY,km,j,i) + bc(m,IBY,km,j,i-1));
      const Real bz = 0.5*(b3f(m,k,j,i) + b3f(m,k,j,i-1));
      const Real bsq = SQR(bx) + SQR(by) + SQR(bz) + 1.0e-300;
      const Real jdb = jx*bx + jy*by + jz*bz;
      const Real de = eta_a*(jy - jdb*by/bsq);
      e2(m,k,j,i) += de;
      if (two_d) e2(m,ke+1,j,i) += de;
    }
    // ---- x3 edge (k-centre, j-face, i-face): k <= ke
    if (k <= ke) {
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j,i-1) + eta(k,j-1,i) + eta(k,j-1,i-1));
      const Real jx = 0.25*(jc1(m,k,j,ib) + jc1(m,k,j,ia)
                          + jc1(m,kp,j,ib) + jc1(m,kp,j,ia));
      const Real jy = 0.5*(j2pair(k,i) + j2pair(kp,i));
      const Real jz = jc3(m,k,j,i);
      const Real bx = 0.5*(b1f(m,k,j,i) + b1f(m,k,j-1,i));
      const Real by = 0.5*(b2f(m,k,j,i) + b2f(m,k,j,i-1));
      const Real bz = 0.25*(bc(m,IBZ,k,j,i) + bc(m,IBZ,k,j,i-1)
                          + bc(m,IBZ,k,j-1,i) + bc(m,IBZ,k,j-1,i-1));
      const Real bsq = SQR(bx) + SQR(by) + SQR(bz) + 1.0e-300;
      const Real jdb = jx*bx + jy*by + jz*bz;
      e3(m,k,j,i) += eta_a*(jz - jdb*bz/bsq);
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::AddEMFAmbipolarGnomonic()
//! \brief cubed sphere.  The angular basis is NOT orthogonal: e_xi.e_eta = c (c = cos of
//! the face angle, s = sin), metric g = [[1, c], [c, 1]].  The three kinds of stored
//! component relate as
//!   face-normal   N_a = V.nhat_a = s V^a          (b0.x2f/x3f, jnorm.x2e/x3e)
//!   covariant     V_a = V.e_a = V^a + c V^b       (the edge EMFs, bcc0 IBY = B_xi)
//!   contravariant V^a = N_a/s
//! and bcc0 IBZ = B.nhat_eta = N_eta (frame {rhat, e_xi, nhat_eta}, orthonormal).  At
//! every edge both J and B are formed in contravariant AND covariant components at the
//! edge, so  J.B = J_r B_r + J^a B_a,  B^2 = B_r^2 + B^a B_a,  and the stored EMF is the
//! covariant component along the edge:  E_a += eta_A (J_a - (J.B) B_a/B^2).
//! J in the face-normal frame comes from the Ohmic two-pass curl (jnorm, filled by
//! AddEMFGnomonicResist just before; its radial-face weights on a stretched grid
//! included), plus the radial-edge layer at i = is-1 and ie+1 the x2/x3 edge averages
//! need.  The averages follow the Cartesian stencils; the halo supplies the seam
//! neighbours in the local panel frame, as for the Ohmic EMF and Poynting flux.

void Resistivity::AddEMFAmbipolarGnomonic(const DvceFaceFld4D<Real> &b0,
                                          const DvceArray5D<Real> &bc,
                                          DvceEdgeFld4D<Real> &efld) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &cxi = pmy_pack->pcoord->cos_face_xi;
  auto &sxi = pmy_pack->pcoord->sin_face_xi;
  auto &cet = pmy_pack->pcoord->cos_face_eta;
  auto &set = pmy_pack->pcoord->sin_face_eta;
  auto &ccl = pmy_pack->pcoord->cos_cell;
  auto &scl = pmy_pack->pcoord->sin_cell;
  auto &dxf2 = pmy_pack->pcoord->dxface.x2f;
  auto &dxf3 = pmy_pack->pcoord->dxface.x3f;
  auto &ae1 = pmy_pack->pcoord->areaedge.x1e;
  auto b = b0;
  auto jn1 = jnorm.x1e, jn2 = jnorm.x2e, jn3 = jnorm.x3e;

  // the radial-edge current one cell beyond the active range in x1 (the formula of
  // res_cs_j1, which has no radial derivative)
  par_for("ad_cs_j1", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, 0, 1,
  KOKKOS_LAMBDA(int m, int k, int j, int n) {
    const int i = (n == 0) ? is-1 : ie+1;
    auto bxi = [&](int kk, int jj) {   // B.e_xi at the x2 face (kk, jj, i)
      const Real b3a = 0.25*(b.x3f(m,kk,jj,i) + b.x3f(m,kk+1,jj,i)
                           + b.x3f(m,kk,jj-1,i) + b.x3f(m,kk+1,jj-1,i));
      return (b.x2f(m,kk,jj,i) + cxi(m,kk,jj)*b3a)/sxi(m,kk,jj);
    };
    auto bet = [&](int kk, int jj) {   // B.e_eta at the x3 face (kk, jj, i)
      const Real b2a = 0.25*(b.x2f(m,kk,jj,i) + b.x2f(m,kk,jj+1,i)
                           + b.x2f(m,kk-1,jj,i) + b.x2f(m,kk-1,jj+1,i));
      return (b.x3f(m,kk,jj,i) + cet(m,kk,jj)*b2a)/set(m,kk,jj);
    };
    jn1(m,k,j,i) = (dxf3(m,k,j,i)*bet(k,j) - dxf3(m,k,j-1,i)*bet(k,j-1)
                  - dxf2(m,k,j,i)*bxi(k,j) + dxf2(m,k-1,j,i)*bxi(k-1,j))/ae1(m,k,j,i);
  });

  auto coef_ = ad_coef;
  const Real cap = max_eta_ad;
  auto e1 = efld.x1e, e2 = efld.x2e, e3 = efld.x3e;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  // r1mspeed-1102: the cell eta_A, evaluated once per cell (the edge loop below reads
  // each cell's value 12 times); the same expression, so bitwise
  if (ad_eta_c.extent_int(0) != ad_coef.extent_int(0) ||
      ad_eta_c.extent_int(3) != ad_coef.extent_int(3)) {
    ad_eta_c = DvceArray4D<Real>("ad_eta_c", ad_coef.extent_int(0), ad_coef.extent_int(1),
                                 ad_coef.extent_int(2), ad_coef.extent_int(3));
  }
  auto etc_ = ad_eta_c;
  par_for("ad_cs_eta", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int kk, const int jj, const int ii) {
    const Real b2 = SQR(bc(m,IBX,kk,jj,ii)) + SQR(bc(m,IBY,kk,jj,ii))
                  + SQR(bc(m,IBZ,kk,jj,ii));
    etc_(m,kk,jj,ii) = fmin(coef_(m,kk,jj,ii)*b2, cap);
  });
  par_for("ad_cs_emf", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    auto eta = [&](int kk, int jj, int ii) {
      return etc_(m,kk,jj,ii);
    };
    // physical radial boundary face: the radial-edge J average uses the active column
    // (as the Cartesian/sp kernel; the Ohmic curl switches stencil there too)
    const bool rin = (mbbcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::block);
    const bool rout = (mbbcs.d_view(m,BoundaryFace::outer_x1) != BoundaryFlag::block);
    const int ia = (i == is && rin) ? i : i-1;
    const int ib = (i == ie+1 && rout) ? i-1 : i;
    // ---- x1 edge (r centre, xi face j, eta face k): along rhat
    if (i <= ie) {
      const Real c = 0.5*(cxi(m,k,j) + cxi(m,k-1,j));
      const Real sn = sqrt(1.0 - c*c);
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j-1,i) + eta(k-1,j,i) + eta(k-1,j-1,i));
      const Real br = 0.25*(bc(m,IBX,k,j,i) + bc(m,IBX,k,j-1,i)
                          + bc(m,IBX,k-1,j,i) + bc(m,IBX,k-1,j-1,i));
      const Real bu1 = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k-1,j,i))/sn;   // B^xi
      const Real bu2 = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k,j-1,i))/sn;   // B^eta
      const Real bd1 = bu1 + c*bu2, bd2 = c*bu1 + bu2;
      const Real jr = jn1(m,k,j,i);
      const Real ju1 = 0.25*(jn2(m,k,j,i) + jn2(m,k,j,i+1)
                           + jn2(m,k,j-1,i) + jn2(m,k,j-1,i+1))/sn;
      const Real ju2 = 0.25*(jn3(m,k,j,i) + jn3(m,k,j,i+1)
                           + jn3(m,k-1,j,i) + jn3(m,k-1,j,i+1))/sn;
      const Real bsq = br*br + bu1*bd1 + bu2*bd2 + 1.0e-300;
      const Real jdb = jr*br + ju1*bd1 + ju2*bd2;
      e1(m,k,j,i) += eta_a*(jr - jdb*br/bsq);
    }
    // ---- x2 edge (r face i, xi centre j, eta face k): along e_xi
    if (j <= je) {
      const Real c = cet(m,k,j), sn = set(m,k,j);
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j,i-1) + eta(k-1,j,i) + eta(k-1,j,i-1));
      const Real br = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k-1,j,i));
      const Real bu2 = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k,j,i-1))/sn;   // B^eta
      const Real bd1 = 0.25*(bc(m,IBY,k,j,i) + bc(m,IBY,k,j,i-1)     // B_xi
                           + bc(m,IBY,k-1,j,i) + bc(m,IBY,k-1,j,i-1));
      const Real bu1 = bd1 - c*bu2, bd2 = c*bu1 + bu2;
      const Real jr = 0.25*(jn1(m,k,j,ib) + jn1(m,k,j,ia)
                          + jn1(m,k,j+1,ib) + jn1(m,k,j+1,ia));
      const Real ju1 = jn2(m,k,j,i)/sn;
      const Real ju2 = 0.25*(jn3(m,k,j,i) + jn3(m,k,j+1,i)
                           + jn3(m,k-1,j,i) + jn3(m,k-1,j+1,i))/sn;
      const Real jd1 = ju1 + c*ju2;
      const Real bsq = br*br + bu1*bd1 + bu2*bd2 + 1.0e-300;
      const Real jdb = jr*br + ju1*bd1 + ju2*bd2;
      e2(m,k,j,i) += eta_a*(jd1 - jdb*bd1/bsq);
    }
    // ---- x3 edge (r face i, xi face j, eta centre k): along e_eta
    if (k <= ke) {
      const Real c = cxi(m,k,j), sn = sxi(m,k,j);
      const Real eta_a = 0.25*(eta(k,j,i) + eta(k,j-1,i) + eta(k,j,i-1) + eta(k,j-1,i-1));
      const Real br = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j-1,i));
      const Real bu1 = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j,i-1))/sn;   // B^xi
      // B.e_eta per cell = c_cell B_xi + s_cell N_eta, then averaged
      auto bde = [&](int jj, int ii) {
        return ccl(m,k,jj)*bc(m,IBY,k,jj,ii) + scl(m,k,jj)*bc(m,IBZ,k,jj,ii);
      };
      const Real bd2 = 0.25*(bde(j,i) + bde(j,i-1) + bde(j-1,i) + bde(j-1,i-1));
      const Real bu2 = bd2 - c*bu1, bd1 = bu1 + c*bu2;
      const Real jr = 0.25*(jn1(m,k,j,ib) + jn1(m,k,j,ia)
                          + jn1(m,k+1,j,ib) + jn1(m,k+1,j,ia));
      const Real ju2 = jn3(m,k,j,i)/sn;
      const Real ju1 = 0.25*(jn2(m,k,j,i) + jn2(m,k+1,j,i)
                           + jn2(m,k,j-1,i) + jn2(m,k+1,j-1,i))/sn;
      const Real jd2 = c*ju1 + ju2;
      const Real bsq = br*br + bu1*bd1 + bu2*bd2 + 1.0e-300;
      const Real jdb = jr*br + ju1*bd1 + ju2*bd2;
      e3(m,k,j,i) += eta_a*(jd2 - jdb*bd2/bsq);
    }
  });
  // ucthshare-1102 cs_seam_mhd = shared_g | shared_o: the seam edges take the TWO-CHART
  // MEAN of the Ohmic + ambipolar EMF (MHD::CsSeamResist, mean), single-valued without
  // the edge-EMF averaging (the Ohmic call skipped its mean so it is taken once, here)
  if (pmy_pack->pmhd != nullptr && pmy_pack->pmhd->cs_seam_mhd >= 3) {
    pmy_pack->pmhd->CsSeamResist(b0, eta_b, jnorm, efld, true);
  }
  return;
}
