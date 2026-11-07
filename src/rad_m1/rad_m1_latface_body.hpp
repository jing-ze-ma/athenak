//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_latface_body.hpp
//! \brief share-1031: the ONE body of RadiationM1::ImplicitLatFaceCoef (spherical
//! polar / Cartesian, M1_LFC_CS 0) and ImplicitLatFaceCoefCs (cubed sphere, M1_LFC_CS 1),
//! included twice from rad_m1_implicit.cpp (the rad_m1_impl_wb_kernel.hpp pattern).  The
//! cs-only parts (cs-s4-1017: the open seams m1bcs, the seam-face mirror-pair distance
//! cs_dxf_eff, the lateral ghost ring of the half-range cell data from the vgd_hr ghosts
//! instead of the per-direction halo exchanges, fs refused) sit in `#if M1_LFC_CS`
//! blocks, so that the sp instantiation is token for token the former
//! ImplicitLatFaceCoef (the same expression trees under nvcc --fmad=true: bitwise).
//! M1_LFC_SFX is the kernel-name suffix ("" or "cs").  No include guard: included once per variant.

  const bool khl_ = impl_kt_lat_harm;   // faceextrap-1106 M1KtLat
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto fw2 = ifw2;
  auto fw3 = ifw3;
  auto &mbsize = pmy_pack->pmb->mb_size;
#if M1_LFC_CS
  auto &mbbcs = m1bcs;   // cs-s4-1017: the seams are open
  auto cdxf = M1Dxf();
  // cs-s4-1017: a seam face takes the mirror-pair distance r angm (cs_dxf_eff)
  auto cdxe = cs_dxf_eff;
  auto cse_ = cs_seam.d_view;
#else
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  auto cdxf = M1Dxf();
#endif
  const bool sph = sph_geom;
  const bool edd = eddington;
  const bool blend = (impl_flux == M1_IFLUX_BLEND);
  const int bkind = impl_blend, bfm = impl_blend_fmode;
  const Real btau0 = impl_blend_tau0;
  const Real bflo = impl_blend_flo, bfhi = impl_blend_fhi;
  const Real ch = chat;
  const Real balph = impl_blend_alpha, br0 = impl_blend_r0;
  const Real bx0 = impl_blend_xthin, dtl = dt_sub;
  const int bxmode = impl_blend_xthin_mode;
  const Real bxwmin = impl_blend_xthin_wmin, bxr0 = impl_blend_xthin_r0;
  if (impl_beam_hr) {
    // hrup-1009: per direction d, h+_d -> S1 and h-_d -> S3 of every active cell
    // (M1HrCell), both through the transverse halo, then the faces (M1HrFace)
    const bool hgd = vgd_on && sph && (vgd_hr.extent_int(0) > 0);
    const int hsrc = vet_sc ? 1 : (hgd ? 2 : (impl_hr_model ? 0 : -1));
    auto hvc_ = vet_cell;
    auto hgh_ = vgd_hr;
    // thinsw-1010 (implicit_thin_switch): the faces are built at the first pass of the
    // solve only (every input is fixed for it); the cell thinness g rides in RES, the
    // formal-solution flux H_d in TA, E^n in TB and J in TC (free at this point)
    const int tmode = impl_thsw;
    const bool thsw = (tmode != M1_THSW_NONE);
    if (thsw && !first) {return;}
    // gdthin-1011: vet_gd (src 2) with C1, eps, J in vgd_hr takes the correction too
    const bool tgx = (hsrc == 2) && (vgd_hr.extent_int(1) >= M1_VGD_HRNX);
    const int tcorr = (hsrc == 1 || tgx) ? impl_thsw_corr : M1_THCR_NONE;
    const bool tfs = (impl_thsw == M1_THSW_FS || impl_thsw == M1_THSW_FSKN ||
                    impl_thsw == M1_THSW_FSQ || impl_thsw == M1_THSW_FSI ||
                    impl_thsw == M1_THSW_FSE);
    const Real tfsa = thsw_fsa, tfsj = thsw_fsj, th0 = thsw_h0, th1 = thsw_h1;
    const Real tkn0 = thsw_kn0, tknp = thsw_knp;
    const Real tqp = (impl_thsw == M1_THSW_FSQ) ? thsw_qp
                     : ((impl_thsw == M1_THSW_FSE) ? -1.0 : 0.0);
    const Real tqi = thsw_qp;
    // thinfix-1011; unify-1012: also vet_gd (src 2) with J, eps in vgd_hr (tgx)
    const bool tjs = thsw && thsw_js && (hsrc == 1 || (tgx && (thsw_jsgd & 2)));
    const bool tjse = tjs && thsw_jse;
    const Real tjl = thsw_jslim;
    const bool tjm = thsw_jsmean;   // jsfix-1012
    // cbfix-1011; unify-sc-gd-1012 (#4): vet_gd on implicit_hr_face_gd = all
    const int hsfb = (impl_hr_sf && (hsrc == 1 || (hsrc == 2 && (impl_hr_sfgd & 2))))
                     ? iw_hsf : -1;
    const bool hcjh = thsw_hcjh;   // unify-1012
    const bool hxs = impl_hr_xs && (hsrc == 2) && (vgd_hr.extent_int(1) >= M1_VGD_HRNX);
    const Real xcfl = impl_hr_xs_cfl;
#if M1_LFC_CS
    // cs-s4-1017: the per-cell beam data of the LATERAL GHOST RING too, from the ghosts
    // of vgd_hr (VetGdCsHrPrep: the seam quad map, this block's orientation) and of E^n,
    // instead of exchanging the per-direction S1/S3/RES/TA.. (a scalar exchange cannot
    // turn a neighbour panel's x3 data into this block's x2 data, nor swap h+/h- on a
    // reversed seam); |H|/J from the cs component of vgd_hr (cs metric)
    const int cfb = vgd_hr.extent_int(1) - 1;
    ImplicitHaloExchange(1, M1_IW_EN);
#endif
    for (int d = 1; d <= (trans_x3 ? 2 : 1); ++d) {
#if M1_LFC_CS
      const int jl0 = (d == 1) ? js - 1 : js, ju0 = (d == 1) ? je + 1 : je;
      const int kl0 = (d == 2) ? ks - 1 : ks, ku0 = (d == 2) ? ke + 1 : ke;
      par_for("m1_impl_lfchrcs", DevExeSpace(), 0, nmb1, kl0, ku0, jl0, ju0, is, ie,
#else
      par_for("m1_impl_lfchr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
#endif
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real hp, hm, fb;
        M1HrCell(hsrc, iw_, hvc_, hgh_, m, k, j, i, d, 0.0, hp, hm, fb);
#if M1_LFC_CS
        if (hsrc == 2) {fb = hgh_(m,cfb,k,j,i);}
#endif
        if (hxs) {   // sptests-1011 P2
          M1HrExcess(iw_, hgh_, m, k, j, i, ch*dtl/mbsize.d_view(m).dx1, xcfl, hp, hm);
        }
        iw_(m,M1_IW_S1,k,j,i) = hp;
        iw_(m,M1_IW_S3,k,j,i) = hm;
        if (thsw) {
          Real hd, jo;
          if (tmode == M1_THSW_FSI && hsrc == 1) {
            M1ThinHr(tqi, iw_(m,M1_IW_EN,k,j,i), fmax(hvc_(m,M1_VET_J,k,j,i), 0.0), fb,
                     th0, th1, hp, hm);
            iw_(m,M1_IW_S1,k,j,i) = hp;
            iw_(m,M1_IW_S3,k,j,i) = hm;
          }
          if (hsrc == 2) {   // gdthin-1011
            fb = M1ThinCellGd(hgh_, m, k, j, i, d, fb, tfs, tqp, iw_(m,M1_IW_EN,k,j,i),
                              hd, jo);
          } else {
          fb = M1ThinCell(hsrc, hvc_, m, k, j, i, d, fb, tfs, tfsa, tfsj, tqp,
                          iw_(m,M1_IW_EN,k,j,i), hd, jo);
          }
          iw_(m,M1_IW_TA,k,j,i) = hd;
          iw_(m,M1_IW_TB,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
          iw_(m,M1_IW_TC,k,j,i) = jo;
          if (tcorr == M1_THCR_HCS || tjse) {
            iw_(m,M1_IW_TR,k,j,i) = 1.0 - ((hsrc == 2) ? hgh_(m,M1_VGD_HREP,k,j,i)
                                                       : hvc_(m,M1_VET_EPS,k,j,i));
          }
        }
        // beam_kn: the light-front flag (M1HrJFlag) rides in RES as fb + 2 (fb <= 1)
        if (bxmode == 2) {fb += 2.0*M1HrJFlag(hsrc, iw_, hvc_, m, k, j, i);}
        iw_(m,M1_IW_RES,k,j,i) = fb;
      });
#if !M1_LFC_CS
      ImplicitHaloExchange(1, M1_IW_S1);
      ImplicitHaloExchange(1, M1_IW_S3);
      ImplicitHaloExchange(1, M1_IW_RES);
      if (thsw) {
        ImplicitHaloExchange(1, M1_IW_TA);
        ImplicitHaloExchange(1, M1_IW_TB);
        if (tcorr != M1_THCR_NONE || tjs) {ImplicitHaloExchange(1, M1_IW_TC);}
        if (tcorr == M1_THCR_HCS || tjse) {ImplicitHaloExchange(1, M1_IW_TR);}
      }
#endif
      auto fw = (d == 1) ? fw2 : fw3;
      // coreL-1012: on the sp wedge the transverse faces only with hcface_gd = all
      const bool thct = thsw_hctau && (!sph || thsw_hctau_lat);
      const int kup = (d == 2) ? ke + 1 : ke, jup = (d == 1) ? je + 1 : je;
      const int inx = (d == 1) ? BoundaryFace::inner_x2 : BoundaryFace::inner_x3;
      const int onx = (d == 1) ? BoundaryFace::outer_x2 : BoundaryFace::outer_x3;
      par_for("m1_impl_lfchf" M1_LFC_SFX, DevExeSpace(), 0, nmb1, ks, kup, js, jup, is,
              ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        BoundaryFlag blo = mbbcs.d_view(m,inx);
        BoundaryFlag bhi = mbbcs.d_view(m,onx);
        bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
        bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
        const int c = (d == 1) ? j : k;
        const int cs = (d == 1) ? js : ks, ce = (d == 1) ? je : ke;
        if ((c == cs && plo) || (c == ce + 1 && phi)) {
          fw(m,M1_IFW_AL,k,j,i) = 0.0;
          fw(m,M1_IFW_HCL,k,j,i) = 0.0;
          fw(m,M1_IFW_HCR,k,j,i) = 0.0;
          if (thsw) {fw(m,M1_IFW_DG,k,j,i) = 0.0;}
          return;
        }
        const int km = (d == 2) ? k - 1 : k, jm = (d == 1) ? j - 1 : j;
#if M1_LFC_CS
        Real dx = (d == 1) ? cdxf.x2f(m,k,j,i) : cdxf.x3f(m,k,j,i);
        // cs-s4-1017: a seam face: the mirror-pair distance
        const bool smf = (c == cs && cse_(m,(d == 1) ? 0 : 2) != 0) ||
                         (c == ce + 1 && cse_(m,(d == 1) ? 1 : 3) != 0);
        if (smf) {dx = (d == 1) ? cdxe.x2f(m,k,j,i) : cdxe.x3f(m,k,j,i);}
#else
        Real dx = (d == 1) ? mbsize.d_view(m).dx2 : mbsize.d_view(m).dx3;
        if (sph) {dx = (d == 1) ? cdxf.x2f(m,k,j,i) : cdxf.x3f(m,k,j,i);}
#endif
        const Real tauf = M1KtLat(iw_(m,M1_IW_KT,km,jm,i), iw_(m,M1_IW_KT,k,j,i), khl_)*dx;
        Real al, hl, hr;
        if (thsw) {   // thinsw-1010
          Real dgt;
          Real hpf = iw_(m,M1_IW_S1,km,jm,i), hmf = iw_(m,M1_IW_S3,k,j,i);
          const Real wjs = !tjs ? 0.0
                           : (tjse ? 1.0 - fmax(iw_(m,M1_IW_TR,km,jm,i),
                                                iw_(m,M1_IW_TR,k,j,i)) : 1.0);
          // coreL-1012 implicit_thin_hcface = tau: optical-depth weights (equal
          // half-cell lengths on a transverse face)
          Real hwl = -1.0, hwr = 0.0;
          if (thct) {
            const Real tl = fmax(iw_(m,M1_IW_KT,km,jm,i), 0.0);
            const Real tr = fmax(iw_(m,M1_IW_KT,k,j,i), 0.0);
            hwl = (tl + tr > 0.0) ? tr/(tl + tr) : 0.5;
            hwr = 1.0 - hwl;
          }
          // jsfix-1012: the js / jse face shape takes the same face interpolation
          if (tjs) {
            M1ThinJs(iw_(m,M1_IW_S1,km,jm,i), iw_(m,M1_IW_S3,km,jm,i),
                     iw_(m,M1_IW_S1,k,j,i), iw_(m,M1_IW_S3,k,j,i),
                     iw_(m,M1_IW_TC,km,jm,i), iw_(m,M1_IW_TC,k,j,i), tjl,
                     wjs, hpf, hmf, tjm ? -1.0 : hwl, hwr);
          }
          if (hsfb >= 0) {   // cbfix-1011
            M1HrSf(iw_, hsfb, m, km, jm, i, k, j, i, d, wjs, wjs, hpf, hmf);
          }
          Real gfl = iw_(m,M1_IW_RES,km,jm,i), gfr = iw_(m,M1_IW_RES,k,j,i);
          if (hsfb >= 0) {   // cbfix-1011
            gfl = M1HrSfFb(iw_, hsfb, m, km, jm, i, gfl);
            gfr = M1HrSfFb(iw_, hsfb, m, k, j, i, gfr);
          }
          // unify-1012 (implicit_thin_hcsf = jh): the hc/hcs correction of a face whose
          // ratios are the sf ratios over the donor scale Jh takes Jh as the cell's
          // formal-solution scale (M1HrSf: E = Jh at the steady fixed point)
          Real jcl = iw_(m,M1_IW_TC,km,jm,i), jcr = iw_(m,M1_IW_TC,k,j,i);
          if (hsfb >= 0 && hcjh) {
            jcl = fmax(jcl, iw_(m,hsfb+M1_HSF_JH,km,jm,i));
            jcr = fmax(jcr, iw_(m,hsfb+M1_HSF_JH,k,j,i));
          }
          M1ThinFace(tmode, tcorr, hpf, hmf, gfl, gfr, tauf,
                     iw_(m,M1_IW_TB,km,jm,i), iw_(m,M1_IW_TB,k,j,i),
                     iw_(m,M1_IW_TA,km,jm,i), iw_(m,M1_IW_TA,k,j,i),
                     jcl, jcr,
                     (tcorr == M1_THCR_HCS) ? fmin(iw_(m,M1_IW_TR,km,jm,i),
                                                   iw_(m,M1_IW_TR,k,j,i)) : 1.0,
                     th0, th1, tkn0, tknp, ch, al, hl, hr, dgt, hwl, hwr);
          fw(m,M1_IFW_AL,k,j,i) = al;
          fw(m,M1_IFW_HCL,k,j,i) = hl;
          fw(m,M1_IFW_HCR,k,j,i) = hr;
          fw(m,M1_IFW_DG,k,j,i) = dgt;
          return;
        }
        Real fbl = iw_(m,M1_IW_RES,km,jm,i), fbr = iw_(m,M1_IW_RES,k,j,i), jfl = 0.0;
        if (fbl > 1.5) {fbl -= 2.0; jfl = 1.0;}
        if (fbr > 1.5) {fbr -= 2.0; jfl = 1.0;}
        Real hpf = iw_(m,M1_IW_S1,km,jm,i), hmf = iw_(m,M1_IW_S3,k,j,i);
        if (hsfb >= 0) {   // cbfix-1011
          M1HrSf(iw_, hsfb, m, km, jm, i, k, j, i, d, 0.0, 0.0, hpf, hmf);
          fbl = M1HrSfFb(iw_, hsfb, m, km, jm, i, fbl);
          fbr = M1HrSfFb(iw_, hsfb, m, k, j, i, fbr);
        }
        M1HrFace(hpf, hmf, fbl, fbr, tauf, blend,
                 bkind, bfm, btau0, bflo, bfhi, balph,
                 ch, iw_(m,M1_IW_EP,km,jm,i), iw_(m,M1_IW_EP,k,j,i), br0,
                 iw_(m,M1_IW_S3,km,jm,i), iw_(m,M1_IW_S1,k,j,i), ch*dtl/dx, bx0, bxmode,
                 bxwmin, bxr0, jfl, al, hl, hr);
        fw(m,M1_IFW_AL,k,j,i) = al;
        fw(m,M1_IFW_HCL,k,j,i) = hl;
        fw(m,M1_IFW_HCR,k,j,i) = hr;
      });
    }
    return;
  }
  // implicit_flux_beam = fs: |f| = |H|/J of the vet_sc formal solution, into the Thomas
  // scratch M1_IW_S1 (free at this point of the pass) and through the transverse halo
  const bool bfs = impl_beam_fs;
  // vet_gd (sp): the beam vector of the twin-free vet_gd intensities, J and H in the
  // local (r, theta, phi) basis: |f| = |H|/J -> S1, n_theta -> S3, n_phi -> RES (all
  // three free at this point of the pass), below the vet_gd cut (vlat_icut) or where
  // J <= 0 the lagged closure's f_Lev(chi) and n.  Exchanged through the halo.
#if M1_LFC_CS
  const bool bgd = false;   // cs-s4-1017: fs is refused on the cubed sphere (init)
#else
  const bool bgd = bfs && vgd_on && sph;
#endif
  if (bgd) {
    auto vi_ = vgd_i;
    auto dir_ = vgd_dir;
    const int nd = vgd_n;
    const int og = vgd_w - indcs.ng;
    const int ilo = is + vlat_icut;
    const bool sth_ = vgd_sth;   // accsp-1102: stretched theta / full sphere node table
    auto thn_ = vgd_thn;
    par_for("m1_impl_lfcgd", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real fm = M1FLev(iw_(m,M1_IW_WCHI,k,j,i));
      Real n2 = iw_(m,M1_IW_N2,k,j,i), n3 = iw_(m,M1_IW_N3,k,j,i);
      if (i >= ilo) {
        const Real th = sth_ ? thn_(m,j+og)
                             : mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
        const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
        const Real st = sin(th), ct = cos(th), sp = sin(ph), cp = cos(ph);
        Real jm = 0.0, hr = 0.0, ht = 0.0, hp = 0.0;
        for (int d = 0; d < nd; ++d) {
          const Real wi = dir_(d,3)*vi_(m,d,k+og,j+og,i);
          jm += wi;
          hr += wi*(dir_(d,0)*st*cp + dir_(d,1)*st*sp + dir_(d,2)*ct);
          ht += wi*(dir_(d,0)*ct*cp + dir_(d,1)*ct*sp - dir_(d,2)*st);
          hp += wi*(-dir_(d,0)*sp + dir_(d,1)*cp);
        }
        const Real hm = sqrt(hr*hr + ht*ht + hp*hp);
        if (jm > 0.0 && hm > 0.0) {
          fm = fmin(hm/jm, 1.0);
          n2 = ht/hm;
          n3 = hp/hm;
        }
      }
      iw_(m,M1_IW_S1,k,j,i) = fm;
      iw_(m,M1_IW_S3,k,j,i) = n2;
      iw_(m,M1_IW_RES,k,j,i) = n3;
    });
    ImplicitHaloExchange(1, M1_IW_S1);
    ImplicitHaloExchange(1, M1_IW_S3);
    ImplicitHaloExchange(1, M1_IW_RES);
  } else if (bfs) {
    auto vc_ = vet_cell;
    par_for("m1_impl_lfcfs" M1_LFC_SFX, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jj = vc_(m,M1_VET_J,k,j,i);
      const Real h1 = vc_(m,M1_VET_H1,k,j,i), h2 = vc_(m,M1_VET_H1+1,k,j,i);
      const Real h3 = vc_(m,M1_VET_H1+2,k,j,i);
      iw_(m,M1_IW_S1,k,j,i) = (jj > 0.0) ? fmin(sqrt(h1*h1 + h2*h2 + h3*h3)/jj, 1.0)
                                         : 0.0;
    });
    ImplicitHaloExchange(1, M1_IW_S1);
  }
  const int cn2 = bgd ? M1_IW_S3 : M1_IW_N2;
  const int cn3 = bgd ? M1_IW_RES : M1_IW_N3;
  auto fmag = [=] KOKKOS_FUNCTION (const int m, const int k, const int j, const int i) {
    return bfs ? iw_(m,M1_IW_S1,k,j,i) : M1FLev(iw_(m,M1_IW_WCHI,k,j,i));
  };
  par_for("m1_impl_lfc2" M1_LFC_SFX, DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((j == js && plo) || (j == je+1 && phi)) {
      fw2(m,M1_IFW_AL,k,j,i) = 0.0;
      fw2(m,M1_IFW_HCL,k,j,i) = 0.0;
      fw2(m,M1_IFW_HCR,k,j,i) = 0.0;
      return;
    }
    const int jm = j - 1;
#if M1_LFC_CS
    Real dx = cdxf.x2f(m,k,j,i);
    if ((j == js && cse_(m,0) != 0) || (j == je+1 && cse_(m,1) != 0)) {
      dx = cdxe.x2f(m,k,j,i);
    }
#else
    Real dx = mbsize.d_view(m).dx2;
    if (sph) {dx = cdxf.x2f(m,k,j,i);}
#endif
    const Real tauf = M1KtLat(iw_(m,M1_IW_KT,k,jm,i), iw_(m,M1_IW_KT,k,j,i), khl_)*dx;
    Real al, hl, hr;
    M1BeamFace(fmag(m,k,jm,i), iw_(m,cn2,k,jm,i),
               fmag(m,k,j,i), iw_(m,cn2,k,j,i), tauf, edd, blend, bkind,
               bfm, btau0, bflo, bfhi, ch, al, hl, hr);
    fw2(m,M1_IFW_AL,k,j,i) = al;
    fw2(m,M1_IFW_HCL,k,j,i) = hl;
    fw2(m,M1_IFW_HCR,k,j,i) = hr;
  });
  if (std::getenv("M1_BLDBG") != nullptr) {   // DEBUG (blendall-1009): face statistics
    auto h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ifw2);
    auto h1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ifw);
    auto hw = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
    Real s[3][4] = {{0.0}};
    int n = 0;
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        for (int c = 0; c < 3; ++c) {
          s[0][c] += h2(0,c,ks,j,i);
          s[1][c] += h1(0,c,ks,j,i);
        }
        s[2][0] += hw(0,M1_IW_WCHI,ks,j,i);
        s[2][1] += fabs(hw(0,M1_IW_N1,ks,j,i));
        s[2][2] += fabs(hw(0,M1_IW_N2,ks,j,i));
        ++n;
      }
    }
    std::cout << "BLDBG x2 <AL,HCL,HCR> " << s[0][0]/n << " " << s[0][1]/n << " "
              << s[0][2]/n << "  x1 " << s[1][0]/n << " " << s[1][1]/n << " "
              << s[1][2]/n << "  <chi,|n1|,|n2|> " << s[2][0]/n << " " << s[2][1]/n
              << " " << s[2][2]/n << std::endl;
  }
  if (!trans_x3) return;
  par_for("m1_impl_lfc3" M1_LFC_SFX, DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
    BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
    bool plo = (blo != BoundaryFlag::block) && (blo != BoundaryFlag::periodic);
    bool phi = (bhi != BoundaryFlag::block) && (bhi != BoundaryFlag::periodic);
    if ((k == ks && plo) || (k == ke+1 && phi)) {
      fw3(m,M1_IFW_AL,k,j,i) = 0.0;
      fw3(m,M1_IFW_HCL,k,j,i) = 0.0;
      fw3(m,M1_IFW_HCR,k,j,i) = 0.0;
      return;
    }
    const int km = k - 1;
#if M1_LFC_CS
    Real dx = cdxf.x3f(m,k,j,i);
    if ((k == ks && cse_(m,2) != 0) || (k == ke+1 && cse_(m,3) != 0)) {
      dx = cdxe.x3f(m,k,j,i);
    }
#else
    Real dx = mbsize.d_view(m).dx3;
    if (sph) {dx = cdxf.x3f(m,k,j,i);}
#endif
    const Real tauf = M1KtLat(iw_(m,M1_IW_KT,km,j,i), iw_(m,M1_IW_KT,k,j,i), khl_)*dx;
    Real al, hl, hr;
    M1BeamFace(fmag(m,km,j,i), iw_(m,cn3,km,j,i),
               fmag(m,k,j,i), iw_(m,cn3,k,j,i), tauf, edd, blend, bkind,
               bfm, btau0, bflo, bfhi, ch, al, hl, hr);
    fw3(m,M1_IFW_AL,k,j,i) = al;
    fw3(m,M1_IFW_HCL,k,j,i) = hl;
    fw3(m,M1_IFW_HCR,k,j,i) = hr;
  });
