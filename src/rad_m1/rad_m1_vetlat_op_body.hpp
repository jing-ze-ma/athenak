//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetlat_op_body.hpp
//! \brief share3-1031: the ONE body of RadiationM1::VetLatOp (spherical polar, M1_VLAT_CS 0)
//! and VetLatOpCs (cubed sphere, M1_VLAT_CS 1), included twice from rad_m1_vetlat.cpp.  The
//! shared part is the face loop: the efix rows, the face flags, theta_f of the face
//! opacity (x1: M1KtX1 with the kirchhoff / wdist options), 1 - AL of the blends.  The
//! basis-specific lateral terms sit in `#if M1_VLAT_CS` blocks: sp M1SphLat (+ M1SphTan
//! under vet_gd_tangential) with the sp areas; cs M1CsLat0 on the x1 faces and M1CsLatF
//! with the face's planar normal and the canonical seam areas on the x2/x3 faces, the
//! open seams m1bcs.  Each variant is token for token its former function (bitwise).
//! M1_VLAT_SFX is the kernel-name suffix ("" or "cs").  No include guard.

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto tt_ = VlatTT();   // meshfb-1106
#if M1_VLAT_CS
  auto mbbcs = m1bcs.d_view;
#else
  auto mbbcs = pmy_pack->pmb->mb_bcs.d_view;
#endif
  auto pos_ = part_pos.d_view;
  const int nblkx1 = part_nblk;
  const bool thrd = trans_x3;
  const int bclo = ibc_x1min, bchi = ibc_x1max;
  const Real cl = c_light, ch = chat, dt = dt_sub;
  const bool fwd = sph_geom && impl_face_wdist;
  auto cx1v = pmy_pack->pcoord->x1v;
#if !M1_VLAT_CS
  auto cx2v = pmy_pack->pcoord->x2v;
  auto cx2t = VtanX2t();   // tensor-sp-1108: theta of the data (M1SphTan)
  auto cx3v = pmy_pack->pcoord->x3v;
#endif
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto cvol = pmy_pack->pcoord->volume;
  auto carea = pmy_pack->pcoord->area;
#if M1_VLAT_CS
  auto car = cs_area_eff;
  auto fn2 = cs_fn2;
  auto fn3 = cs_fn3;
  // (the seam foot -> midpoint shift of M1CsLat0 stays lagged: corr = false here)
  auto cse_ = cs_seam.d_view;
  auto g2_ = csg2;
  auto g3_ = csg3;
  // cs-tan-1106: Q in the operator (vet_gd_tan_operator), the frozen terms of M1CsTanFaceOp
  const bool tanq = vgd_cstan && vgd_tanop;
  auto cdxf_ = M1Dxf();
  const auto tg_ = M1MakeCsTanGeo(cs_qf, cdxf_, pmy_pack->pcoord->sin_face_xi,
                                  pmy_pack->pcoord->cos_face_xi,
                                  pmy_pack->pcoord->sin_face_eta,
                                  pmy_pack->pcoord->cos_face_eta, csg2, csg3,
                                  cs_seam.d_view, cs_fn2, cs_fn3, cs_cn2, cs_cn3, cx1v,
                                  js, je, ks, ke);
#endif
  const int cx = xc, cy = yc;
  const Real sg = sgn;
  const int c0 = M1_TT_LAT0;
#if !M1_VLAT_CS
  const bool tanop = vgd_on && vgd_tan && vgd_tanop;
#endif
  // sp-blend-1008: under implicit_flux = berthon | blend the x1 face keeps 1 - AL of the
  // central face flux, the lateral term included (as the row and the face update)
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);
  auto ifw_ = ifw;
  // blendall-1009: the lateral faces keep 1 - AL of their central flux too
  const bool blt = blat_on;
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  par_for("m1_vlat_op" M1_VLAT_SFX, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int ipos = pos_(m);
    const bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
    if ((i == is && botb && bclo == M1_IBC_EFIX) ||
        (i == ie && topb && bchi == M1_IBC_EFIX)) {
      return;
    }
    BoundaryFlag q1 = mbbcs(m,BoundaryFace::inner_x1);
    BoundaryFlag q2 = mbbcs(m,BoundaryFace::outer_x1);
    BoundaryFlag q3 = mbbcs(m,BoundaryFace::inner_x2);
    BoundaryFlag q4 = mbbcs(m,BoundaryFace::outer_x2);
    BoundaryFlag q5 = mbbcs(m,BoundaryFace::inner_x3);
    BoundaryFlag q6 = mbbcs(m,BoundaryFace::outer_x3);
    const bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
    const bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
    const bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
    const bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
#if M1_VLAT_CS
    int il = is, iu = ie;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
#else
    int il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke;
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {il = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {iu = ie+1;}
    if (!p2lo) {jl = js-1;}
    if (!p2hi) {ju = je+1;}
    if (thrd && !p3lo) {kl = ks-1;}
    if (thrd && !p3hi) {ku = ke+1;}
#endif
    const Real cr = ch/cl, kk = ch*cl*dt;
    const Real iv = dt/cvol(m,k,j,i);
    Real y = 0.0;
#if M1_VLAT_CS
    auto lat0 = [&](const int iq) {
      return M1CsLat0(iw_, tt_, c0, cvol, car, fn2, fn3, m, k, j, iq, thrd, cx, false,
                      cse_, g2_, g3_, js, je, ks, ke);
    };
    auto tanf = [&](const int e, const int fk, const int fj) {
      Real tv = 0.0;
      M1CsTanFaceOp(tt_, c0, tg_, e, m, fk, fj, i,
                    [&](const int qk, const int qj, const Real c) {
                      tv += c*iw_(m,cx,qk,qj,i);
                    });
      return tv;
    };
#else
    auto lat = [&](const int d, const int kq, const int jq, const int iq) {
      Real v = M1SphLat(iw_, tt_, c0, cx1v, cx2v, cx3v, m, d, kq, jq, iq, thrd, il, iu,
                        jl, ju, kl, ku, cx);
      if (tanop && d > 0) {   // vet_gd: the tangential cross terms in the operator too
        v += M1SphTan(iw_, tt_, c0, cx1v, cx2v, cx2t, cx3v, m, d, kq, jq, iq, thrd, jl,
                      ju, kl, ku, cx);
      }
      return v;
    };
#endif
    // the two x1 faces
    if (i < ie || !topb) {
      const int ip = i + 1;
      const Real ktf = M1KtX1(iw_, kir_, m, k, j, i, ip, cx1f, fwd);
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (aphll) {th *= 1.0 - ifw_(m,M1_IFW_AL,k,j,i+1);}
#if M1_VLAT_CS
      y -= carea.x1f(m,k,j,i+1)*iv*cr*th*kk*0.5*(lat0(i) + lat0(ip));
#else
      y -= carea.x1f(m,k,j,i+1)*iv*cr*th*kk*0.5*(lat(0,k,j,i) + lat(0,k,j,ip));
#endif
    }
    if (i > is || !botb) {
      const int im = i - 1;
      const Real ktf = M1KtX1(iw_, kir_, m, k, j, im, i, cx1f, fwd);
      Real th = 1.0/(1.0 + ch*dt*ktf);
      if (aphll) {th *= 1.0 - ifw_(m,M1_IFW_AL,k,j,i);}
#if M1_VLAT_CS
      y += carea.x1f(m,k,j,i)*iv*cr*th*kk*0.5*(lat0(im) + lat0(i));
#else
      y += carea.x1f(m,k,j,i)*iv*cr*th*kk*0.5*(lat(0,k,j,im) + lat(0,k,j,i));
#endif
    }
    // the two x2 faces (cs: the face's own planar normal for both cells)
    if (!(j == je && p2hi)) {
      Real th = 1.0/(1.0 + ch*dt*0.5*(iw_(m,M1_IW_KT,k,j,i)
                                      + iw_(m,M1_IW_KT,k,j+1,i)));
      if (blt) {th *= 1.0 - bw2_(m,M1_IFW_AL,k,j+1,i);}
#if M1_VLAT_CS
      const Real nf[3] = {fn2(m,j+1,0), fn2(m,j+1,1), fn2(m,j+1,2)};
      y -= car.x2f(m,k,j+1,i)*iv*cr*th*kk*0.5
           *(M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j, i, il, iu, cx)
             + M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j+1, i, il, iu, cx));
      if (tanq) {y -= car.x2f(m,k,j+1,i)*iv*cr*th*kk*tanf(1, k, j+1);}
#else
      y -= carea.x2f(m,k,j+1,i)*iv*cr*th*kk*0.5*(lat(1,k,j,i) + lat(1,k,j+1,i));
#endif
    }
    if (!(j == js && p2lo)) {
      Real th = 1.0/(1.0 + ch*dt*0.5*(iw_(m,M1_IW_KT,k,j-1,i)
                                      + iw_(m,M1_IW_KT,k,j,i)));
      if (blt) {th *= 1.0 - bw2_(m,M1_IFW_AL,k,j,i);}
#if M1_VLAT_CS
      const Real nf[3] = {fn2(m,j,0), fn2(m,j,1), fn2(m,j,2)};
      y += car.x2f(m,k,j,i)*iv*cr*th*kk*0.5
           *(M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j-1, i, il, iu, cx)
             + M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j, i, il, iu, cx));
      if (tanq) {y += car.x2f(m,k,j,i)*iv*cr*th*kk*tanf(1, k, j);}
#else
      y += carea.x2f(m,k,j,i)*iv*cr*th*kk*0.5*(lat(1,k,j-1,i) + lat(1,k,j,i));
#endif
    }
    // the two x3 faces
    if (thrd) {
      if (!(k == ke && p3hi)) {
        Real th = 1.0/(1.0 + ch*dt*0.5*(iw_(m,M1_IW_KT,k,j,i)
                                        + iw_(m,M1_IW_KT,k+1,j,i)));
        if (blt) {th *= 1.0 - bw3_(m,M1_IFW_AL,k+1,j,i);}
#if M1_VLAT_CS
        const Real nf[3] = {fn3(m,k+1,0), fn3(m,k+1,1), fn3(m,k+1,2)};
        y -= car.x3f(m,k+1,j,i)*iv*cr*th*kk*0.5
             *(M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j, i, il, iu, cx)
               + M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k+1, j, i, il, iu, cx));
        if (tanq) {y -= car.x3f(m,k+1,j,i)*iv*cr*th*kk*tanf(2, k+1, j);}
#else
        y -= carea.x3f(m,k+1,j,i)*iv*cr*th*kk*0.5*(lat(2,k,j,i) + lat(2,k+1,j,i));
#endif
      }
      if (!(k == ks && p3lo)) {
        Real th = 1.0/(1.0 + ch*dt*0.5*(iw_(m,M1_IW_KT,k-1,j,i)
                                        + iw_(m,M1_IW_KT,k,j,i)));
        if (blt) {th *= 1.0 - bw3_(m,M1_IFW_AL,k,j,i);}
#if M1_VLAT_CS
        const Real nf[3] = {fn3(m,k,0), fn3(m,k,1), fn3(m,k,2)};
        y += car.x3f(m,k,j,i)*iv*cr*th*kk*0.5
             *(M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k-1, j, i, il, iu, cx)
               + M1CsLatF(iw_, tt_, c0, cx1v, nf, m, k, j, i, il, iu, cx));
        if (tanq) {y += car.x3f(m,k,j,i)*iv*cr*th*kk*tanf(2, k, j);}
#else
        y += carea.x3f(m,k,j,i)*iv*cr*th*kk*0.5*(lat(2,k-1,j,i) + lat(2,k,j,i));
#endif
      }
    }
    iw_(m,cy,k,j,i) += sg*y;
  });
