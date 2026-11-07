//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetlat_st_body.hpp
//! \brief share3-1031: the ONE body of RadiationM1::VetLatStencilAdd (spherical polar, M1_VLAT_CS 0)
//! and VetLatStencilAddCs (cubed sphere, M1_VLAT_CS 1), included twice from rad_m1_vetlat.cpp.  The
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
  auto st_ = ost;
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
  auto cx3v = pmy_pack->pcoord->x3v;
#endif
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto cvol = pmy_pack->pcoord->volume;
  auto carea = pmy_pack->pcoord->area;
#if M1_VLAT_CS
  auto car = cs_area_eff;
  auto fn2 = cs_fn2;
  auto fn3 = cs_fn3;
#endif
  const int c0 = M1_TT_LAT0;
#if !M1_VLAT_CS
  const bool tanop = vgd_on && vgd_tan && vgd_tanop;
#endif
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);   // sp-blend-1008, as VetLatOp
  auto ifw_ = ifw;
  const bool blt = blat_on;   // blendall-1009, as VetLatOp
  auto bw2_ = ifw2;
  auto bw3_ = ifw3;
  par_for("m1_vlat_stencil" M1_VLAT_SFX, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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
    int lo[3] = {is, js, ks}, hi[3] = {ie, je, ke};
    if ((q1 == BoundaryFlag::block) || (q1 == BoundaryFlag::periodic)) {lo[0] = is-1;}
    if ((q2 == BoundaryFlag::block) || (q2 == BoundaryFlag::periodic)) {hi[0] = ie+1;}
    if (!p2lo) {lo[1] = js-1;}
    if (!p2hi) {hi[1] = je+1;}
    if (thrd && !p3lo) {lo[2] = ks-1;}
    if (thrd && !p3hi) {hi[2] = ke+1;}
#endif
    const Real cr = ch/cl, kk = ch*cl*dt;
    const Real iv = dt/cvol(m,k,j,i);
#if M1_VLAT_CS
    auto add = [&](const int pi, const int pj, const int pk, const Real c) {
      st_(m,VlatStIdx(pi - i, pj - j, pk - k),k,j,i) += c;
    };
    auto vdot = [&](const int qk, const int qj, const int qi, const Real *nf) {
      return tt_(m,c0+1,qk,qj,qi)*nf[0] + tt_(m,c0+2,qk,qj,qi)*nf[1]
             + tt_(m,c0+3,qk,qj,qi)*nf[2];
    };
    // w M1CsLat0(x; cell (qk, j, qi)) distributed (x1-face cells: k, j fixed)
    auto lat0c = [&](const int qi, const Real w) {
      const Real wv = w/cvol(m,k,j,qi);
      {
        const Real nh[3] = {fn2(m,j+1,0), fn2(m,j+1,1), fn2(m,j+1,2)};
        const Real nl[3] = {fn2(m,j,0), fn2(m,j,1), fn2(m,j,2)};
        const Real ah = 0.5*wv*car.x2f(m,k,j+1,qi), al = 0.5*wv*car.x2f(m,k,j,qi);
        add(qi, j, k, ah*vdot(k, j, qi, nh) - al*vdot(k, j, qi, nl));
        add(qi, j+1, k, ah*vdot(k, j+1, qi, nh));
        add(qi, j-1, k, -al*vdot(k, j-1, qi, nl));
      }
      if (thrd) {
        const Real nh[3] = {fn3(m,k+1,0), fn3(m,k+1,1), fn3(m,k+1,2)};
        const Real nl[3] = {fn3(m,k,0), fn3(m,k,1), fn3(m,k,2)};
        const Real ah = 0.5*wv*car.x3f(m,k+1,j,qi), al = 0.5*wv*car.x3f(m,k,j,qi);
        add(qi, j, k, ah*vdot(k, j, qi, nh) - al*vdot(k, j, qi, nl));
        add(qi, j, k+1, ah*vdot(k+1, j, qi, nh));
        add(qi, j, k-1, -al*vdot(k-1, j, qi, nl));
      }
    };
    // w M1CsLatF(x; nf, cell (qk, qj, i)) distributed
    auto latfc = [&](const int qk, const int qj, const Real *nf, const Real w) {
      const int ia = (i + 1 <= iu) ? i + 1 : i;
      const int ib = (i - 1 >= il) ? i - 1 : i;
      if (ia == ib) {return;}
      const Real r = cx1v(m,i), ra = cx1v(m,ia), rb = cx1v(m,ib);
      const Real f = w/(r*r*r*(ra - rb));
      add(ia, qj, qk, f*ra*ra*ra*vdot(qk, qj, ia, nf));
      add(ib, qj, qk, -f*rb*rb*rb*vdot(qk, qj, ib, nf));
    };
#else
    // x coefficient c at the cell p (absolute indices) of this row
    auto add = [&](const int pi, const int pj, const int pk, const Real c) {
      st_(m,VlatStIdx(pi - i, pj - j, pk - k),k,j,i) += c;
    };
    // w M1SphLat(x; q, d): the frozen coefficients, distributed (the clamps of M1SphLat)
    auto latc = [&](const int d, const int qi, const int qj, const int qk, const Real w) {
      const Real r = cx1v(m,qi);
      if (d == 0) {
        const Real sn = sin(cx2v(m,qj));
        const Real ct = cos(cx2v(m,qj))/sn;
        const int ja = (qj + 1 <= hi[1]) ? qj + 1 : qj;
        const int jb = (qj - 1 >= lo[1]) ? qj - 1 : qj;
        if (ja != jb) {
          const Real f = w/(r*(cx2v(m,ja) - cx2v(m,jb)));
          add(qi, ja, qk, f*tt_(m,c0+1,qk,ja,qi));
          add(qi, jb, qk, -f*tt_(m,c0+1,qk,jb,qi));
        }
        add(qi, qj, qk, w*ct*tt_(m,c0+1,qk,qj,qi)/r);
        if (thrd) {
          const int ka = (qk + 1 <= hi[2]) ? qk + 1 : qk;
          const int kb = (qk - 1 >= lo[2]) ? qk - 1 : qk;
          if (ka != kb) {
            const Real f = w/(r*sn*(cx3v(m,ka) - cx3v(m,kb)));
            add(qi, qj, ka, f*tt_(m,c0+2,ka,qj,qi));
            add(qi, qj, kb, -f*tt_(m,c0+2,kb,qj,qi));
          }
        }
      } else {
        const int ia = (qi + 1 <= hi[0]) ? qi + 1 : qi;
        const int ib = (qi - 1 >= lo[0]) ? qi - 1 : qi;
        if (ia != ib) {
          const Real ra = cx1v(m,ia), rb = cx1v(m,ib);
          const Real f = w/(r*r*r*(ra - rb));
          add(ia, qj, qk, f*ra*ra*ra*tt_(m,c0+d,qk,qj,ia));
          add(ib, qj, qk, -f*rb*rb*rb*tt_(m,c0+d,qk,qj,ib));
        }
        if (tanop) {
          // M1SphTan (vet_gd), the same clamps: d = 1: (1/(r sin)) d_ph (D_tp E) +
          // 2 cot a E / r; d = 2: (1/(r sin^2)) d_th (sin^2 D_tp E)
          const Real sn = sin(cx2v(m,qj));
          const Real ct = cos(cx2v(m,qj))/sn;
          if (d == 1) {
            const int ka = (thrd && qk + 1 <= hi[2]) ? qk + 1 : qk;
            const int kb = (thrd && qk - 1 >= lo[2]) ? qk - 1 : qk;
            if (ka != kb) {
              const Real f = w/(r*sn*(cx3v(m,ka) - cx3v(m,kb)));
              add(qi, qj, ka, f*tt_(m,c0+4,ka,qj,qi));
              add(qi, qj, kb, -f*tt_(m,c0+4,kb,qj,qi));
            }
            add(qi, qj, qk, w*2.0*ct*tt_(m,c0+3,qk,qj,qi)/r);
          } else {
            const int ja = (qj + 1 <= hi[1]) ? qj + 1 : qj;
            const int jb = (qj - 1 >= lo[1]) ? qj - 1 : qj;
            if (ja != jb) {
              const Real sa = sin(cx2v(m,ja)), sb = sin(cx2v(m,jb));
              const Real f = w/(r*sn*sn*(cx2v(m,ja) - cx2v(m,jb)));
              add(qi, ja, qk, f*sa*sa*tt_(m,c0+4,qk,ja,qi));
              add(qi, jb, qk, -f*sb*sb*tt_(m,c0+4,qk,jb,qi));
            }
          }
        }
      }
    };
#endif
    // the faces: upper y -= a od, lower y += a od, od = 0.5 [lat(c) + lat(nb)]
    if (i < ie || !topb) {
      const Real ktf = M1KtX1(iw_, kir_, m, k, j, i, i+1, cx1f, fwd);
      Real w = -0.5*carea.x1f(m,k,j,i+1)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (aphll) {w *= 1.0 - ifw_(m,M1_IFW_AL,k,j,i+1);}
#if M1_VLAT_CS
      lat0c(i, w);
      lat0c(i+1, w);
#else
      latc(0, i, j, k, w);
      latc(0, i+1, j, k, w);
#endif
    }
    if (i > is || !botb) {
      const Real ktf = M1KtX1(iw_, kir_, m, k, j, i-1, i, cx1f, fwd);
      Real w = 0.5*carea.x1f(m,k,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (aphll) {w *= 1.0 - ifw_(m,M1_IFW_AL,k,j,i);}
#if M1_VLAT_CS
      lat0c(i-1, w);
      lat0c(i, w);
#else
      latc(0, i-1, j, k, w);
      latc(0, i, j, k, w);
#endif
    }
    if (!(j == je && p2hi)) {
      const Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j+1,i));
#if M1_VLAT_CS
      Real w = -0.5*car.x2f(m,k,j+1,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (blt) {w *= 1.0 - bw2_(m,M1_IFW_AL,k,j+1,i);}
      const Real nf[3] = {fn2(m,j+1,0), fn2(m,j+1,1), fn2(m,j+1,2)};
      latfc(k, j, nf, w);
      latfc(k, j+1, nf, w);
#else
      Real w = -0.5*carea.x2f(m,k,j+1,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (blt) {w *= 1.0 - bw2_(m,M1_IFW_AL,k,j+1,i);}
      latc(1, i, j, k, w);
      latc(1, i, j+1, k, w);
#endif
    }
    if (!(j == js && p2lo)) {
      const Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + iw_(m,M1_IW_KT,k,j,i));
#if M1_VLAT_CS
      Real w = 0.5*car.x2f(m,k,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (blt) {w *= 1.0 - bw2_(m,M1_IFW_AL,k,j,i);}
      const Real nf[3] = {fn2(m,j,0), fn2(m,j,1), fn2(m,j,2)};
      latfc(k, j-1, nf, w);
      latfc(k, j, nf, w);
#else
      Real w = 0.5*carea.x2f(m,k,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
      if (blt) {w *= 1.0 - bw2_(m,M1_IFW_AL,k,j,i);}
      latc(1, i, j-1, k, w);
      latc(1, i, j, k, w);
#endif
    }
    if (thrd) {
      if (!(k == ke && p3hi)) {
        const Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k+1,j,i));
#if M1_VLAT_CS
        Real w = -0.5*car.x3f(m,k+1,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
        if (blt) {w *= 1.0 - bw3_(m,M1_IFW_AL,k+1,j,i);}
        const Real nf[3] = {fn3(m,k+1,0), fn3(m,k+1,1), fn3(m,k+1,2)};
        latfc(k, j, nf, w);
        latfc(k+1, j, nf, w);
#else
        Real w = -0.5*carea.x3f(m,k+1,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
        if (blt) {w *= 1.0 - bw3_(m,M1_IFW_AL,k+1,j,i);}
        latc(2, i, j, k, w);
        latc(2, i, j, k+1, w);
#endif
      }
      if (!(k == ks && p3lo)) {
        const Real ktf = 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + iw_(m,M1_IW_KT,k,j,i));
#if M1_VLAT_CS
        Real w = 0.5*car.x3f(m,k,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
        if (blt) {w *= 1.0 - bw3_(m,M1_IFW_AL,k,j,i);}
        const Real nf[3] = {fn3(m,k,0), fn3(m,k,1), fn3(m,k,2)};
        latfc(k-1, j, nf, w);
        latfc(k, j, nf, w);
#else
        Real w = 0.5*carea.x3f(m,k,j,i)*iv*cr*kk/(1.0 + ch*dt*ktf);
        if (blt) {w *= 1.0 - bw3_(m,M1_IFW_AL,k,j,i);}
        latc(2, i, j, k-1, w);
        latc(2, i, j, k, w);
#endif
      }
    }
  });
