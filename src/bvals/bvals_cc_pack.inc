//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file bvals_cc_pack.inc
//! \brief Body of the PackAndSendCC team kernel, including the cubed-sphere seam
//! transform and along-seam resample. With BVCC_PACK_KSEAM / BVCC_PACK_KWB defined
//! (CUDA) kseam_ and kwb_ are compile-time constants. Needs BVCC_PACK_NAME defined.
#if defined(BVCC_PACK_KSEAM)
  auto BVCC_PACK_NAME = KOKKOS_LAMBDA(const TeamMember_t &tmember) BV_CUDA_INLINE {
    constexpr bool kseam_ = BVCC_PACK_KSEAM;
    constexpr bool kwb_ = BVCC_PACK_KWB;
#else
  auto BVCC_PACK_NAME = KOKKOS_LAMBDA(const TeamMember_t &tmember, const auto kseam_,
                             const auto kwb_) {
#endif
    const int lr_ = tmember.league_rank()/nsp_;
    const int sp_ = tmember.league_rank() - lr_*nsp_;
    const int m = (lr_)/(nnghbr*nvar);
    const int n = (lr_ - m*(nnghbr*nvar))/nvar;
    const int vl = (lr_ - m*(nnghbr*nvar) - n*nvar);
#if defined(KOKKOS_ENABLE_CUDA)
    const int v = BvVarOf(vmask_, vl);   // (bvals.hpp)
#else
    const int v = vlst_.v[vl];
#endif

    // only load buffers when neighbor exists. CUBED-SPHERE CUBE VERTEX: skipped on
    // both sides; FillPanelCornersCC overwrites exactly this corner block. See bvals.hpp.
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n)) &&
        !(skipd_ && IsX2X3DiagSlot(n))) {
      // if neighbor is at coarser level, use coar indices to pack buffer
      int il, iu, jl, ju, kl, ku;
      if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        il = sbuf[n].icoar[0].bis;
        iu = sbuf[n].icoar[0].bie;
        jl = sbuf[n].icoar[0].bjs;
        ju = sbuf[n].icoar[0].bje;
        kl = sbuf[n].icoar[0].bks;
        ku = sbuf[n].icoar[0].bke;
      // if neighbor is at same level, use same indices to pack buffer
      } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
        il = sbuf[n].isame[0].bis;
        iu = sbuf[n].isame[0].bie;
        jl = sbuf[n].isame[0].bjs;
        ju = sbuf[n].isame[0].bje;
        kl = sbuf[n].isame[0].bks;
        ku = sbuf[n].isame[0].bke;
      // if neighbor is at finer level, use fine indices to pack buffer
      } else {
        il = sbuf[n].ifine[0].bis;
        iu = sbuf[n].ifine[0].bie;
        jl = sbuf[n].ifine[0].bjs;
        ju = sbuf[n].ifine[0].bje;
        kl = sbuf[n].ifine[0].bks;
        ku = sbuf[n].ifine[0].bke;
      }
      int ni = iu - il + 1;
      int nj = ju - jl + 1;
      int nk = ku - kl + 1;
      int nkj  = nk*nj;
      // bvsplit: this team's chunk of the (k, j) range
      const int q0_ = static_cast<int>((static_cast<int64_t>(nkj)*sp_)/nsp_);
      const int q1_ = static_cast<int>((static_cast<int64_t>(nkj)*(sp_ + 1))/nsp_);
      // the window narrows the LOOP only: ni and the index map stay as they were
      int ilp = il, iup = iu;
      // same-level x2/x3 halos only: coarse/fine ranges differ between the two sides
      // and the cubed-sphere corner fill reads the full strips
      if (wlo_ >= 0 && !use_cs && nghbr.d_view(m,n).lev == mblev.d_view(m) &&
          il <= mbis_ && iu >= mbie_) {
        ilp = (il > wlo_) ? il : wlo_;
        iup = (iu < whi_) ? iu : whi_;
      }

      // indices of recv'ing (destination) MB and buffer: MB IDs are stored sequentially
      // in MeshBlockPacks, so array index equals (target_id - first_id)
      int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
      int dn = nghbr.d_view(m,n).dest;

      const bool do_cs = use_cs &&
                         (nghbr.d_view(m,n).panel != mbpanel.d_view(m));
      const bool do_pole = use_pole &&
                           (nghbr.d_view(m,n).polar > 0);

      if (static_cast<bool>(kseam_) && (do_cs || do_pole)) {
        int aj = 1, bj = 0;
        int ak = 1, bk = 0;
        int signvar = 1;
        int sj = 1, sk = nj;
        int vv = v;
        bool cs_xform = false;
        // cs- the FACE-NORMAL pair transform of cs_perm_pairs (see below) and
        // the half-range quad of cs_hrq0
        bool cs_fxform = false;
        bool cs_hrquad = false;
          // component role of slot v: member of a tangential vector pair?
        bool vec_ = (v == IVY) || (v == IVZ);
        int va_ = IVY, vb_ = IVZ;
        if (vtab_) {
          const int r_ = vrole(v);
          vec_ = (r_ >= 0);
          if (vec_) {
            if ((r_ & 1) == 0) {
              va_ = v; vb_ = r_ >> 1;
            } else {
              va_ = r_ >> 1; vb_ = v;
            }
          }
        }
          // a cross-panel neighbour at a coarser level is served from the restricted
          // (coarse) array
        const bool cs_coar = (nghbr.d_view(m,n).lev < mblev.d_view(m));
        int cs_srcpanel = 0, cs_dstpanel = 0;
        // 0 = no along-seam resample; 2 = x2-face seam (resample in k);
        // 3 = x3-face seam (resample in j). See the note on seamval below.
        int cs_seam = 0;

        if (do_cs) {
          const auto ngh = nghbr.d_view(m,n);
          const int my_panel = mbpanel.d_view(m);
          PanelBoundaries pb;
          pb = GetPanelBoundary(my_panel, ngh.panel);

          // x1 is radial and no seam crosses it: i and IVX pass unchanged; the
          // tangential pair is transformed to the destination basis
          cs_xform = vec_;
          cs_dstpanel = ngh.panel;
          cs_srcpanel = my_panel;

          // pb.rev_a/rev_b are expressed in the DESTINATION panel's axes, so undo the
          // swap to get the reversal that applies to this block's own j and k.
          int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
          int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
          if (pb.swap_ax == 1) {
          // transpose the two tangential axes in the buffer across a swap seam, so the
          // receiver's generic unpack lands the data on its own j,k
            sj = nk;
            sk = 1;
          }
          if (rev_a_preswap) {
            aj = -1;
            bj = jl + ju;
          }
          if (rev_b_preswap) {
            ak = -1;
            bk = kl + ku;
          }
          // cs_perm_pairs: the tangential pairs are face-normal components (implicit M1
          // work array), transformed with cubed_sphere::TransformFieldToDstNormals
          if (cperm_ && vec_) {
            cs_xform = false;
            cs_fxform = true;
          }
          // cs- the half-range quad (cs_hrq0, see bvals.hpp)
          if (hrq_ >= 0 && v >= hrq_ && v < hrq_ + 4) {
            cs_xform = false;
            cs_fxform = false;
            cs_hrquad = true;
            const int ad = (v - hrq_)/2, pm = (v - hrq_) % 2;
            const bool sw = (pb.swap_ax == 1);
            const int as = sw ? (1 - ad) : ad;
            const int rva = sw ? rev_b_preswap : rev_a_preswap;
            const int rvb = sw ? rev_a_preswap : rev_b_preswap;
            const bool rv = ((ad == 0) ? rva : rvb) != 0;
            vv = hrq_ + 2*as + (rv ? (1 - pm) : pm);
            signvar = rv ? -1 : 1;
          }

          // along-seam direction from the buffer slot (nghbr_index.hpp): x2-face type
          // buffers resample in x3 (cs_seam = 2), x3-face type in x2 (cs_seam = 3)
          if (n >= 8 && n < 24) {
            cs_seam = 2;
          } else if (n >= 24 && n < 40) {
            cs_seam = 3;
          } else if (n >= 40 && n < 56) {
            // x2x3 edge (40-47) and corner (48-55) buffers are ghost in both tangential
            // directions; the along-seam axis is taken from the neighbour table
            int fj, fk;
            if (n < 48) {
              const int q = n - 40;
              fj = (q/2) & 1;
              fk = (q/4) & 1;
            } else {
              const int q = n - 48;
              fj = (q >> 1) & 1;
              fk = (q >> 2) & 1;
            }
            const int nface2 = (fj == 0) ? 8 : 12;         // the flanking x2 face
            const int nface3 = (fk == 0) ? 24 : 28;        // the flanking x3 face
            const bool s2 = (nghbr.d_view(m,nface2).gid >= 0) &&
                            (nghbr.d_view(m,nface2).panel != my_panel);
            const bool s3 = (nghbr.d_view(m,nface3).gid >= 0) &&
                            (nghbr.d_view(m,nface3).panel != my_panel);
            if (s2 && !s3) {
              cs_seam = 2;
            } else if (s3 && !s2) {
              cs_seam = 3;
            }
          }
          // the resample needs 3 cells along the seam
          const int seam_extent = cs_coar
              ? ((cs_seam == 2) ? cs_indcs.cnx3 : cs_indcs.cnx2)
              : ((cs_seam == 2) ? cs_indcs.nx3 : cs_indcs.nx2);
          if (cs_seam != 0 && seam_extent < 3) { cs_seam = 0; }
          if (no_rs_cc || nors_) { cs_seam = 0; }
          if (cs_fxform || cs_hrquad) { cs_seam = 0; }
        } else if (do_pole) {
          aj = -1;
          bj = jl + ju;
          if (vec_) signvar = -1;
        }

          // seam geometry depends on (xi,eta) only and is computed once per (k,j),
          // outside the radial loop
        const int js_ = cs_coar ? cs_indcs.cjs : cs_indcs.js;
        const int ks_ = cs_coar ? cs_indcs.cks : cs_indcs.ks;
        const int nx2_ = cs_coar ? cs_indcs.cnx2 : cs_indcs.nx2;
        const int nx3_ = cs_coar ? cs_indcs.cnx3 : cs_indcs.nx3;
        const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
        const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;

        // Middle loop over k,j
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          int kk = ak*k + bk;
          int jj = aj*j + bj;

          // along-seam resample: the two charts share the seam-normal coordinate but not
          // the seam-parallel one; a source cell at (n, a) lies at seam-parallel angle
          // atan(tan(a)/tan|n|) of the destination chart, so ghost values are
          // interpolated along the seam (3-point stencil, or linear with cs_lin_resample)
          int kst[3] = {kk, kk, kk};
          int jst[3] = {jj, jj, jj};
          int nst = 1;
          Real wm = 1.0, w0 = 0.0, wp = 0.0;
          // is the stencil interpolating (b is clamped to the source's active range)?
          bool cs_interp = true;
          Real wpos = 0.0;    // wb_seam: the target's along-seam index position
          int winw = 0;       // and the inward seam-normal step in the source chart
          if (cs_seam != 0) {
            const Real xi  = 0.25*M_PI*CellCenterX(jj-js_, nx2_, x2mn, x2mx);
            const Real eta = 0.25*M_PI*CellCenterX(kk-ks_, nx3_, x3mn, x3mx);
            Real ang, nrm, dang;
            int sc, blo, bhi;
            if (cs_seam == 2) {
              ang = eta; nrm = xi;
              dang = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3_);
              sc = kk;
              blo = ks_; bhi = ks_ + nx3_ - 3;
            } else {
              ang = xi; nrm = eta;
              dang = 0.25*M_PI*(x2mx - x2mn)/static_cast<Real>(nx2_);
              sc = jj;
              blo = js_; bhi = js_ + nx2_ - 3;
            }
          // <mesh>/cs_seam_blockend = extend (second pass): at an internal block end the
          // stencil may extend into the source block's ghost cells
            if (ext_ && !cs_coar) {
              const Real w_ = (cs_seam == 2) ? (x3mx - x3mn) : (x2mx - x2mn);
              const Real x0_ = (cs_seam == 2) ? x3mn : x2mn;
              const int nbp_ = static_cast<int>(floor(2.0/w_ + 0.5));
              const int ib_ = static_cast<int>(floor((x0_ + 1.0)/w_ + 0.5));
              if (ib_ > 0) {blo -= ngx_;}
              if (ib_ < nbp_ - 1) {bhi += ngx_;}
            }
          // the stencil is clamped to the active range of the source block (not of the
          // buffer), since it reads the array directly
            const Real wp_ = (cs_seam == 2) ? (x3mx - x3mn) : (x2mx - x2mn);
            const Real xp_ = (cs_seam == 2) ? x3mn : x2mn;
            const int spn_ = static_cast<int>(floor((xp_ + 1.0)/wp_ + 0.5))
                             *((cs_seam == 2) ? nx3_ : nx2_);
            const Real pos = cs_seam::SeamPos(sc, spn_,
                                              (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang);
            wpos = pos;
            winw = (nrm > 0.0) ? -1 : 1;
            int b = cs_seam::SeamRound(pos) - 1;   // tie-safe
            b = (b < blo) ? blo : ((b > bhi) ? bhi : b);
            const Real u = pos - static_cast<Real>(b + 1);
            cs_interp = (u >= -1.0) && (u <= 1.0);
            wm = 0.5*u*(u - 1.0);
            w0 = 1.0 - u*u;
            wp = 0.5*u*(u + 1.0);
            if (rs_lin && cs_interp) {
              // <mesh>/cs_seam_resample = linear: the two bracketing donors, weights >= 0
              wm = (u < 0.0) ? -u : 0.0;
              wp = (u > 0.0) ? u : 0.0;
              w0 = 1.0 - wm - wp;
            }
            nst = 3;
            for (int s=0; s<3; ++s) {
              if (cs_seam == 2) {
                kst[s] = b + s;
              } else {
                jst[s] = b + s;
              }
            }
          }

          // tangent-basis geometry of each stencil cell (the two charts have different
          // tangent bases at the same point)
          cubed_sphere::SeamXform xf[3];
          if (cs_xform) {
            for (int s=0; s<nst; ++s) {
              const Real xis = 0.25*M_PI*CellCenterX(jst[s]-js_, nx2_, x2mn, x2mx);
              const Real etas = 0.25*M_PI*CellCenterX(kst[s]-ks_, nx3_, x3mn, x3mx);
              cubed_sphere::SeamXformAt(cs_srcpanel, cs_dstpanel, xis, etas, xf[s]);
            }
          }
          // cs- the face-normal pair (cs_perm_pairs), one source cell (nst = 1)
          cubed_sphere::SeamFieldXform fxf;
          if (cs_fxform) {
            const Real xis = 0.25*M_PI*CellCenterX(jst[0]-js_, nx2_, x2mn, x2mx);
            const Real etas = 0.25*M_PI*CellCenterX(kst[0]-ks_, nx3_, x3mn, x3mx);
            cubed_sphere::SeamFieldXformAt(cs_srcpanel, cs_dstpanel, xis, etas, fxf);
          }

          // value of variable v at stencil cell s, radial index i (tangential pairs
          // transformed to the destination basis)
          auto sval = [&](const int s, const int i) BV_CUDA_INLINE {
            const int kq = kst[s], jq = jst[s];
            if (cs_fxform) {
              Real f2o, f3o;
              const Real fy_ = cs_coar ? ca(m,va_,kq,jq,i) : a(m,va_,kq,jq,i);
              const Real fz_ = cs_coar ? ca(m,vb_,kq,jq,i) : a(m,vb_,kq,jq,i);
              cubed_sphere::ApplyFieldXform(fxf, fy_, fz_, f2o, f3o);
              return (v == va_) ? f2o : f3o;
            }
            if (!cs_xform) {
              return (cs_coar ? ca(m,vv,kq,jq,i) : a(m,vv,kq,jq,i))*signvar;
            }
            Real m2o, m3o;
            const Real my_ = cs_coar ? ca(m,va_,kq,jq,i) : a(m,va_,kq,jq,i);
            const Real mz_ = cs_coar ? ca(m,vb_,kq,jq,i) : a(m,vb_,kq,jq,i);
            cubed_sphere::ApplyMomentumXform(xf[s], my_, mz_, m2o, m3o);
            return (v == va_) ? m2o : m3o;
          };

          // monotonicity limit: where the stencil data are monotone the resampled value
          // is clamped to their range (prevents overshoot at shocks); smooth extrema are
          // left unlimited
          cubed_sphere::SeamXform xfc;
          const bool rgon = (rg_ > 0.0) && (cs_seam != 0) && (v >= 0);
          if (rgon && cs_xform) {
            const Real xic = 0.25*M_PI*CellCenterX(jj-js_, nx2_, x2mn, x2mx);
            const Real etac = 0.25*M_PI*CellCenterX(kk-ks_, nx3_, x3mn, x3mx);
            cubed_sphere::SeamXformAt(cs_srcpanel, cs_dstpanel, xic, etac, xfc);
          }
          auto rgbad = [&](const int i) BV_CUDA_INLINE {
            const Real d0 = cs_coar ? ca(m,0,kst[0],jst[0],i) : a(m,0,kst[0],jst[0],i);
            const Real d1 = cs_coar ? ca(m,0,kst[1],jst[1],i) : a(m,0,kst[1],jst[1],i);
            const Real d2 = cs_coar ? ca(m,0,kst[2],jst[2],i) : a(m,0,kst[2],jst[2],i);
            Real dq = wm*d0 + w0*d1 + wp*d2;
            if (cs_interp) {
              dq = fmin(fmax(d0, fmax(d1, d2)), fmax(fmin(d0, fmin(d1, d2)), dq));
            }
            const Real dc = cs_coar ? ca(m,0,kk,jj,i) : a(m,0,kk,jj,i);
            return (dq > rg_*dc) || (dc > rg_*dq);
          };
          auto plainval = [&](const int i) BV_CUDA_INLINE {
            if (!cs_xform) {
              return (cs_coar ? ca(m,vv,kk,jj,i) : a(m,vv,kk,jj,i))*signvar;
            }
            Real m2o, m3o;
            const Real my_ = cs_coar ? ca(m,va_,kk,jj,i) : a(m,va_,kk,jj,i);
            const Real mz_ = cs_coar ? ca(m,vb_,kk,jj,i) : a(m,vb_,kk,jj,i);
            cubed_sphere::ApplyMomentumXform(xfc, my_, mz_, m2o, m3o);
            return (v == va_) ? m2o : m3o;
          };
          // wb_seam: the state of stencil cell s walked to the ghost's potential
          // (density for v = 0, total energy for v = wbsien_); see bvals.hpp
          const bool wbsv = static_cast<bool>(kwb_) && wbs_ && (cs_seam != 0) &&
                            !cs_coar && (v == 0 || v == wbsien_);
          Real lw5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
          int c5 = 0;
          if (wbsv) {
            c5 = static_cast<int>(floor(wpos + 0.5));
            const Real x = wpos - static_cast<Real>(c5);
            for (int q = 0; q < 5; ++q) {
              Real l = 1.0;
              for (int r = 0; r < 5; ++r) {
                if (r != q) l *= (x - static_cast<Real>(r - 2))/static_cast<Real>(q - r);
              }
              lw5[q] = l;
            }
          }
          auto wbval = [&](const int s, const int i) BV_CUDA_INLINE {
            const int kq = kst[s], jq = jst[s];
            Real phig = 0.0;
            const int nA = (cs_seam == 2) ? wbsphi_.extent_int(1) : wbsphi_.extent_int(2);
            for (int q = 0; q < 5; ++q) {
              int ix = c5 - 2 + q;
              ix = (ix < 0) ? 0 : ((ix > nA - 1) ? nA - 1 : ix);
              phig += lw5[q]*((cs_seam == 2) ? wbsphi_(m,ix,jj,i) : wbsphi_(m,kk,ix,i));
            }
            auto prim = [&](const int kc, const int jc, Real &d, Real &e, Real &ke,
                            Real &ph) {
              d = a(m,0,kc,jc,i);
              const Real m1 = a(m,1,kc,jc,i), m2 = a(m,2,kc,jc,i), m3 = a(m,3,kc,jc,i);
              const Real xs = 0.25*M_PI*CellCenterX(jc-js_, nx2_, x2mn, x2mx);
              const Real es = 0.25*M_PI*CellCenterX(kc-ks_, nx3_, x3mn, x3mx);
              const Real cc = -sin(xs)*sin(es);
              ke = 0.5*(m1*m1 + (m2*m2 - 2.0*cc*m2*m3 + m3*m3)/(1.0 - cc*cc))/d;
              ph = wbsphi_(m,kc,jc,i);
              e = a(m,wbsien_,kc,jc,i) - ke - (wbsetg_ ? d*ph : 0.0);
              if (wbsmhd_) {
                const Real b1 = wbsbcc_(m,0,kc,jc,i), b2 = wbsbcc_(m,1,kc,jc,i);
                const Real b3 = wbsbcc_(m,2,kc,jc,i);
                ke += 0.5*(b1*b1 + b2*b2 + b3*b3);   // carried as part of E like KE
                e -= 0.5*(b1*b1 + b2*b2 + b3*b3);
              }
            };
            Real d, e, ke, ph;
            prim(kq, jq, d, e, ke, ph);
            Real dd = d, ee = e, tt = -1.0;
            if (wbsopt_ == 3) {
              const int ki = (cs_seam == 3) ? kq + winw : kq;
              const int ji = (cs_seam == 2) ? jq + winw : jq;
              Real di, ei, kei, phi_;
              prim(ki, ji, di, ei, kei, phi_);
              const Real ts = wbseos_.Temperature(d, e), ti = wbseos_.Temperature(di, ei);
              const Real dph = ph - phi_;
              const Real al = (fabs(dph) > 1.0e-12*(fabs(ph) + 1.0e-300)) ? (ts - ti)/dph : 0.0;
              WBAdvance(wbseos_, 3, d, e, phig - ph, dd, ee, tt, ts, al, ts);
            } else {
              WBAdvance(wbseos_, 1, d, e, phig - ph, dd, ee, tt);
            }
            if (v == 0) return dd;
            return ee + ke + (wbsetg_ ? dd*phig : 0.0);
          };
          auto seamval = [&](const int i) BV_CUDA_INLINE {
            if (cs_seam == 0) return sval(0,i);
            if (wbsv) {
              const Real sv0 = wbval(0,i), sv1 = wbval(1,i), sv2 = wbval(2,i);
              const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
              if (!cs_interp) return qq;
              return fmin(fmax(sv0, fmax(sv1, sv2)), fmax(fmin(sv0, fmin(sv1, sv2)), qq));
            }
            if (rgon && rgbad(i)) return plainval(i);
            const Real sv0 = sval(0,i);
            const Real sv1 = sval(1,i);
            const Real sv2 = sval(2,i);
            const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
            if (!cs_interp) return qq;
            const Real lo0 = fmin(sv0, fmin(sv1, sv2));
            const Real hi0 = fmax(sv0, fmax(sv1, sv2));
            return fmin(hi0, fmax(lo0, qq));
          };

          // Inner (vector) loop over i
          // copy directly into recv buffer if MeshBlocks on same rank

          // seamval reads u0 or coarse_u0 according to cs_coar, so one expression
          // serves a neighbour at ANY level.
          if (nghbr.d_view(m,n).rank == my_rank) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*vl);
              rbuf[dn].vars(dm, index) = val;
            });

          // else copy into send buffer for MPI communication below

          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*vl);
              sbuf[n].vars(m,index) = val;
            });
          }
        });
      } else {   // normal boundary exchange
      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_), [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // Inner (vector) loop over i
        // copy directly into recv buffer if MeshBlocks on same rank

        if (nghbr.d_view(m,n).rank == my_rank) {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = ca(m,v,k,j,i);
            });
          }

        // else copy into send buffer for MPI communication below

        } else {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = ca(m,v,k,j,i);
            });
          }
        }
      });
      }  // end if-do-cs/do-pole block
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  };
