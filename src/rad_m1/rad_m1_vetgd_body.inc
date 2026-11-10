//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_body.inc
//! \brief wallbug-1027: the shell kernel body of RadiationM1::VetGdSweep (moved verbatim
//! from rad_m1_vetgd.cpp), included twice there: M1_VGD_BODY_MIR 0 (the default kernel,
//! the code as before) and 1 (vet_gd_theta_wall = mirror: theta ghosts read their mirror
//! row).  Two kernels, so that the mirror read costs the default nothing (compiled into
//! the one kernel it cost +3.6 % GPU time on the accretor deck, 12414551 vs 12414756).
//! Needs M1_VGD_BODY_NAME (the lambda's name) and M1_VGD_BODY_MIR defined.

      auto M1_VGD_BODY_NAME = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                            const int d, const int i, const int l,
                                            const bool vpa) {
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        // the face midpoints (uniform in index space, as the bilinear reads assume)
        // accsp-1102: a stretched theta / the full sphere: the row's node (sph2.inc)
        const Real th = sth_ ? thn_(m,j)
                             : lg_.Th(mbsize, m, j, js);
        // accpolar-1103 vet_gd_band_global: phi from the global index (layout exact)
        const Real ph = gphb_ ? phcb_(m, gpwb_ + k - ks)
                              : lg_.Ph(mbsize, m, k, ks);
        // vgdfuse-1010: sin / cos of the cell's (theta, phi) from the per-cell table (the
        // same expressions, evaluated once per cell instead of once per direction)
        const Real st = tg_(m,k-ks,j-js,0), ct = tg_(m,k-ks,j-js,1);
        const Real sp = tg_(m,k-ks,j-js,2), cp = tg_(m,k-ks,j-js,3);
        const Real mr = nx*st*cp + ny*st*sp + nz*ct;
        if (inw == (mr >= 0.0)) {return;}
        // vet_gd_twin_fuse: the same ray for the shell-mean (twin) field in the same
        // thread; one launch and one halo per shell for both sweeps
        auto solve = [&](const DvceArray5D<Real> &csx, const VgdIView &vix,
                         const bool twn) -> Real {
          // vetgdup-1110 vet_gd_upwind: a read of another cell / direction than the
          // stencil's own (a wall ghost: periodic partner, mirror row, mapped or
          // interpolated directions) shifted along the shell profile of ii from the read
          // value's own n . r_hat (mc) to the ghost's (mt); the symmetric part of the
          // ghost then is that of the ghost's own position and direction.  Opt-in
          // (vet_gd_upwind_wall); on-rank wall partners only (a remote partner's band
          // copy is read unshifted, so the result then depends on the decomposition)
          const bool upw = upwl_ && (upf_ == 1 || upf_ == 2) && !twn;
          auto ush = [&](const int ii, const Real mt, const Real mc) -> Real {
            if (!upw || mt == mc || ((mt < 0.0) != (mc < 0.0))) {return 0.0;}
            Real ft, fc, f1, f2;
            if (!VgdUpProf(upp_, ii, upk_, mt < 0.0, mt, ft, f1, f2)) {return 0.0;}
            if (!VgdUpProf(upp_, ii, upk_, mc < 0.0, mc, fc, f1, f2)) {return 0.0;}
            return ft - fc;
          };
          // mt: the target cosine (the ghost's own; used when upw)
          auto rd0t = [&](const int kk, const int jj, const int ii, const int dq,
                          const Real mt) -> Real {
            const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
            const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
            if (ek == 0 && ej == 0) {   // (a mirror row of the MIR kernel: shifted)
              const Real v0 = vix(m,dq,kk,jj,ii);
              return upw ? (v0 + ush(ii, mt, mr_(m,kk,jj,dq))) : v0;
            }
            const int oo = 3*(ek + 1) + (ej + 1);
            const int sl = 8*m + ((oo < 4) ? oo : (oo - 1));
            const int nl = hl_(sl);
            if (nl < 0) {return vix(m,dq,kk,jj,ii);}
            // accsp-1102: an on-rank pole partner: its rows reversed about the pole
            const int jn = (hp_(sl) != 0) ? VgdPoleRow(jj, wb, nx2b) : (jj - ej*nx2b);
            if (wint && wl_(m,kk,jj) != 0) {
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                const int dt = m3_(m,kk,jj,dq,t);
                const Real vt = vix(nl,dt,kk - ek*nx3b,jn,ii);
                v += w3_(m,kk,jj,dq,t)*(upw ? (vt + ush(ii, mt,
                                                       mr_(nl,kk - ek*nx3b,jn,dt)))
                                            : vt);
              }
              return v;
            }
            const int dd = (wl_(m,kk,jj) != 0) ? mp_(m,kk,jj,dq) : dq;
            const Real vd = vix(nl,dd,kk - ek*nx3b,jn,ii);
            // (a plain neighbour ghost is the same cell: no shift)
            return (upw && wl_(m,kk,jj) != 0) ?
                   (vd + ush(ii, mt, mr_(nl,kk - ek*nx3b,jn,dd))) : vd;
          };
          auto rd0 = [&](const int kk, const int jj, const int ii, const int dq) -> Real {
            return rd0t(kk, jj, ii, dq, upw ? mr_(m,kk,jj,dq) : 0.0);
          };
          auto rd = [&](const int kk, const int jj, const int ii) -> Real {
            // wallbug-1027: a theta ghost of vet_gd_theta_wall = mirror (wall flag 2)
            // reads its mirror row (then as any read, e.g. a phi wall) through its own
            // wall tables; no extra capture (the closure stays below the launch limit)
#if M1_VGD_BODY_MIR
            if ((jj < wb || jj >= wb + nx2b) && wl_(m,kk,jj) >= 2) {
              const bool hi = (jj >= wb + nx2b);
              // vet_gd_theta_wall = vacuum (flag 3): no intensity beyond the wall;
              // extrap (4): the edge row's, same direction (zero gradient)
              // (main field only: the twin keeps the periodic read, rd0 below)
              if (wl_(m,kk,jj) == 3) {
                if (!twn) {return 0.0;}
                return rd0(kk,jj,ii,d);
              }
              if (wl_(m,kk,jj) == 4) {
                // the edge value only when it is of the ghost's branch (the one the
                // stencil selected it for): a turned ray reads its own shell, whose
                // other branch is being written by this launch; else the periodic read
                const int je = hi ? (wb + nx2b - 1) : wb;
                if (!twn && ((mr_(m,kk,je,d) < 0.0) == (mr_(m,kk,jj,d) < 0.0))) {
                  return rd0(kk,je,ii,d);
                }
                return rd0(kk,jj,ii,d);
              }
              const int g = hi ? (jj - wb - nx2b) : (wb - 1 - jj);
              const int jm = hi ? (wb + nx2b - 1 - g) : (wb + g);
              const Real mtg = upw ? mr_(m,kk,jj,d) : 0.0;   // vetgdup-1110
              if (!wint) {return rd0t(kk,jm,ii,mp_(m,kk,jj,d),mtg);}
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                const Real wt = w3_(m,kk,jj,d,t);
                if (t == 0 || wt != 0.0) {v += wt*rd0t(kk,jm,ii,m3_(m,kk,jj,d,t),mtg);}
              }
              return v;
            }
#endif
            return rd0(kk,jj,ii,d);
          };
          const Real r = cx1v(m,i);
          const Real zr = r*fabs(mr);
          const Real p2 = fmax(r*r - zr*zr, 0.0);
          const Real ch0 = exp(csx(m,0,k,j,i));
          const Real s0 = exp(csx(m,1,k,j,i));
          Real iv;
          if (inw && l == n1 - 1) {
            // from the vacuum top face (the end cell's S extrapolated, as vet_col_lat),
            // or with the top bath (vgdbath-1011) its incident intensity attenuated
            const Real rtop = cx1f(m,ie+1);
            const Real ds = sqrt(fmax(rtop*rtop - p2, 0.0)) - zr;
            iv = VgdTopInflow(ch0, ds, s0, exp(csx(m,1,k,j,i-1)), tbh, itop);
          } else if (!inw && zr < 1.0e-12*r) {
            iv = s0;                     // tangent at the cell centre (measure zero)
          } else {
            // the upwind shell, the chord and its own-cell part (rad_m1_vgdray.hpp);
            // typ 1 upwind shell, 2 inner face, 3 own shell
            int typ, iu;
            Real ds, dso;
            VgdRayChord(inw, l, lcut, m, i, zr, p2, cx1v, cx1f, spl, typ, iu, ds, dso);
            // the upwind point, its lateral cell and bilinear fractions (index space)
            const Real xu = r*st*cp - ds*nx;
            const Real yu = r*st*sp - ds*ny;
            const Real zu = r*ct - ds*nz;
            const Real ru = sqrt(xu*xu + yu*yu + zu*zu);
            Real thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
            Real dph = atan2(yu, xu) - ph;
            dph -= twopi*floor((dph + M_PI)/twopi);
            // accsp-1102: behind the polar axis -> the extended chart; the fractional
            // row from the node table (rad_m1_vgdtheta.hpp)
            VgdPoleChart(pol_, th, thu, dph);
            // layinv-1105: the fractional upwind index at its GLOBAL magnitude (fjg, fkg;
            // sj_ = sk_ = 0 without <mesh>/layout_invariant: the old fj, fk exactly); the
            // stencil base and the bilinear fractions come from it, so they are the same
            // numbers on every MeshBlock layout
            const int sj_ = lg_.Sj(m, js), sk_ = lg_.Sk(m, ks);
            const Real fjg = sth_ ? VgdThFj(thn_, m, j, thu, c2w_, sj_)
                                  : (j + sj_) + (thu - th)/lg_.Dx2(mbsize, m);
            const Real fkg = (k + sk_) + dph/lg_.Dx3(mbsize, m);
            const Real fj = fjg - sj_, fk = fkg - sk_;
            int j0 = static_cast<int>(floor(fjg)) - sj_;
            int k0 = static_cast<int>(floor(fkg)) - sk_;
            Real uj = fjg - floor(fjg), uk = fkg - floor(fkg);
            bool clp = false;
            // accsp-1102 polar RING (sph2.inc): both stencil rows inside a polar zone ->
            // the corners are (global phi K, ring code RB0 + 2 (p nzr + q) + across-axis)
            // read from the gathered ring buffers at their true phi (no band, no clamp)
            bool rng = false;
            int rk0 = 0, rk1 = 0, rj0 = 0, rj1 = 0;
            Real ruk = 0.0;
            if (rbon_) {
              const int J0 = lxy_(m,0)*nx2b + (j0 - wb);
              int zp = -1;
              if (J0 >= -nzr_ && J0 + 1 <= nzr_ - 1) {
                zp = 0;
              } else if (J0 >= nt2_ - nzr_ && J0 + 1 <= nt2_ + nzr_ - 1) {
                zp = 1;
              }
              if (zp >= 0) {
                rng = true;
                auto code = [&](const int J) -> int {
                  int q, sh;
                  if (zp == 0) {
                    sh = (J < 0) ? 1 : 0;
                    q = sh ? (-1 - J) : J;
                  } else {
                    sh = (J >= nt2_) ? 1 : 0;
                    q = sh ? (J - nt2_) : (nt2_ - 1 - J);
                  }
                  return VGD_RB0 + 2*(zp*nzr_ + q) + sh;
                };
                rj0 = code(J0);
                rj1 = code(J0 + 1);
                const Real kg = lxy_(m,1)*nx3b + (fk - wb);
                const Real kf = floor(kg);
                ruk = kg - kf;
                int kq = static_cast<int>(kf) % nx3t_;
                if (kq < 0) {kq += nx3t_;}
                rk0 = kq;
                rk1 = (kq + 1) % nx3t_;
              }
            }
            // the full sphere: reads stop at the shell's halo depth
            // (the polar rings' long reads; clamped and counted like any clamp)
            int jlo_ = jlo, jhi_ = jhi, klo_ = klo, khi_ = khi;
            if (pol_) {
              // the depth the halo of the read data filled: the inward band of shell iu
              // (read by the inward pass and by a turned ray), else its outward band
              // (every side, on-rank or not: the sweep is then independent of the
              // decomposition, as the band construction makes it on a wedge)
              const int wu = wsd_((inw || typ == 3) ? 0 : 1, iu);
              klo_ = wb - wu;
              khi_ = wb + nx3b + wu - 2;
              jlo_ = wb - wu;
              jhi_ = wb + nx2b + wu - 2;
              // accpolar-1103 vet_gd_band_global: the CELL's window (the same cells for
              // every layout; inside the block's band since wb <= k, j < wb + nx)
              if (pol_ == 2) {
                klo_ = k - wu;
                khi_ = k + wu - 1;
                jlo_ = j - wu;
                jhi_ = j + wu - 1;
              }
            }
            if (!rng) {
              if (j0 < jlo_) {j0 = jlo_; uj = 0.0; clp = true;}
              if (j0 > jhi_) {j0 = jhi_; uj = 1.0; clp = true;}
              if (k0 < klo_) {k0 = klo_; uk = 0.0; clp = true;}
              if (k0 > khi_) {k0 = khi_; uk = 1.0; clp = true;}
            }
            if (clp) {Kokkos::atomic_add(&cnt_(0), 1.0);}
            int j1 = j0 + 1, k1 = k0 + 1;
            if (rng) {
              j0 = rj0; j1 = rj1; k0 = rk0; k1 = rk1; uk = ruk;
            }
            // share3-1031: the stencil and the field reads of the shared upwind core
            const VgdSt4 st4{k0, k1, j0, j1, uj, uk};
            // accsp-1102: a ring code (jj >= VGD_RB0) reads the ring buffers
            auto rkp = [&](const int kk, const int jj) -> int {
              return (kk + ((jj - VGD_RB0) & 1)*(nx3t_/2)) % nx3t_;
            };
            auto fC = [&](const int c, const int kk, const int jj, const int ii) -> Real {
              if (jj >= VGD_RB0) {
                return rbc_(((static_cast<size_t>(ii)*2 + c)*2*nzr_ + ((jj - VGD_RB0) >> 1))
                            *nx3t_ + rkp(kk,jj));
              }
              return csx(m,c,kk,jj,ii);
            };
            auto fIu = [&](const int kk, const int jj) -> Real {
              if (jj >= VGD_RB0) {
                return rbi_(((static_cast<size_t>(iu)*2*nzr_ + ((jj - VGD_RB0) >> 1))*nx3t_
                             + rkp(kk,jj))*n + d);
              }
              return rd(kk,jj,iu);
            };
            // the branch n . r_hat of a stencil cell (ring: at its true position)
            auto mrq = [&](const int kk, const int jj) -> Real {
              if (jj >= VGD_RB0) {
                const int pq = (jj - VGD_RB0) >> 1;
                const Real tq = rbth_(pq/nzr_, pq % nzr_);
                const Real pk = mp0 + (rkp(kk,jj) + 0.5)*(mpl/nx3t_);
                return nx*sin(tq)*cos(pk) + ny*sin(tq)*sin(pk) + nz*cos(tq);
              }
              return mr_(m,kk,jj,d);
            };
            bool bdone = false;
#if M1_VGD_BODY_MIR
            // wallbug-1027 vet_gd_theta_wall = vacuum: the upwind point lies beyond a
            // physical vacuum theta wall (the mesh's theta limits; flag 3 on the band row
            // nearest the point).  The ray starts at its wall crossing (the cone theta =
            // theta_wall, exact) with I = 0, as from the vacuum top face; chi, S there
            // from the row next to the wall (linear in phi index, linear in r between
            // the shells i and iu as vet_gd_band_exit)
            const bool vlo = (thu < mt0), vhi = (thu > mt0 + mtl);
            const int jrw = static_cast<int>(fmin(fmax(floor(fjg + 0.5) - sj_, Real(jlo)),
                                                  Real(jhi + 1)));
            if (!twn && typ != 2 && (vlo || vhi) && wl_(m,k,jrw) == 3) {
              const bool hi = vhi;
              const Real thw = hi ? (mt0 + mtl) : mt0;
              // the row next to the wall (the block's own edge row with one MeshBlock in
              // theta), band index
              const int jw = static_cast<int>(fmin(fmax(
                lg_.Jw(mbsize, m, js, thw)
                - (hi ? 1 : 0), Real(jlo)), Real(jhi + 1)));
              const Real cw = cos(thw), c2w = cw*cw;
              const Real x0 = r*st*cp, y0 = r*st*sp, z0 = r*ct;
              // |x0 - t n| on the cone: (z0 - t nz)^2 = c2w |x0 - t n|^2, nappe sign(cw)
              const Real xn = x0*nx + y0*ny + z0*nz;
              const Real qa = nz*nz - c2w, qb = 2.0*(c2w*xn - z0*nz);
              const Real qc = z0*z0 - c2w*r*r;
              Real tc = -1.0;
              auto tryt = [&](const Real t) {
                if (!(t > 0.0) || t > ds*(1.0 + 1.0e-9)) {return;}
                if (fabs(cw) > 1.0e-12 && (z0 - t*nz)*cw < 0.0) {return;}
                if (tc < 0.0 || t < tc) {tc = t;}
              };
              if (fabs(qa) > 1.0e-14) {
                const Real dsc = qb*qb - 4.0*qa*qc;
                if (dsc >= 0.0) {
                  const Real sq = sqrt(dsc);
                  const Real qq = -0.5*(qb + ((qb < 0.0) ? -sq : sq));
                  tryt(qq/qa);
                  if (qq != 0.0) {tryt(qc/qq);}
                }
              } else if (qb != 0.0) {
                tryt(-qc/qb);
              }
              if (tc < 0.0) {   // (round-off) the index-linear estimate
                const Real fw = jw + (hi ? 0.5 : -0.5);
                tc = ds*fmin(fmax((fw - j)/(fj - j), 0.0), 1.0);
              }
              const Real dsx = fmin(tc, ds);
              const Real xx = x0 - dsx*nx, yx = y0 - dsx*ny, zx = z0 - dsx*nz;
              const Real rx = sqrt(xx*xx + yx*yx + zx*zx);
              Real dpx = atan2(yx, xx) - ph;
              dpx -= twopi*floor((dpx + M_PI)/twopi);
              const Real fkx = fmin(fmax(k + dpx/lg_.Dx3(mbsize, m), Real(klo)),
                                    Real(khi + 1));
              int kx0 = static_cast<int>(floor(fkx));
              kx0 = (kx0 > khi) ? khi : kx0;
              const Real ukx = fkx - kx0;
              const Real wr = (iu != i) ?
                              fmin(fmax((rx - r)/(cx1v(m,iu) - r), 0.0), 1.0) : 0.0;
              auto linw = [&](const int c, const int ii) {
                return VgdLerp(csx(m,c,kx0,jw,ii), csx(m,c,kx0+1,jw,ii), ukx);
              };
              const Real chx = exp((1.0 - wr)*linw(0,i) + wr*linw(0,iu));
              const Real sx = exp((1.0 - wr)*linw(1,i) + wr*linw(1,iu));
              if (spl) {
                const Real dxo = fmin(dso, dsx);
                iv = fmax(VetSegSplit(0.0, sx, s0, chx*(dsx - dxo), ch0*dxo, spe), 0.0);
              } else {
                Real ex, w0, wu;
                VetSegW(0.5*(chx + ch0)*dsx, ex, w0, wu);
                iv = fmax(wu*sx + w0*s0, 0.0);
              }
              bdone = true;
            }
            if (!bdone && clp && bandx && typ == 1) {
#else
            if (clp && bandx && typ == 1) {
#endif
              // vet_gd_band_exit: the upwind point lies beyond the ghost band.  Instead
              // of reading the clamped (wrong) column, stop the segment where the ray
              // leaves the band (index linear in the path fraction s) and take I, chi, S
              // there, bilinear on the band-edge ghost columns of the two shells i and iu
              // (the lagged inflow of the last sweep) and linear in r between them
              const Real gj = fj - j, gk = fk - k;
              Real s = 1.0;
              if (fj < jlo) {s = fmin(s, (jlo - j)/gj);}
              if (fj > jhi + 1) {s = fmin(s, (jhi + 1 - j)/gj);}
              if (fk < klo) {s = fmin(s, (klo - k)/gk);}
              if (fk > khi + 1) {s = fmin(s, (khi + 1 - k)/gk);}
              s = fmax(s*(1.0 - 1.0e-9), 0.0);
              const Real dsx = s*ds;
              const Real xx = r*st*cp - dsx*nx, yx = r*st*sp - dsx*ny, zx = r*ct - dsx*nz;
              const Real rx = sqrt(xx*xx + yx*yx + zx*zx);
              const Real thx = acos(fmin(fmax(zx/rx, -1.0), 1.0));
              Real dpx = atan2(yx, xx) - ph;
              dpx -= twopi*floor((dpx + M_PI)/twopi);
              const Real fjx = fmin(fmax(j + (thx - th)/lg_.Dx2(mbsize, m), Real(jlo)),
                                    Real(jhi + 1));
              const Real fkx = fmin(fmax(k + dpx/lg_.Dx3(mbsize, m), Real(klo)),
                                    Real(khi + 1));
              int jx0 = static_cast<int>(floor(fjx)), kx0 = static_cast<int>(floor(fkx));
              jx0 = (jx0 > jhi) ? jhi : jx0;
              kx0 = (kx0 > khi) ? khi : kx0;
              const Real ujx = fjx - jx0, ukx = fkx - kx0;
              const Real wr = fmin(fmax((rx - r)/(cx1v(m,iu) - r), 0.0), 1.0);
              const bool want_in = inw;
              Real wsum = 0.0, isum = 0.0;
              for (int q = 0; q < 4; ++q) {
                const int kk = (q < 2) ? kx0 : (kx0 + 1);
                const int jj = (q % 2 == 0) ? jx0 : (jx0 + 1);
                const Real wq = ((q < 2) ? (1.0 - ukx) : ukx)
                                *((q % 2 == 0) ? (1.0 - ujx) : ujx);
                const Real mq = mr_(m,kk,jj,d);
                if (((want_in && mq < 0.0) || (!want_in && mq >= 0.0)) && wq > 0.0) {
                  wsum += wq;
                  isum += wq*((1.0 - wr)*rd(kk,jj,i) + wr*rd(kk,jj,iu));
                }
              }
              if (wsum > 0.0) {
                auto linx = [&](const int c, const int ii) {
                  return VgdLerp(VgdLerp(csx(m,c,kx0,jx0,ii), csx(m,c,kx0,jx0+1,ii), ujx),
                                 VgdLerp(csx(m,c,kx0+1,jx0,ii), csx(m,c,kx0+1,jx0+1,ii),
                                         ujx), ukx);
                };
                const Real chx = exp((1.0 - wr)*linx(0,i) + wr*linx(0,iu));
                const Real sx = exp((1.0 - wr)*linx(1,i) + wr*linx(1,iu));
                if (spl) {
                  const Real dxo = fmin(dso, dsx);
                  iv = fmax(VetSegSplit(isum/wsum, sx, s0, chx*(dsx - dxo),
                                        ch0*dxo, spe), 0.0);
                } else {
                  Real ex, w0, wu;
                  VetSegW(0.5*(chx + ch0)*dsx, ex, w0, wu);
                  iv = fmax((isum/wsum)*ex + wu*sx + w0*s0, 0.0);
                }
                bdone = true;
              }
            }
            if (bdone) {
              // done: the shortened segment from the band edge
            } else if (typ == 2 && emi) {
              // sptests-1011 vet_gd_inner = emit: the inner boundary's emission at the
              // upwind point of the inner face, carried to the cell with the face source
              Real ie0 = lhx ? lhxi_(m,k-ks,j-js) : iemit;   // linhao-prod-1104 per column
              if (bm) {
                ie0 = 0.0;
                const int bq = (d == bdir0) ? 0 : ((nbm > 1 && d == bdir1) ? 1 : -1);
                if (bq >= 0) {
                  const Real bib = (bq == 0) ? bib0 : bib1;
                  if (twn) {
                    ie0 = bib*((bq == 0) ? bfp0 : bfp1);   // the twin: shell-mean inflow
                  } else {
                    Real tu = thu - mt0, pu = ph + dph - mp0;
                    tu -= mtl*floor(tu/mtl);
                    pu -= mpl*floor(pu/mpl);
                    tu += mt0;
                    pu += mp0;
                    const bool in = (bq == 0) ?
                      (tu >= b0t0 && tu <= b0t1 && pu >= b0p0 && pu <= b0p1) :
                      (tu >= b1t0 && tu <= b1t1 && pu >= b1p0 && pu <= b1p1);
                    if (in) {ie0 = bib;}
                  }
                }
              }
              iv = VgdInnerFace(st4, fC, i, true, ie0, bfc, cx1f(m,i), cx1v(m,i),
                                cx1v(m,i+1), p2, ch0, ds, s0);
            } else if (typ == 2) {
              // the diffusion intensity at the inner face of the first shell,
              // S + mu (-dS/dtau_r); vet_bc_diff = face (gdface-1011): S of the FACE,
              // 1.5 S_a - 0.5 S_b (the S the half-cell carry below starts from), not of
              // the first shell's cell (O(dtau_cell) low isotropic part otherwise)
              iv = VgdInnerFace(st4, fC, i, false, 0.0, bfc, cx1f(m,i), cx1v(m,i),
                                cx1v(m,i+1), p2, ch0, ds, s0);
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
                const Real mq = mrq(kk,jj);
                ok[q] = want_in ? (mq < 0.0) : (mq >= 0.0);
                all = all && ok[q];
              }
              // vetgdup-1110 vet_gd_upwind (rad_m1_vgdray.hpp VgdUpCorr): each stencil
              // cell's intensity shifted along the shell's mean angular profile from its
              // own n . r_hat to the upwind point's (main field only; off: no shift)
              Real cq_[4] = {0.0, 0.0, 0.0, 0.0};
              bool upo = false;
              if ((upf_ == 1 || upf_ == 2) && !twn) {
                const Real muu = (nx*xu + ny*yu + nz*zu)/ru;
                auto mqe = [&](const int e) -> Real {return mrq(st4.K(e), st4.J(e));};
                upo = VgdUpCorr(upp_, iu, upk_, want_in, upf_, muu, st4, mqe, ok, all,
                                cq_);
              }
              auto fIc = [&](const int kk, const int jj) -> Real {
                if (!upo) {return fIu(kk,jj);}
                int e = 0;
                while (e < 3 && !(st4.K(e) == kk && st4.J(e) == jj)) {++e;}
                return fIu(kk,jj) + cq_[e];
              };
              // none (-1): the upwind source (below)
              const Real ivu = VgdUpwindI(st4, fIc, ok, all);
              auto lin = [&](const int c) {
                return VgdLerp(VgdLerp(fC(c,k0,j0,iu), fC(c,k0,j1,iu), uj),
                               VgdLerp(fC(c,k1,j0,iu), fC(c,k1,j1,iu), uj), uk);
              };
              const Real chu = exp(lin(0));
              const Real su = exp(lin(1));
              const Real rfi = cx1f(m,i);
              // open-od-1013 vet_path = cells (vet_sc's VetPathStep, rad_m1_vetseg.hpp):
              // a ray whose segment jumps more than one lateral cell (|fj - j| or
              // |fk - k| > 1: dr tan(psi)/(r dtheta) > 1) sees chi and S only at its two
              // ends and misses a lateral material edge narrower than the jump.  Where
              // the cells the step passes (the foot stencil, the cell that holds each
              // node, the own cell) have a chi contrast > kVetCornerRatio (or an S
              // contrast > kVetCornerRatio among those not thin over the step, chi ds
              // > 0.01), the step is subdivided at the NODES where it crosses a theta-row
              // or phi-layer centre (index linear in the path fraction, as
              // vet_gd_band_exit): each node takes chi, S trilinear in ln (theta, phi
              // index; radius between the own and upwind shells at the node's true r),
              // or at a material edge (corner contrast > kVetCornerRatio) the cell that
              // holds it (own shell below the radial face dso, upwind above).  Each
              // sub-segment is the code's own segment (vet_seg split / jump), split at
              // the radial face where that lies inside it, at its midpoint otherwise;
              // the first one from the foot (vet_interp corner / auto as below).  A
              // smooth step (no contrast) keeps the single segment: bitwise.  Not taken
              // on a clamped foot or on the turned ray that dips below its inner face
              // (vet_seg = jump branch below: piecewise per cell in radius already).
              bool pdone = false;
              if (vpa && !clp && !(spl && spe == 2 && typ == 3 && l > lcut &&
                                   p2 < rfi*rfi)) {
                const Real gj = fjg - (j + sj_), gk = fkg - (k + sk_);   // layinv-1105
                const Real a2 = fabs(gj), a3 = fabs(gk);
                const int n2 = (a2 > 1.0) ? (static_cast<int>(ceil(a2)) - 1) : 0;
                const int n3 = (a3 > 1.0) ? (static_cast<int>(ceil(a3)) - 1) : 0;
                if (n2 + n3 > 0) {
                  const int g2 = (gj < 0.0) ? -1 : 1, g3 = (gk < 0.0) ? -1 : 1;
                  const Real rup = cx1v(m,iu);
                  const bool own3 = (typ != 1);
                  // the node: path fraction sn from the foot (0) to the cell (1), its
                  // lateral fractional indices; returns ln chi, ln S (edge: holder cell)
                  auto node = [&](const Real sn, const Real yj, const Real yk, Real &lcn,
                                  Real &lsn, Real &lch, Real &lsh) {
                    const Real q = (1.0 - sn)*ds;   // the path distance from the cell
                    const int jn = static_cast<int>(floor(yj)) - sj_;
                    const int kn = static_cast<int>(floor(yk)) - sk_;
                    const Real aj = yj - floor(yj), ak = yk - floor(yk);
                    Real wr = 0.0;
                    if (!own3) {
                      const Real rq = sqrt(fmax(r*r - 2.0*q*r*mr + q*q, 0.0));
                      wr = fmin(fmax((rq - r)/(rup - r), 0.0), 1.0);
                    }
                    const int ih = (own3 || q < dso) ? i : iu;
                    const int jh = jn + ((aj < 0.5) ? 0 : 1);
                    const int kh = kn + ((ak < 0.5) ? 0 : 1);
                    lch = csx(m,0,kh,jh,ih);
                    lsh = csx(m,1,kh,jh,ih);
                    Real lmn = 1.0e300, lmx = -1.0e300, ws = 0.0, lc = 0.0, ls = 0.0;
                    for (int c = 0; c < 2; ++c) {
                      const Real wc = (c == 0) ? (1.0 - wr) : wr;
                      if (!(wc > 0.0)) {continue;}
                      const int ic = (c == 0) ? i : iu;
                      for (int b = 0; b < 2; ++b) {
                        const Real wk = (b == 0) ? (1.0 - ak) : ak;
                        if (!(wk > 0.0)) {continue;}
                        for (int a = 0; a < 2; ++a) {
                          const Real wj = (a == 0) ? (1.0 - aj) : aj;
                          if (!(wj > 0.0)) {continue;}
                          const Real w = wc*wk*wj;
                          const Real cl = csx(m,0,kn + b,jn + a,ic);
                          lmn = fmin(lmn, cl);
                          lmx = fmax(lmx, cl);
                          ws += w;
                          lc += w*cl;
                          ls += w*csx(m,1,kn + b,jn + a,ic);
                        }
                      }
                    }
                    if ((lmx - lmn) > log(kVetCornerRatio)) {
                      lcn = lch;
                      lsn = lsh;
                    } else {
                      lcn = lc/ws;
                      lsn = ls/ws;
                    }
                  };
                  // the next node after (t2, t3) in path order from the foot
                  auto next = [&](int &t2, int &t3, Real &sn, Real &yj, Real &yk) {
                    const Real s2 = (t2 < n2) ? (1.0 - (n2 - t2)/a2) : 2.0;
                    const Real s3 = (t3 < n3) ? (1.0 - (n3 - t3)/a3) : 2.0;
                    if (s2 <= s3) {
                      sn = s2;
                      // layinv-1105: yj, yk at the global magnitude (node() and the
                      // gate subtract sj_, sk_ from their integer parts)
                      yj = (j + sj_) + g2*(n2 - t2);
                      yk = (k + sk_) + gk*(1.0 - sn);
                      ++t2;
                      if (s3 == s2) {++t3;}   // a row-layer corner: one node
                    } else {
                      sn = s3;
                      yk = (k + sk_) + g3*(n3 - t3);
                      yj = (j + sj_) + gj*(1.0 - sn);
                      ++t3;
                    }
                  };
                  // the gate (ln values): foot stencil, node holder cells, own cell
                  const Real ltn = log(1.0e-2/fmax(ds, 1.0e-300));
                  Real cmn = csx(m,0,k,j,i), cmx = cmn, smn = 1.0e300, smx = -1.0e300;
                  if (cmn > ltn) {smn = csx(m,1,k,j,i); smx = smn;}
                  for (int q = 0; q < 4; ++q) {
                    const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                    *((q % 2 == 0) ? (1.0 - uj) : uj);
                    if (!(wq > 0.0)) {continue;}
                    const int kk = (q < 2) ? k0 : k1, jj = (q % 2 == 0) ? j0 : j1;
                    const Real cl = csx(m,0,kk,jj,iu);
                    cmn = fmin(cmn, cl);
                    cmx = fmax(cmx, cl);
                    if (cl > ltn) {
                      smn = fmin(smn, csx(m,1,kk,jj,iu));
                      smx = fmax(smx, csx(m,1,kk,jj,iu));
                    }
                  }
                  {
                    int t2 = 0, t3 = 0;
                    while (t2 < n2 || t3 < n3) {
                      Real sn, yj, yk;
                      next(t2, t3, sn, yj, yk);
                      const Real q = (1.0 - sn)*ds;
                      const int ih = (own3 || q < dso) ? i : iu;
                      const int jn = static_cast<int>(floor(yj)) - sj_;   // layinv-1105
                      const int kn = static_cast<int>(floor(yk)) - sk_;
                      const int jh = jn + ((yj - floor(yj) < 0.5) ? 0 : 1);
                      const int kh = kn + ((yk - floor(yk) < 0.5) ? 0 : 1);
                      const Real cl = csx(m,0,kh,jh,ih);
                      cmn = fmin(cmn, cl);
                      cmx = fmax(cmx, cl);
                      if (cl > ltn) {
                        smn = fmin(smn, csx(m,1,kh,jh,ih));
                        smx = fmax(smx, csx(m,1,kh,jh,ih));
                      }
                    }
                  }
                  const Real lr = log(kVetCornerRatio);
                  if ((cmx - cmn) > lr || (smx - smn) > lr) {
                    // one sub-segment from (I, chi, S) at path distance qa to (chi, S)
                    // at qb < qa: split at the radial face dso if inside, else midway
                    auto seg = [&](const Real ia, const Real ca, const Real sa,
                                   const Real qa, const Real cb, const Real sb,
                                   const Real qb) -> Real {
                      const Real f = (typ == 1 && spl && dso < qa && dso > qb) ? dso
                                     : 0.5*(qa + qb);
                      if (spl) {return VetSegSplit(ia, sa, sb, ca*(qa - f), cb*(f - qb),
                                                   spe);}
                      Real ex, w0, wu;
                      VetSegW(0.5*(ca + cb)*(qa - qb), ex, w0, wu);
                      return ia*ex + wu*sa + w0*sb;
                    };
                    int t2 = 0, t3 = 0;
                    Real sn, yj, yk, lcn, lsn, lch, lsh;
                    next(t2, t3, sn, yj, yk);
                    node(sn, yj, yk, lcn, lsn, lch, lsh);
                    Real cb = exp(lcn), sb = exp(lsn), qb = (1.0 - sn)*ds;
                    Real ivp;
                    if (vim == 2 ||
                        (vim == 3 && VgdCornerEdge(csx, m, k0, k1, j0, j1, iu, uk, uj))) {
                      Real ws = 0.0, vs = 0.0;
                      for (int q = 0; q < 4; ++q) {
                        const int kk = (q < 2) ? k0 : k1;
                        const int jj = (q % 2 == 0) ? j0 : j1;
                        const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                        *((q % 2 == 0) ? (1.0 - uj) : uj);
                        if (!(wq > 0.0) || !(ok[q] || ivu < 0.0)) {continue;}
                        const Real cq = exp(csx(m,0,kk,jj,iu)), sq = exp(csx(m,1,kk,jj,iu));
                        const Real iq = (ivu < 0.0) ? sq : fIc(kk,jj);
                        ws += wq;
                        vs += wq*seg(iq, cq, sq, ds, cb, sb, qb);
                      }
                      ivp = vs/ws;
                    } else {
                      ivp = seg((ivu < 0.0) ? su : ivu, chu, su, ds, cb, sb, qb);
                    }
                    while (t2 < n2 || t3 < n3) {
                      const Real ca = cb, sa = sb, qa = qb;
                      next(t2, t3, sn, yj, yk);
                      node(sn, yj, yk, lcn, lsn, lch, lsh);
                      cb = exp(lcn);
                      sb = exp(lsn);
                      qb = (1.0 - sn)*ds;
                      ivp = seg(ivp, ca, sa, qa, cb, sb, qb);
                    }
                    iv = fmax(seg(ivp, cb, sb, qb, ch0, s0, 0.0), 0.0);
                    pdone = true;
                    Kokkos::atomic_add(&cnt_(1), 1.0);
                  }
                }
              }
              if (!pdone) {
                // the closing segment (vet_seg = jump below the inner face, vet_interp =
                // corner, vet_seg = split, the plain segment: rad_m1_vgdray.hpp)
                auto edgef = [&]() {
                  auto csr = [&](const int, const int c, const int kk, const int jj,
                                 const int ii) -> Real {return fC(c,kk,jj,ii);};
                  return VgdCornerEdge(csr, m, k0, k1, j0, j1, iu, uk, uj);
                };
                iv = VgdSegFinish(st4, fIc, fC, edgef, ok, iu, k, j, i, ivu, chu, su, ch0,
                                  s0, ds, dso, spl, spe, typ, l, lcut, p2, rfi, zr, vim);
              }
            }
          }
          return iv;
        };
        vi_(m,d,k,j,i) = solve(cs_, vi_, tws);
        if (two) {vt_(m,d,k,j,i) = solve(ct_, vt_, true);}
      };
