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
        const Real th = mbsize.d_view(m).x2min + (j - js + 0.5)*mbsize.d_view(m).dx2;
        const Real ph = mbsize.d_view(m).x3min + (k - ks + 0.5)*mbsize.d_view(m).dx3;
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
          auto rd0 = [&](const int kk, const int jj, const int ii, const int dq) -> Real {
            const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
            const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
            if (ek == 0 && ej == 0) {return vix(m,dq,kk,jj,ii);}
            const int oo = 3*(ek + 1) + (ej + 1);
            const int nl = hl_(8*m + ((oo < 4) ? oo : (oo - 1)));
            if (nl < 0) {return vix(m,dq,kk,jj,ii);}
            if (wint && wl_(m,kk,jj) != 0) {
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                v += w3_(m,kk,jj,dq,t)*vix(nl,m3_(m,kk,jj,dq,t),kk - ek*nx3b,
                                           jj - ej*nx2b,ii);
              }
              return v;
            }
            const int dd = (wl_(m,kk,jj) != 0) ? mp_(m,kk,jj,dq) : dq;
            return vix(nl,dd,kk - ek*nx3b,jj - ej*nx2b,ii);
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
              if (!wint) {return rd0(kk,jm,ii,mp_(m,kk,jj,d));}
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                const Real wt = w3_(m,kk,jj,d,t);
                if (t == 0 || wt != 0.0) {v += wt*rd0(kk,jm,ii,m3_(m,kk,jj,d,t));}
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
            const Real stp = fmax(1.5*s0 - 0.5*exp(csx(m,1,k,j,i-1)), 1.0e-300);
            Real ex, w0, wu;
            VgdW(ch0*ds, ex, w0, wu);
            if (tbh) {
              iv = itop*ex + wu*stp + w0*s0;
            } else {
              iv = wu*stp + w0*s0;
            }
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
            // segfix-1011 (vet_seg = split): the path length from the cell centre back
            // to the radial face it shares with the upwind shell (typ 3: half the chord)
            Real dso = 0.5*ds;
            if (spl && typ == 1) {
              const Real rf = inw ? cx1f(m,i+1) : cx1f(m,i);
              dso = inw ? (sqrt(fmax(rf*rf - p2, 0.0)) - zr)
                        : (zr - sqrt(fmax(rf*rf - p2, 0.0)));
              dso = fmin(fmax(dso, 0.0), ds);
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
            bool bdone = false;
#if M1_VGD_BODY_MIR
            // wallbug-1027 vet_gd_theta_wall = vacuum: the upwind point lies beyond a
            // physical vacuum theta wall (the mesh's theta limits; flag 3 on the band row
            // nearest the point).  The ray starts at its wall crossing (the cone theta =
            // theta_wall, exact) with I = 0, as from the vacuum top face; chi, S there
            // from the row next to the wall (linear in phi index, linear in r between
            // the shells i and iu as vet_gd_band_exit)
            const bool vlo = (thu < mt0), vhi = (thu > mt0 + mtl);
            const int jrw = static_cast<int>(fmin(fmax(floor(fj + 0.5), Real(jlo)),
                                                  Real(jhi + 1)));
            if (!twn && typ != 2 && (vlo || vhi) && wl_(m,k,jrw) == 3) {
              const bool hi = vhi;
              const Real thw = hi ? (mt0 + mtl) : mt0;
              // the row next to the wall (the block's own edge row with one MeshBlock in
              // theta), band index
              const int jw = static_cast<int>(fmin(fmax(
                js + floor((thw - mbsize.d_view(m).x2min)/mbsize.d_view(m).dx2 + 0.5)
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
              const Real fkx = fmin(fmax(k + dpx/mbsize.d_view(m).dx3, Real(klo)),
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
                VgdW(0.5*(chx + ch0)*dsx, ex, w0, wu);
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
              const Real fjx = fmin(fmax(j + (thx - th)/mbsize.d_view(m).dx2, Real(jlo)),
                                    Real(jhi + 1));
              const Real fkx = fmin(fmax(k + dpx/mbsize.d_view(m).dx3, Real(klo)),
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
                  VgdW(0.5*(chx + ch0)*dsx, ex, w0, wu);
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
              Real ie0 = iemit;
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
              const int i1 = i + 1;
              auto sbt = [&](const int kk, const int jj) {
                return fmax(1.5*exp(csx(m,1,kk,jj,i)) - 0.5*exp(csx(m,1,kk,jj,i1)),
                            1.0e-300);
              };
              Real ex, w0, wu;
              VgdW(ch0*ds, ex, w0, wu);
              const Real qv = VgdLerp(VgdLerp(sbt(k0,j0), sbt(k0,j1), uj),
                                      VgdLerp(sbt(k1,j0), sbt(k1,j1), uj), uk);
              iv = ie0*ex + wu*qv + w0*s0;
            } else if (typ == 2) {
              // the diffusion intensity at the inner face of the first shell,
              // S + mu (-dS/dtau_r); vet_bc_diff = face (gdface-1011): S of the FACE,
              // 1.5 S_a - 0.5 S_b (the S the half-cell carry below starts from), not of
              // the first shell's cell (O(dtau_cell) low isotropic part otherwise)
              const Real rf = cx1f(m,i);
              const Real muf = sqrt(fmax(1.0 - p2/(rf*rf), 0.0));
              const int i1 = i + 1;
              auto ib = [&](const int kk, const int jj) {
                const Real sa = exp(csx(m,1,kk,jj,i)), sb = exp(csx(m,1,kk,jj,i1));
                const Real g0 = -(sb - sa)/(cx1v(m,i1) - cx1v(m,i))/exp(csx(m,0,kk,jj,i));
                return (bfc ? (1.5*sa - 0.5*sb) : sa) + g0*muf;
              };
              auto sbt = [&](const int kk, const int jj) {
                return fmax(1.5*exp(csx(m,1,kk,jj,i)) - 0.5*exp(csx(m,1,kk,jj,i1)),
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
                ivu = VgdLerp(VgdLerp(rd(k0,j0,iu), rd(k0,j1,iu), uj),
                              VgdLerp(rd(k1,j0,iu), rd(k1,j1,iu), uj), uk);
              } else {
                Real ws = 0.0, vs = 0.0;
                for (int q = 0; q < 4; ++q) {
                  const int kk = (q < 2) ? k0 : k1;
                  const int jj = (q % 2 == 0) ? j0 : j1;
                  const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                  *((q % 2 == 0) ? (1.0 - uj) : uj);
                  if (ok[q] && wq > 0.0) {ws += wq; vs += wq*rd(kk,jj,iu);}
                }
                ivu = (ws > 0.0) ? vs/ws : -1.0;   // none: the upwind source (below)
              }
              auto lin = [&](const int c) {
                return VgdLerp(VgdLerp(csx(m,c,k0,j0,iu), csx(m,c,k0,j1,iu), uj),
                               VgdLerp(csx(m,c,k1,j0,iu), csx(m,c,k1,j1,iu), uj), uk);
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
                const Real gj = fj - j, gk = fk - k;
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
                    const int jn = static_cast<int>(floor(yj));
                    const int kn = static_cast<int>(floor(yk));
                    const Real aj = yj - jn, ak = yk - kn;
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
                      yj = j + g2*(n2 - t2);
                      yk = k + gk*(1.0 - sn);
                      ++t2;
                      if (s3 == s2) {++t3;}   // a row-layer corner: one node
                    } else {
                      sn = s3;
                      yk = k + g3*(n3 - t3);
                      yj = j + gj*(1.0 - sn);
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
                      const int jn = static_cast<int>(floor(yj));
                      const int kn = static_cast<int>(floor(yk));
                      const int jh = jn + ((yj - jn < 0.5) ? 0 : 1);
                      const int kh = kn + ((yk - kn < 0.5) ? 0 : 1);
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
                      VgdW(0.5*(ca + cb)*(qa - qb), ex, w0, wu);
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
                        const Real iq = (ivu < 0.0) ? sq : rd(kk,jj,iu);
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
              if (pdone) {
                // done: the subdivided step (vet_path = cells)
              } else if (spl && spe == 2 && typ == 3 && l > lcut && p2 < rfi*rfi) {
                // coreL-1012 vet_seg = jump: the turned ray dips below the own inner face;
                // its chord there takes chi, S of the cell below (log-mean of the own
                // column and the upwind point), the two own parts their own: piecewise
                // constant per cell (the limb of an emitting body inside the half cell)
                auto lind = [&](const int c) {
                  const Real a = VgdLerp(VgdLerp(csx(m,c,k0,j0,i-1), csx(m,c,k0,j1,i-1), uj),
                                         VgdLerp(csx(m,c,k1,j0,i-1), csx(m,c,k1,j1,i-1),
                                                 uj), uk);
                  return exp(0.5*(a + csx(m,c,k,j,i-1)));
                };
                const Real zf = sqrt(fmax(rfi*rfi - p2, 0.0));
                const Real e1 = exp(-chu*(zr - zf)), e2 = exp(-lind(0)*2.0*zf);
                const Real e3 = exp(-ch0*(zr - zf));
                const Real ia = ((ivu < 0.0) ? su : ivu)*e1 + su*(1.0 - e1);
                iv = fmax((ia*e2 + lind(1)*(1.0 - e2))*e3 + s0*(1.0 - e3), 0.0);
              } else if (vim == 2 ||
                         (vim == 3 && VgdCornerEdge(csx, m, k0, k1, j0, j1, iu, uk, uj))) {
                // unify-sc-gd-1012 (#6) vet_interp = corner (VetUpCorner of vet_sc): one
                // segment per stencil corner from its own (I, chi, S), the outcome
                // interpolated (renormalised over the corners of the right branch, as the
                // intensity above; none: from the corner's own source)
                Real ws = 0.0, vs = 0.0;
                for (int q = 0; q < 4; ++q) {
                  const int kk = (q < 2) ? k0 : k1;
                  const int jj = (q % 2 == 0) ? j0 : j1;
                  const Real wq = ((q < 2) ? (1.0 - uk) : uk)
                                  *((q % 2 == 0) ? (1.0 - uj) : uj);
                  if (!(wq > 0.0) || !(ok[q] || ivu < 0.0)) {continue;}
                  const Real cq = exp(csx(m,0,kk,jj,iu)), sq = exp(csx(m,1,kk,jj,iu));
                  const Real iq = (ivu < 0.0) ? sq : rd(kk,jj,iu);
                  Real vq;
                  if (spl) {
                    vq = VetSegSplit(iq, sq, s0, cq*(ds - dso), ch0*dso, spe);
                  } else {
                    Real ex, w0, wu;
                    VgdW(0.5*(cq + ch0)*ds, ex, w0, wu);
                    vq = iq*ex + wu*sq + w0*s0;
                  }
                  ws += wq;
                  vs += wq*vq;
                }
                iv = fmax(vs/ws, 0.0);
              } else if (spl) {
                iv = fmax(VetSegSplit((ivu < 0.0) ? su : ivu, su, s0, chu*(ds - dso),
                                      ch0*dso, spe), 0.0);
              } else {
                Real ex, w0, wu;
                VgdW(0.5*(chu + ch0)*ds, ex, w0, wu);
                iv = fmax(((ivu < 0.0) ? su : ivu)*ex + wu*su + w0*s0, 0.0);
              }
            }
          }
          return iv;
        };
        vi_(m,d,k,j,i) = solve(cs_, vi_, tws);
        if (two) {vt_(m,d,k,j,i) = solve(ct_, vt_, true);}
      };
