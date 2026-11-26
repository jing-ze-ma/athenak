//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_sweep_body.inc
//! \brief shell kernel body of RadiationVet::ScSweep, included twice:
//! VET_VGD_BODY_MIR 0 (default) and 1 (vet_gd_theta_wall = mirror), so the mirror read
//! costs the default kernel nothing. Needs VET_VGD_BODY_NAME and VET_VGD_BODY_MIR.
//! Returns true for an item deferred by dfr (bool or std::false_type): its stencil
//! reaches the band, which an exchange still in flight fills (interior / edge split,
//! RadiationVet::vgd_mid); the last argument is unused.  With VET_VGD_SWEEP_SLAB
//! (sc_sph_sweep_shell.inc) the Cartesian slab's geometry, x1 ends and face scalars;
//! with VET_VGD_SWEEP_FS the spherical-polar rays' r-face intensities (ScSpFaceSum).

// nvcc rejects generic __host__ __device__ lambdas: on CUDA the arguments are plain bools
#if defined(KOKKOS_ENABLE_CUDA) && defined(VET_VGD_BODY_GX0)
      auto VET_VGD_BODY_NAME = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                            const int d, const int i, const int l,
                                            const bool dfr, const bool) -> bool {
#else
#if defined(KOKKOS_ENABLE_CUDA)
#define VET_VGD_TT bool
#else
#define VET_VGD_TT auto
#endif
      auto VET_VGD_BODY_NAME = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                            const int d, const int i, const int l,
                                            const VET_VGD_TT dfr,
                                            const VET_VGD_TT) -> bool {
#undef VET_VGD_TT
#endif
        // several MeshBlocks in x1: first shell / top face of the innermost / outermost
        const int lcb = (x1f_(m) & 1) ? lcut : -1;
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
#if VET_VGD_SWEEP_SLAB
        const Real mr = nx;   // slab: r_hat = x1_hat
#else
        // face midpoints (uniform in index space); stretched theta / full sphere: node
        const Real th = sth_ ? thn_(m,j)
                             : lg_.Th(mbsize, m, j, js);
        // vet_gd_band_global: phi from the global index
        const Real ph = gphb_ ? phcb_(m, gpwb_ + k - ks)
                              : lg_.Ph(mbsize, m, k, ks);
        // sin / cos of the cell's (theta, phi), from the per-cell table
        const Real st = tg_(m,k-ks,j-js,0), ct = tg_(m,k-ks,j-js,1);
        const Real sp = tg_(m,k-ks,j-js,2), cp = tg_(m,k-ks,j-js,3);
        const Real mr = nx*st*cp + ny*st*sp + nz*ct;
#endif
        if (inw == (mr >= 0.0)) {return false;}
        bool dfd = false;
#if VET_VGD_SWEEP_SLAB
        // closure_consistency (sfs_): the ray's x1 face-sum scalars (entry, crossing,
        // exit; the x2 / x3 ones in ScSlabFaceSum), as the Cartesian sweep's
        Real gq[3] = {0.0, 0.0, 0.0};
#endif
#if VET_VGD_SWEEP_FS
        // closure_consistency (ScSpFaceSum): the intensity through the ray's r face
        // (ScSpFsSlot; the exit's is stored where it is made)
        Real gq[1] = {0.0};
#endif
        // the ray of direction d through cell (k, j, i) on source csx, intensities vix
        auto solve = [&](const DvceArray5D<Real> &csx, const ScIView &vix) -> Real {
          auto rd0 = [&](int kk, const int jj, const int ii, int dq) -> Real {
            // axisymmetric (nx3 = 1): a phi band column is the own column with the
            // direction rotated by its phi offset
            if (!raw_ && nx3b == 1 && kk != wb) {
              dq = mp_(m,kk,wb,dq);
              kk = wb;
            }
            const int ek = (kk < wb) ? -1 : ((kk >= wb + nx3b) ? 1 : 0);
            const int ej = (jj < wb) ? -1 : ((jj >= wb + nx2b) ? 1 : 0);
            if (ek == 0 && ej == 0) {   // MIR kernel: a mirror row is shifted
              return vix(m,dq,kk,jj,ii);
            }
            const int oo = 3*(ek + 1) + (ej + 1);
            const int sl = 8*m + ((oo < 4) ? oo : (oo - 1));
            const int ow = raw_ ? own_(m,kk,jj,0) : 0;
            const int nl = raw_ ? ((ow < 0) ? m : ow) : hl_(sl);
            if (nl < 0) {
              // band copy (the unpack applied the wall map)
              return vix(m,dq,kk,jj,ii);
            }
            // an on-rank pole partner: its rows reversed about the pole
            const int jn = raw_ ? ((ow < 0) ? jj : own_(m,kk,jj,2))
                                : ((hp_(sl) != 0) ? ScPoleRow(jj, wb, nx2b)
                                                  : (jj - ej*nx2b));
            const int kn = raw_ ? ((ow < 0) ? kk : own_(m,kk,jj,1)) : (kk - ek*nx3b);
            if (wint && wl_(m,kk,jj) != 0) {
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                v += w3_(m,kk,jj,dq,t)*vix(nl,m3_(m,kk,jj,dq,t),kn,jn,ii);
              }
              return v;
            }
            const int dd = (wl_(m,kk,jj) != 0) ? mp_(m,kk,jj,dq) : dq;
            return vix(nl,dd,kn,jn,ii);
          };
          // rdq: direction dq at a stencil cell; rd: the thread's own direction d
          auto rdq = [&](const int kk, const int jj, const int ii, const int dq) -> Real {
            // theta_wall = mirror (flag 2): a theta ghost reads its mirror row (no extra
            // capture: the closure stays below the launch limit)
#if VET_VGD_BODY_MIR
            // slab: beyond a vacuum face
            if (VET_VGD_SWEEP_SLAB && wl_(m,kk,jj) >= 8) {return 0.0;}
            if ((jj < wb || jj >= wb + nx2b) && wl_(m,kk,jj) >= 2) {
              const bool hi = (jj >= wb + nx2b);
              // vacuum (flag 3): zero; extrap (4): the edge row's (zero gradient)
              if (wl_(m,kk,jj) == 3) {return 0.0;}
              if (wl_(m,kk,jj) == 4) {
                // edge value only if of the ghost's branch (the other branch of the
                // own shell is being written); else the periodic read
                const int je = hi ? (wb + nx2b - 1) : wb;
                if ((mr_(m,kk,je,dq) < 0.0) == (mr_(m,kk,jj,dq) < 0.0)) {
                  return rd0(kk,je,ii,dq);
                }
                return rd0(kk,jj,ii,dq);
              }
              const int g = hi ? (jj - wb - nx2b) : (wb - 1 - jj);
              const int jm = ScMirRow(g, hi, wb, nx2b);
              if (!wint) {return rd0(kk,jm,ii,mp_(m,kk,jj,dq));}
              Real v = 0.0;
              for (int t = 0; t < 3; ++t) {
                const Real wt = w3_(m,kk,jj,dq,t);
                if (t == 0 || wt != 0.0) {v += wt*rd0(kk,jm,ii,m3_(m,kk,jj,dq,t));}
              }
              return v;
            }
#endif
            return rd0(kk,jj,ii,dq);
          };
          auto rd = [&](const int kk, const int jj, const int ii) -> Real {
            return rdq(kk, jj, ii, d);
          };
          const Real r = cx1v(m,i);
          const Real zr = r*fabs(mr);
          const Real p2 = fmax(r*r - zr*zr, 0.0);
          const Real ch0 = exp(csx(m,0,k,j,i));
          const Real s0 = exp(csx(m,1,k,j,i));
          Real iv;
#if VET_VGD_SWEEP_SLAB
          // x1 periodic: the first pass enters with the face S (psrc), the later ones
          // from the periodic image of the other end (img)
          const bool psrc = (x1p_ == 1) && x1si_;
          const bool img = (x1p_ == 2) && (inw ? (l == n1 - 1) : (l == lcb));
          // the boundary beams entering the inner face in direction d (own cell)
          auto bmv = [&]() -> Real {
            Real b = 0.0;
            for (int q = 0; q < bmn_; ++q) {
              const int kb = (bmi_.extent_int(2) > 1) ? k - ks : 0;
              if (d == bmr_[q]) {b += bmi_(m,q,kb,j - js);}
            }
            return b;
          };
#else
          constexpr bool psrc = false, img = false;
#endif
          if (inw && l == n1 - 1 && (x1f_(m) & 2) && !img) {
            // top face: extrapolated S (vacuum), or the bath / mirrored (ScTopMirror)
            // inflow attenuated
            const Real rtop = cx1f(m,ie+1);
#if VET_VGD_SWEEP_SLAB
            const Real ds = (rtop - r)/fabs(mr);
#else
            const Real ds = sqrt(fmax(rtop*rtop - p2, 0.0)) - zr;
#endif
            const bool rfl = (x1f_(m) & 4);
            const Real sbl = exp(csx(m,1,k,j,i-1));
            const Real tin = rfl ? vix(m,d,k,j,i+1)
                                 : (psrc ? fmax(1.5*s0 - 0.5*sbl, 1.0e-300) : itop);
            iv = ScTopInflow(ch0, ds, s0, sbl, tbh || rfl || psrc, tin);
#if VET_VGD_SWEEP_SLAB
            // (a mirrored face carries no entry, as the Cartesian sweep)
            if (sfs_ && !rfl) {gq[0] = dir_(d,3)*nx*((tbh || psrc) ? tin : 0.0);}
#endif
#if VET_VGD_SWEEP_FS
            gq[SPFS_XV] = (tbh || rfl) ? tin : 0.0;   // (entry through the top face)
#endif
          } else if (!VET_VGD_SWEEP_SLAB && !inw && zr < 1.0e-12*r) {
            iv = s0;                     // tangent at the cell centre (measure zero)
#if VET_VGD_SWEEP_FS
            gq[SPFS_XV] = (l > lcut) ? vix(m,d,k,j,i-1) : (emi ? iemit : s0);
#endif
          } else {
            // upwind shell, chord, own-cell part; typ 1 upwind, 2 inner face, 3 own shell
            int typ, iu;
            Real ds, dso;
#if VET_VGD_SWEEP_SLAB
            {
              // slab: the upwind cell-centre plane, or the inner face from the first
              // shell; chords dx1/|n_1|
              typ = (!inw && l == lcb) ? 2 : 1;
              iu = inw ? (i + 1) : ((typ == 2) ? i : (i - 1));
              const Real amr = fabs(mr);
              ds = fabs(((typ == 2) ? cx1f(m,i) : cx1v(m,iu)) - r)/amr;
              dso = (typ == 2) ? 0.5*ds
                               : fmin(fabs(cx1f(m,inw ? (i + 1) : i) - r)/amr, ds);
              if (img) {
                // the image of the other end's cell centre, one period away
                typ = 1;
                iu = inw ? is : ie;
                dso = fabs(cx1f(m,inw ? (ie + 1) : is) - r)/amr;
                ds = dso + (inw ? (cx1v(m,is) - cx1f(m,is))
                                : (cx1f(m,ie+1) - cx1v(m,ie)))/amr;
              }
            }
#else
            ScRayChord(inw, l, lcb, m, i, zr, p2, cx1v, cx1f, spl, typ, iu, ds, dso);
#endif
            // upwind point, its lateral cell and bilinear fractions (index space);
            // fractional upwind index at its global magnitude (sj_ = sk_ = 0 without
            // layout_invariant): stencil base and fractions are layout independent
            const int sj_ = lg_.Sj(m, js), sk_ = lg_.Sk(m, ks);
#if VET_VGD_SWEEP_SLAB
            // slab: a constant offset per direction (2-D: none in x3); r_hat = x1_hat
            // (the inner-face mirror)
            const Real xu = 1.0, yu = 0.0, zu = 0.0, ru = 1.0;
            const Real fjg = (j + sj_) - ds*ny/lg_.Dx2(mbsize, m);
            const Real fkg = (k + sk_) - ((nx3b > 1) ? ds*nz/lg_.Dx3(mbsize, m) : 0.0);
#else
            Real xu, yu, zu, ru, thu, fjg, fkg;
            {
              xu = r*st*cp - ds*nx;
              yu = r*st*sp - ds*ny;
              zu = r*ct - ds*nz;
              ru = sqrt(xu*xu + yu*yu + zu*zu);
              thu = acos(fmin(fmax(zu/ru, -1.0), 1.0));
              Real dph = atan2(yu, xu) - ph;
              dph -= twopi*floor((dph + M_PI)/twopi);
              // behind the pole -> extended chart; fractional row from the node table
              ScPoleChart(pol_, th, thu, dph);
              fjg = sth_ ? ScThFj(thn_, m, j, thu, c2w_, sj_)
                         : (j + sj_) + (thu - th)/lg_.Dx2(mbsize, m);
              fkg = (k + sk_) + dph/lg_.Dx3(mbsize, m);
            }
#endif
            const Real fj = fjg - sj_, fk = fkg - sk_;
            int j0 = static_cast<int>(floor(fjg)) - sj_;
            int k0 = static_cast<int>(floor(fkg)) - sk_;
            Real uj = fjg - floor(fjg), uk = fkg - floor(fkg);
            bool clp = false;
            // full sphere: reads stop at the shell's halo depth (clamped and counted)
            int jlo_ = jlo, jhi_ = jhi, klo_ = klo, khi_ = khi;
            if (pol_) {
              // filled halo depth: inward band of shell iu (inward pass, turned ray),
              // else outward band; on every side: decomposition independent
              const int wu = wsd_((inw || typ == 3) ? 0 : 1, iu);
              klo_ = wb - wu;
              khi_ = wb + nx3b + wu - 2;
              jlo_ = wb - wu;
              jhi_ = wb + nx2b + wu - 2;
              // vet_gd_band_global: the cell's window (layout independent; inside the
              // block's band since wb <= k, j < wb + nx)
              if (pol_ == 2) {
                klo_ = k - wu;
                khi_ = k + wu - 1;
                jlo_ = j - wu;
                jhi_ = j + wu - 1;
              }
            }
            if (j0 < jlo_) {j0 = jlo_; uj = 0.0; clp = true;}
            if (j0 > jhi_) {j0 = jhi_; uj = 1.0; clp = true;}
            if (k0 < klo_) {k0 = klo_; uk = 0.0; clp = true;}
            if (k0 > khi_) {k0 = khi_; uk = 1.0; clp = true;}
            if (clp) {Kokkos::atomic_add(&cnt_(0), 1.0);}
            int j1 = j0 + 1, k1 = (VET_VGD_SWEEP_SLAB && nx3b == 1) ? k0 : (k0 + 1);
            // interior launch: an item whose stencil reaches the band waits for it
            if (dfr && (k0 - xr_ < wb || k1 + xr_ >= wb + nx3b || j0 - xr_ < wb ||
                        j1 + xr_ >= wb + nx2b)) {
              dfd = true;
              return 0.0;
            }
            // stencil and field reads of the shared upwind core
            const ScSt4 st4{k0, k1, j0, j1, uj, uk};
            auto fC = [&](const int c, int kk, int jj, const int ii) -> Real {
#if VET_VGD_BODY_MIR
              // slab, beyond a vacuum x2 (x3) face (flag 8 + 1 (x2) + 2 (x3)): the edge
              // cell's chi and S
              const int wf = VET_VGD_SWEEP_SLAB ? wl_(m,kk,jj) : 0;
              if (wf >= 8 && (wf & 1)) {jj = (jj < wb) ? wb : (wb + nx2b - 1);}
              if (wf >= 8 && (wf & 2)) {kk = (kk < wb) ? wb : (wb + nx3b - 1);}
#endif
              return csx(m,c,kk,jj,ii);
            };
            // active stencil cells straight from shell iu's slab (the element rd reads)
            int uc3;
            const int64_t ub = vix.Base(m, iu, uc3);
            auto fIu = [&](const int kk, const int jj) -> Real {
              if (kk >= wb && kk < wb + nx3b && jj >= wb && jj < wb + nx2b) {
                return vix.d(ub + d + static_cast<int64_t>(vix.n)*(kk + uc3*jj));
              }
              return rd(kk,jj,iu);
            };
            // branch n . r_hat of a stencil cell
            auto mrq = [&](const int kk, const int jj) -> Real {return mr_(m,kk,jj,d);};
            bool bdone = false;
#if VET_VGD_BODY_MIR && !VET_VGD_SWEEP_SLAB
            // theta_wall = vacuum, upwind point beyond the wall (flag 3): the ray
            // starts at its exact crossing of the cone theta = theta_wall with I = 0;
            // chi, S from the row next to the wall (linear in phi index and r)
            const bool vlo = (thu < mt0), vhi = (thu > mt0 + mtl);
            const int jrw = static_cast<int>(fmin(fmax(floor(fjg + 0.5) - sj_, Real(jlo)),
                                                  Real(jhi + 1)));
            if (typ != 2 && (vlo || vhi) && wl_(m,k,jrw) == 3) {
              const bool hi = vhi;
              const Real thw = hi ? (mt0 + mtl) : mt0;
              // band index of the row next to the wall
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
              if (tc < 0.0) {   // round-off: the index-linear estimate
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
                return ScLerp(csx(m,c,kx0,jw,ii), csx(m,c,kx0+1,jw,ii), ukx);
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
#endif
            if (bdone) {
              // segment from the theta wall
#if VET_VGD_SWEEP_FS
              gq[SPFS_XV] = (typ == 1) ? vix(m,d,k,j,iu)
                                       : ((l > lcut) ? vix(m,d,k,j,i-1) : s0);
#endif
            } else if (typ == 2 && (emi || psrc)) {
              // vet_gd_inner = emit: boundary emission at the inner-face point; x1
              // periodic first pass: the face S
#if VET_VGD_SWEEP_SLAB
              const Real ie0 = psrc ? fmax(1.5*s0 - 0.5*exp(csx(m,1,k,j,i+1)), 0.0)
                                    : iemit;
              Real fc = 0.0;
              iv = ScInnerFace(st4, fC, i, true, ie0 + bmv(), bfc, cx1f(m,i), cx1v(m,i),
                               cx1v(m,i+1), p2, ch0, ds, s0, -1.0, 0.0,
                               sfs_ ? &fc : nullptr);
              gq[0] = dir_(d,3)*nx*fc;
#else
              iv = ScInnerFace(st4, fC, i, true, iemit, bfc, cx1f(m,i), cx1v(m,i),
                               cx1v(m,i+1), p2, ch0, ds, s0, -1.0, 0.0,
                               VET_VGD_FSP(SPFS_XV));
#endif
            } else if (typ == 2 && imir) {
              // reflecting inner face: first-shell inward intensity in the nearest
              // specular direction (inward stencil cells) times the flux ratio
              // (ScInnerMirror)
              const Real ux = xu/ru, uy = yu/ru, uz = zu/ru;
              const Real mu = nx*ux + ny*uy + nz*uz;
              int dm = 0;
              Real best = -2.0;
              for (int q = 0; q < n; ++q) {   // nearest inward one, as the ratio
                const Real c = dir_(q,0)*(nx - 2.0*mu*ux) + dir_(q,1)*(ny - 2.0*mu*uy) +
                               dir_(q,2)*(nz - 2.0*mu*uz);
                const Real mq = dir_(q,0)*ux + dir_(q,1)*uy + dir_(q,2)*uz;
                if (mq < 0.0 && c > best + 1.0e-12) {best = c; dm = q;}
              }
              bool ok[4];
              bool all = true;
              for (int q = 0; q < 4; ++q) {
                ok[q] = (mr_(m,(q < 2) ? k0 : k1,(q % 2 == 0) ? j0 : j1,dm) < 0.0);
                all = all && ok[q];
              }
              const Real im = ScUpwindI(st4, [&](const int kk, const int jj) {
                return rdq(kk,jj,i,dm);
              }, ok, all);
              iv = ScInnerFace(st4, fC, i, true, (im < 0.0) ? s0 : isc_(m,k-ks,j-js)*im,
                               bfc, cx1f(m,i), cx1v(m,i), cx1v(m,i+1), p2, ch0, ds, s0,
                               -1.0, 0.0, VET_VGD_FSP(SPFS_XV));
            } else if (typ == 2) {
              // diffusion intensity S + mu (-dS/dtau_r) at the inner face; vet_bc_diff =
              // face: face S = 1.5 S_a - 0.5 S_b (cell S is O(dtau_cell) low)
#if VET_VGD_SWEEP_SLAB
              Real fc = 0.0;
              iv = ScInnerFace(st4, fC, i, false, 0.0, bfc, cx1f(m,i), cx1v(m,i),
                               cx1v(m,i+1), p2, ch0, ds, s0, fabs(mr), bmv(),
                               sfs_ ? &fc : nullptr);
              gq[0] = dir_(d,3)*nx*fc;
#else
              iv = ScInnerFace(st4, fC, i, false, 0.0, bfc, cx1f(m,i), cx1v(m,i),
                               cx1v(m,i+1), p2, ch0, ds, s0, -1.0, 0.0,
                               VET_VGD_FSP(SPFS_XV));
#endif
            } else {
              // one intensity per (cell, direction), of the pass n_d . r_hat selects;
              // stencil cells of the other pass are dropped, weights renormalised
              const bool want_in = inw || (typ == 3);
              bool ok[4];
              bool all = true;
              // active stencil far enough from this cell's tangent: every corner's
              // branch is the cell's own (|n.r_hat' - n.r_hat| <= |dtheta| + |dphi| on
              // a uniform grid; the slab's n.r_hat is the same everywhere), so the
              // tests below would give exactly these values
              bool sure = false;
              if (!sth_ && !pol_ && k0 >= wb && k1 < wb + nx3b && j0 >= wb &&
                  j1 < wb + nx2b) {
                const Real mo = mr_(m,k,j,d);
                const int djm = (j1 - j > j - j0) ? (j1 - j) : (j - j0);
                const int dkm = (k1 - k > k - k0) ? (k1 - k) : (k - k0);
                const Real bd = VET_VGD_SWEEP_SLAB ? 0.0 : (fabs(djm*lg_.Dx2(mbsize, m)) +
                                              fabs(dkm*lg_.Dx3(mbsize, m)))*(1.0 + 1.0e-9)
                                             + 1.0e-12;
                sure = fabs(mo) > bd;
                if (sure) {
                  all = want_in ? (mo < 0.0) : (mo >= 0.0);
                  for (int q = 0; q < 4; ++q) {ok[q] = all;}
                }
              }
              for (int q = 0; !sure && q < 4; ++q) {
                const int kk = (q < 2) ? k0 : k1;
                const int jj = (q % 2 == 0) ? j0 : j1;
                const Real mq = mrq(kk,jj);
                ok[q] = want_in ? (mq < 0.0) : (mq >= 0.0);
                all = all && ok[q];
              }
              // none (-1): start from the upwind source
              const Real ivu = ScUpwindI(st4, fIu, ok, all);
              auto lin = [&](const int c) {
                return ScLerp(ScLerp(fC(c,k0,j0,iu), fC(c,k0,j1,iu), uj),
                              ScLerp(fC(c,k1,j0,iu), fC(c,k1,j1,iu), uj), uk);
              };
              const Real chu = exp(lin(0));
              const Real su = exp(lin(1));
              const Real rfi = cx1f(m,i);
              // closing segment
              auto edgef = [&]() {
                auto csr = [&](const int, const int c, const int kk, const int jj,
                               const int ii) -> Real {return fC(c,kk,jj,ii);};
                return ScCornerEdge(csr, m, k0, k1, j0, j1, iu, uk, uj);
              };
              // jump plateau test: the shell one further upwind (in the block)
              const int iu2 = 2*iu - i;
              const int iuu = (typ == 1 && iu2 >= is && iu2 <= ie) ? iu2 : iu;
              const int id2 = 2*i - iu;
              const int idd = (typ == 1 && id2 >= is && id2 <= ie) ? id2 : i;
#if VET_VGD_SWEEP_SLAB
              Real fc = 0.0;
              iv = ScSegFinish(st4, fIu, fC, edgef, ok, iu, k, j, i, ivu, chu, su, ch0,
                               s0, ds, dso, spl, spe, typ, l, lcb, p2, rfi, zr, vim, iuu,
                               idd, sfs_ ? &fc : nullptr);
              // x1: the intensity where the segment crosses the face
              gq[1] = dir_(d,3)*nx*fmax(fc, 0.0);
#else
              iv = ScSegFinish(st4, fIu, fC, edgef, ok, iu, k, j, i, ivu, chu, su, ch0,
                               s0, ds, dso, spl, spe, typ, l, lcb, p2, rfi, zr, vim, iuu,
                               idd, VET_VGD_FSP(SPFS_XV));
#endif
#if VET_VGD_SWEEP_FS
              // a turned ray crosses no r face of its own: its r face carries the cell
              // below (the inner face's rule on the first shell)
              if (typ == 3) {
                gq[SPFS_XV] = (l > lcut) ? vix(m,d,k,j,i-1) : (emi ? iemit : s0);
              }
#endif
            }
          }
#if VET_VGD_SWEEP_SLAB
          // exit through an x1 boundary face (not periodic): the half cell at chi_0, S_0
          if (sfs_ && !x1pr_ && (inw ? ((x1f_(m) & 1) && l == lcb)
                                     : ((x1f_(m) & 2) && l == n1 - 1))) {
            const Real ex = exp(-ch0*fabs(cx1f(m,inw ? i : (i + 1)) - r)/fabs(mr));
            gq[2] = dir_(d,3)*nx*(iv*ex + (1.0 - ex)*s0);
          }
#endif
#if VET_VGD_SWEEP_FS
          // leaving through the inner (not past its tangent, else -1) or top face: the
          // chord from the centre at chi_0, S_0
          if (inw ? ((x1f_(m) & 1) && l == lcb) : ((x1f_(m) & 2) && l == n1 - 1)) {
            const Real rb = cx1f(m,inw ? i : (i + 1));
            Real xe = -1.0;
            if (!inw || p2 < rb*rb) {
              const Real t = inw ? (zr - sqrt(rb*rb - p2))
                                 : (sqrt(fmax(rb*rb - p2, 0.0)) - zr);
              const Real ex = exp(-ch0*fmax(t, 0.0));
              xe = iv*ex + (1.0 - ex)*s0;
            }
            rgs_((l & 1)*nmbt_ + m,SPFS_XE*n + d,k-ks,j-js) = xe;
          }
#endif
          return iv;
        };
        const Real v = solve(cs_, vi_);
        if (dfd) {return true;}
        vi_(m,d,k,j,i) = v;
#if VET_VGD_SWEEP_SLAB
        for (int q = 0; q < 3 && sfs_; ++q) {
          rgs_((l & 1)*nmbt_ + m,q*n + d,k-ks,j-js) = gq[q];
        }
#endif
#if VET_VGD_SWEEP_FS
        rgs_((l & 1)*nmbt_ + m,SPFS_XV*n + d,k-ks,j-js) = fmax(gq[SPFS_XV], 0.0);
#endif
        return false;
      };
