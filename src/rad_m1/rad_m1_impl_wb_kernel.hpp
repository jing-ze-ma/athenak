//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
// fragment included twice by rad_m1_implicit.cpp, inside the m1_impl_wb lambda:
// M1_RCP 0 = the default kernel (byte-identical to the code before
// implicit_realisable_coupling), M1_RCP 1 = the kernel with the key on.
// cs-open-1018: M1_CSW 1 (the cubed-sphere launches only) adds the cs blocks of the
// rad_inertia and conservative paths (`#if M1_CSW`); with M1_CSW undefined (the sp and
// Cartesian launches) the preprocessed kernel is exactly the one before.
// NOLINT(build/header_guard)
    // pgrst implicit_hold: a held cell keeps the Dirichlet E (= its old vector); its F
    // and gas belong to the problem generator; no stage slope
    if (hwo_ && hwm_(m,k,j,i) > 0.0) {
      u0_(m,M1_E,k,j,i) = iw_(m,M1_IW_EN,k,j,i);
      if (t2k) {
        for (int n = 0; n < M1_T2_F1; ++n) {kk_(m,n,k,j,i) = 0.0;}
        if (gcs && coupling_ && dbgf) {
          ckw_(m,0,k,j,i) = 0.0;
          ckw_(m,1,k,j,i) = 0.0;
          ckw_(m,2,k,j,i) = 0.0;
        }
      }
      return;
    }
    Real ep = iw_(m,M1_IW_EP,k,j,i);
    // F3 (implicit_pos_floor_solve): the iterate is floored only at 1e-6 e_floor in the
    // T-solve / accept; the energy that residual raise created over a solved S2 below it
    // (gas units, per volume; rare) is charged like the e_floor raise below
    Real s2def = 0.0;
    if (pfs) {
      const Real s2 = iw_(m,M1_IW_S2,k,j,i);
      if (!(ep > 1.0e-6*efl) && s2 < ep) {
        s2def = (cl/ch)*(1.0 + iw_(m,M1_IW_SRCB,k,j,i))*(ep - s2);
      }
    }
    // an M1_IBC_EFIX end cell: its row was replaced by E' = EN (the Dirichlet value),
    // and that is the E it keeps -- the work term below must not move it.  It used to:
    // the replaced row has no transverse coupling, so nothing damped the x2 structure
    // the per-step work kicks left in the end column, and at the OUTFLOW end of the
    // radiative shock (T7) it grew into an x2 mode (E, F2, gas v2) until the end cell
    // hit the E floor (tests_m1/runs_4d_efix).
    const int ipw = pos_.d_view(m);
    const bool efc = !cyclic && ((i == is && ipw == 0 && bclo == M1_IBC_EFIX) ||
                                 (i == ie && ipw == nblkx1-1 && bchi == M1_IBC_EFIX));
    Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
    Real fp1 = 0.5*(fl + fr) + iw_(m,M1_IW_ADV,k,j,i)*ep;
    // MILESTONE 3b phase B: the derived cell-centred transverse fluxes, the face means
    // plus the enthalpy flux, exactly as in x1.
    Real fp2 = 0.0, fp3 = 0.0;
    Real g2l = 0.0, g2r = 0.0, g3l = 0.0, g3r = 0.0;
    if (trans) {
      g2l = f2_(m,k,j,i);
      g2r = f2_(m,k,j+1,i);
      fp2 = 0.5*(g2l + g2r) + iw_(m,M1_IW_A2,k,j,i)*ep;
      if (thrd) {
        g3l = f3_(m,k,j,i);
        g3r = f3_(m,k+1,j,i);
        fp3 = 0.5*(g3l + g3r) + iw_(m,M1_IW_A3,k,j,i)*ep;
      }
    }
    if (csw && trans && thrd) {
      const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
      const Real a = fp2, b = fp3;
      fp2 = (a + c*b)*si;
      fp3 = (b + c*a)*si;
    }

    Real work = 0.0, dm1 = 0.0, dmref = 0.0, eg = 0.0, ekin = 0.0, egrv = 0.0;
    Real dm2 = 0.0, dm3 = 0.0;
    Real go1 = 0.0, go2 = 0.0, go3 = 0.0;   // consrc-1011: F^old of the cell
    Real dd = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0;
    if (have_hydro) {
      dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      v1 = uh(m,IM1,k,j,i)*idd;
      v2 = uh(m,IM2,k,j,i)*idd;
      v3 = uh(m,IM3,k,j,i)*idd;
      ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                  SQR(uh(m,IM3,k,j,i)))*idd;
      if (csw) {
        const Real c = cclw(m,k,j), id2 = 1.0/(1.0 - c*c);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        v2 = (m2 - c*m3)*idd*id2;   // contravariant
        v3 = (m3 - c*m2)*idd*id2;
        ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + (m2*m2 + m3*m3 - 2.0*c*m2*m3)*id2)*idd;
      }
      egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      if (mhd) egrv += emag_(m,k,j,i);
      eg = iw_(m,M1_IW_EGN,k,j,i);
      // (a) the energy the RADIATION gained from the gas over the step.  It is taken
      // from the ASSEMBLED row, q = SRCR - SRCB E', and not from rho kappa_P a T'^4:
      // the two differ by the Picard remainder of the linearisation, and only the first
      // is the amount the solved E actually received.  Setting the gas energy from it
      // makes e_gas + (c/chat) E change by the face fluxes and the work term ALONE, to
      // round-off -- measured: the T'^4 form drifted 2.8e-11 over 2000 steps of T5, this
      // one 0.  At convergence the two agree, so T' stays the consistent temperature.
      if (coupling_ && dbgh) {
        Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*ep;
        eg -= (cl/ch)*qq;
      }
      // (b) MOMENTUM.  Each x1 face hands dt (rho k_t)_f F0'_f/c to the gas, half to
      // each of its two cells (a physical boundary face gives all of it to its one
      // interior cell), which is exactly what the implicit face source removed from the
      // radiation: sum_cells dm = sum_faces dt (rho k_t)_f F0'_f/c.
      if (coupling_ && dbgf && !gcs) {
        Real ktl, ktr;
        Real wl = 0.5, wr = 0.5;
        int ipos = pos_.d_view(m);
        if (i == is && ipos == 0 && !cyclic) {
          ktl = iw_(m,M1_IW_KT,k,j,i);
          wl = bmhalf ? 0.5 : 1.0;
        } else {
          int im = (cyclic && i == is) ? ie : (i-1);
          ktl = M1KtX1(iw_, kir_, m, k, j, im, i, cx1f, fwd);
        }
        if (i == ie && ipos == nblkx1-1 && !cyclic) {
          ktr = iw_(m,M1_IW_KT,k,j,i);
          wr = bmhalf ? 0.5 : 1.0;
        } else {
          int ip = (cyclic && i == ie) ? is : (i+1);
          ktr = M1KtX1(iw_, kir_, m, k, j, i, ip, cx1f, fwd);
        }
        dm1 = (dt/cl)*(wl*ktl*fl + wr*ktr*fr);
        dmref = fref ? (dt*dd*aref_(m,k,j,i)) : 0.0;
        if (trans && dbgft) {
          // the same rule per transverse direction: each face hands
          // dt (rho k_t)_f F0_f/c to the gas, half to each of its two cells, and a
          // PHYSICAL boundary face (where F0 is zero anyway) half as well under
          // implicit_bmom_half.
          BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
          BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
          bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
          bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
          Real ktc = iw_(m,M1_IW_KT,k,j,i);
          Real kl2 = (j == js && p2lo) ? ktc
                     : 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + ktc);
          Real kr2 = (j == je && p2hi) ? ktc
                     : 0.5*(ktc + iw_(m,M1_IW_KT,k,j+1,i));
          Real u2 = ((j == js && p2lo) && !bmhalf) ? 1.0 : 0.5;
          Real w2 = ((j == je && p2hi) && !bmhalf) ? 1.0 : 0.5;
          dm2 = (dt/cl)*(u2*kl2*g2l + w2*kr2*g2r);
          if (thrd) {
            BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
            BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
            bool p3lo = (q5 != BoundaryFlag::block) && (q5 != BoundaryFlag::periodic);
            bool p3hi = (q6 != BoundaryFlag::block) && (q6 != BoundaryFlag::periodic);
            Real kl3 = (k == ks && p3lo) ? ktc
                       : 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + ktc);
            Real kr3 = (k == ke && p3hi) ? ktc
                       : 0.5*(ktc + iw_(m,M1_IW_KT,k+1,j,i));
            Real u3 = ((k == ks && p3lo) && !bmhalf) ? 1.0 : 0.5;
            Real w3 = ((k == ke && p3hi) && !bmhalf) ? 1.0 : 0.5;
            dm3 = (dt/cl)*(u3*kl3*g3l + w3*kr3*g3r);
          }
          if (csw && thrd) {
            const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
            const Real a = dm2, b = dm3;
            dm2 = (a + c*b)*si;   // covariant, as the hydro momentum
            dm3 = (b + c*a)*si;
          }
        }
        // moveq-1011 (implicit_rad_inertia): the face equation also removes
        // chat dt v_f g0_f from the radiation flux; the gas gets dt v_f g0_f/c, half per
        // interior face (as the k F0 deposit; a physical face has no face equation).
        // Exactly 0 where v = 0 (static gas: bitwise)
        if (rinr) {
          int ipq = pos_.d_view(m);
          if (!(i == is && ipq == 0 && !cyclic)) {
            const int im = (cyclic && i == is) ? ie : (i-1);
            Real g = M1FaceAvgX1(iw_(m,M1_IW_V1,k,j,im), iw_(m,M1_IW_V1,k,j,i), cx1f, m,
                                 im, i, fwd)
                     *M1FaceAvgX1(iw_(m,M1_IW_G0,k,j,im), iw_(m,M1_IW_G0,k,j,i), cx1f, m,
                                  im, i, fwd);
            if (aphll) {g *= (1.0 - ifw_(m,M1_IFW_AL,k,j,i));}
            if (g != 0.0) {dm1 += 0.5*(dt/cl)*g;}
          }
          if (!(i == ie && ipq == nblkx1-1 && !cyclic)) {
            const int ip = (cyclic && i == ie) ? is : (i+1);
            Real g = M1FaceAvgX1(iw_(m,M1_IW_V1,k,j,i), iw_(m,M1_IW_V1,k,j,ip), cx1f, m,
                                 i, ip, fwd)
                     *M1FaceAvgX1(iw_(m,M1_IW_G0,k,j,i), iw_(m,M1_IW_G0,k,j,ip), cx1f, m,
                                  i, ip, fwd);
            if (aphll) {g *= (1.0 - ifw_(m,M1_IFW_AL,k,j,i+1));}
            if (g != 0.0) {dm1 += 0.5*(dt/cl)*g;}
          }
          if (trans && dbgft && !csw) {
            // x2 (and x3): the arithmetic face means of the transverse face equation
            auto vg = [&](const int c, const int k0, const int j0, const int k1,
                          const int j1) {
              return 0.125*(dt/cl)*(iw_(m,c,k0,j0,i) + iw_(m,c,k1,j1,i))
                     *(iw_(m,M1_IW_G0,k0,j0,i) + iw_(m,M1_IW_G0,k1,j1,i));
            };
            BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
            BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
            const bool p2lo = (q3 != BoundaryFlag::block) &&
                              (q3 != BoundaryFlag::periodic);
            const bool p2hi = (q4 != BoundaryFlag::block) &&
                              (q4 != BoundaryFlag::periodic);
            Real g = (j == js && p2lo) ? 0.0 : vg(M1_IW_V2, k, j-1, k, j);
            if (!(j == je && p2hi)) {g += vg(M1_IW_V2, k, j, k, j+1);}
            if (g != 0.0) {dm2 += g;}
            if (thrd) {
              BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
              BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
              const bool p3lo = (q5 != BoundaryFlag::block) &&
                                (q5 != BoundaryFlag::periodic);
              const bool p3hi = (q6 != BoundaryFlag::block) &&
                                (q6 != BoundaryFlag::periodic);
              g = (k == ks && p3lo) ? 0.0 : vg(M1_IW_V3, k-1, j, k, j);
              if (!(k == ke && p3hi)) {g += vg(M1_IW_V3, k, j, k+1, j);}
              if (g != 0.0) {dm3 += g;}
            }
          }
#if M1_CSW
          // cs-open-1018: the same transverse deposit on the cubed sphere.  The face
          // means are those of the sp rule (face-normal V of the cell and its neighbour,
          // a seam ghost re-expressed on this block's face normals; m1bcs: a seam is an
          // interior face), i.e. the face-normal momentum of the face equation's own
          // v_f g0_f term, mapped to the covariant dm2, dm3 as the k F0 deposit above
          if (trans && dbgft && thrd) {
            auto vgc = [&](const int c, const int k0, const int j0, const int k1,
                           const int j1) {
              return 0.125*(dt/cl)*(iw_(m,c,k0,j0,i) + iw_(m,c,k1,j1,i))
                     *(iw_(m,M1_IW_G0,k0,j0,i) + iw_(m,M1_IW_G0,k1,j1,i));
            };
            const BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
            const BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
            const BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
            const BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
            const bool p2lo = (q3 != BoundaryFlag::block) &&
                              (q3 != BoundaryFlag::periodic);
            const bool p2hi = (q4 != BoundaryFlag::block) &&
                              (q4 != BoundaryFlag::periodic);
            const bool p3lo = (q5 != BoundaryFlag::block) &&
                              (q5 != BoundaryFlag::periodic);
            const bool p3hi = (q6 != BoundaryFlag::block) &&
                              (q6 != BoundaryFlag::periodic);
            Real ga = (j == js && p2lo) ? 0.0 : vgc(M1_IW_V2, k, j-1, k, j);
            if (!(j == je && p2hi)) {ga += vgc(M1_IW_V2, k, j, k, j+1);}
            Real gb = (k == ks && p3lo) ? 0.0 : vgc(M1_IW_V3, k-1, j, k, j);
            if (!(k == ke && p3hi)) {gb += vgc(M1_IW_V3, k, j, k+1, j);}
            if (ga != 0.0 || gb != 0.0) {
              const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
              dm2 += (ga + c*gb)*si;
              dm3 += (gb + c*ga)*si;
            }
          }
#endif
        }
#if M1_RCP
        // implicit_realisable_coupling: the deposit of the realisable cell flux, scaled
        // by s = min(1, c E/|F_iter|) (the factor the write-back clip applies to F); the
        // work below follows dm, so the radiation loses exactly what the gas gains.  On
        // the cubed sphere fp2, fp3 and dm2, dm3 are covariant here: |F| is the clip's
        // metric norm (M1FluxNormCs) and the covariant dm takes the same scalar factor
        {
          const Real fq = csw ? M1FluxNormCs(cclw(m,k,j), fp1, fp2, fp3)
                              : sqrt(fp1*fp1 + fp2*fp2 + fp3*fp3);
          const Real fc = cl*fmax(ep, efl);
          if (fq > fc) {
            const Real sc = fc/fq;
            dm1 *= sc;
            dm2 *= sc;
            dm3 *= sc;
          }
        }
#endif
        if (feedback) {
          Real idg = 1.0/fmax(dd, 1.0e-300);
          Real w1 = (uh(m,IM1,k,j,i) + dm1)*idg;
          work = 0.5*(v1 + w1)*dm1;
          if (trans && dbgft) {
            Real w2n = (uh(m,IM2,k,j,i) + dm2)*idg;
            work += 0.5*(v2 + w2n)*dm2;
            if (thrd) {
              Real w3n = (uh(m,IM3,k,j,i) + dm3)*idg;
              work += 0.5*(v3 + w3n)*dm3;
            }
            if (csw && thrd) {
              // v^a dm_a with the contravariant v^a before and after the kick
              const Real c = cclw(m,k,j), id2 = 1.0/(1.0 - c*c);
              const Real n2 = uh(m,IM2,k,j,i) + dm2, n3 = uh(m,IM3,k,j,i) + dm3;
              const Real u2 = (n2 - c*n3)*idg*id2, u3 = (n3 - c*n2)*idg*id2;
              work = 0.5*(v1 + w1)*dm1 + 0.5*(v2 + u2)*dm2 + 0.5*(v3 + u3)*dm3;
            }
          }
          // moveq-1011 (implicit_rad_inertia): the work moves (chat/c) W of energy into
          // E after the solve, and the cell's enthalpy flux a E carries its momentum
          // a W/c^2: the gas pays it (else the next solve, whose E^n includes W, counts
          // it as created).  W depends on the kick: a few fixed-point passes (the
          // factor is v a/c^2).  Off where the row already carried the work (vstage)
          const Real ax = iw_(m,M1_IW_ADV,k,j,i);
          const Real ay = (trans && dbgft) ? iw_(m,M1_IW_A2,k,j,i) : 0.0;
          const Real az = (trans && dbgft && thrd) ? iw_(m,M1_IW_A3,k,j,i) : 0.0;
          if (rinr && !rwk && !efc && !csw && (ax != 0.0 || ay != 0.0 || az != 0.0)) {
            const Real ic2 = 1.0/(cl*cl);
            Real wq = work;
            Real x1 = dm1, x2 = dm2, x3 = dm3;
            for (int q = 0; q < 3; ++q) {
              x1 = dm1 + ax*wq*ic2;
              x2 = dm2 + ay*wq*ic2;
              x3 = dm3 + az*wq*ic2;
              wq = 0.5*(v1 + (uh(m,IM1,k,j,i) + x1)*idg)*x1;
              if (trans && dbgft) {
                wq += 0.5*(v2 + (uh(m,IM2,k,j,i) + x2)*idg)*x2;
                if (thrd) {wq += 0.5*(v3 + (uh(m,IM3,k,j,i) + x3)*idg)*x3;}
              }
            }
            dm1 = x1;
            dm2 = x2;
            dm3 = x3;
            work = wq;
            // the cell flux takes the a E of the E it ends on
            fp1 -= ax*(ch/cl)*work;
            fp2 -= ay*(ch/cl)*work;
            fp3 -= az*(ch/cl)*work;
          }
#if M1_CSW
          // cs-open-1018: the same fixed point on the cubed sphere.  (ay, az) are the
          // face-normal components of the enthalpy-flux coefficient a (A = a E), so its
          // momentum a W/c^2 enters the covariant dm and the covariant cell F through
          // the face-normal -> covariant map; the kinetic energy takes the contravariant
          // velocity (as kew of the conservative path)
          if (rinr && !rwk && !efc && thrd && (ax != 0.0 || ay != 0.0 || az != 0.0)) {
            const Real ic2 = 1.0/(cl*cl);
            const Real c = cclw(m,k,j), sn = csnw(m,k,j), id2 = 1.0/(1.0 - c*c);
            const Real by = (ay + c*az)/sn, bz = (az + c*ay)/sn;
            Real wq = work;
            Real x1 = dm1, x2 = dm2, x3 = dm3;
            for (int q = 0; q < 3; ++q) {
              x1 = dm1 + ax*wq*ic2;
              x2 = dm2 + by*wq*ic2;
              x3 = dm3 + bz*wq*ic2;
              wq = 0.5*(v1 + (uh(m,IM1,k,j,i) + x1)*idg)*x1;
              const Real n2 = uh(m,IM2,k,j,i) + x2, n3 = uh(m,IM3,k,j,i) + x3;
              const Real u2 = (n2 - c*n3)*idg*id2, u3 = (n3 - c*n2)*idg*id2;
              wq += 0.5*(v2 + u2)*x2 + 0.5*(v3 + u3)*x3;
            }
            dm1 = x1;
            dm2 = x2;
            dm3 = x3;
            work = wq;
            fp1 -= ax*(ch/cl)*work;
            fp2 -= by*(ch/cl)*work;
            fp3 -= bz*(ch/cl)*work;
          }
#endif
          if (!efc) {
            ep -= (ch/cl)*work;
          }
        } else if (nfbw) {
          // moveq-1011 (implicit_nofb_work): gas_feedback = false is a gas of infinite
          // inertia, not a gas at rest: the radiation still does the work v.dm on it
          // (the lab-frame v.G of the E equation), at the fixed v.  Without it a pulse
          // diffusing in a moving medium is advected at (4/3) v (enthalpy flux alone)
          Real wnf = v1*dm1;
          if (trans && dbgft) {
            wnf += v2*dm2;
            if (thrd) {wnf += v3*dm3;}
          }
          if (!efc && wnf != 0.0 && !nfwr) {
            ep -= (ch/cl)*wnf;
          }
        }
      }
      if (coupling_ && dbgf && gcs) {
        // consrc-1011, implicit_gas_src = conservative (gas_feedback = true): the gas
        // takes minus the radiation's own source of the solve,
        //   dm = -(F'_c - F^old_c)/(c chat) - (1/2) sum_{interior faces} Z_f
        //        + sum_{physical x1 faces} [(F0'_b - F0old_b)/(2 c chat) + w dt k F0'/c],
        // F_c the STORED cell flux (u0; F^old of a hesdirk2 stage + its increment) and
        // Z_f the face's transport part (face kernels).  A physical face has no face
        // equation: its F0 change is inflow, not a source.  The work W (exact kinetic
        // energy change of dm) leaves E as before; F' depends on dm through the clip
        // (and the vimp re-forming), so a few fixed-point passes, the last dm being the
        // one of the stored F' (exact conservation)
        const Real icc = 1.0/(cl*ch);
        const int ipq = pos_.d_view(m);
        Real x1 = 0.0, x2 = 0.0, x3 = 0.0;
        // a physical face: its F0 change is inflow, and it keeps the matched drag
        // deposit of the boundary half cell (dt k F0_b/c, a boundary momentum source)
        const Real wbq = bmhalf ? 0.5 : 1.0;
        if (i == is && ipq == 0 && !cyclic) {
          x1 += 0.5*(fl - f0n_(m,k,j,i))*icc + wbq*(dt/cl)*iw_(m,M1_IW_KT,k,j,i)*fl;
        } else {
          x1 -= 0.5*gz1w_(m,k,j,i);
        }
        if (i == ie && ipq == nblkx1-1 && !cyclic) {
          x1 += 0.5*(fr - f0n_(m,k,j,i+1))*icc + wbq*(dt/cl)*iw_(m,M1_IW_KT,k,j,i)*fr;
        } else {
          x1 -= 0.5*gz1w_(m,k,j,i+1);
        }
        const bool tq = trans && dbgft;
        if (tq) {
          BoundaryFlag q3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
          BoundaryFlag q4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
          const bool p2lo = (q3 != BoundaryFlag::block) && (q3 != BoundaryFlag::periodic);
          const bool p2hi = (q4 != BoundaryFlag::block) && (q4 != BoundaryFlag::periodic);
          if (!(j == js && p2lo)) {x2 -= 0.5*gz2w_(m,k,j,i);}
          if (!(j == je && p2hi)) {x2 -= 0.5*gz2w_(m,k,j+1,i);}
          if (thrd) {
            BoundaryFlag q5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
            BoundaryFlag q6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
            const bool p3lo = (q5 != BoundaryFlag::block) &&
                              (q5 != BoundaryFlag::periodic);
            const bool p3hi = (q6 != BoundaryFlag::block) &&
                              (q6 != BoundaryFlag::periodic);
            if (!(k == ks && p3lo)) {x3 -= 0.5*gz3w_(m,k,j,i);}
            if (!(k == ke && p3hi)) {x3 -= 0.5*gz3w_(m,k+1,j,i);}
          }
          if (csw && thrd) {
            const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
            const Real a = x2, b = x3;
            x2 = (a + c*b)*si;   // covariant, as the hydro momentum and u0 F
            x3 = (b + c*a)*si;
          }
        }
        // F^old: the stored cell flux of the old vector
        go1 = u0_(m,M1_F1,k,j,i);
        go2 = u0_(m,M1_F2,k,j,i);
        go3 = u0_(m,M1_F3,k,j,i);
        if (gfst) {
          // the first solve of the run: the problem generator's F is not (yet) the
          // radiation's state (F0 faces + a E), whose old value is taken instead (the
          // difference is a cold-start mismatch, not a source)
          const Real eo = iw_(m,M1_IW_EN,k,j,i);
          go1 = 0.5*(f0n_(m,k,j,i) + f0n_(m,k,j,i+1)) + iw_(m,M1_IW_ADV,k,j,i)*eo;
          if (trans) {
            Real a = 0.5*(f2nw_(m,k,j,i) + f2nw_(m,k,j+1,i)) + iw_(m,M1_IW_A2,k,j,i)*eo;
            Real b = thrd ? (0.5*(f3nw_(m,k,j,i) + f3nw_(m,k+1,j,i))
                             + iw_(m,M1_IW_A3,k,j,i)*eo) : 0.0;
            if (csw && thrd) {
              const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
              const Real a0 = a;
              a = (a0 + c*b)*si;
              b = (b + c*a0)*si;
            }
            go2 = a;
            go3 = b;
          }
        }
        if (t2st) {
          if (csok) {
            go1 += ciw_(m,0,k,j,i);
            go2 += ciw_(m,1,k,j,i);
            go3 += ciw_(m,2,k,j,i);
          } else {
            // no cell-F increment (restart, multirate): the face means of the F0
            // increments plus a times the E increment
            const Real tie = t2i_(m,M1_T2_E,k,j,i);
            go1 += 0.5*(t2i_(m,M1_T2_F1,k,j,i) + t2i_(m,M1_T2_F1,k,j,i+1))
                   + iw_(m,M1_IW_ADV,k,j,i)*tie;
            if (trans) {
              go2 += 0.5*(t2i_(m,M1_T2_F1+1,k,j,i) + t2i_(m,M1_T2_F1+1,k,j+1,i))
                     + iw_(m,M1_IW_A2,k,j,i)*tie;
              if (thrd) {
                go3 += 0.5*(t2i_(m,M1_T2_F1+2,k,j,i) + t2i_(m,M1_T2_F1+2,k+1,j,i))
                       + iw_(m,M1_IW_A3,k,j,i)*tie;
              }
            }
          }
        }
#if M1_CSW
        // cs-open-1018: the stage increment without a cell-F increment is the face-normal
        // face mean + A tie: its transverse part to the covariant basis of u0 F
        if (t2st && !csok && trans && thrd) {
          const Real tie = t2i_(m,M1_T2_E,k,j,i);
          const Real a = 0.5*(t2i_(m,M1_T2_F1+1,k,j,i) + t2i_(m,M1_T2_F1+1,k,j+1,i))
                         + iw_(m,M1_IW_A2,k,j,i)*tie;
          const Real b = 0.5*(t2i_(m,M1_T2_F1+2,k,j,i) + t2i_(m,M1_T2_F1+2,k+1,j,i))
                         + iw_(m,M1_IW_A3,k,j,i)*tie;
          const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
          go2 += (a + c*b)*si - a;
          go3 += (b + c*a)*si - b;
        }
#endif
        const Real idg = 1.0/fmax(dd, 1.0e-300);
        // the exact kinetic-energy change of the kick (the matched path's work)
        auto kew = [&](const Real a1, const Real a2, const Real a3) {
          Real w = 0.5*(v1 + (uh(m,IM1,k,j,i) + a1)*idg)*a1;
          if (tq) {
            if (csw && thrd) {
              const Real c = cclw(m,k,j), id2 = 1.0/(1.0 - c*c);
              const Real n2 = uh(m,IM2,k,j,i) + a2, n3 = uh(m,IM3,k,j,i) + a3;
              const Real u2 = (n2 - c*n3)*idg*id2, u3 = (n3 - c*n2)*idg*id2;
              w += 0.5*(v2 + u2)*a2 + 0.5*(v3 + u3)*a3;
            } else {
              w += 0.5*(v2 + (uh(m,IM2,k,j,i) + a2)*idg)*a2;
              if (thrd) {w += 0.5*(v3 + (uh(m,IM3,k,j,i) + a3)*idg)*a3;}
            }
          }
          return w;
        };
        dmref = fref ? (dt*dd*aref_(m,k,j,i)) : 0.0;
        // gsdrag-1021 implicit_cs2_dacap: the drag of the matched rule (dt k_f F0'_f/c,
        // half per face, as above) and the bound of the conservative excess over it,
        // |dm - dm_drag| <= 2 K |v| max(E^old, E')/(c chat), |v| >= vmin c: the
        // radiation-inertia share a real enthalpy flux can hand over
        Real gd1 = 0.0, gd2 = 0.0, gd3 = 0.0, gbd = 0.0;
        if (gcap) {
          const Real ktc = iw_(m,M1_IW_KT,k,j,i);
          Real ktl = ktc, ktr = ktc, wl = 0.5, wr = 0.5;
          if (i == is && ipq == 0 && !cyclic) {
            wl = bmhalf ? 0.5 : 1.0;
          } else {
            const int im = (cyclic && i == is) ? ie : (i-1);
            ktl = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,im), ktc, cx1f, m, im, i, fwd);
          }
          if (i == ie && ipq == nblkx1-1 && !cyclic) {
            wr = bmhalf ? 0.5 : 1.0;
          } else {
            const int ip = (cyclic && i == ie) ? is : (i+1);
            ktr = M1FaceAvgX1(ktc, iw_(m,M1_IW_KT,k,j,ip), cx1f, m, i, ip, fwd);
          }
          gd1 = (dt/cl)*(wl*ktl*fl + wr*ktr*fr);
          if (tq) {
            BoundaryFlag b3 = mbbcs.d_view(m,BoundaryFace::inner_x2);
            BoundaryFlag b4 = mbbcs.d_view(m,BoundaryFace::outer_x2);
            const bool l2 = (b3 != BoundaryFlag::block) && (b3 != BoundaryFlag::periodic);
            const bool h2 = (b4 != BoundaryFlag::block) && (b4 != BoundaryFlag::periodic);
            const Real kl2 = (j == js && l2) ? ktc : 0.5*(iw_(m,M1_IW_KT,k,j-1,i) + ktc);
            const Real kr2 = (j == je && h2) ? ktc : 0.5*(ktc + iw_(m,M1_IW_KT,k,j+1,i));
            const Real u2 = ((j == js && l2) && !bmhalf) ? 1.0 : 0.5;
            const Real w2 = ((j == je && h2) && !bmhalf) ? 1.0 : 0.5;
            gd2 = (dt/cl)*(u2*kl2*g2l + w2*kr2*g2r);
            if (thrd) {
              BoundaryFlag b5 = mbbcs.d_view(m,BoundaryFace::inner_x3);
              BoundaryFlag b6 = mbbcs.d_view(m,BoundaryFace::outer_x3);
              const bool l3 = (b5 != BoundaryFlag::block) &&
                              (b5 != BoundaryFlag::periodic);
              const bool h3 = (b6 != BoundaryFlag::block) &&
                              (b6 != BoundaryFlag::periodic);
              const Real kl3 = (k == ks && l3) ? ktc
                                               : 0.5*(iw_(m,M1_IW_KT,k-1,j,i) + ktc);
              const Real kr3 = (k == ke && h3) ? ktc
                                               : 0.5*(ktc + iw_(m,M1_IW_KT,k+1,j,i));
              const Real u3 = ((k == ks && l3) && !bmhalf) ? 1.0 : 0.5;
              const Real w3 = ((k == ke && h3) && !bmhalf) ? 1.0 : 0.5;
              gd3 = (dt/cl)*(u3*kl3*g3l + w3*kr3*g3r);
            }
            if (csw && thrd) {
              const Real c = cclw(m,k,j), si = 1.0/csnw(m,k,j);
              const Real a = gd2, b = gd3;
              gd2 = (a + c*b)*si;
              gd3 = (b + c*a)*si;
            }
          }
#if M1_RCP
          {
            // implicit_realisable_coupling: the matched drag of the realisable cell flux
            // (the factor min(1, c E/|F_iter|) of the matched deposit above)
            const Real fq = csw ? M1FluxNormCs(cclw(m,k,j), fp1, fp2, fp3)
                                : sqrt(fp1*fp1 + fp2*fp2 + fp3*fp3);
            const Real fc = cl*fmax(ep, efl);
            if (fq > fc) {
              const Real sc = fc/fq;
              gd1 *= sc;
              gd2 *= sc;
              gd3 *= sc;
            }
          }
#endif
          const Real vq = fmax(sqrt(v1*v1 + v2*v2 + v3*v3), gvmw);
          gbd = gkw*vq*fmax(fmax(iw_(m,M1_IW_EN,k,j,i), ep), 0.0)*icc;
        }
        const Real r1 = fp1, r2 = fp2, r3 = fp3;   // the cell flux before the kick
        Real q1 = fp1, q2 = fp2, q3 = fp3;
        // cs2clip-1020 (implicit_cs2_clip = rad): the gas is charged with the cell flux
        // BEFORE the realisability limit (u*); the limit's change stays with the radiation
        // (as on the matched path, which limits after the deposit).  gas: the limited q.
#if defined(M1_C2CR) && M1_C2CR
        Real u1 = fp1, u2 = fp2, u3 = fp3;
#endif
        // pgrst-vimp-1021 implicit_cs2_wbsolve: the fixed point of the passes below in
        // closed form (vstage, not on the cubed sphere): per active component a,
        // dm_a = C_a - om_a icc e_a with e = F0 mean + E (I + D) w, w = (p + dm)/rho
        // -> (I + G) dm = R, G_ab = om_a icc E M_ab / rho, M = I + D
        const bool wcf = wbs && vfx && !csw;
        if (wcf) {
          const int na = tq ? (thrd ? 3 : 2) : 1;
          const Real es = iw_(m,M1_IW_EP,k,j,i);
          const Real chi = iw_(m,M1_IW_WCHI,k,j,i);
          const Real nn1 = iw_(m,M1_IW_N1,k,j,i), nn2 = iw_(m,M1_IW_N2,k,j,i);
          const Real nn3 = iw_(m,M1_IW_N3,k,j,i);
          Real dg[3] = {M1EddDiag(chi,nn1), M1EddDiag(chi,nn2), M1EddDiag(chi,nn3)};
          Real dof[3] = {M1EddOff(chi,nn1,nn2), M1EddOff(chi,nn1,nn3),
                         M1EddOff(chi,nn2,nn3)};
          if (dfull) {
            for (int a = 0; a < 3; ++a) {dg[a] = vd_(m,M1_VET_D11+a,k,j,i);}
            for (int a = 0; a < 3; ++a) {dof[a] = vd_(m,M1_VET_D11+3+a,k,j,i);}
          }
          Real M[3][3];
          M[0][0] = 1.0 + dg[0]; M[1][1] = 1.0 + dg[1]; M[2][2] = 1.0 + dg[2];
          M[0][1] = M[1][0] = dof[0]; M[0][2] = M[2][0] = dof[1];
          M[1][2] = M[2][1] = dof[2];
          const Real fm[3] = {0.5*(fl + fr), 0.5*(g2l + g2r), 0.5*(g3l + g3r)};
          const Real xa[3] = {x1, x2, x3}, goa[3] = {go1, go2, go3};
          Real om[3], cc[3], w0[3];
          w0[0] = (uh(m,IM1,k,j,i) - dmref)*idg;
          w0[1] = tq ? uh(m,IM2,k,j,i)*idg : v2;
          w0[2] = (tq && thrd) ? uh(m,IM3,k,j,i)*idg : v3;
          for (int a = 0; a < 3; ++a) {
            om[a] = c2h ? iw_(m,c2a+6+a,k,j,i) : 1.0;
            cc[a] = xa[a] + goa[a]*icc;
            if (c2h) {cc[a] -= (1.0 - om[a])*(fm[a] + iw_(m,c2a+a,k,j,i))*icc;}
          }
          Real A[3][3], R[3];
          for (int a = 0; a < na; ++a) {
            Real mw = 0.0;
            for (int b = 0; b < 3; ++b) {mw += M[a][b]*w0[b];}
            R[a] = cc[a] - om[a]*icc*(fm[a] + es*mw);
            for (int b = 0; b < na; ++b) {
              A[a][b] = ((a == b) ? 1.0 : 0.0) + om[a]*icc*es*M[a][b]*idg;
            }
          }
          Real dm[3] = {0.0, 0.0, 0.0};
          if (na == 1) {
            dm[0] = R[0]/A[0][0];
          } else if (na == 2) {
            const Real det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
            dm[0] = (R[0]*A[1][1] - A[0][1]*R[1])/det;
            dm[1] = (A[0][0]*R[1] - A[1][0]*R[0])/det;
          } else {
            const Real det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
                           - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
                           + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
            for (int c = 0; c < 3; ++c) {
              Real B[3][3];
              for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {B[a][b] = (b == c) ? R[a] : A[a][b];}
              }
              dm[c] = (B[0][0]*(B[1][1]*B[2][2] - B[1][2]*B[2][1])
                     - B[0][1]*(B[1][0]*B[2][2] - B[1][2]*B[2][0])
                     + B[0][2]*(B[1][0]*B[2][1] - B[1][1]*B[2][0]))/det;
            }
          }
          Real wv[3] = {w0[0], w0[1], w0[2]};
          for (int a = 0; a < na; ++a) {wv[a] += dm[a]*idg;}
          Real ev[3];
          for (int a = 0; a < 3; ++a) {
            ev[a] = fm[a] + es*(M[a][0]*wv[0] + M[a][1]*wv[1] + M[a][2]*wv[2]);
          }
          dm1 = dm[0]; dm2 = dm[1]; dm3 = dm[2];
          q1 = ev[0];
          q2 = trans ? ev[1] : r2;
          q3 = (trans && thrd) ? ev[2] : r3;
          const Real wq = kew(dm1, dm2, dm3);
          Real eq = efc ? ep : (ep - (ch/cl)*wq);
          M1ApplyLimits(cl, efl, eq, q1, q2, q3);
#if !(defined(M1_C2CR) && M1_C2CR)
          // the gas charged with the limited flux (implicit_cs2_clip = gas): one pass
          {
            const Real qa[3] = {q1, q2, q3};
            Real dmn[3] = {0.0, 0.0, 0.0};
            for (int a = 0; a < na; ++a) {dmn[a] = cc[a] - om[a]*icc*qa[a];}
            dm1 = dmn[0]; dm2 = dmn[1]; dm3 = dmn[2];
          }
#endif
        }
        for (int q = 0; q < 4 && !wcf; ++q) {
#if defined(M1_C2CR) && M1_C2CR
          const Real e1 = u1, e2 = u2, e3 = u3;
#else
          const Real e1 = q1, e2 = q2, e3 = q3;
#endif
          dm1 = x1 - (e1 - go1)*icc;
          dm2 = tq ? (x2 - (e2 - go2)*icc) : 0.0;
          dm3 = (tq && thrd) ? (x3 - (e3 - go3)*icc) : 0.0;
#if M1_CSW
          // cs-open-1018: the same drag rule on the cubed sphere.  A^old (iw_cs2 + 1..2,
          // ImplicitCs2DeltaA) and the face means are face-normal, q the covariant cell
          // F: the transverse dA is formed face-normal, weighted per axis, and mapped
          // to the covariant dm
          if (c2h) {
            dm1 += (1.0 - iw_(m,c2a+6,k,j,i))
                   *(e1 - 0.5*(fl + fr) - iw_(m,c2a,k,j,i))*icc;
            if (tq && thrd) {
              const Real c = cclw(m,k,j), sn = csnw(m,k,j);
              const Real qa = (e2 - c*e3)/sn, qb = (e3 - c*e2)/sn;
              const Real ta = (1.0 - iw_(m,c2a+7,k,j,i))
                              *(qa - 0.5*(g2l + g2r) - iw_(m,c2a+1,k,j,i))*icc;
              const Real tb = (1.0 - iw_(m,c2a+8,k,j,i))
                              *(qb - 0.5*(g3l + g3r) - iw_(m,c2a+2,k,j,i))*icc;
              dm2 += (ta + c*tb)/sn;
              dm3 += (tb + c*ta)/sn;
            }
          }
#else
          if (c2h) {
            // consfix3 (implicit_cs2_hll = drag): the AP-HLL part of dA stays with the
            // radiation (no face equation there): the gas pays om_bar dA, not dA
            dm1 += (1.0 - iw_(m,c2a+6,k,j,i))
                   *(e1 - 0.5*(fl + fr) - iw_(m,c2a,k,j,i))*icc;
            if (tq) {
              dm2 += (1.0 - iw_(m,c2a+7,k,j,i))
                     *(e2 - 0.5*(g2l + g2r) - iw_(m,c2a+1,k,j,i))*icc;
              if (thrd) {
                dm3 += (1.0 - iw_(m,c2a+8,k,j,i))
                       *(e3 - 0.5*(g3l + g3r) - iw_(m,c2a+2,k,j,i))*icc;
              }
            }
          }
#endif
          if (gcap) {
            // gsdrag-1021: the conservative excess over the drag, capped (the rest of
            // it stays with the radiation)
            const Real x1q = dm1 - gd1, x2q = dm2 - gd2, x3q = dm3 - gd3;
            const Real xn = sqrt(x1q*x1q + x2q*x2q + x3q*x3q);
            if (xn > gbd) {
              const Real sc = gbd/xn;
              dm1 = gd1 + sc*x1q;
              dm2 = gd2 + sc*x2q;
              dm3 = gd3 + sc*x3q;
            }
          }
          if (q == 3) break;
          const Real wq = kew(dm1, dm2, dm3);
          Real eq = efc ? ep : (ep - (ch/cl)*wq);
          q1 = r1;  q2 = r2;  q3 = r3;
          if (vfx && !csw) {
            // m1-sp-order2b: the cell F at the velocity the kick gives the gas
            const Real w1 = (uh(m,IM1,k,j,i) + dm1 - dmref)*idg;
            Real w2 = v2, w3 = v3;
            if (tq) {
              w2 = (uh(m,IM2,k,j,i) + dm2)*idg;
              if (thrd) {w3 = (uh(m,IM3,k,j,i) + dm3)*idg;}
            }
            const Real chi = iw_(m,M1_IW_WCHI,k,j,i);
            const Real n1 = iw_(m,M1_IW_N1,k,j,i), n2 = iw_(m,M1_IW_N2,k,j,i);
            const Real n3 = iw_(m,M1_IW_N3,k,j,i);
            Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2);
            Real d33 = M1EddDiag(chi,n3);
            Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
            Real d23 = M1EddOff(chi,n2,n3);
            if (dfull) {
              d11 = vd_(m,M1_VET_D11,k,j,i);
              d22 = vd_(m,M1_VET_D11+1,k,j,i);
              d33 = vd_(m,M1_VET_D11+2,k,j,i);
              d12 = vd_(m,M1_VET_D11+3,k,j,i);
              d13 = vd_(m,M1_VET_D11+4,k,j,i);
              d23 = vd_(m,M1_VET_D11+5,k,j,i);
            }
            const Real es = iw_(m,M1_IW_EP,k,j,i);
            q1 = 0.5*(fl + fr) + (w1 + (w1*d11 + w2*d12 + w3*d13))*es;
            if (trans) {
              q2 = 0.5*(g2l + g2r) + (w2 + (w1*d12 + w2*d22 + w3*d23))*es;
              if (thrd) {
                q3 = 0.5*(g3l + g3r) + (w3 + (w1*d13 + w2*d23 + w3*d33))*es;
              }
            }
          }
#if M1_CSW
          // cs-open-1018: the same re-forming on the cubed sphere (as the m1_impl_wb_csvfx
          // kernel of the matched path): face-normal velocities s v^a of the kicked
          // covariant momentum, the face-normal A, then the covariant cell F
          if (vfx && trans && thrd) {
            const Real c = cclw(m,k,j), sn = csnw(m,k,j), s2 = 1.0 - c*c;
            const Real w1 = (uh(m,IM1,k,j,i) + dm1 - dmref)*idg;
            const Real p2 = uh(m,IM2,k,j,i) + (tq ? dm2 : 0.0);
            const Real p3 = uh(m,IM3,k,j,i) + (tq ? dm3 : 0.0);
            const Real u2 = sn*(p2 - c*p3)*idg/s2, u3 = sn*(p3 - c*p2)*idg/s2;
            const Real chi = iw_(m,M1_IW_WCHI,k,j,i);
            const Real n1 = iw_(m,M1_IW_N1,k,j,i), n2 = iw_(m,M1_IW_N2,k,j,i);
            const Real n3 = iw_(m,M1_IW_N3,k,j,i);
            Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2);
            Real d33 = M1EddDiag(chi,n3);
            Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
            Real d23 = M1EddOff(chi,n2,n3);
            if (dfull) {
              d11 = vd_(m,M1_VET_D11,k,j,i);
              d22 = vd_(m,M1_VET_D11+1,k,j,i);
              d33 = vd_(m,M1_VET_D11+2,k,j,i);
              d12 = vd_(m,M1_VET_D11+3,k,j,i);
              d13 = vd_(m,M1_VET_D11+4,k,j,i);
              d23 = vd_(m,M1_VET_D11+5,k,j,i);
            }
            const Real es = iw_(m,M1_IW_EP,k,j,i);
            q1 = 0.5*(fl + fr) + (w1 + (w1*d11 + u2*d12 + u3*d13))*es;
            const Real a2 = 0.5*(g2l + g2r) + (u2 + (w1*d12 + u2*d22 + u3*d23))*es;
            const Real a3 = 0.5*(g3l + g3r) + (u3 + (w1*d13 + u2*d23 + u3*d33))*es;
            q2 = (a2 + c*a3)/sn;
            q3 = (a3 + c*a2)/sn;
          }
#endif
#if defined(M1_C2CR) && M1_C2CR
          u1 = q1;  u2 = q2;  u3 = q3;   // cs2clip-1020: the pre-limit cell flux
#endif
          if (csw) {
            M1ApplyLimitsCs(cl, efl, cclw(m,k,j), eq, q1, q2, q3);
          } else {
            M1ApplyLimits(cl, efl, eq, q1, q2, q3);
          }
        }
        fp1 = q1;
        fp2 = q2;
        fp3 = q3;
        work = kew(dm1, dm2, dm3);
        if (!efc) {
          ep -= (ch/cl)*work;
        }
      }
    }

    // m1-sp-order2b: a hesdirk2 stage solve under implicit_vimp.  The derived cell flux
    // F = F0 + a E above takes a = v + v.D of the STAGE-START velocity, lagged by the
    // radiative kick of the stage (the solved E already carries the implicit v'): an
    // O(dt) error in the cell F (rsw_u_T: F order 1.0, E, T, gas 1.9-2.0).  It is
    // re-formed with the velocity the write-back below gives the gas.
    if (vfx && have_hydro && feedback && !gcs) {
      Real idg = 1.0/fmax(dd, 1.0e-300);
      Real w1 = (uh(m,IM1,k,j,i) + dm1 - dmref)*idg;
      Real w2 = v2, w3 = v3;
      if (trans && dbgft) {
        w2 = (uh(m,IM2,k,j,i) + dm2)*idg;
        if (thrd) {
          w3 = (uh(m,IM3,k,j,i) + dm3)*idg;
        }
      }
      Real chi = iw_(m,M1_IW_WCHI,k,j,i);
      Real n1 = iw_(m,M1_IW_N1,k,j,i), n2 = iw_(m,M1_IW_N2,k,j,i);
      Real n3 = iw_(m,M1_IW_N3,k,j,i);
      Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2), d33 = M1EddDiag(chi,n3);
      Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
      Real d23 = M1EddOff(chi,n2,n3);
      if (dfull) {
        d11 = vd_(m,M1_VET_D11,k,j,i);
        d22 = vd_(m,M1_VET_D11+1,k,j,i);
        d33 = vd_(m,M1_VET_D11+2,k,j,i);
        d12 = vd_(m,M1_VET_D11+3,k,j,i);
        d13 = vd_(m,M1_VET_D11+4,k,j,i);
        d23 = vd_(m,M1_VET_D11+5,k,j,i);
      }
      // E of the solve (as the stage-start form above): the rows carry the work of
      // the last pass (ImplicitWorkRow), so it is the stage's E to the tolerance
      const Real es = iw_(m,M1_IW_EP,k,j,i);
      fp1 = 0.5*(fl + fr) + (w1 + (w1*d11 + w2*d12 + w3*d13))*es;
      if (trans) {
        fp2 = 0.5*(g2l + g2r) + (w2 + (w1*d12 + w2*d22 + w3*d23))*es;
        if (thrd) {
          fp3 = 0.5*(g3l + g3r) + (w3 + (w1*d13 + w2*d23 + w3*d33))*es;
        }
      }
    }
    // m1-positivity.  Both moves keep e_gas + (c/chat) E of the cell exactly; eg is the
    // gas internal energy the write-back sets (the kinetic part is ekin + work).
    //  implicit_pos_floor: an E below e_floor is raised to it with the energy taken from
    //    the gas (down to its limit below); what the gas cannot cover is counted.
    //  implicit_pos_gas: a gas eint below e_min = max(e(rho, tfloor),
    //    implicit_pos_gas_frac x max(e_gas old, 0)) is raised to it with energy taken
    //    from the radiation (down to e_floor).
    if (pany && !efc) {
      const Real vol = psph ? pvol(m,k,j,i) : (mbsize.d_view(m).dx1*mbsize.d_view(m).dx2*
                                               mbsize.d_view(m).dx3);
      Real emin = 0.0;
      if (have_hydro) {
        emin = pgf*fmax(iw_(m,M1_IW_EGN,k,j,i), 0.0);
        if (peos.tfloor > 0.0) {
          Real te, pp, cr, ct, tcv;
          peos.ThermoAt(dd, peos.tfloor, te, pp, cr, ct, tcv);
          emin = fmax(emin, te);
        }
      }
      if (pflr && (!(ep > efl) || s2def > 0.0)) {
        const Real need = (cl/ch)*fmax(efl - ep, 0.0) + s2def;   // gas units
        const Real avail = (have_hydro && feedback) ? fmax(eg - emin, 0.0) : 0.0;
        const Real take = fmin(need, avail);
        if (have_hydro && feedback) {eg -= take;}
        ep = fmax(ep, efl);
        Kokkos::atomic_add(&pc_(M1_POS_FLR), 1.0);
        Kokkos::atomic_add(&pc_(M1_POS_FLR_DE), take*vol);
        Kokkos::atomic_add(&pc_(M1_POS_FLR_UN), (need - take)*vol);
        if (pfs) {Kokkos::atomic_add(&pc_(M1_POS_FLRS), s2def*vol);}
      }
      if (pgas && feedback && eg < emin) {
        const Real need = emin - eg;                     // gas units
        const Real take = fmin(need, fmax((cl/ch)*(ep - efl), 0.0));
        eg += take;
        ep -= (ch/cl)*take;
        Kokkos::atomic_add(&pc_(M1_POS_GAS), 1.0);
        Kokkos::atomic_add(&pc_(M1_POS_GAS_DE), take*vol);
      }
    }
    {   // diagnostic counters only: the flux scaled back to |F| = c E by M1ApplyLimits
      const Real fq = sqrt(fp1*fp1 + fp2*fp2 + fp3*fp3);
      if (ep > efl && fq > cl*ep) {
        Kokkos::atomic_add(&pc_(M1_POS_FCLIP), 1.0);
        Kokkos::atomic_add(&pc_(M1_POS_FCLIPM), fq/(cl*ep) - 1.0);
      }
    }
    const Real cb1 = fp1, cb2 = fp2, cb3 = fp3;
    if (csw) {
      M1ApplyLimitsCs(cl, efl, cclw(m,k,j), ep, fp1, fp2, fp3);
    } else {
      M1ApplyLimits(cl, efl, ep, fp1, fp2, fp3);
    }
    if (gcs && coupling_ && dbgf && (fp1 != cb1 || fp2 != cb2 || fp3 != cb3)) {
      // consrc-1011: a clip the positivity moves made bind after the kick (rare): the
      // gas takes the clipped momentum too, its kinetic-energy change from eint (the
      // total gas energy is unchanged)
      const Real icc = 1.0/(cl*ch), idg = 1.0/fmax(dd, 1.0e-300);
      const bool tq = trans && dbgft;
      const Real a1 = dm1 - (fp1 - cb1)*icc;
      const Real a2 = tq ? (dm2 - (fp2 - cb2)*icc) : dm2;
      const Real a3 = (tq && thrd) ? (dm3 - (fp3 - cb3)*icc) : dm3;
      Real wn = 0.5*(v1 + (uh(m,IM1,k,j,i) + a1)*idg)*a1;
      if (tq) {
        if (csw && thrd) {
          const Real c = cclw(m,k,j), id2 = 1.0/(1.0 - c*c);
          const Real n2 = uh(m,IM2,k,j,i) + a2, n3 = uh(m,IM3,k,j,i) + a3;
          const Real u2 = (n2 - c*n3)*idg*id2, u3 = (n3 - c*n2)*idg*id2;
          wn += 0.5*(v2 + u2)*a2 + 0.5*(v3 + u3)*a3;
        } else {
          wn += 0.5*(v2 + (uh(m,IM2,k,j,i) + a2)*idg)*a2;
          if (thrd) {wn += 0.5*(v3 + (uh(m,IM3,k,j,i) + a3)*idg)*a3;}
        }
      }
      eg -= (wn - work);
      work = wn;
      dm1 = a1;
      dm2 = a2;
      dm3 = a3;
    }
    // hesdirk2: the slope of this solve, K = (Y - old vector)/dt_solve
    if (t2k) {
      const Real fk = 1.0/dt;
      kk_(m,M1_T2_E,k,j,i) = (ep - iw_(m,M1_IW_EN,k,j,i))*fk;
      if (eso) {kk_(m,M1_T2_E,k,j,i) += dtes*fk*es_(m,k,j,i);}
      bool gk = have_hydro && feedback;
      kk_(m,M1_T2_M1,k,j,i) = gk ? ((dm1 - dmref)*fk) : 0.0;
      kk_(m,M1_T2_M1+1,k,j,i) = (gk && trans && dbgft) ? (dm2*fk) : 0.0;
      kk_(m,M1_T2_M1+2,k,j,i) = (gk && trans && dbgft && thrd) ? (dm3*fk) : 0.0;
      kk_(m,M1_T2_EN,k,j,i) = gk ? ((eg + ekin + egrv + work - uh(m,IEN,k,j,i))*fk)
                                 : 0.0;
      if (gcs && coupling_ && dbgf) {
        // consrc-1011: the slope of the stored cell F, which the next stage's old
        // vector takes (Time2FormStage -> csinc)
        ckw_(m,0,k,j,i) = (fp1 - go1)*fk;
        ckw_(m,1,k,j,i) = (fp2 - go2)*fk;
        ckw_(m,2,k,j,i) = (fp3 - go3)*fk;
      }
    }
    u0_(m,M1_E,k,j,i) = ep;
    u0_(m,M1_F1,k,j,i) = fp1;
    u0_(m,M1_F2,k,j,i) = fp2;
    u0_(m,M1_F3,k,j,i) = fp3;
    if (have_hydro && feedback) {
      uh(m,IM1,k,j,i) = uh(m,IM1,k,j,i) + dm1 - dmref;
      if (trans && dbgft) {
        uh(m,IM2,k,j,i) = uh(m,IM2,k,j,i) + dm2;
        if (thrd) {uh(m,IM3,k,j,i) = uh(m,IM3,k,j,i) + dm3;}
      }
      uh(m,IEN,k,j,i) = eg + ekin + egrv + work;
    }
