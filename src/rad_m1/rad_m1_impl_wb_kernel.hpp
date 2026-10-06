//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
// fragment included twice by rad_m1_implicit.cpp, inside the m1_impl_wb lambda:
// M1_RCP 0 = the default kernel (byte-identical to the code before
// implicit_realisable_coupling), M1_RCP 1 = the kernel with the key on.
// NOLINT(build/header_guard)
    Real ep = iw_(m,M1_IW_EP,k,j,i);
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
      if (coupling_ && dbgf) {
        Real ktl, ktr;
        Real wl = 0.5, wr = 0.5;
        int ipos = pos_.d_view(m);
        if (i == is && ipos == 0 && !cyclic) {
          ktl = iw_(m,M1_IW_KT,k,j,i);
          wl = bmhalf ? 0.5 : 1.0;
        } else {
          int im = (cyclic && i == is) ? ie : (i-1);
          ktl = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,im), iw_(m,M1_IW_KT,k,j,i), cx1f, m, im, i,
                            fwd);
        }
        if (i == ie && ipos == nblkx1-1 && !cyclic) {
          ktr = iw_(m,M1_IW_KT,k,j,i);
          wr = bmhalf ? 0.5 : 1.0;
        } else {
          int ip = (cyclic && i == ie) ? is : (i+1);
          ktr = M1FaceAvgX1(iw_(m,M1_IW_KT,k,j,i), iw_(m,M1_IW_KT,k,j,ip), cx1f, m, i, ip,
                            fwd);
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
          if (!efc) {
            ep -= (ch/cl)*work;
          }
        }
      }
    }

    // m1-sp-order2b: a hesdirk2 stage solve under implicit_vimp.  The derived cell flux
    // F = F0 + a E above takes a = v + v.D of the STAGE-START velocity, lagged by the
    // radiative kick of the stage (the solved E already carries the implicit v'): an
    // O(dt) error in the cell F (rsw_u_T: F order 1.0, E, T, gas 1.9-2.0).  It is
    // re-formed with the velocity the write-back below gives the gas.
    if (vfx && have_hydro && feedback) {
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
      if (pflr && !(ep > efl)) {
        const Real need = (cl/ch)*(efl - ep);           // gas units
        const Real avail = (have_hydro && feedback) ? fmax(eg - emin, 0.0) : 0.0;
        const Real take = fmin(need, avail);
        if (have_hydro && feedback) {eg -= take;}
        ep = efl;
        Kokkos::atomic_add(&pc_(M1_POS_FLR), 1.0);
        Kokkos::atomic_add(&pc_(M1_POS_FLR_DE), take*vol);
        Kokkos::atomic_add(&pc_(M1_POS_FLR_UN), (need - take)*vol);
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
    if (csw) {
      M1ApplyLimitsCs(cl, efl, cclw(m,k,j), ep, fp1, fp2, fp3);
    } else {
      M1ApplyLimits(cl, efl, ep, fp1, fp2, fp3);
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
