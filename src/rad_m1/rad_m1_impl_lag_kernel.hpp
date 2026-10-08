//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
// fragment included twice by rad_m1_implicit.cpp, inside the m1_impl_lag lambda:
// M1_RCP 0 = the default kernel (byte-identical to the code before
// implicit_realisable_coupling), M1_RCP 1 = the kernel with the key on.
// NOLINT(build/header_guard)
      Real e = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real f1 = iw_(m,M1_IW_F1,k,j,i);
      Real rf = f1/(cl*e);
      if (rf > 1.0) {rf = 1.0;}
      if (rf < -1.0) {rf = -1.0;}
      Real chi = edd ? (1.0/3.0) : M1Chi(fabs(rf), chk);
      Real v1 = iw_(m,M1_IW_V1,k,j,i);
      Real de0;
      if (!trans) {
        iw_(m,M1_IW_WCHI,k,j,i) = chi;
        iw_(m,M1_IW_ADV,k,j,i) = v1*(1.0 + chi);
        // E0 - E with F_2 = F_3 = 0: P_11 = chi E, P_22 = P_33 = (1-chi) E/2
        Real b1 = v1/cl;
#if M1_RCP
        // implicit_realisable_coupling: |F| <= c E in beta.F
        const Real f1r = fmin(fmax(f1, -cl*e), cl*e);
        de0 = ovc ? (-2.0*b1*f1r/cl) : (b1*b1*e - 2.0*b1*f1r/cl + b1*b1*chi*e);
#else
        de0 = ovc ? (-2.0*b1*f1/cl) : (b1*b1*e - 2.0*b1*f1/cl + b1*b1*chi*e);
#endif
      } else {
        // MILESTONE 3b phase B: the closure of the MULTI-DIMENSIONAL solve.  The reduced
        // flux is the MAGNITUDE |F|/(c E) and the Eddington tensor is built around the
        // unit flux direction n, D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b.  In
        // 1-D this is the branch above with n = (+-1,0,0), which is why the 1-D mesh is
        // handled there (gate G3) and nothing about implicit_x1 changes.
        Real f2c = iw_(m,M1_IW_F2,k,j,i);
        Real f3c = iw_(m,M1_IW_F3,k,j,i);
        Real fm = sqrt(f1*f1 + f2c*f2c + f3c*f3c);
        Real rfm = fm/(cl*e);
        if (rfm > 1.0) {rfm = 1.0;}
        chi = edd ? (1.0/3.0) : M1Chi(rfm, chk);
        Real ifm = 1.0/fmax(fm, 1.0e-300);
        Real n1 = f1*ifm, n2 = f2c*ifm, n3 = f3c*ifm;
        Real v2 = iw_(m,M1_IW_V2,k,j,i), v3 = iw_(m,M1_IW_V3,k,j,i);
        // MILESTONE 3b phase D: how fast the LAGGED closure is allowed to move.
        //   implicit_closure_lag = step  freezes (chi, n) at the start-of-step state for
        //     the whole step (the Eddington tensor is then explicit in time, as in a VET
        //     code that reuses the previous step's tensor), leaving ONE linear solve plus
        //     the temperature nonlinearity per step;
        //   implicit_closure_relax = w   under-relaxes them between passes, optionally
        //     only where the cell is optically thin (theta > 1/2), which is where the
        //     closure feeds back on the solve through (c dt/dx)^2.
        if (ttau) {
          // column optical depth to the top of the block at j-1, j, j+1 (cell centres)
          Real dx1 = mbsize.d_view(m).dx1;
          Real dx2 = mbsize.d_view(m).dx2;
          Real tm = 0.5*iw_(m,M1_IW_KT,k,j-1,i)*dx1;
          Real tc = 0.5*iw_(m,M1_IW_KT,k,j,i)*dx1;
          Real tq = 0.5*iw_(m,M1_IW_KT,k,j+1,i)*dx1;
          for (int ii = i+1; ii <= ie; ++ii) {
            tm += iw_(m,M1_IW_KT,k,j-1,ii)*dx1;
            tc += iw_(m,M1_IW_KT,k,j,ii)*dx1;
            tq += iw_(m,M1_IW_KT,k,j+1,ii)*dx1;
          }
          // exact grey plane-parallel K/J = (tau + q_inf)/(3 (tau + q(tau))), q = Hopf
          Real qh = 0.710446 - 0.133054*exp(-3.4488*tc);
          chi = (tc + 0.710446)/(3.0*(tc + qh));
          Real g1 = -iw_(m,M1_IW_KT,k,j,i);
          Real g2 = (tq - tm)/(2.0*dx2);
          Real ign = 1.0/fmax(sqrt(g1*g1 + g2*g2), 1.0e-300);
          n1 = -g1*ign;
          n2 = -g2*ign;
          n3 = 0.0;
        } else if (ttilt) {
          Real dx2 = mbsize.d_view(m).dx2;
          Real x2v = mbsize.d_view(m).x2min + (static_cast<Real>(j - js) + 0.5)*dx2;
          Real al = talp*sin(2.0*M_PI*(x2v - tx2min)/tx2len);
          Real ca = cos(al), sa = sin(al);
          Real r1 = ca*n1 - sa*n2, r2 = sa*n1 + ca*n2;
          n1 = r1;
          n2 = r2;
        }
        if (vetsc) {
          chi = vc_(m,M1_VET_CHI,k,j,i);
          n1 = vc_(m,M1_VET_N1,k,j,i);
          n2 = vc_(m,M1_VET_N1+1,k,j,i);
          n3 = vc_(m,M1_VET_N1+2,k,j,i);
        } else if (tkeep) {
          chi = iw_(m,M1_IW_WCHI,k,j,i);
          n1 = iw_(m,M1_IW_N1,k,j,i);
          n2 = iw_(m,M1_IW_N2,k,j,i);
          n3 = iw_(m,M1_IW_N3,k,j,i);
        } else if (dofreeze) {
          chi = iw_(m,M1_IW_WCHI,k,j,i);
          n1 = iw_(m,M1_IW_N1,k,j,i);
          n2 = iw_(m,M1_IW_N2,k,j,i);
          n3 = iw_(m,M1_IW_N3,k,j,i);
        } else if (dorel && (!crthin || (ch*dt*iw_(m,M1_IW_KT,k,j,i) < 1.0))) {
          Real w1 = 1.0 - crw;
          chi = w1*iw_(m,M1_IW_WCHI,k,j,i) + crw*chi;
          n1 = w1*iw_(m,M1_IW_N1,k,j,i) + crw*n1;
          n2 = w1*iw_(m,M1_IW_N2,k,j,i) + crw*n2;
          n3 = w1*iw_(m,M1_IW_N3,k,j,i) + crw*n3;
          Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
          if (nn > 0.0) {
            Real inn = 1.0/nn;
            n1 *= inn;
            n2 *= inn;
            n3 *= inn;
          }
        }
        if (tauc) {
          chi = tt_(m,0,k,j,i);
          n1 = tt_(m,1,k,j,i);
          n2 = tt_(m,2,k,j,i);
          n3 = tt_(m,3,k,j,i);
        }
        // implicit_closure_thin_relax (tests_m1/runs_5c_thinstab): the lagged closure's
        // dependence on F is explicit, and at c dt >> dx a flux perturbation comes back
        // amplified by g ~ max(chi', b/f)/tau_c each step (tau_c the cell optical depth
        // along n).  Relaxing (chi, n) from the previous step's values with
        // w = 2/(1 + G^2), G >= |g|, makes the step map contract (|1 - w + w g| < 1 for
        // Re g < 1) to the SAME fixed point; cells with G <= 1 are left untouched.
        if (ctr && !edd && !vetsc && !tkeep && !tauc) {
          if (ctri) {
            Real dx1 = mbsize.d_view(m).dx1;
            if (ctstr) {dx1 = ctx1f(m,i+1) - ctx1f(m,i);}
            Real kt = iw_(m,M1_IW_KT,k,j,i);
            Real tc = kt*dx1;
            if (!ctsph) {
              Real dx2 = mbsize.d_view(m).dx2;
              Real dx3 = mbsize.d_view(m).dx3;
              tc = kt/fmax(fabs(n1)/dx1 + fabs(n2)/dx2 + fabs(n3)/dx3, 1.0e-300);
            }
            Real fh = fmin(fmax(rfm, 1.0e-3), 0.999);
            Real cp = (M1Chi(fh + 1.0e-3, chk) - M1Chi(fh - 1.0e-3, chk))/2.0e-3;
            Real bf = 0.5*(3.0*chi - 1.0)/fh;
            Real gg = ctc*fmax(cp, bf)/fmax(tc, 1.0e-300)
                      + fh*cp/fmax(fmin(chi, 1.0 - chi), 1.0e-3);
            Real w = 2.0/(1.0 + gg*gg);
            if (w < 1.0) {
              Real w1 = 1.0 - w;
              chi = w1*cm_(m,0,k,j,i) + w*chi;
              n1 = w1*cm_(m,1,k,j,i) + w*n1;
              n2 = w1*cm_(m,2,k,j,i) + w*n2;
              n3 = w1*cm_(m,3,k,j,i) + w*n3;
              Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
              if (nn > 0.0) {
                Real inn = 1.0/nn;
                n1 *= inn;
                n2 *= inn;
                n3 *= inn;
              }
            }
          }
          cm_(m,0,k,j,i) = chi;
          cm_(m,1,k,j,i) = n1;
          cm_(m,2,k,j,i) = n2;
          cm_(m,3,k,j,i) = n3;
        }
        iw_(m,M1_IW_WCHI,k,j,i) = chi;
        iw_(m,M1_IW_N1,k,j,i) = n1;
        iw_(m,M1_IW_N2,k,j,i) = n2;
        iw_(m,M1_IW_N3,k,j,i) = n3;
        // a_d = v_d + (v.D)_d, so that the enthalpy flux A_d = v_d E + (v.P)_d = a_d E
        Real d11 = M1EddDiag(chi,n1), d22 = M1EddDiag(chi,n2), d33 = M1EddDiag(chi,n3);
        Real d12 = M1EddOff(chi,n1,n2), d13 = M1EddOff(chi,n1,n3);
        Real d23 = M1EddOff(chi,n2,n3);
        if (dfull) {
          // vet_tensor = full: the guarded K/J of the formal solution, all six components
          d11 = vd_(m,M1_VET_D11,k,j,i);
          d22 = vd_(m,M1_VET_D11+1,k,j,i);
          d33 = vd_(m,M1_VET_D11+2,k,j,i);
          d12 = vd_(m,M1_VET_D11+3,k,j,i);
          d13 = vd_(m,M1_VET_D11+4,k,j,i);
          d23 = vd_(m,M1_VET_D11+5,k,j,i);
        }
        iw_(m,M1_IW_ADV,k,j,i) = v1 + (v1*d11 + v2*d12 + v3*d13);
        iw_(m,M1_IW_A2,k,j,i) = v2 + (v1*d12 + v2*d22 + v3*d23);
        iw_(m,M1_IW_A3,k,j,i) = v3 + (v1*d13 + v2*d23 + v3*d33);
        if (t2dav) {
          // hesdirk2: the part of a carried by the old vector's velocity increment
          Real idg = 1.0/fmax(uh(m,IDN,k,j,i), 1.0e-300);
          Real w1 = t2i_(m,M1_T2_M1,k,j,i)*idg, w2 = t2i_(m,M1_T2_M1+1,k,j,i)*idg;
          Real w3 = t2i_(m,M1_T2_M1+2,k,j,i)*idg;
          iw_(m,t2da,k,j,i) = w1 + (w1*d11 + w2*d12 + w3*d13);
          iw_(m,t2da+1,k,j,i) = w2 + (w1*d12 + w2*d22 + w3*d23);
          iw_(m,t2da+2,k,j,i) = w3 + (w1*d13 + w2*d23 + w3*d33);
        }
        // E0 - E to O(beta^2), with the full pressure tensor.  m1-sp-order2b
        // (time2_vstage): in a hesdirk2 stage under implicit_vimp at the iterate's
        // velocity, v_old + dv^k (the previous pass's write-back increment), not the
        // stage-start one, which lags by the stage's radiative kick (O(dt))
        Real u1 = v1, u2 = v2, u3 = v3;
        if (vdv) {
          u1 += iw_(m,ivd,k,j,i);
          u2 += iw_(m,ivd+1,k,j,i);
          u3 += iw_(m,ivd+2,k,j,i);
        }
        Real b1 = u1/cl, b2 = u2/cl, b3 = u3/cl;
        Real bf = (b1*f1 + b2*f2c + b3*f3c)/cl;
#if M1_RCP
        // implicit_realisable_coupling: the realisable flux, |F| <= c E, in beta.F.  On
        // the cubed sphere f2c, f3c are face-normal: |F| is the clip's metric norm of the
        // covariant pair (a + c b)/s, (b + c a)/s (the write-back's map)
        Real fmr = fm;
        if (rcs) {
          const Real c = rccl(m,k,j), si = 1.0/rcsn(m,k,j);
          fmr = M1FluxNormCs(c, f1, (f2c + c*f3c)*si, (f3c + c*f2c)*si);
        }
        if (fmr > cl*e) {bf *= cl*e/fmr;}
#endif
        Real bpb = (b1*b1*d11 + b2*b2*d22 + b3*b3*d33
                    + 2.0*(b1*b2*d12 + b1*b3*d13 + b2*b3*d23))*e;
        Real b2sq = b1*b1 + b2*b2 + b3*b3;
        de0 = ovc ? (-2.0*bf) : (b2sq*e - 2.0*bf + bpb);
      }
      iw_(m,M1_IW_DE0,k,j,i) = de0;
      Real rkev = opac_(m,M1_OP_E,k,j,i);
      Real rkpv = opac_(m,M1_OP_P,k,j,i);
      Real tp = iw_(m,M1_IW_TP,k,j,i);
      Real t2 = tp*tp;
      iw_(m,M1_IW_G0,k,j,i) = rkev*(e + de0) - rkpv*ar*t2*t2;
      // The COMOVING reduced flux of the iterate, which is what the HLL part of the
      // ap_hll flux lags (the HLL acts on F0 only; the enthalpy flux A is added back
      // upwinded, exactly as in the explicit advective split).  The LAB f above still
      // drives the closure chi, as it does in the explicit scheme.
      //
      // It is NOT F0_cell/(c E) with F0_cell the arithmetic mean of the two faces.  That
      // is design risk R4 and it is fatal here: in free streaming the upwind face flux is
      // c E_{i-1}, so the cell mean is c (E_{i-1}+E_i)/2 and the derived f is
      // (1 + E_{i-1}/E_i)/2, i.e. 0.5 rather than 1 on the steep side of a pulse.  The
      // wave speeds then reopen to +-c/sqrt(3), the HLL flux turns CENTRED, and the
      // I6 pulse is flattened to its box mean in one crossing (amplitude ratio 0.0014).
      // Each FACE flux is therefore normalised by the E of the cell it comes FROM, which
      // is exactly 1 for an upwind free-streaming face, and the cell value is the mean of
      // the two face ratios.  On a cold start (f0x1 is zero-initialised and the problem
      // generator's state lives in u0) the cell flux is used instead.
      Real r0;
      Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
      if (rdgx) {
        fl -= clch*ifw_(m,M1_IFW_DG,k,j,i);
        fr -= clch*ifw_(m,M1_IFW_DG,k,j,i+1);
      }
      if (fabs(fl) + fabs(fr) > 0.0) {
        int iml = (i > is) ? (i-1)
                  : (cyclic ? ie : ((pos_.d_view(m) > 0) ? (is-1) : is));
        int ipr = (i < ie) ? (i+1)
                  : (cyclic ? is : ((pos_.d_view(m) < nblkx1-1) ? (ie+1) : ie));
        Real eul = fmax((fl > 0.0) ? iw_(m,M1_IW_EP,k,j,iml) : e, efl);
        Real eur = fmax((fr > 0.0) ? e : iw_(m,M1_IW_EP,k,j,ipr), efl);
        r0 = 0.5*(fl/(cl*eul) + fr/(cl*eur));
      } else {
        r0 = (f1 - iw_(m,M1_IW_ADV,k,j,i)*e)/(cl*e);
      }
      if (r0 > 1.0) {r0 = 1.0;}
      if (r0 < -1.0) {r0 = -1.0;}
      iw_(m,M1_IW_RF0,k,j,i) = r0;
