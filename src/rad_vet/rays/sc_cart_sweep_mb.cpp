//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_sweep_mb.cpp
//! \brief Cartesian short characteristics on several MeshBlocks: the same sweep, exact.
//!
//!  called from: rays/sc_cart_sweep.cpp
//!  calls into: rad_vet_imex.cpp, rays/sc_cart_setup.cpp

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_cart_exchange.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScCartSweepMB
//! \brief the short-characteristics sweep over several MeshBlocks, exact.
//! Launch L does global layer L for the up rays and nx1g-1-L for the down rays in
//! whichever block holds them.  Extinction and source go into banded planes with their
//! band and x1 face layers exchanged first; after each launch the new intensity plane's
//! band is exchanged, and at an x1 block face the upwind plane moves to the next block
//! (VetX1Move).  Bit-identical to the single-block sweep on one MeshBlock.

void RadiationVet::ScCartSweepMB() {
  // vet_source_noesrc: source and boundary intensity read the physical E^n = u0
  const bool noes_ = VsrcNoes();
  auto esv_ = noes_ ? u0 : iw;
  const int esc_ = noes_ ? VET_E : VET_IW_EN;
  ScCartMBState &st = *vet_mbs;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2;
  const int nx3 = trans_x3 ? indcs.nx3 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const bool hrs_ = impl_beam_hr;   // the half-range sums (cell kernel)
  const bool nls_ = thsw_c1;   // sum w |I - S| (VETRAY_NL)
  const int nray = vet_nray;
  const int nsb = st.nsb, nrl = st.nrl, nhl = st.nhl;
  const int nx1g = st.nx1g, nbx1 = st.nbx1;
  // w2, w3 here = the band (plane offset of the active cells)
  const int w2 = st.b2, w3 = st.b3;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto cs_ = st.csw;
  auto ip_ = st.ipl;
  auto lx_ = st.lx1;
  auto ang_ = vet_ang;
  auto cnt_ = vray_cnt;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real efl = e_floor;
  const Real bblo = VetBathLo();   // vet_bc_inner
  const Real bbhi = (vet_bc_bath && ibc_x1max == VET_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool bfc = vet_bc_face;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  const bool plm_ = true;

  // (a) banded extinction/source planes, layer i - is + 2 (moments zeroed by the caller)
  par_for("vet_vet_mb_csw", DevExeSpace(), 0, nmb-1, is, ie, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int i, const int k, const int j) {
    cs_(m,i-is+2,0,k-ks+w3,j-js+w2) = vc_(m,VETRAY_CHX,k,j,i);
    cs_(m,i-is+2,1,k-ks+w3,j-js+w2) = vc_(m,VETRAY_SRC,k,j,i);
  });
  // the neighbours' two face layers, nearer first (with nx1 = 1 the farther one is the
  // neighbour's own face layer from its neighbour)
  VetX1Move(st, st.csw, 0, -2, 1, nx1 + 1, 0, 2, nsb);
  VetX1Move(st, st.csw, 0, -2, 0, nx1, 0, 2, nsb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 2, 2, 0, 2, nsb);
  VetX1Move(st, st.csw, 1, -2, nx1 + 3, 3, 0, 2, nsb);
  VetBandExchange(st, st.csw, 0, nx1 + 4, 0, 2, nsb, nx2, nx3);

  const int nring = st.nring;
  // vet_x1_periodic: npass sweeps, each pass entering with the exit plane of the one
  // before; moments from the last pass.  closure_consistency: face sums of a neighbour's
  // cell go to spl (VetX1Move after)
  const int npass = vet_x1per ? vet_x1npass : 1, ngl = npass*nx1g;
  const bool x1src = vet_x1per && vet_x1src;
  const bool kps_ = true, f2s_ = cln_consist, x1per = vet_x1per;
  const int bmn_ = ibm_nr;
  const int bmr0 = ibm_ray[0], bmr1 = ibm_ray[1], bmr2 = ibm_ray[2], bmr3 = ibm_ray[3];
  auto bmi_ = ibm_i;
  auto sp_ = st.spl;
  if (f2s_) Kokkos::deep_copy(st.spl, 0.0);
  for (int gl = 0; gl < ngl; ++gl) {
    const int l = gl%nx1g;
    const bool wrap = (gl >= nx1g), accum = (gl >= ngl - nx1g);
    // ring slot of this launch's plane and of the upwind one
    const int pw = gl%nring, pr = (gl + nring - 1)%nring;
    if (gl > 0 && (l % nx1) == 0) {
      // entering the next block of the stack: up rays from below, down from above
      const int q = l/nx1;
      VetX1Move(st, st.ipl, 0, q, pr, pr, 0, nhl, nsb);
      VetX1Move(st, st.ipl, 1, nbx1 - 1 - q, pr, pr, nhl, nrl - nhl, nsb);
    }
    par_for("vet_vet_mb_sweep", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const int lb = lx_.d_view(m);
      const int lu = l - lb*nx1;
      const int ld = (nx1g - 1 - l) - lb*nx1;
      const bool ua = (lu >= 0 && lu < nx1), da = (ld >= 0 && ld < nx1);
      if (!ua && !da) return;
      const Real dx1 = mbsize.d_view(m).dx1;
      const Real dx2 = mbsize.d_view(m).dx2;
      const Real dx3 = mbsize.d_view(m).dx3;
      const int iu = is + lu, id = is + ld;
      const int jc = j - js + w2, kc = k - ks + w3;
      Real acu[10], acd[10];
      for (int n = 0; n < 10; ++n) {
        acu[n] = 0.0;
        acd[n] = 0.0;
      }
      Real hpu[3] = {0.0, 0.0, 0.0}, hpd[3] = {0.0, 0.0, 0.0};
      Real nlu = 0.0, nld = 0.0;
      Real kpu[VETRAY_NKP], kpd[VETRAY_NKP];
      for (int n = 0; n < VETRAY_NKP; ++n) {kpu[n] = 0.0; kpd[n] = 0.0;}
      // face fluxes into the cells iu, id (F1, FX) and their upwind cells (F2, FX, F3)
      Real f2u[VETRAY_NFS] = {0.0}, f2d[VETRAY_NFS] = {0.0};
      Real fnu[VETRAY_NFS] = {0.0}, fnd[VETRAY_NFS] = {0.0};
      for (int r = 0; r < nray; ++r) {
        const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
        const bool up = (m1 > 0.0);
        if (up ? !ua : !da) continue;
        const int i = up ? iu : id;
        const Real am1 = fabs(m1);
        const Real c0 = vc_(m,VETRAY_CHX,k,j,i);
        const Real s0 = vc_(m,VETRAY_SRC,k,j,i);
        Real iv;
        auto add4 = [&](Real *fx, const Real g) {
          fx[0] += g;
          fx[1] += g*m1;
          fx[2] += g*m2;
          fx[3] += g*m3;
        };
        Real ifc = 0.0;   // closure_consistency: the intensity at the face crossing
        if (l == 0 && !wrap) {
          // physical x1 face: as the single-block sweep
          const int iin = up ? (i + 1) : (i - 1);
          const Real s1 = vc_(m,VETRAY_SRC,k,j,iin);
          const Real su = fmax(1.5*s0 - 0.5*s1, 0.0);
          Real ib = up ? 0.0 : bbhi;
          if (x1src) {
            ib = su;
          } else if (up) {
            if (bblo >= 0.0) {
              ib = bblo;
            } else {
              Real e = esv_(m,esc_,k,j,i);
              if (bfc) e += 0.5*(s0 - s1);   // vet_bc_diff = face
              e = fmax(e, efl);
              Real fo = iw_(m,VET_IW_F1,k,j,i)*m1 + iw_(m,VET_IW_F2,k,j,i)*m2
                        + iw_(m,VET_IW_F3,k,j,i)*m3;
              ib = fmax(e + 3.0*fo/cl, 0.0);
            }
            const int q = (r == bmr0) ? 0 : ((r == bmr1) ? 1
                          : ((r == bmr2) ? 2 : ((r == bmr3) ? 3 : -1)));
            if (q >= 0 && q < bmn_) {ib += bmi_(m,q,k-ks,j-js);}
          }
          if (f2s_) add4((up ? f2u : f2d) + 4, wr*m1*ib);
          const Real dtau = 0.5*c0*dx1/am1;
          Real ex, w0, wu;
          VetSegW(dtau, ex, w0, wu);
          iv = ib*ex + wu*su + w0*s0;
        } else {
          // upwind plane: banded layer i -+ 1 (1 and nx1+2 = the neighbour's face)
          const int li = up ? (i - is + 1) : (i - is + 3);
          const Real sh2 = -dx1*m2/(am1*dx2);
          const Real fl2 = floor(sh2);
          const Real a2 = sh2 - fl2;
          const int o2 = static_cast<int>(fl2);
          const int ja = jc + o2;
          const int jb = ja + 1;
          int ka = kc, kb = kc, o3 = 0;
          Real a3 = 0.0;
          const Real sh3 = thrd ? -dx1*m3/(am1*dx3) : 0.0;
          if (thrd) {
            const Real fl3 = floor(sh3);
            a3 = sh3 - fl3;
            o3 = static_cast<int>(fl3);
            ka = kc + o3;
            kb = ka + 1;
          }
          Real dv = 0.0, fl[2] = {0.0, 0.0}, dq[4] = {0.0, 0.0, 0.0, 0.0};
          if (plm_) {
            dv = VetPlmFoot([&](const int dk, const int dj) {
              return ip_(m,pr,r,ka + dk,ja + dj);
            }, a2, a3, thrd, fl, dq);
          }
          if (f2s_) {
            // upwind plane's intensity carried across the low x2 (x3) face by the foot
            // shift (ScCartSweep)
            auto g2 = [&](const int d) {
              if (!thrd) return ip_(m,pr,r,kc,jc + d);
              return 0.5*(ip_(m,pr,r,kc,jc + d) + (1.0 - a3)*ip_(m,pr,r,ka,jc + d)
                          + a3*ip_(m,pr,r,kb,jc + d));
            };
            auto g3 = [&](const int d) {
              return 0.5*(ip_(m,pr,r,kc + d,jc) + (1.0 - a2)*ip_(m,pr,r,kc + d,ja)
                          + a2*ip_(m,pr,r,kc + d,jb));
            };
            auto phi = [&](const Real sh, const auto &g) {
              const Real sa = fabs(sh), fa = sa - floor(sa);
              const int na = static_cast<int>(floor(sa)), sg = (sh < 0.0) ? 1 : -1;
              Real p = fa*g((sg > 0) ? -1 - na : na);
              for (int q = 0; q < na; ++q) p += g((sg > 0) ? -1 - q : q);
              return sg*p;
            };
            add4(up ? fnu : fnd, wr*am1*dx2/dx1*(phi(sh2, g2) + fl[0]));
            if (thrd) {
              add4((up ? fnu : fnd) + 12, wr*am1*dx3/dx1*(phi(sh3, g3) + fl[1]));
            }
          }
          const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
          const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
          const Real iup_v = waa*ip_(m,pr,r,ka,ja) + wab*ip_(m,pr,r,ka,jb)
                             + wba*ip_(m,pr,r,kb,ja) + wbb*ip_(m,pr,r,kb,jb) + dv;
          const Real cu = waa*cs_(m,li,0,ka,ja) + wab*cs_(m,li,0,ka,jb)
                          + wba*cs_(m,li,0,kb,ja) + wbb*cs_(m,li,0,kb,jb);
          const Real su = waa*cs_(m,li,1,ka,ja) + wab*cs_(m,li,1,ka,jb)
                          + wba*cs_(m,li,1,kb,ja) + wbb*cs_(m,li,1,kb,jb);
          // jump: the corners' S one x1 layer further upwind (0 beyond a physical x1
          // face: q = 1)
          const bool juu = spl && spe == 2;
          const int lu2 = up ? (li - 1) : (li + 1);
          auto s2 = [&](const int kq, const int jq) {
            return juu ? cs_(m,lu2,1,kq,jq) : 0.0;
          };
          const Real squ[4] = {s2(ka,ja), s2(ka,jb), s2(kb,ja), s2(kb,jb)};
          const Real *suq = juu ? squ : nullptr;
          // and the own cell's S one x1 layer downstream
          const Real sdd = juu ? cs_(m,up ? (li + 2) : (li - 2),1,kc,jc) : -1.0;
          if (vim != 0) {   // vet_interp: log | corner
            const Real wq[4] = {waa, wab, wba, wbb};
            const Real iq[4] = {ip_(m,pr,r,ka,ja) + dq[0], ip_(m,pr,r,ka,jb) + dq[1],
                                ip_(m,pr,r,kb,ja) + dq[2], ip_(m,pr,r,kb,jb) + dq[3]};
            const Real cq[4] = {cs_(m,li,0,ka,ja), cs_(m,li,0,ka,jb), cs_(m,li,0,kb,ja),
                                cs_(m,li,0,kb,jb)};
            const Real sq[4] = {cs_(m,li,1,ka,ja), cs_(m,li,1,ka,jb), cs_(m,li,1,kb,ja),
                                cs_(m,li,1,kb,jb)};
            iv = VetUpCorner(vim, wq, iq, cq, sq, s0, c0, dx1/am1, spl, spe,
                             f2s_ ? &ifc : nullptr, suq, sdd);
          } else if (spl) {
            // split: the two halves at the x1 face
            iv = VetSegSplit(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1, spe,
                             f2s_ ? &ifc : nullptr,
                             juu ? (waa*squ[0] + wab*squ[1] + wba*squ[2] + wbb*squ[3])
                                    : -1.0, sdd);
          } else {
            const Real dtau = 0.5*(cu + c0)*dx1/am1;
            Real ex, w0, wu;
            VetSegW(dtau, ex, w0, wu);
            iv = iup_v*ex + wu*su + w0*s0;
            // face crossing at the face (0.5 chi_u hl from u), not at half the tau
            if (f2s_) ifc = VetSegFaceI(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1,
                                        0);
          }
          if (f2s_) add4(up ? f2u + 8 : fnd + 8, wr*m1*fmax(ifc, 0.0));
        }
        if (iv < 0.0) Kokkos::atomic_add(&cnt_(0), 1.0);   // (none by construction)
        iv = fmax(iv, 0.0);
        ip_(m,pw,r,kc,jc) = iv;
        if (f2s_ && !x1per && l == nx1g - 1) {
          // ray leaving through the x1 boundary face: the half cell at chi_0, S_0
          const Real ex = exp(-0.5*c0*dx1/am1);
          add4((up ? f2u : f2d) + 4, wr*m1*(iv*ex + (1.0 - ex)*s0));
        }
        Real *ac = up ? acu : acd;
        const Real a = wr*iv;
        ac[0] += a;
        ac[1] += a*m1*m1;
        ac[2] += a*m2*m2;
        ac[3] += a*m3*m3;
        ac[4] += a*m1*m2;
        ac[5] += a*m1*m3;
        ac[6] += a*m2*m3;
        ac[7] += a*m1;
        ac[8] += a*m2;
        ac[9] += a*m3;
        if (hrs_) {   // the half-range sums, as ScCartBuild
          Real *hp = up ? hpu : hpd;
          hp[0] += a*fmax(m1, 0.0);
          hp[1] += a*fmax(m2, 0.0);
          hp[2] += a*fmax(m3, 0.0);
        }
        if (kps_) {
          Real *kp = up ? kpu : kpd;
          const Real mu[3] = {m1, m2, m3};
          for (int d = 0; d < 3; ++d) {
            if (!(mu[d] > 0.0)) continue;
            for (int b = 0; b < 3; ++b) kp[3*d+b] += a*mu[b]*mu[d];
          }
        }
        if (nls_) {
          if (up) {nlu += wr*fabs(iv - s0);} else {nld += wr*fabs(iv - s0);}
        }
      }
      if (!accum) return;
      for (int n = 0; n < 10; ++n) {
        if (ua) {
          vc_(m,VETRAY_J+n,k,j,iu) += acu[n];
        }
        if (da) {
          vc_(m,VETRAY_J+n,k,j,id) += acd[n];
        }
      }
      if (hrs_) {
        for (int n = 0; n < 3; ++n) {
          if (ua) {vc_(m,VETRAY_HP1+n,k,j,iu) += hpu[n];}
          if (da) {vc_(m,VETRAY_HP1+n,k,j,id) += hpd[n];}
        }
      }
      if (nls_) {
        if (ua) {vc_(m,VETRAY_NL,k,j,iu) += nlu;}
        if (da) {vc_(m,VETRAY_NL,k,j,id) += nld;}
      }
      if (kps_) {
        for (int n = 0; n < VETRAY_NKP; ++n) {
          if (ua) {vc_(m,VETRAY_KP+n,k,j,iu) += kpu[n];}
          if (da) {vc_(m,VETRAY_KP+n,k,j,id) += kpd[n];}
        }
      }
      if (f2s_) {
        for (int n = 0; n < VETRAY_NFS; ++n) {
          if (ua) {
            vc_(m,VETRAY_F2+n,k,j,iu) += f2u[n];
            if (lu > 0) {
              vc_(m,VETRAY_F2+n,k,j,iu-1) += fnu[n];
            } else {
              sp_(m,0,n,kc,jc) += fnu[n];
            }
          }
          if (da) {
            vc_(m,VETRAY_F2+n,k,j,id) += f2d[n];
            if (ld < nx1 - 1) {
              vc_(m,VETRAY_F2+n,k,j,id+1) += fnd[n];
            } else {
              sp_(m,1,n,kc,jc) += fnd[n];
            }
          }
        }
      }
    });
    if (gl < ngl - 1) {
      for (int p = 1; p < st.npb; ++p) {
        VetPlanePost(st, st.ipl, pw, nrl, nsb, nx2, nx3, p);
        VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
      }
      VetPlanePost(st, st.ipl, pw, nrl, nsb, nx2, nx3, st.npb);
      VetPlaneWait(st, st.ipl, nrl, nsb, nx2, nx3);
    }
  }
  if (f2s_) {
    // face sums of the x1 neighbours' boundary cells
    VetX1Move(st, st.spl, 1, -2, 2, 0, 0, VETRAY_NFS, nsb);
    VetX1Move(st, st.spl, 0, -2, 3, 1, 0, VETRAY_NFS, nsb);
    par_for("vet_vet_mb_spill", DevExeSpace(), 0, nmb-1, 0, VETRAY_NFS-1, ks, ke, js,
            je,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j) {
      vc_(m,VETRAY_F2+n,k,j,ie) += sp_(m,2,n,k-ks+w3,j-js+w2);
      vc_(m,VETRAY_F2+n,k,j,is) += sp_(m,3,n,k-ks+w3,j-js+w2);
    });
  }
}

}  // namespace radvet
