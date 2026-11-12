//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_cs_shell.inc
//! \brief d6fast4-1112: the cs shell kernel body of RadiationM1::VetGdCsSweep (moved
//! verbatim from rad_m1_vetgd_cs.cpp), included there once as shell_ and, on CUDA, a
//! second time as shellf_ with M1_VGDCS_SHELL_FAST defined: upf_ = 0 and two = false
//! are then compile-time constants inside, so nvcc drops the upwind-profile and
//! fused-twin branches (the same arithmetic on every path taken).
//! Needs M1_VGDCS_SHELL_NAME defined.
      auto M1_VGDCS_SHELL_NAME = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                               const int d, const bool cnt) {
#if defined(M1_VGDCS_SHELL_FAST)
        constexpr int upf_ = 0;
        constexpr bool two = false;
#endif
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        const Real mr = mr_(m,k,j,d);
        if (inw == (mr >= 0.0)) {return;}
        const Real rx = cr_(m,k,j,0), ry = cr_(m,k,j,1), rz = cr_(m,k,j,2);
        const Real r = cx1v(m,i);
        const Real zr = r*fabs(mr);
        const Real p2 = fmax(r*r - zr*zr, 0.0);
        // cstwin-1031: the geometry (branch, upwind shell, chord, stencil) of the cell
        // and direction, then the value of field f (main; the fused twin's from vgd_cst
        // / vgd_itw) by the same code as the unfused sweep
        // br: 0 the top shell's inflow, 1 the radial ray at the outward turn, 2 upwind
        const int br = (inw && l == n1 - 1) ? 0 : ((!inw && zr < 1.0e-12*r) ? 1 : 2);
        int typ = 1;
        int iu = i + 1;
        Real ds = 0.0, dso = 0.0;
        Real muu = 0.0;   // vetgdup-1110: n . r_hat at the upwind point
        VgdCsSt st;
        if (br == 0) {
          const Real rtop = cx1f(m,ie+1);
          ds = sqrt(fmax(rtop*rtop - p2, 0.0)) - zr;
        } else if (br == 2) {
          // the upwind shell, the chord and its own-cell part (rad_m1_vgdray.hpp)
          VgdRayChord(inw, l, lcut, m, i, zr, p2, cx1v, cx1f, spl, typ, iu, ds, dso);
          // the upwind point (Cartesian), its direction and the cs stencil
          Real xu[3] = {r*rx - ds*nx, r*ry - ds*ny, r*rz - ds*nz};
          const Real ru = sqrt(xu[0]*xu[0] + xu[1]*xu[1] + xu[2]*xu[2]);
          xu[0] /= ru;
          xu[1] /= ru;
          xu[2] /= ru;
          muu = nx*xu[0] + ny*xu[1] + nz*xu[2];
          const bool sok = VgdCsStencil(xu, bl_(m,0), bl_(m,1), bl_(m,2), wb, c2b, c3b,
                                        nc, h, cm_, nb_, st);
          // (counted per field swept, as the unfused twin's own sweep)
          if (!sok && cnt) {Kokkos::atomic_add(&cnt_(0), two ? 2.0 : 1.0);}
          if (st.n == 0) {   // nothing in the band (counted): the own column
            st.n = 1;
            st.k[0] = k;
            st.j[0] = j;
            st.w[0] = 1.0;
          }
        }
        bool ok[4];
        bool all = true;
        if (br == 2 && typ != 2) {
          const bool want_in = inw || (typ == 3);
          for (int e = 0; e < st.n; ++e) {
            const Real mq = mr_(m,st.k[e],st.j[e],d);
            ok[e] = want_in ? (mq < 0.0) : (mq >= 0.0);
            all = all && ok[e];
          }
        }
        for (int f = 0; f < (two ? 2 : 1); ++f) {
        const DvceArray5D<Real> &cv = (f == 0) ? cs_ : ct_;   // ln chi, ln S of field f
        const VgdIView &va = (f == 0) ? vi_ : vt_;            // its intensities
        const Real ch0 = exp(cv(m,0,k,j,i));
        const Real s0 = exp(cv(m,1,k,j,i));
        Real iv;
        // share3-1031: the stencil and the field reads of the shared upwind core
        const VgdCsA acc{st};
        auto fC = [&](const int c, const int kk, const int jj, const int ii) -> Real {
          return cv(m,c,kk,jj,ii);
        };
        if (br == 0) {
          iv = VgdTopInflow(ch0, ds, s0, exp(cv(m,1,k,j,i-1)), tbh, itop);
        } else if (br == 1) {
          iv = s0;
        } else {
          if (typ == 2) {
            iv = VgdInnerFace(acc, fC, i, emi, iemit, bfc, cx1f(m,i), cx1v(m,i),
                              cx1v(m,i+1), p2, ch0, ds, s0, izf);
          } else {
            // the intensity of a stencil cell: its owner's interior when on this rank
            // (stream order: shell iu is complete; the other branch of shell i, which
            // this launch writes, is masked), else the band copy of the remote value
            auto rdi = [&](const int kk, const int jj) -> Real {
              if (ccl && cnt) {
                const int x = rx_(m,kk,jj);
                if (x >= 0) {
                  const int64_t t = static_cast<int64_t>(x)*per_ + d*n1 + (iu - is);
                  if (((msk_(t >> 4) >> (2*(t & 15))) & 3u) == 0) {
                    Kokkos::atomic_add(&miss_(0), 1);
                  }
                }
              }
              return va(ow_(m,kk,jj,0),d,ow_(m,kk,jj,1),ow_(m,kk,jj,2),iu);
            };
            // vetgdup-1110 vet_gd_upwind (rad_m1_vgdray.hpp VgdUpCorr): the stencil
            // cells' intensities shifted along the shell's mean angular profile (main
            // field only; off: no shift)
            Real cq_[4] = {0.0, 0.0, 0.0, 0.0};
            bool upo = false;
            if (upf_ == 1 || upf_ == 2) {   // main and fused twin, each its own profile
              auto mqe = [&](const int e) -> Real {return mr_(m,st.k[e],st.j[e],d);};
              auto iqe = [&](const int e) -> Real {return rdi(st.k[e], st.j[e]);};
              upo = VgdUpCorr((f == 0) ? upp_ : upp2_, iu, upk_, inw || (typ == 3), upf_,
                              upt_, muu, acc,
                              mqe, iqe, ok, all, cq_);
            }
            auto rdc = [&](const int kk, const int jj) -> Real {
              if (!upo) {return rdi(kk,jj);}
              int e = 0;
              while (e < st.n - 1 && !(st.k[e] == kk && st.j[e] == jj)) {++e;}
              return rdi(kk,jj) + cq_[e];
            };
            const Real ivu = VgdUpwindI(acc, rdc, ok, all);
            auto lin = [&](const int c, const int ii) {
              Real s = 0.0;
              for (int e = 0; e < st.n; ++e) {s += st.w[e]*cv(m,c,st.k[e],st.j[e],ii);}
              return s;
            };
            const Real chu = exp(lin(0, iu));
            const Real su = exp(lin(1, iu));
            const Real rfi = cx1f(m,i);
            bool edge = false;
            if (vim == 3) {
              Real lmn = 1.0e300, lmx = -1.0e300;
              for (int e = 0; e < st.n; ++e) {
                if (!(st.w[e] > 0.0)) {continue;}
                const Real lc = cv(m,0,st.k[e],st.j[e],iu);
                lmn = fmin(lmn, lc);
                lmx = fmax(lmx, lc);
              }
              edge = (lmx - lmn) > lcr;
            }
            auto edgef = [&]() {return edge;};
            // the closing segment (vet_seg = jump below the inner face, vet_interp =
            // corner, vet_seg = split, the plain segment: rad_m1_vgdray.hpp)
            iv = VgdSegFinish(acc, rdc, fC, edgef, ok, iu, k, j, i, ivu, chu, su, ch0, s0,
                              ds, dso, spl, spe, typ, l, lcut, p2, rfi, zr, vim);
          }
        }
        va(m,d,k,j,i) = iv;
        }   // f
      };
