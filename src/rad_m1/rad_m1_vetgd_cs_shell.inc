//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_vetgd_cs_shell.inc
//! \brief the cs shell kernel body of RadiationM1::VetGdCsSweep, included
//! there as shell_ (the current set) and shellx_ (the next set of a fused cross-fade).
//! Needs M1_VGDCS_SHELL_NAME defined.
      auto M1_VGDCS_SHELL_NAME = KOKKOS_LAMBDA(const int m, const int k, const int j,
                                               const int d, const bool cnt) {
        const Real nx = dir_(d,0), ny = dir_(d,1), nz = dir_(d,2);
        const Real mr = mr_(m,k,j,d);
        if (inw == (mr >= 0.0)) {return;}
        const Real rx = cr_(m,k,j,0), ry = cr_(m,k,j,1), rz = cr_(m,k,j,2);
        const Real r = cx1v(m,i);
        const Real zr = r*fabs(mr);
        const Real p2 = fmax(r*r - zr*zr, 0.0);
        // the geometry (branch, upwind shell, chord, stencil) of the cell and direction,
        // then its value
        // br: 0 the top shell's inflow, 1 the radial ray at the outward turn, 2 upwind
        const int br = (inw && l == n1 - 1) ? 0 : ((!inw && zr < 1.0e-12*r) ? 1 : 2);
        int typ = 1;
        int iu = i + 1;
        Real ds = 0.0, dso = 0.0;
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
          const bool sok = VgdCsStencil(xu, bl_(m,0), bl_(m,1), bl_(m,2), wb, c2b, c3b,
                                        nc, h, cm_, nb_, st);
          if (!sok && cnt) {Kokkos::atomic_add(&cnt_(0), 1.0);}
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
        const DvceArray5D<Real> &cv = cs_;   // ln chi, ln S
        const VgdIView &va = vi_;            // the intensities
        const Real ch0 = exp(cv(m,0,k,j,i));
        const Real s0 = exp(cv(m,1,k,j,i));
        Real iv;
        // the stencil and the field reads of the shared upwind core
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
            const Real ivu = VgdUpwindI(acc, rdi, ok, all);
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
            // vet_seg = jump's plateau test: the shell one further upwind (in the block)
            const int iu2 = 2*iu - i;
            const int iuu = (typ == 1 && iu2 >= is && iu2 <= ie) ? iu2 : iu;
            const int id2 = 2*i - iu;
            const int idd = (typ == 1 && id2 >= is && id2 <= ie) ? id2 : i;
            iv = VgdSegFinish(acc, rdi, fC, edgef, ok, iu, k, j, i, ivu, chu, su, ch0, s0,
                              ds, dso, spl, spe, typ, l, lcut, p2, rfi, zr, vim, iuu,
                              idd);
          }
        }
        va(m,d,k,j,i) = iv;
      };
