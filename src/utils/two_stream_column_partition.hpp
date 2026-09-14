#ifndef UTILS_TWO_STREAM_COLUMN_PARTITION_HPP_
#define UTILS_TWO_STREAM_COLUMN_PARTITION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_column_partition.hpp
//! \brief THE SAME EXACT IMPLICIT COLUMN SOLVE as two_stream_column_implicit.hpp, with
//! the column's serial recurrences PARTITIONED over a Kokkos thread team
//! (<problem>/rt_impl_solver = pcr).  Same linear system, same Newton loop, same clamps;
//! only the ORDER of the eliminations differs, so the two agree to round-off.
//!
//! WHY.  Mode 3's cost is not arithmetic, it is STARVATION.  One thread per column is
//! ~1.2e4 threads per GPU for a production FeCZ box, each running a 336-step serial
//! block recurrence: a 300-CU MI300A gets ~40 threads per CU and every one of them is
//! latency-bound on a dependent chain.  Nothing about the system forces that chain --
//! a block-tridiagonal solve parallelises -- so the fix is to give each column a TEAM.
//!
//! THE ALGORITHM, in five lines (P = team size = number of segments):
//!   1. cut [icut, ie] into P contiguous segments, one per thread;
//!   2. each thread Thomas-eliminates its OWN segment with the incoming unknown carried
//!      symbolically, i.e. x_i = d_i - G_i x_{i+1} - H_i x_L (H is the extra 5x3 spike);
//!   3. a backward pass inside the segment expresses its FIRST cell as
//!      x_{i0} = p - Q y_s - R y_{s-1} in the two segment-boundary unknowns;
//!   4. the boundary cells then satisfy a block-tridiagonal system of P 5x5 rows -- the
//!      REDUCED system -- which one thread solves by the same block Thomas;
//!   5. every thread back-substitutes its own segment in parallel.
//! The two formal-solution sweeps of step 4a are partitioned the same way: each thread
//! runs its segment with a ZERO incoming intensity and also accumulates the homogeneous
//! factor prod t^2, one thread composes the P segment entry values, and each thread then
//! adds h*D_entry to its cells.  Exact for a linear recurrence, and it costs one extra
//! multiply-add per cell instead of a second pass of the (expensive) source evaluation.
//!
//! WHAT IS NOT BITWISE.  The clamp of step 4c (|db| <= 3b, >= -0.75b) propagates INTO the
//! next cell's back-substitution.  Inside a segment and along the reduced chain that is
//! reproduced exactly, but a clamp firing in the interior of segment s cannot reach
//! segment s-1's boundary value the way it does in the serial sweep.  So the two paths
//! agree to round-off whenever nclamp = 0 (the converged regime; nclamp is reported every
//! cycle) and differ at O(clamp) when it is not.  The Newton FIXED POINT is the same
//! either way: at convergence db -> 0 and no clamp fires.
//!
//! NOT IMPLEMENTED HERE: the rt_outer_verbose per-cell assembly dump (rt_col3_it /
//! rt_col3_flux / rt_col3_cell).  It is a cycle-0 debugging aid; run rt_impl_solver =
//! thomas to get it.

#include <math.h>

#include "athena.hpp"
#include "utils/two_stream_column_implicit.hpp"

namespace two_stream_rt {

// the per-cell workspace of the PARTITIONED path: the 35 slots of the serial solver,
// plus the spike H (5x3) and the two homogeneous sweep factors (one per angle)
#define RTCOL3_NWP 54
// the per-SEGMENT workspace, see the slot list in RTCol3TeamSolve
#define RTCOL3_NRD 96

//----------------------------------------------------------------------------------------
//! \fn void RTCol3TeamSolve
//! \brief one column on one team.  Every thread walks the same code; the work is split
//! by SEGMENT (thread s owns cells [i0(s), i1(s)]), which keeps each cell's data on one
//! thread through the whole Newton pass.

KOKKOS_INLINE_FUNCTION
void RTCol3TeamSolve(const RTCol3 &c, const TeamMember_t &tm,
                     const int m, const int k, const int j) {
  const int ie = c.ie;
  const int ic = c.icut(m,k,j);
  if (ic > ie) return;                        // nothing radiative in this column
  const int nc = ie - ic + 1;
  const int nsg = (c.nseg < nc) ? c.nseg : nc;

  // per-cell slots: the serial solver's layout, then H and the sweep factors
  const int G0 = 0, DP = 15, EE = 20, CI = 22, CO = 24, BB = 26, DD = 27, UU = 29,
            EX = 31, SA = 32, ES = 34, HH = 35, HD = 50, HU = 52;
  // per-segment slots of c.rd(m, s*RTCOL3_NRD + slot, k, j)
  const int PP = 0, QQ = 5, RR = 30, GR = 45, DR = 70, DE = 75, UE = 77, YY = 79,
            YR = 84, BUD = 85, BSC = 86, RTM = 87, UBM = 88, RHS = 89, SRS = 90,
            NCL = 91;
  const int nq = c.nq;
  const Real sopi = c.sigma/M_PI;
  // the three column components of the reduced (x[0], x[1], x[4]) coupling
  const int rc3[3] = {0, 1, 4};

  // ---- 1. the frozen per-cell layer coefficients, and b^0 --------------------------
  Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
    const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
    for (int i=i0; i<=i1; ++i) {
      const Real h = c.Ht(m,k,j,i);
      for (int q=0; q<nq; ++q) {
        const Real x = h/c.mu[q];
        const Real e0 = -expm1(-x);
        c.wk(m,EE+q,i,k,j) = e0;
        c.wk(m,CI+q,i,k,j) = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
        c.wk(m,CO+q,i,k,j) = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
      }
      c.wk(m,BB,i,k,j) = c.Bb(m,0,i,k,j);
      c.wk(m,ES,i,k,j) = c.Ei(m,k,j,i);
    }
  });
  tm.team_barrier();
  Real emax = 0.0;
  if (c.norm == 1) {
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
    [&](const int i, Real &mx) {
      const Real e = c.wk(m,ES,i,k,j);
      if (e > mx) mx = e;
    }, Kokkos::Max<Real>(emax));
  }
  const Real eoff = c.norm_eps*emax;

  // ---- 2. the frozen boundary data, formed redundantly on every thread --------------
  Real dbdtau = 0.0;
  if (!c.cut_legacy && ic + 1 <= ie) {
    const Real dtc = c.Ht(m,k,j,ic) + c.Ht(m,k,j,ic+1);
    if (dtc > 0.0) dbdtau = (c.Bb(m,0,ic,k,j) - c.Bb(m,0,ic+1,k,j))/dtc;
  }
  // problem/rt_bottom_flux: the cut IS the bottom wall and carries the imposed internal
  // flux, so its gradient is set by that flux, not by the column's own two deepest cells.
  if (c.bot_flux > 0.0) dbdtau = 3.0*c.bot_flux/(4.0*M_PI);
  const Real cutc = dbdtau*c.Ht(m,k,j,ic);
  Real Dtop[2], Ucut[2];
  for (int q=0; q<2; ++q) {
    Dtop[q] = c.dtop(m,q,k,j);
    Ucut[q] = cutc + (c.int_at_cut ? c.Iint : 0.0) + c.mu[q]*dbdtau;
  }

  // ---- 3. the tau-blend handover, frozen from the entry-state sweep -----------------
  Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
    const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
    for (int i=i0; i<=i1; ++i) {
      Real ex = 0.0;
      if (c.taublend) {
        const Real dxi = c.Dx(m,k,j,i);
        const Real wb = 0.5*(c.wblend(m,k,j,i) + c.wblend(m,k,j,i+1));
        const Real ft = c.Fb(m,0,i+1,k,j), fb = c.Fb(m,0,i,k,j);
        if (c.direct) {
          ex = wb*c.Src(m,0,i,k,j)
             + (c.wblend(m,k,j,i+1)*ft - c.wblend(m,k,j,i)*fb)/dxi;
        } else {
          ex = -((1.0 - c.wblend(m,k,j,i+1))*ft - (1.0 - c.wblend(m,k,j,i))*fb)/dxi
             + (1.0 - wb)*(ft - fb)/dxi;
        }
      }
      c.wk(m,EX,i,k,j) = ex + c.Qb(m,0,i,k,j);
    }
  });

  // ---- 4. Newton ---------------------------------------------------------------------
  int nit = 0;
  Real dbmax = 0.0;
  Real rfin = 0.0;
  for (int it=0; it<c.maxit; ++it) {
    nit = it + 1;
    tm.team_barrier();
    // ---- 4a. the formal solution, each segment with a ZERO incoming intensity --------
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      Real L[2] = {0.0, 0.0}, hg[2] = {1.0, 1.0};
      for (int i=i1; i>=i0; --i) {
        Real sl, su, sfu, sfd;
        c.SourceVals(m, k, j, i, ic, cutc, sl, su, sfu, sfd);
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.wk(m,EE+q,i,k,j), t = 1.0 - E;
          const Real ci = c.wk(m,CI+q,i,k,j), co = c.wk(m,CO+q,i,k,j);
          const Real P = ci*sfu + co*sl;
          const Real Q = ci*su + co*sfd;
          acc += c.wf[q]*W*(E*(1.0 + t)*L[q] + E*P - (P + Q));
          L[q] = t*t*L[q] + t*P + Q;
          hg[q] = t*t*hg[q];
          c.wk(m,DD+q,i,k,j) = L[q];
          c.wk(m,HD+q,i,k,j) = hg[q];
        }
        c.wk(m,SA,i,k,j) = acc;
      }
      Real Lu[2] = {0.0, 0.0}, hu[2] = {1.0, 1.0};
      for (int i=i0; i<=i1; ++i) {
        Real sl, su, sfu, sfd;
        c.SourceVals(m, k, j, i, ic, cutc, sl, su, sfu, sfd);
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.wk(m,EE+q,i,k,j), t = 1.0 - E;
          const Real ci = c.wk(m,CI+q,i,k,j), co = c.wk(m,CO+q,i,k,j);
          const Real Pu = ci*sfd + co*su;
          const Real Qu = ci*sl + co*sfu;
          acc += c.wf[q]*W*(E*(1.0 + t)*Lu[q] + E*Pu - (Pu + Qu));
          Lu[q] = t*t*Lu[q] + t*Pu + Qu;
          hu[q] = t*t*hu[q];
          c.wk(m,UU+q,i,k,j) = Lu[q];
          c.wk(m,HU+q,i,k,j) = hu[q];
        }
        c.wk(m,SA,i,k,j) += acc;
      }
    });
    tm.team_barrier();
    // the P segment entry intensities: a serial scan over SEGMENTS, not over cells
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      Real de[2], ue[2];
      for (int q=0; q<nq; ++q) de[q] = Dtop[q];
      for (int s=nsg-1; s>=0; --s) {
        const int i0 = ic + (nc*s)/nsg;
        for (int q=0; q<nq; ++q) c.rd(m,s*RTCOL3_NRD+DE+q,k,j) = de[q];
        for (int q=0; q<nq; ++q) {
          de[q] = c.wk(m,DD+q,i0,k,j) + c.wk(m,HD+q,i0,k,j)*de[q];
        }
      }
      for (int q=0; q<nq; ++q) ue[q] = c.wk(m,BB,ic,k,j) + Ucut[q];
      for (int s=0; s<nsg; ++s) {
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        for (int q=0; q<nq; ++q) c.rd(m,s*RTCOL3_NRD+UE+q,k,j) = ue[q];
        for (int q=0; q<nq; ++q) {
          ue[q] = c.wk(m,UU+q,i1,k,j) + c.wk(m,HU+q,i1,k,j)*ue[q];
        }
      }
    });
    tm.team_barrier();
    // and the correction each cell owes to its segment's entry value
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      Real de[2], ue[2];
      for (int q=0; q<nq; ++q) {
        de[q] = c.rd(m,s*RTCOL3_NRD+DE+q,k,j);
        ue[q] = c.rd(m,s*RTCOL3_NRD+UE+q,k,j);
      }
      for (int i=i0; i<=i1; ++i) {
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.wk(m,EE+q,i,k,j), t = 1.0 - E;
          const Real hd = (i == i1) ? 1.0 : c.wk(m,HD+q,i+1,k,j);
          const Real hup = (i == i0) ? 1.0 : c.wk(m,HU+q,i-1,k,j);
          acc += c.wf[q]*W*E*(1.0 + t)*(hd*de[q] + hup*ue[q]);
          c.wk(m,DD+q,i,k,j) += c.wk(m,HD+q,i,k,j)*de[q];
          c.wk(m,UU+q,i,k,j) += c.wk(m,HU+q,i,k,j)*ue[q];
        }
        c.wk(m,SA,i,k,j) += acc;
      }
    });
    tm.team_barrier();
    // ---- 4a'. the handover, re-formed from THIS column's OWN flux ------------------
    if (c.ex_iter && c.taublend && c.direct) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
        Real f3lo = 0.0, f3hi = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real ulo = (i == ic) ? (c.wk(m,BB,ic,k,j) + Ucut[q])
                                     : c.wk(m,UU+q,i-1,k,j);
          const Real dhi = (i == ie) ? Dtop[q] : c.wk(m,DD+q,i+1,k,j);
          f3lo += c.wf[q]*(ulo - c.wk(m,DD+q,i,k,j));
          f3hi += c.wf[q]*(c.wk(m,UU+q,i,k,j) - dhi);
        }
        const Real wlo = c.wblend(m,k,j,i), whi = c.wblend(m,k,j,i+1);
        c.wk(m,EX,i,k,j) = 0.5*(wlo + whi)*c.wk(m,SA,i,k,j)
                         + (whi*f3hi - wlo*f3lo)/c.Dx(m,k,j,i) + c.Qb(m,0,i,k,j);
      });
      tm.team_barrier();
    }

    // ---- 4b. the PARTITIONED forward elimination ------------------------------------
    Real rmax = 0.0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nsg),
    [&](const int s, Real &rmx) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      Real Gp[5][3], Hp[5][3], dp[5];
      for (int i=i0; i<=i1; ++i) {
        Real A3[5][3], Bm[5][5], C3[5][3], rv[5], AH[5][3];
        Real rsc = 1.0;
        c.BuildRow(m, k, j, i, ic, cutc, it, A3, Bm, C3, rv, rsc);
        const Real den = rsc + eoff;
        const Real rr = (den > 0.0) ? fabs(rv[4])/den : 0.0;
        if (rr > rmx) rmx = rr;
        Real sg = 1.0;
        if (i > i0) {
          sg = -1.0;
          for (int r=0; r<5; ++r) {
            for (int cc=0; cc<3; ++cc) {
              Real sa = 0.0, sh = 0.0;
              for (int t=0; t<3; ++t) {
                sa += A3[r][t]*Gp[t+2][cc];
                sh += A3[r][t]*Hp[t+2][cc];
              }
              Bm[r][rc3[cc]] -= sa;
              AH[r][cc] = sh;
            }
            Real s2 = 0.0;
            for (int t=0; t<3; ++t) s2 += A3[r][t]*dp[t+2];
            rv[r] -= s2;
          }
        } else {
          for (int r=0; r<5; ++r) {
            for (int cc=0; cc<3; ++cc) AH[r][cc] = A3[r][cc];
          }
        }
        Real Bi[5][5];
        if (!RTCol3Inv5(Bm, Bi)) {
          Kokkos::atomic_add(&c.stat(6), 1.0);
          rmx = HUGE_VAL;
          return;
        }
        for (int r=0; r<5; ++r) {
          Real sd = 0.0;
          for (int cc=0; cc<5; ++cc) sd += Bi[r][cc]*rv[cc];
          dp[r] = sd;
          c.wk(m,DP+r,i,k,j) = sd;
          for (int cc=0; cc<3; ++cc) {
            Real g = 0.0, h = 0.0;
            for (int t=0; t<5; ++t) {
              g += Bi[r][t]*C3[t][cc];
              h += Bi[r][t]*AH[t][cc];
            }
            Gp[r][cc] = g;
            Hp[r][cc] = sg*h;
            c.wk(m,G0+3*r+cc,i,k,j) = g;
            c.wk(m,HH+3*r+cc,i,k,j) = sg*h;
          }
        }
      }
      // ---- 4b'. the segment's FIRST cell in the two boundary unknowns ---------------
      // x_i = p_i - Q_i y_s - R_i y_{s-1}, started at i1 from x_{i1} = y_s.  With one
      // segment there is no reduced system to feed and this is pure waste, so skip it:
      // nsg = 1 is then exactly the serial block Thomas (the host-backend path).
      if (nsg < 2) return;
      Real pq[5], Qm[5][5], Rm[5][3];
      for (int r=0; r<5; ++r) {
        pq[r] = 0.0;
        for (int cc=0; cc<5; ++cc) Qm[r][cc] = (r == cc) ? -1.0 : 0.0;
        for (int cc=0; cc<3; ++cc) Rm[r][cc] = 0.0;
      }
      for (int i=i1-1; i>=i0; --i) {
        Real pn[5], Qn[5][5], Rn[5][3];
        for (int r=0; r<5; ++r) {
          Real sp = c.wk(m,DP+r,i,k,j);
          for (int cc=0; cc<3; ++cc) sp -= c.wk(m,G0+3*r+cc,i,k,j)*pq[rc3[cc]];
          pn[r] = sp;
          for (int col=0; col<5; ++col) {
            Real sq = 0.0;
            for (int cc=0; cc<3; ++cc) sq += c.wk(m,G0+3*r+cc,i,k,j)*Qm[rc3[cc]][col];
            Qn[r][col] = -sq;
          }
          for (int col=0; col<3; ++col) {
            Real sr = c.wk(m,HH+3*r+col,i,k,j);
            for (int cc=0; cc<3; ++cc) sr -= c.wk(m,G0+3*r+cc,i,k,j)*Rm[rc3[cc]][col];
            Rn[r][col] = sr;
          }
        }
        for (int r=0; r<5; ++r) {
          pq[r] = pn[r];
          for (int col=0; col<5; ++col) Qm[r][col] = Qn[r][col];
          for (int col=0; col<3; ++col) Rm[r][col] = Rn[r][col];
        }
      }
      const int b0 = s*RTCOL3_NRD;
      for (int r=0; r<5; ++r) {
        c.rd(m,b0+PP+r,k,j) = pq[r];
        for (int col=0; col<5; ++col) c.rd(m,b0+QQ+5*r+col,k,j) = Qm[r][col];
        for (int col=0; col<3; ++col) c.rd(m,b0+RR+3*r+col,k,j) = Rm[r][col];
      }
    }, Kokkos::Max<Real>(rmax));
    rfin = rmax;
    if (!(rmax < HUGE_VAL)) return;     // a singular block: leave the column alone
    if (rmax < c.tol) break;

    // ---- 4c. the REDUCED system, and the clamped update -----------------------------
    tm.team_barrier();
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      for (int s=0; s<nsg; ++s) {
        const int b0 = s*RTCOL3_NRD;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real Bh[5][5], Ch[5][5], rh[5];
        for (int r=0; r<5; ++r) {
          rh[r] = c.wk(m,DP+r,i1,k,j);
          for (int col=0; col<5; ++col) {
            Bh[r][col] = (r == col) ? 1.0 : 0.0;
            Ch[r][col] = 0.0;
          }
        }
        if (s + 1 < nsg) {
          const int b1 = (s+1)*RTCOL3_NRD;
          for (int r=0; r<5; ++r) {
            for (int cc=0; cc<3; ++cc) {
              const Real g = c.wk(m,G0+3*r+cc,i1,k,j);
              if (g == 0.0) continue;
              rh[r] -= g*c.rd(m,b1+PP+rc3[cc],k,j);
              for (int col=0; col<5; ++col) {
                Ch[r][col] -= g*c.rd(m,b1+QQ+5*rc3[cc]+col,k,j);
              }
              // R's three columns are components 2,3,4 of the left unknown, which
              // for segment s+1 IS y_s -- so this lands in columns 2..4 of B, not 0,1,4
              for (int col=0; col<3; ++col) {
                Bh[r][col+2] -= g*c.rd(m,b1+RR+3*rc3[cc]+col,k,j);
              }
            }
          }
        }
        if (s > 0) {
          const int bm1 = (s-1)*RTCOL3_NRD;
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) {
              Real sa = 0.0;
              for (int t=0; t<3; ++t) {
                sa += c.wk(m,HH+3*r+t,i1,k,j)*c.rd(m,bm1+GR+5*(t+2)+col,k,j);
              }
              Bh[r][col] -= sa;
            }
            Real s2 = 0.0;
            for (int t=0; t<3; ++t) {
              s2 += c.wk(m,HH+3*r+t,i1,k,j)*c.rd(m,bm1+DR+t+2,k,j);
            }
            rh[r] -= s2;
          }
        }
        Real Bi[5][5];
        if (!RTCol3Inv5(Bh, Bi)) {
          for (int r=0; r<5; ++r) {
            c.rd(m,b0+DR+r,k,j) = 0.0;
            for (int col=0; col<5; ++col) c.rd(m,b0+GR+5*r+col,k,j) = 0.0;
          }
          continue;
        }
        for (int r=0; r<5; ++r) {
          Real sd = 0.0;
          for (int t=0; t<5; ++t) sd += Bi[r][t]*rh[t];
          c.rd(m,b0+DR+r,k,j) = sd;
          for (int col=0; col<5; ++col) {
            Real g = 0.0;
            for (int t=0; t<5; ++t) g += Bi[r][t]*Ch[t][col];
            c.rd(m,b0+GR+5*r+col,k,j) = g;
          }
        }
      }
      // the reduced back-substitution, with the SAME clamp, propagated the same way
      Real ynext[5];
      for (int r=0; r<5; ++r) ynext[r] = 0.0;
      for (int s=nsg-1; s>=0; --s) {
        const int b0 = s*RTCOL3_NRD;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real y[5];
        for (int r=0; r<5; ++r) {
          Real sy = c.rd(m,b0+DR+r,k,j);
          if (s + 1 < nsg) {
            for (int col=0; col<5; ++col) sy -= c.rd(m,b0+GR+5*r+col,k,j)*ynext[col];
          }
          y[r] = sy;
        }
        c.rd(m,b0+YR,k,j) = y[4];
        const Real b = c.wk(m,BB,i1,k,j);
        if (b > 0.0) {
          if (y[4] > 3.0*b) {
            y[4] = 3.0*b;
          } else if (y[4] < -0.75*b) {
            y[4] = -0.75*b;
          }
        }
        for (int r=0; r<5; ++r) {
          c.rd(m,b0+YY+r,k,j) = y[r];
          ynext[r] = y[r];
        }
      }
    });
    tm.team_barrier();
    // ---- 4c'. the segments, in parallel ---------------------------------------------
    Real dbm = 0.0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nsg),
    [&](const int s, Real &dmx) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      const int b0 = s*RTCOL3_NRD;
      int ncl = 0;
      Real yL[3] = {0.0, 0.0, 0.0};
      if (s > 0) {
        const int bm1 = (s-1)*RTCOL3_NRD;
        for (int r=0; r<3; ++r) yL[r] = c.rd(m,bm1+YY+r+2,k,j);
      }
      Real ynext[5];
      for (int r=0; r<5; ++r) ynext[r] = c.rd(m,b0+YY+r,k,j);
      // the boundary cell itself: the clamp is re-applied to the RAW increment, which is
      // the same arithmetic the reduced pass did, so the two agree bit for bit
      {
        const Real b = c.wk(m,BB,i1,k,j);
        Real db = c.rd(m,b0+YR,k,j);
        if (b > 0.0) {
          if (db > 3.0*b) {
            db = 3.0*b;
            ++ncl;
          } else if (db < -0.75*b) {
            db = -0.75*b;
            ++ncl;
          }
          const Real rel = fabs(db)/b;
          if (rel > dmx) dmx = rel;
          c.wk(m,BB,i1,k,j) = b + db;
        }
      }
      for (int i=i1-1; i>=i0; --i) {
        Real y[5];
        for (int r=0; r<5; ++r) {
          Real sy = c.wk(m,DP+r,i,k,j);
          for (int cc=0; cc<3; ++cc) {
            sy -= c.wk(m,G0+3*r+cc,i,k,j)*ynext[rc3[cc]]
                + c.wk(m,HH+3*r+cc,i,k,j)*yL[cc];
          }
          y[r] = sy;
        }
        const Real b = c.wk(m,BB,i,k,j);
        Real db = y[4];
        if (b > 0.0) {
          if (db > 3.0*b) {
            db = 3.0*b;
            ++ncl;
          } else if (db < -0.75*b) {
            db = -0.75*b;
            ++ncl;
          }
          const Real rel = fabs(db)/b;
          if (rel > dmx) dmx = rel;
          c.wk(m,BB,i,k,j) = b + db;
        }
        y[4] = db;
        for (int r=0; r<5; ++r) ynext[r] = y[r];
      }
      c.rd(m,b0+NCL,k,j) = static_cast<Real>(ncl);
    }, Kokkos::Max<Real>(dbm));
    if (dbm > dbmax) dbmax = dbm;
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      Real ncl = 0.0;
      for (int s=0; s<nsg; ++s) ncl += c.rd(m,s*RTCOL3_NRD+NCL,k,j);
      if (ncl > 0.0) Kokkos::atomic_add(&c.stat(3), ncl);
    });
    if (c.dstop && dbm < c.tol) break;
    if (dbm < 1.0e-14) break;
  }

  // ---- 5. apply, and the column energy budget --------------------------------------
  tm.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
    const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
    const int b0 = s*RTCOL3_NRD;
    Real budget = 0.0, bscale = 0.0, rtmax = 0.0, ubmax = 0.0;
    Real rhsum = 0.0, srsum = 0.0;
    for (int i=i0; i<=i1; ++i) {
      const Real dxi = c.Dx(m,k,j,i);
      const Real wb = c.taublend
          ? (1.0 - 0.5*(c.wblend(m,k,j,i) + c.wblend(m,k,j,i+1))) : 1.0;
      rhsum += c.bdt*(wb*c.wk(m,SA,i,k,j) + c.wk(m,EX,i,k,j))*dxi;
      srsum += c.wk(m,SA,i,k,j)*dxi;
      const Real tk = c.Tg(m,k,j,i);
      if (!(tk > 0.0)) continue;
      const Real rho = c.Rho(m,k,j,i);
      const Real es = c.Ei(m,k,j,i);
      const Real b = c.wk(m,BB,i,k,j);
      if (!(b > 0.0)) continue;
      const Real tnew = sqrt(sqrt(b/sopi));
      const Real enew = c.EFromT(rho, tnew);
      if (!(enew > 0.0) || !isfinite(enew)) {
        Kokkos::atomic_add(&c.stat(7), 1.0);
        continue;
      }
      const Real eref = c.EFromT(rho, tk);
      if (es > 0.0 && eref > 0.0) {
        const Real rt = fabs(eref - es)/es;
        if (rt > rtmax) rtmax = rt;
      }
      Real de = enew - es;
      const Real bz = c.Bb(m,0,i,k,j);
      if (bz > 0.0) {
        const Real rb = fabs(b/bz - 1.0);
        if (rb > ubmax) ubmax = rb;
      }
      if (!((es + de) > 0.0)) {
        de = -0.999*es;
        Kokkos::atomic_add(&c.stat(8), 1.0);
      }
      c.u0(m,IEN,k,j,i) += de;
      budget += de*dxi;
      bscale += fabs(de)*dxi;
    }
    c.rd(m,b0+BUD,k,j) = budget;
    c.rd(m,b0+BSC,k,j) = bscale;
    c.rd(m,b0+RTM,k,j) = rtmax;
    c.rd(m,b0+UBM,k,j) = ubmax;
    c.rd(m,b0+RHS,k,j) = rhsum;
    c.rd(m,b0+SRS,k,j) = srsum;
  });
  tm.team_barrier();
  Kokkos::single(Kokkos::PerTeam(tm), [&]() {
    Real budget = 0.0, bscale = 0.0, rtmax = 0.0, ubmax = 0.0;
    Real rhsum = 0.0, srsum = 0.0;
    for (int s=0; s<nsg; ++s) {
      const int b0 = s*RTCOL3_NRD;
      budget += c.rd(m,b0+BUD,k,j);
      bscale += c.rd(m,b0+BSC,k,j);
      rhsum += c.rd(m,b0+RHS,k,j);
      srsum += c.rd(m,b0+SRS,k,j);
      const Real rt = c.rd(m,b0+RTM,k,j);
      if (rt > rtmax) rtmax = rt;
      const Real ub = c.rd(m,b0+UBM,k,j);
      if (ub > ubmax) ubmax = ub;
    }
    Real fnet = 0.0;
    for (int q=0; q<nq; ++q) {
      const Real ftop = c.wf[q]*(c.wk(m,UU+q,ie,k,j) - Dtop[q]);
      const Real fcut = c.wf[q]*((c.wk(m,BB,ic,k,j) + Ucut[q]) - c.wk(m,DD+q,ic,k,j));
      fnet += fcut - ftop;
    }
    Real ftop3 = 0.0;
    for (int q=0; q<nq; ++q) ftop3 += c.wf[q]*(c.wk(m,UU+q,ie,k,j) - Dtop[q]);
    Kokkos::atomic_add(&c.stat(12), rhsum);
    Kokkos::atomic_add(&c.stat(13), srsum);
    Kokkos::atomic_add(&c.stat(14), fnet);
    Kokkos::atomic_add(&c.stat(15), fabs(srsum) + fabs(fnet));
    Kokkos::atomic_add(&c.stat(16), ftop3);
    Kokkos::atomic_add(&c.stat(17), c.Fb(m,0,ie+1,k,j));
    Kokkos::atomic_add(&c.stat(0), static_cast<Real>(nit));
    Kokkos::atomic_add(&c.stat(1), 1.0);
    Kokkos::atomic_max(&c.stat(2), static_cast<Real>(nit));
    Kokkos::atomic_max(&c.stat(4), dbmax);
    Kokkos::atomic_max(&c.stat(5), rtmax);
    Kokkos::atomic_max(&c.stat(18), ubmax);
    Kokkos::atomic_add(&c.stat(9), budget);
    Kokkos::atomic_add(&c.stat(10), bscale);
    Kokkos::atomic_max(&c.stat(19), rfin);
    Kokkos::atomic_add(&c.stat(20), rfin);
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RTCol3TeamLaunch
//! \brief one TEAM per column.  The team size IS the number of segments, so it is what
//! sets both the parallelism and the size of the reduced system; <problem>/rt_impl_nseg
//! asks for it and Kokkos clamps it to what the backend allows (1 on Serial, where the
//! partition degenerates to exactly the serial block Thomas).

inline void RTCol3TeamLaunch(const RTCol3 &c, const int nmb1, const int ks, const int ke,
                             const int js, const int je) {
  RTCol3 cc = c;
  const int nk = ke - ks + 1, nj = je - js + 1;
  const int nteam = (nmb1 + 1)*nk*nj;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nteam, cc.nseg);
  Kokkos::parallel_for("rt_col3_p", policy, KOKKOS_LAMBDA(TeamMember_t tm) {
    const int q = tm.league_rank();
    const int m = q/(nk*nj);
    const int k = ks + (q - m*nk*nj)/nj;
    const int j = js + (q - m*nk*nj)%nj;
    RTCol3TeamSolve(cc, tm, m, k, j);
  });
}

}  // namespace two_stream_rt
#endif  // UTILS_TWO_STREAM_COLUMN_PARTITION_HPP_
