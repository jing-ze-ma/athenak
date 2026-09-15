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
//! THE HYBRID (<problem>/rt_col3_hybrid_tau) IS IMPLEMENTED HERE.  On the same equal
//! partition, a segment is DEEP if its top cell lies below the interface isp (every cell
//! a scalar diffusion row: one unknown, one divide, no 5x5 inverse), THIN if its bottom
//! cell lies at or above isp, and the one segment that straddles isp carries the SPLICE
//! -- the 5x1 interface column into the thin row and the 1x3 interface row on the deep
//! top face -- inside its own lane's recurrence, exactly as the serial hybrid forms it.
//! The reduced system over the segment boundaries then has MIXED block sizes (1 for a
//! deep boundary, 5 for a thin one); it is carried as 5x5 rows with the deep ones PADDED
//! by a unit diagonal on the four dead components.  That padding is exact -- those
//! components' right-hand sides and off-diagonal rows are identically zero, so they
//! solve to zero -- and it leaves BOTH reduced solvers (the serial block Thomas and the
//! PCR of rt_impl_redpar) untouched.  It wastes one 5x5 inverse per deep segment, nsg of
//! them against nc cells, while every deep CELL drops from a 5x5 assembly plus inverse
//! to a single divide.
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
// ... and the same, with room for the PARALLEL CYCLIC REDUCTION of the reduced system
// (<problem>/rt_impl_redpar): two buffers of the block row (A, B, C, r = 80 slots) plus
// the diagonal block's inverse.  Allocated only when the switch is on.
#define RTCOL3_NRDP 281

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
  // FL aliases DD: a DEEP cell has no downward intensity unknown and the slot carries
  // that cell's LOWER-face diffusion flux instead.
  const int FL = 27;
  // row 4 -- the b row -- is the only live row of a deep cell's G, H and d
  const int GD = G0 + 12, HHD = HH + 12, DPD = DP + 4;
  // per-segment slots of c.rd(m, k, j, s*nrd + slot)
  const int PP = 0, QQ = 5, RR = 30, GR = 45, DR = 70, DE = 75, UE = 77, YY = 79,
            YR = 84, BUD = 85, BSC = 86, RTM = 87, UBM = 88, RHS = 89, SRS = 90,
            NCL = 91;
  // the PCR slots of the extended per-segment workspace: A, B, C, r in two buffers
  // (stride PSTR), then the diagonal block's inverse
  const int PA = 96, PSTR = 80, PBI = 256;
  const int nrd = c.redpar ? RTCOL3_NRDP : RTCOL3_NRD;
  // the reduced system is solved by PCR only when lane s can own block row s and the
  // round count is exactly log2(nsg); otherwise the serial reduced solve runs.
  const bool usepcr = c.redpar && (nsg >= 2) && (nsg <= 64)
                      && ((nsg & (nsg - 1)) == 0);
  const int nq = c.nq;
  const Real sopi = c.sigma/M_PI;
  // the three column components of the reduced (x[0], x[1], x[4]) coupling
  const int rc3[3] = {0, 1, 4};

  // ---- 0. THE HYBRID SPLIT (problem/rt_col3_hybrid_tau) ----------------------------
  // The same interface the serial path picks, on the same partition: nsg equal segments
  // over [ic, ie], so a segment is DEEP if i1 < ib (every cell a SCALAR diffusion row,
  // one unknown, one divide), THIN if i0 >= ib (the 5x5 machinery unchanged), and the
  // one segment straddling ib carries the SPLICE -- its deep cells scalar, its thin
  // cells 5x5, and the 5x1 / 1x3 interface pair inside that one lane's own recurrence,
  // exactly as the serial hybrid forms it.
  //
  // THE REDUCED SYSTEM then has MIXED block sizes: a deep segment's boundary unknown is
  // the scalar b, a thin one's the five (D0, D1, U0, U1, b).  Rather than carry two
  // block types through both reduced solvers, the deep rows are PADDED to 5x5 with a
  // unit diagonal on the four dead components: their right-hand sides and their
  // off-diagonal rows are identically zero, so those components solve to exactly zero
  // and the padded system is the mixed system, bit for bit.  Both the serial reduced
  // Thomas and the PCR of rt_impl_redpar therefore run UNCHANGED.  The padding wastes
  // one 5x5 inverse per deep segment -- nsg of them against nc cells -- while the
  // per-cell saving is the whole 5x5 assembly, inverse and back-substitution.
  RTCol3Hyb hb;
  hb.isp = ic;
  Real kflx = 0.0, wsum = 0.0;
  for (int q=0; q<nq; ++q) {
    kflx += 2.0*c.wf[q]*c.mu[q];
    wsum += c.wf[q];
  }
  {
    const int isp0 = c.Interface(m, k, j, ic);
    if (isp0 > ic) {
      const Real hm = c.Ht(m,k,j,isp0-1), h0 = c.Ht(m,k,j,isp0);
      if (hm + h0 > 0.0) {
        hb.on = true;
        hb.isp = isp0;
        hb.dtc = hm + h0;
        hb.g = h0/hb.dtc;
        for (int q=0; q<nq; ++q) hb.am[q] = hb.g + c.mu[q]/hb.dtc;
      }
    }
  }
  const int ib = hb.isp;              // the bottom cell of the TWO-STREAM segment

  // ---- 1. the frozen per-cell layer coefficients, and b^0 --------------------------
  Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
    const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
    for (int i=i0; i<=i1; ++i) {
      const Real h = c.Ht(m,k,j,i);
      for (int q=0; q<nq; ++q) {
        const Real x = h/c.mu[q];
        const Real e0 = -expm1(-x);
        c.Wk<true>(m,EE+q,i,k,j) = e0;
        c.Wk<true>(m,CI+q,i,k,j) = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
        c.Wk<true>(m,CO+q,i,k,j) = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0 - SQR(x)/6.0);
      }
      c.Wk<true>(m,BB,i,k,j) = RTCol3WarmStart(c, m, k, j, i, c.Bb(m,0,i,k,j));
      c.Wk<true>(m,ES,i,k,j) = c.Ei(m,k,j,i);
    }
  });
  tm.team_barrier();
  Real emax = 0.0;
  if (c.norm == 1) {
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
    [&](const int i, Real &mx) {
      const Real e = c.Wk<true>(m,ES,i,k,j);
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
  // the deep segment's BOTTOM face flux, the two-stream cut flux in its deep limit
  const Real fbot = kflx*dbdtau + (c.int_at_cut ? wsum*c.Iint : 0.0);
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
      c.Wk<true>(m,EX,i,k,j) = ex + c.Qb(m,0,i,k,j);
    }
  });

  // ---- 4. Newton ---------------------------------------------------------------------
  int nit = 0;
  Real dbmax = 0.0;
  Real rfin = 0.0;
  // THE INTERFACE STATE, re-formed on every lane from the current b at every pass: the
  // deep limit the thin segment's bottom cell sees in place of the cut's frozen Ucut.
  // Equal to the cut values when the hybrid is off, which keeps that path bitwise.
  Real cutc_i = cutc;
  Real Ucut_i[2] = {Ucut[0], Ucut[1]};
  Real fif = 0.0;                             // the interface net flux
  for (int it=0; it<c.maxit; ++it) {
    nit = it + 1;
    tm.team_barrier();
    if (hb.on) {
      const Real dbd = (c.Wk<true>(m,BB,ib-1,k,j) - c.Wk<true>(m,BB,ib,k,j))/hb.dtc;
      cutc_i = dbd*c.Ht(m,k,j,ib);
      for (int q=0; q<nq; ++q) Ucut_i[q] = cutc_i + c.mu[q]*dbd;
    }
    // ---- 4a. the formal solution, each segment with a ZERO incoming intensity --------
    // problem/rt_impl_ablate bit 1 repeats the SEGMENT SWEEPS, which write absolutely
    // and are therefore idempotent.  The entry-value scan and the correction that follow
    // accumulate, so they are not repeated and their cost is not measured here.
    for (int rep=(c.ablate & 1); rep>=0; --rep) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      if (i1 < ib) return;                    // a wholly DEEP segment: no two-stream here
      const int t0 = (i0 > ib) ? i0 : ib;     // this segment's bottom THIN cell
      Real L[2] = {0.0, 0.0}, hg[2] = {1.0, 1.0};
      for (int i=i1; i>=t0; --i) {
        Real sl, su, sfu, sfd;
        c.SourceVals<true>(m, k, j, i, ib, hb, cutc, sl, su, sfu, sfd);
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.Wk<true>(m,EE+q,i,k,j), t = 1.0 - E;
          const Real ci = c.Wk<true>(m,CI+q,i,k,j), co = c.Wk<true>(m,CO+q,i,k,j);
          const Real P = ci*sfu + co*sl;
          const Real Q = ci*su + co*sfd;
          acc += c.wf[q]*W*(E*(1.0 + t)*L[q] + E*P - (P + Q));
          L[q] = t*t*L[q] + t*P + Q;
          hg[q] = t*t*hg[q];
          c.Wk<true>(m,DD+q,i,k,j) = L[q];
          c.Wk<true>(m,HD+q,i,k,j) = hg[q];
        }
        c.Wk<true>(m,SA,i,k,j) = acc;
      }
      Real Lu[2] = {0.0, 0.0}, hu[2] = {1.0, 1.0};
      for (int i=t0; i<=i1; ++i) {
        Real sl, su, sfu, sfd;
        c.SourceVals<true>(m, k, j, i, ib, hb, cutc, sl, su, sfu, sfd);
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.Wk<true>(m,EE+q,i,k,j), t = 1.0 - E;
          const Real ci = c.Wk<true>(m,CI+q,i,k,j), co = c.Wk<true>(m,CO+q,i,k,j);
          const Real Pu = ci*sfd + co*su;
          const Real Qu = ci*sl + co*sfu;
          acc += c.wf[q]*W*(E*(1.0 + t)*Lu[q] + E*Pu - (Pu + Qu));
          Lu[q] = t*t*Lu[q] + t*Pu + Qu;
          hu[q] = t*t*hu[q];
          c.Wk<true>(m,UU+q,i,k,j) = Lu[q];
          c.Wk<true>(m,HU+q,i,k,j) = hu[q];
        }
        c.Wk<true>(m,SA,i,k,j) += acc;
      }
    });
    tm.team_barrier();
    }
    // the P segment entry intensities: a serial scan over SEGMENTS, not over cells
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      Real de[2], ue[2];
      for (int q=0; q<nq; ++q) de[q] = Dtop[q];
      for (int s=nsg-1; s>=0; --s) {
        const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
        if (i1 < ib) continue;
        const int t0 = (i0 > ib) ? i0 : ib;
        for (int q=0; q<nq; ++q) c.rd(m,k,j,s*nrd+DE+q) = de[q];
        for (int q=0; q<nq; ++q) {
          de[q] = c.Wk<true>(m,DD+q,t0,k,j) + c.Wk<true>(m,HD+q,t0,k,j)*de[q];
        }
      }
      for (int q=0; q<nq; ++q) ue[q] = c.Wk<true>(m,BB,ib,k,j) + Ucut_i[q];
      for (int s=0; s<nsg; ++s) {
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        if (i1 < ib) continue;
        for (int q=0; q<nq; ++q) c.rd(m,k,j,s*nrd+UE+q) = ue[q];
        for (int q=0; q<nq; ++q) {
          ue[q] = c.Wk<true>(m,UU+q,i1,k,j) + c.Wk<true>(m,HU+q,i1,k,j)*ue[q];
        }
      }
    });
    tm.team_barrier();
    // and the correction each cell owes to its segment's entry value
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      if (i1 < ib) return;
      const int t0 = (i0 > ib) ? i0 : ib;
      Real de[2], ue[2];
      for (int q=0; q<nq; ++q) {
        de[q] = c.rd(m,k,j,s*nrd+DE+q);
        ue[q] = c.rd(m,k,j,s*nrd+UE+q);
      }
      for (int i=t0; i<=i1; ++i) {
        const Real W = 1.0/c.Dx(m,k,j,i);
        Real acc = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real E = c.Wk<true>(m,EE+q,i,k,j), t = 1.0 - E;
          const Real hd = (i == i1) ? 1.0 : c.Wk<true>(m,HD+q,i+1,k,j);
          const Real hup = (i == t0) ? 1.0 : c.Wk<true>(m,HU+q,i-1,k,j);
          acc += c.wf[q]*W*E*(1.0 + t)*(hd*de[q] + hup*ue[q]);
          c.Wk<true>(m,DD+q,i,k,j) += c.Wk<true>(m,HD+q,i,k,j)*de[q];
          c.Wk<true>(m,UU+q,i,k,j) += c.Wk<true>(m,HU+q,i,k,j)*ue[q];
        }
        c.Wk<true>(m,SA,i,k,j) += acc;
      }
    });
    tm.team_barrier();
    // ---- 4a''. THE DEEP SEGMENT: the diffusion fluxes and their divergence --------
    // The interface flux is the two-stream's OWN net flux at that face, built from the
    // deep-limit U (implicit in b(isp-1), b(isp)) and the solved D(isp), so the deep
    // cell loses upward exactly what the thin cell gains at its lower face and the
    // column still telescopes across the interface.  Two passes with a barrier between
    // them: F is a face quantity each lane writes for its own cells, and the divergence
    // reads the face above, which for the top cell of a lane belongs to the next lane.
    if (hb.on) {
      fif = 0.0;
      for (int q=0; q<nq; ++q) {
        fif += c.wf[q]*((c.Wk<true>(m,BB,ib,k,j) + Ucut_i[q])
                        - c.Wk<true>(m,DD+q,ib,k,j));
      }
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
        const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
        const int d1 = (i1 < ib-1) ? i1 : (ib-1);
        for (int i=i0; i<=d1; ++i) {
          Real flo = fbot;
          if (i > ic) {
            const Real dtm = c.Ht(m,k,j,i-1) + c.Ht(m,k,j,i);
            flo = (dtm > 0.0) ? kflx*(c.Wk<true>(m,BB,i-1,k,j)
                                      - c.Wk<true>(m,BB,i,k,j))/dtm : 0.0;
          }
          c.Wk<true>(m,FL,i,k,j) = flo;
        }
      });
      tm.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
        const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
        const int d1 = (i1 < ib-1) ? i1 : (ib-1);
        for (int i=i0; i<=d1; ++i) {
          const Real fhi = (i + 1 < ib) ? c.Wk<true>(m,FL,i+1,k,j) : fif;
          c.Wk<true>(m,SA,i,k,j) = (c.Wk<true>(m,FL,i,k,j) - fhi)/c.Dx(m,k,j,i);
        }
      });
      tm.team_barrier();
    }

    // ---- 4a'. the handover, re-formed from THIS column's OWN flux ------------------
    if (c.ex_iter && c.taublend && c.direct) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
        const Real wlo = c.wblend(m,k,j,i), whi = c.wblend(m,k,j,i+1);
        if (i < ib) {                         // a DEEP cell: the diffusion faces
          const Real fhi = (i + 1 < ib) ? c.Wk<true>(m,FL,i+1,k,j) : fif;
          c.Wk<true>(m,EX,i,k,j) = 0.5*(wlo + whi)*c.Wk<true>(m,SA,i,k,j)
              + (whi*fhi - wlo*c.Wk<true>(m,FL,i,k,j))/c.Dx(m,k,j,i) + c.Qb(m,0,i,k,j);
          return;
        }
        Real f3lo = 0.0, f3hi = 0.0;
        for (int q=0; q<nq; ++q) {
          const Real ulo = (i == ib) ? (c.Wk<true>(m,BB,ib,k,j) + Ucut_i[q])
                                     : c.Wk<true>(m,UU+q,i-1,k,j);
          const Real dhi = (i == ie) ? Dtop[q] : c.Wk<true>(m,DD+q,i+1,k,j);
          f3lo += c.wf[q]*(ulo - c.Wk<true>(m,DD+q,i,k,j));
          f3hi += c.wf[q]*(c.Wk<true>(m,UU+q,i,k,j) - dhi);
        }
        c.Wk<true>(m,EX,i,k,j) = 0.5*(wlo + whi)*c.Wk<true>(m,SA,i,k,j)
                         + (whi*f3hi - wlo*f3lo)/c.Dx(m,k,j,i) + c.Qb(m,0,i,k,j);
      });
      tm.team_barrier();
    }

    // ---- 4a''. THE RESIDUAL-ONLY CONVERGENCE TEST (problem/rt_impl_rescheck) --------
    // The same test step 4b makes, made before the factorisation instead of after it, so
    // that a pass which only confirms convergence never assembles or inverts a block.
    // ResidRel is the shared helper both solver paths use, and it reproduces
    // fabs(rv[4])/(rsc + eoff) of BuildRow bit for bit, so the break is bitwise the
    // break the standard path takes.  The partition is the same i0..i1 the sweep used,
    // so each thread reads only cells it has just written; the team reduce broadcasts
    // rpre, which is what makes the collective break legal.
    if (c.rescheck || (c.ablate & 4)) {
      Real rpre = 0.0;
      for (int rep=((c.ablate & 4) ? 1 : 0); rep>=0; --rep) {
        rpre = 0.0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nsg),
        [&](const int s, Real &rmx) {
          const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
          for (int i=i0; i<=i1; ++i) {
            const Real rr = c.ResidRel<true>(m, k, j, i, eoff);
            if (rr > rmx) rmx = rr;
          }
        }, Kokkos::Max<Real>(rpre));
      }
      if (c.rescheck) {
        rfin = rpre;
        if (rpre < c.tol && !c.fixit) break;
      }
    }

    // ---- 4b. the PARTITIONED forward elimination ------------------------------------
    // ablate bit 2 repeats it: it writes absolutely out of b, so it is idempotent.
    Real rmax = 0.0;
    for (int rep=(c.ablate & 2); rep>=0; --rep) {
    rmax = 0.0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nsg),
    [&](const int s, Real &rmx) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      Real Gp[5][3], Hp[5][3], dp[5];
      for (int i=i0; i<=i1; ++i) {
        if (i < ib) {
          // ---- the DEEP cell: a SCALAR row, padded into the b slot ------------------
          // x_i = d_i - g_i . x_{i+1}|(D0,D1,b) - h_i . x_L|(U0,U1,b), all three of them
          // scalars times a 3-vector, and the four dead components of this cell carry
          // nothing: their d, G and H rows are written as zero ONLY at the segment
          // boundary cell, which is the only one the reduced system reads in full.
          Real aa, dd, cup, C5[3], rvd, rscd;
          c.DeepRow<true>(m, k, j, i, ic, hb, kflx, it, aa, dd, cup, C5, rvd, rscd);
          const Real dend = rscd + eoff;
          const Real rrd = (dend > 0.0) ? fabs(rvd)/dend : 0.0;
          if (rrd > rmx) rmx = rrd;
          Real cv[3] = {0.0, 0.0, cup};       // the upper coupling, on (D0, D1, b)
          if (i + 1 >= ib) {
            for (int t=0; t<3; ++t) cv[t] = C5[t];   // the 1x3 INTERFACE row
          }
          Real hv[3] = {0.0, 0.0, 0.0};
          Real piv = dd, rvv = rvd;
          if (i > i0) {
            piv -= aa*Gp[4][2];
            rvv -= aa*dp[4];
            for (int t=0; t<3; ++t) hv[t] = -aa*Hp[4][t];
          } else {
            hv[2] = aa;                       // the spike on the left segment's unknown
          }
          if (!(fabs(piv) > 0.0)) {
            Kokkos::atomic_add(&c.stat(6), 1.0);
            rmx = HUGE_VAL;
            return;
          }
          dp[4] = rvv/piv;
          dp[2] = 0.0;
          dp[3] = 0.0;
          c.Wk<true>(m,DPD,i,k,j) = dp[4];
          for (int t=0; t<3; ++t) {
            Gp[2][t] = 0.0;
            Gp[3][t] = 0.0;
            Hp[2][t] = 0.0;
            Hp[3][t] = 0.0;
            Gp[4][t] = cv[t]/piv;
            Hp[4][t] = hv[t]/piv;
            c.Wk<true>(m,GD+t,i,k,j) = Gp[4][t];
            c.Wk<true>(m,HHD+t,i,k,j) = Hp[4][t];
          }
          if (i == i1) {          // pad the boundary row for the reduced system
            for (int r=0; r<4; ++r) {
              c.Wk<true>(m,DP+r,i,k,j) = 0.0;
              for (int t=0; t<3; ++t) {
                c.Wk<true>(m,G0+3*r+t,i,k,j) = 0.0;
                c.Wk<true>(m,HH+3*r+t,i,k,j) = 0.0;
              }
            }
          }
          continue;
        }
        Real A3[5][3], Bm[5][5], C3[5][3], rv[5], AH[5][3];
        Real rsc = 1.0;
        Real A1z[5];
        c.BuildRow<true>(m, k, j, i, ib, hb, cutc, it, A3, A1z, Bm, C3, rv, rsc);
        if (i == ib && hb.on) {
          // THE SPLICE.  The thin segment's bottom cell couples to the deep unknown
          // b(ib-1) through the 5x1 column A1; in the A3 mapping (the columns are
          // components 2, 3, 4 = U0, U1, b of the cell below) that is column 2, so the
          // shared elimination below needs no other special case.
          for (int r=0; r<5; ++r) {
            A3[r][0] = 0.0;
            A3[r][1] = 0.0;
            A3[r][2] = A1z[r];
          }
        }
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
          c.Wk<true>(m,DP+r,i,k,j) = sd;
          for (int cc=0; cc<3; ++cc) {
            Real g = 0.0, h = 0.0;
            for (int t=0; t<5; ++t) {
              g += Bi[r][t]*C3[t][cc];
              h += Bi[r][t]*AH[t][cc];
            }
            Gp[r][cc] = g;
            Hp[r][cc] = sg*h;
            c.Wk<true>(m,G0+3*r+cc,i,k,j) = g;
            c.Wk<true>(m,HH+3*r+cc,i,k,j) = sg*h;
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
        if (i < ib) {
          // the DEEP cell: only the b row of (p, Q, R) is live.  The four dead rows are
          // zeroed once, when the recurrence first enters the deep part -- either at the
          // interface or at a wholly deep segment's own boundary cell -- and stay zero.
          Real p4 = c.Wk<true>(m,DPD,i,k,j);
          Real Q4[5], R4[3];
          for (int col=0; col<5; ++col) Q4[col] = 0.0;
          for (int col=0; col<3; ++col) R4[col] = c.Wk<true>(m,HHD+col,i,k,j);
          for (int cc=0; cc<3; ++cc) {
            const Real g = c.Wk<true>(m,GD+cc,i,k,j);
            if (g == 0.0) continue;
            p4 -= g*pq[rc3[cc]];
            for (int col=0; col<5; ++col) Q4[col] -= g*Qm[rc3[cc]][col];
            for (int col=0; col<3; ++col) R4[col] -= g*Rm[rc3[cc]][col];
          }
          if (i == ib-1 || i == i1-1) {
            for (int r=0; r<4; ++r) {
              pq[r] = 0.0;
              for (int col=0; col<5; ++col) Qm[r][col] = 0.0;
              for (int col=0; col<3; ++col) Rm[r][col] = 0.0;
            }
          }
          pq[4] = p4;
          for (int col=0; col<5; ++col) Qm[4][col] = Q4[col];
          for (int col=0; col<3; ++col) Rm[4][col] = R4[col];
          continue;
        }
        Real pn[5], Qn[5][5], Rn[5][3];
        for (int r=0; r<5; ++r) {
          Real sp = c.Wk<true>(m,DP+r,i,k,j);
          for (int cc=0; cc<3; ++cc) sp -= c.Wk<true>(m,G0+3*r+cc,i,k,j)*pq[rc3[cc]];
          pn[r] = sp;
          for (int col=0; col<5; ++col) {
            Real sq = 0.0;
            for (int cc=0; cc<3; ++cc) sq += c.Wk<true>(m,G0+3*r+cc,i,k,j)*Qm[rc3[cc]][col];
            Qn[r][col] = -sq;
          }
          for (int col=0; col<3; ++col) {
            Real sr = c.Wk<true>(m,HH+3*r+col,i,k,j);
            for (int cc=0; cc<3; ++cc) sr -= c.Wk<true>(m,G0+3*r+cc,i,k,j)*Rm[rc3[cc]][col];
            Rn[r][col] = sr;
          }
        }
        for (int r=0; r<5; ++r) {
          pq[r] = pn[r];
          for (int col=0; col<5; ++col) Qm[r][col] = Qn[r][col];
          for (int col=0; col<3; ++col) Rm[r][col] = Rn[r][col];
        }
      }
      const int b0 = s*nrd;
      for (int r=0; r<5; ++r) {
        c.rd(m,k,j,b0+PP+r) = pq[r];
        for (int col=0; col<5; ++col) c.rd(m,k,j,b0+QQ+5*r+col) = Qm[r][col];
        for (int col=0; col<3; ++col) c.rd(m,k,j,b0+RR+3*r+col) = Rm[r][col];
      }
    }, Kokkos::Max<Real>(rmax));
    }
    rfin = rmax;
    if (!(rmax < HUGE_VAL)) return;     // a singular block: leave the column alone
    if (rmax < c.tol && !c.fixit) break;

    // ---- 4c. the REDUCED system, and the clamped update -----------------------------
    // The segment boundary unknowns y_s = x_{i1(s)} satisfy, exactly,
    //     A_s y_{s-1} + B_s y_s + C_s y_{s+1} = r_s,   s = 0 .. nsg-1,
    // a NON-periodic block-tridiagonal system of nsg 5x5 rows (the column has a top and
    // a bottom wall, so there is no wrap-around row).  A_s is the spike H of the
    // boundary cell and touches only components 2..4 of y_{s-1}; C_s is that cell's G
    // composed with the next segment's Q; the next segment's R lands in B_s.  The
    // default path eliminates that system by the SERIAL block Thomas below, one lane
    // walking nsg rows while the other nsg-1 idle.
    tm.team_barrier();
    if (usepcr) {
      // <problem>/rt_impl_redpar: PARALLEL CYCLIC REDUCTION instead.  Lane s owns block
      // row s; each of the log2(nsg) rounds eliminates the neighbours at distance d by
      // A'_s = -A_s B_{s-d}^-1 A_{s-d},  C'_s = -C_s B_{s+d}^-1 C_{s+d},
      // B'_s = B_s - A_s B_{s-d}^-1 C_{s-d} - C_s B_{s+d}^-1 A_{s+d}, and the same for
      // r, with out-of-range blocks zero.  After the last round A = C = 0 on every lane
      // and y_s = B_s^-1 r_s directly: there is NO reduced back-substitution left.  The
      // arithmetic is a different elimination ORDER of the same system, so it agrees
      // with the serial solve to round-off (and the clamp, which the serial pass
      // propagates down the reduced chain, is applied per segment here -- the same
      // difference the partitioned path already carries, and zero once converged).
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
        const int b0 = s*nrd, pb = b0 + PA;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real Bh[5][5], Ch[5][5], Ah[5][5], rh[5];
        for (int r=0; r<5; ++r) {
          rh[r] = c.Wk<true>(m,DP+r,i1,k,j);
          for (int col=0; col<5; ++col) {
            Bh[r][col] = (r == col) ? 1.0 : 0.0;
            Ch[r][col] = 0.0;
            Ah[r][col] = 0.0;
          }
        }
        if (s + 1 < nsg) {
          const int b1 = (s+1)*nrd;
          for (int r=0; r<5; ++r) {
            for (int cc=0; cc<3; ++cc) {
              const Real g = c.Wk<true>(m,G0+3*r+cc,i1,k,j);
              if (g == 0.0) continue;
              rh[r] -= g*c.rd(m,k,j,b1+PP+rc3[cc]);
              for (int col=0; col<5; ++col) {
                Ch[r][col] -= g*c.rd(m,k,j,b1+QQ+5*rc3[cc]+col);
              }
              for (int col=0; col<3; ++col) {
                Bh[r][col+2] -= g*c.rd(m,k,j,b1+RR+3*rc3[cc]+col);
              }
            }
          }
        }
        if (s > 0) {
          for (int r=0; r<5; ++r) {
            for (int t=0; t<3; ++t) {
              Ah[r][t+2] = c.Wk<true>(m,HH+3*r+t,i1,k,j);
            }
          }
        }
        for (int r=0; r<5; ++r) {
          for (int col=0; col<5; ++col) {
            c.rd(m,k,j,pb+5*r+col) = Ah[r][col];
            c.rd(m,k,j,pb+25+5*r+col) = Bh[r][col];
            c.rd(m,k,j,pb+50+5*r+col) = Ch[r][col];
          }
          c.rd(m,k,j,pb+75+r) = rh[r];
        }
      });
      int pcur = 0;
      for (int d=1; d<nsg; d*=2) {
        tm.team_barrier();
        // the diagonal block's inverse, one per lane
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
          const int b0 = s*nrd, pb = b0 + PA + pcur*PSTR;
          Real Bm[5][5], Bi[5][5];
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) Bm[r][col] = c.rd(m,k,j,pb+25+5*r+col);
          }
          if (!RTCol3Inv5(Bm, Bi)) {
            for (int r=0; r<5; ++r) {
              for (int col=0; col<5; ++col) Bi[r][col] = 0.0;
            }
          }
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) c.rd(m,k,j,b0+PBI+5*r+col) = Bi[r][col];
          }
        });
        tm.team_barrier();
        // the round itself: read buffer pcur, write buffer 1-pcur, no hazard
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
          const int b0 = s*nrd, pb = b0 + PA + pcur*PSTR;
          const int pn = b0 + PA + (1-pcur)*PSTR;
          Real nA[5][5], nB[5][5], nC[5][5], nr[5], W[5][5];
          for (int r=0; r<5; ++r) {
            nr[r] = c.rd(m,k,j,pb+75+r);
            for (int col=0; col<5; ++col) {
              nA[r][col] = 0.0;
              nC[r][col] = 0.0;
              nB[r][col] = c.rd(m,k,j,pb+25+5*r+col);
            }
          }
          if (s - d >= 0) {
            const int bm = (s-d)*nrd, pm = bm + PA + pcur*PSTR;
            for (int r=0; r<5; ++r) {
              for (int col=0; col<5; ++col) {
                Real w = 0.0;
                for (int t=0; t<5; ++t) {
                  w += c.rd(m,k,j,pb+5*r+t)*c.rd(m,k,j,bm+PBI+5*t+col);
                }
                W[r][col] = w;
              }
            }
            for (int r=0; r<5; ++r) {
              for (int col=0; col<5; ++col) {
                Real sa = 0.0, sc = 0.0;
                for (int t=0; t<5; ++t) {
                  sa += W[r][t]*c.rd(m,k,j,pm+5*t+col);
                  sc += W[r][t]*c.rd(m,k,j,pm+50+5*t+col);
                }
                nA[r][col] = -sa;
                nB[r][col] -= sc;
              }
              Real sr = 0.0;
              for (int t=0; t<5; ++t) sr += W[r][t]*c.rd(m,k,j,pm+75+t);
              nr[r] -= sr;
            }
          }
          if (s + d < nsg) {
            const int bp = (s+d)*nrd, pp = bp + PA + pcur*PSTR;
            for (int r=0; r<5; ++r) {
              for (int col=0; col<5; ++col) {
                Real w = 0.0;
                for (int t=0; t<5; ++t) {
                  w += c.rd(m,k,j,pb+50+5*r+t)*c.rd(m,k,j,bp+PBI+5*t+col);
                }
                W[r][col] = w;
              }
            }
            for (int r=0; r<5; ++r) {
              for (int col=0; col<5; ++col) {
                Real sc = 0.0, sa = 0.0;
                for (int t=0; t<5; ++t) {
                  sc += W[r][t]*c.rd(m,k,j,pp+50+5*t+col);
                  sa += W[r][t]*c.rd(m,k,j,pp+5*t+col);
                }
                nC[r][col] = -sc;
                nB[r][col] -= sa;
              }
              Real sr = 0.0;
              for (int t=0; t<5; ++t) sr += W[r][t]*c.rd(m,k,j,pp+75+t);
              nr[r] -= sr;
            }
          }
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) {
              c.rd(m,k,j,pn+5*r+col) = nA[r][col];
              c.rd(m,k,j,pn+25+5*r+col) = nB[r][col];
              c.rd(m,k,j,pn+50+5*r+col) = nC[r][col];
            }
            c.rd(m,k,j,pn+75+r) = nr[r];
          }
        });
        pcur = 1 - pcur;
      }
      tm.team_barrier();
      // every lane now holds a DECOUPLED 5x5 row: y_s = B_s^-1 r_s, then the clamp
      Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
        const int b0 = s*nrd, pb = b0 + PA + pcur*PSTR;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real Bm[5][5], Bi[5][5], y[5];
        for (int r=0; r<5; ++r) {
          for (int col=0; col<5; ++col) Bm[r][col] = c.rd(m,k,j,pb+25+5*r+col);
        }
        if (!RTCol3Inv5(Bm, Bi)) {
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) Bi[r][col] = 0.0;
          }
        }
        for (int r=0; r<5; ++r) {
          Real sy = 0.0;
          for (int t=0; t<5; ++t) sy += Bi[r][t]*c.rd(m,k,j,pb+75+t);
          y[r] = sy;
        }
        c.rd(m,k,j,b0+YR) = y[4];
        const Real b = c.Wk<true>(m,BB,i1,k,j);
        if (b > 0.0) {
          if (y[4] > 3.0*b) {
            y[4] = 3.0*b;
          } else if (y[4] < -0.75*b) {
            y[4] = -0.75*b;
          }
        }
        for (int r=0; r<5; ++r) c.rd(m,k,j,b0+YY+r) = y[r];
      });
    } else {
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      for (int s=0; s<nsg; ++s) {
        const int b0 = s*nrd;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real Bh[5][5], Ch[5][5], rh[5];
        for (int r=0; r<5; ++r) {
          rh[r] = c.Wk<true>(m,DP+r,i1,k,j);
          for (int col=0; col<5; ++col) {
            Bh[r][col] = (r == col) ? 1.0 : 0.0;
            Ch[r][col] = 0.0;
          }
        }
        if (s + 1 < nsg) {
          const int b1 = (s+1)*nrd;
          for (int r=0; r<5; ++r) {
            for (int cc=0; cc<3; ++cc) {
              const Real g = c.Wk<true>(m,G0+3*r+cc,i1,k,j);
              if (g == 0.0) continue;
              rh[r] -= g*c.rd(m,k,j,b1+PP+rc3[cc]);
              for (int col=0; col<5; ++col) {
                Ch[r][col] -= g*c.rd(m,k,j,b1+QQ+5*rc3[cc]+col);
              }
              // R's three columns are components 2,3,4 of the left unknown, which
              // for segment s+1 IS y_s -- so this lands in columns 2..4 of B, not 0,1,4
              for (int col=0; col<3; ++col) {
                Bh[r][col+2] -= g*c.rd(m,k,j,b1+RR+3*rc3[cc]+col);
              }
            }
          }
        }
        if (s > 0) {
          const int bm1 = (s-1)*nrd;
          for (int r=0; r<5; ++r) {
            for (int col=0; col<5; ++col) {
              Real sa = 0.0;
              for (int t=0; t<3; ++t) {
                sa += c.Wk<true>(m,HH+3*r+t,i1,k,j)*c.rd(m,k,j,bm1+GR+5*(t+2)+col);
              }
              Bh[r][col] -= sa;
            }
            Real s2 = 0.0;
            for (int t=0; t<3; ++t) {
              s2 += c.Wk<true>(m,HH+3*r+t,i1,k,j)*c.rd(m,k,j,bm1+DR+t+2);
            }
            rh[r] -= s2;
          }
        }
        Real Bi[5][5];
        if (!RTCol3Inv5(Bh, Bi)) {
          for (int r=0; r<5; ++r) {
            c.rd(m,k,j,b0+DR+r) = 0.0;
            for (int col=0; col<5; ++col) c.rd(m,k,j,b0+GR+5*r+col) = 0.0;
          }
          continue;
        }
        for (int r=0; r<5; ++r) {
          Real sd = 0.0;
          for (int t=0; t<5; ++t) sd += Bi[r][t]*rh[t];
          c.rd(m,k,j,b0+DR+r) = sd;
          for (int col=0; col<5; ++col) {
            Real g = 0.0;
            for (int t=0; t<5; ++t) g += Bi[r][t]*Ch[t][col];
            c.rd(m,k,j,b0+GR+5*r+col) = g;
          }
        }
      }
      // the reduced back-substitution, with the SAME clamp, propagated the same way
      Real ynext[5];
      for (int r=0; r<5; ++r) ynext[r] = 0.0;
      for (int s=nsg-1; s>=0; --s) {
        const int b0 = s*nrd;
        const int i1 = ic + (nc*(s+1))/nsg - 1;
        Real y[5];
        for (int r=0; r<5; ++r) {
          Real sy = c.rd(m,k,j,b0+DR+r);
          if (s + 1 < nsg) {
            for (int col=0; col<5; ++col) sy -= c.rd(m,k,j,b0+GR+5*r+col)*ynext[col];
          }
          y[r] = sy;
        }
        c.rd(m,k,j,b0+YR) = y[4];
        const Real b = c.Wk<true>(m,BB,i1,k,j);
        if (b > 0.0) {
          if (y[4] > 3.0*b) {
            y[4] = 3.0*b;
          } else if (y[4] < -0.75*b) {
            y[4] = -0.75*b;
          }
        }
        for (int r=0; r<5; ++r) {
          c.rd(m,k,j,b0+YY+r) = y[r];
          ynext[r] = y[r];
        }
      }
    });
    }
    tm.team_barrier();
    // ---- 4c'. the segments, in parallel ---------------------------------------------
    Real dbm = 0.0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, nsg),
    [&](const int s, Real &dmx) {
      const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
      const int b0 = s*nrd;
      int ncl = 0;
      Real yL[3] = {0.0, 0.0, 0.0};
      if (s > 0) {
        const int bm1 = (s-1)*nrd;
        for (int r=0; r<3; ++r) yL[r] = c.rd(m,k,j,bm1+YY+r+2);
      }
      Real ynext[5];
      for (int r=0; r<5; ++r) ynext[r] = c.rd(m,k,j,b0+YY+r);
      // the boundary cell itself: the clamp is re-applied to the RAW increment, which is
      // the same arithmetic the reduced pass did, so the two agree bit for bit
      {
        const Real b = c.Wk<true>(m,BB,i1,k,j);
        Real db = c.rd(m,k,j,b0+YR);
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
          c.Wk<true>(m,BB,i1,k,j) = b + db;
        }
      }
      for (int i=i1-1; i>=i0; --i) {
        if (i < ib) {                         // the DEEP cell: one row, one unknown
          Real sy = c.Wk<true>(m,DPD,i,k,j);
          for (int cc=0; cc<3; ++cc) {
            sy -= c.Wk<true>(m,GD+cc,i,k,j)*ynext[rc3[cc]]
                + c.Wk<true>(m,HHD+cc,i,k,j)*yL[cc];
          }
          const Real bd = c.Wk<true>(m,BB,i,k,j);
          Real dbd = sy;
          if (bd > 0.0) {
            if (dbd > 3.0*bd) {
              dbd = 3.0*bd;
              ++ncl;
            } else if (dbd < -0.75*bd) {
              dbd = -0.75*bd;
              ++ncl;
            }
            const Real rel = fabs(dbd)/bd;
            if (rel > dmx) dmx = rel;
            c.Wk<true>(m,BB,i,k,j) = bd + dbd;
          }
          for (int r=0; r<4; ++r) ynext[r] = 0.0;
          ynext[4] = dbd;
          continue;
        }
        Real y[5];
        for (int r=0; r<5; ++r) {
          Real sy = c.Wk<true>(m,DP+r,i,k,j);
          for (int cc=0; cc<3; ++cc) {
            sy -= c.Wk<true>(m,G0+3*r+cc,i,k,j)*ynext[rc3[cc]]
                + c.Wk<true>(m,HH+3*r+cc,i,k,j)*yL[cc];
          }
          y[r] = sy;
        }
        const Real b = c.Wk<true>(m,BB,i,k,j);
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
          c.Wk<true>(m,BB,i,k,j) = b + db;
        }
        y[4] = db;
        for (int r=0; r<5; ++r) ynext[r] = y[r];
      }
      c.rd(m,k,j,b0+NCL) = static_cast<Real>(ncl);
    }, Kokkos::Max<Real>(dbm));
    if (dbm > dbmax) dbmax = dbm;
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      Real ncl = 0.0;
      for (int s=0; s<nsg; ++s) ncl += c.rd(m,k,j,s*nrd+NCL);
      if (ncl > 0.0) Kokkos::atomic_add(&c.stat(3), ncl);
    });
    if (c.fixit) continue;
    if (c.dstop && dbm < c.tol) break;
    if (dbm < 1.0e-14) break;
  }

  // ---- 5. apply, and the column energy budget --------------------------------------
  tm.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, nsg), [&](const int s) {
    const int i0 = ic + (nc*s)/nsg, i1 = ic + (nc*(s+1))/nsg - 1;
    const int b0 = s*nrd;
    Real budget = 0.0, bscale = 0.0, rtmax = 0.0, ubmax = 0.0;
    Real rhsum = 0.0, srsum = 0.0;
    for (int i=i0; i<=i1; ++i) {
      const Real dxi = c.Dx(m,k,j,i);
      const Real wb = c.taublend
          ? (1.0 - 0.5*(c.wblend(m,k,j,i) + c.wblend(m,k,j,i+1))) : 1.0;
      rhsum += c.bdt*(wb*c.Wk<true>(m,SA,i,k,j) + c.Wk<true>(m,EX,i,k,j))*dxi;
      srsum += c.Wk<true>(m,SA,i,k,j)*dxi;
      RTCol3WarmStore(c, m, k, j, i, c.Wk<true>(m,BB,i,k,j));
      const Real tk = c.Tg(m,k,j,i);
      if (!(tk > 0.0)) continue;
      const Real rho = c.Rho(m,k,j,i);
      const Real es = c.Ei(m,k,j,i);
      const Real b = c.Wk<true>(m,BB,i,k,j);
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
      // problem/rt_src_theta: the same re-centring the serial solver makes
      if (c.theta != 1.0) {
        de = RTCol3ThetaDe(de, c.bdt*(wb*c.Src(m,0,i,k,j) + c.Wk<true>(m,EX,i,k,j)),
                           c.theta);
      }
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
    c.rd(m,k,j,b0+BUD) = budget;
    c.rd(m,k,j,b0+BSC) = bscale;
    c.rd(m,k,j,b0+RTM) = rtmax;
    c.rd(m,k,j,b0+UBM) = ubmax;
    c.rd(m,k,j,b0+RHS) = rhsum;
    c.rd(m,k,j,b0+SRS) = srsum;
  });
  tm.team_barrier();
  Kokkos::single(Kokkos::PerTeam(tm), [&]() {
    Real budget = 0.0, bscale = 0.0, rtmax = 0.0, ubmax = 0.0;
    Real rhsum = 0.0, srsum = 0.0;
    for (int s=0; s<nsg; ++s) {
      const int b0 = s*nrd;
      budget += c.rd(m,k,j,b0+BUD);
      bscale += c.rd(m,k,j,b0+BSC);
      rhsum += c.rd(m,k,j,b0+RHS);
      srsum += c.rd(m,k,j,b0+SRS);
      const Real rt = c.rd(m,k,j,b0+RTM);
      if (rt > rtmax) rtmax = rt;
      const Real ub = c.rd(m,k,j,b0+UBM);
      if (ub > ubmax) ubmax = ub;
    }
    Real ftop3 = 0.0;
    for (int q=0; q<nq; ++q) ftop3 += c.wf[q]*(c.Wk<true>(m,UU+q,ie,k,j) - Dtop[q]);
    Real fnet = 0.0;
    if (hb.on) {
      fnet = c.Wk<true>(m,FL,ic,k,j) - ftop3;   // the imposed deep bottom flux
    } else {
      for (int q=0; q<nq; ++q) {
        fnet += c.wf[q]*((c.Wk<true>(m,BB,ic,k,j) + Ucut[q])
                         - c.Wk<true>(m,DD+q,ic,k,j));
      }
      fnet -= ftop3;
    }
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
  // problem/rt_col3_skip_sweep: the entry sweep did not run, so write this solve's own
  // converged face flux into Fb.  Same definition as the serial path, one thread per
  // face; the workspace is global, so reading i-1 across a segment boundary is legal
  // here (the whole column has converged and the barrier above has passed).
  if (c.wrflux) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
      if (i < ib) {                           // a DEEP face carries the diffusion flux
        c.Fb(m,0,i,k,j) = c.Wk<true>(m,FL,i,k,j);
        return;
      }
      Real f3lo = 0.0;
      for (int q=0; q<nq; ++q) {
        const Real ulo = (i == ib) ? (c.Wk<true>(m,BB,ib,k,j) + Ucut_i[q])
                                   : c.Wk<true>(m,UU+q,i-1,k,j);
        f3lo += c.wf[q]*(ulo - c.Wk<true>(m,DD+q,i,k,j));
      }
      c.Fb(m,0,i,k,j) = f3lo;
      if (i == ie) {
        Real f3hi = 0.0;
        for (int q=0; q<nq; ++q) {
          f3hi += c.wf[q]*(c.Wk<true>(m,UU+q,ie,k,j) - Dtop[q]);
        }
        c.Fb(m,0,ie+1,k,j) = f3hi;
      }
    });
    tm.team_barrier();
  }
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
