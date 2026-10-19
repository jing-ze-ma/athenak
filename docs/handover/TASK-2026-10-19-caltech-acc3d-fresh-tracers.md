# TASK 2026-10-19 raven -> caltech: 3-D accretor FRESH run with tracers on H200 (USER: "move the fresh accretor to caltech")

USER 10-19 (~21:00 CEST): start the 3-D accretor FRESH (t=0) with the newest accretor binary + passive tracers, on Caltech.
Raven's fresh copy (run12) sits ~17 h in the queue. This supersedes all earlier accretor TASKs (v10/v11/v11am): if any
Caltech acc3d chain is still queued or running, scancel it and KEEP its run dir/rst/dumps. Reply with t reached + GPU-h.

What's new vs v11 (all gated on Raven A100, BITWISE where stated):
- accwall-1019 4d2d5628: `problem/theta_wall_closed = true`. The hse theta walls are now exactly closed (in v11 the walls leaked ~1.8x the stream into the
  star). Plus the Ebnd/Eirr booking fix (energy budget closes to 1e-3 of |Ebnd|, was 23 %).
- acctrc 4960ac49: `<hydro>/nscalars = 2`: s0 = stream gas (nozzle), s1 = boundary inflow; diagnostics in ryper.trc.hst;
  hydro_w bins gain s_00, s_01 (+40 % disk per dump).
- am_transport = true (exact L_z, since v11am).

1. Build: origin/acctrc-1019b @ 4960ac49 (= ab758e94 + accwall 4d2d5628 + tracer commit). H200, your nofma recipe as for v11. Report md5.
2. Deck: docs/handover/acc3d-prod-1018/acc3d_prod_v11amwtr.athinput.in (md5 of the Raven original 9fca5608). Replace @FILES@ with
   your copy of this directory; the column + both tables are unchanged (5b813886 / d81b7a46 / 99d668f6). Same layout as your v11 (2 H200, 8 MeshBlocks).
3. Quick check on H200 (USER rule), as for v11:
   - static phase to t ~0.002: C and C2 must match bitwise;
   - restart C from a mid rst (R): R must match C;
   - compare hst every cycle, incl. ryper.trc.hst + final bins.
   If DIFF: report, don't start.
4. Production: your acc3d_chain.sh, FRESH from t=0 (static phase to t 0.05 with problem/stream=false, then the deck as is). Use the longest
   links you have, as many links as v11 had queued. STOP only on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc.
5. Report (NOTE caltech->raven):
   - md5, check results, link ids + start estimate;
   - after the stream starts: s/cycle, dt, retries, ryper.trc.hst s0 mass vs the nozzle rate (254.84 per code time);
   - the hst + every 4th hydro_w dump path, so Raven can rsync them for the page.
