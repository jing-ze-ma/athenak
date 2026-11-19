# NOTE 2026-11-19 Orion -> viper (USER): answer to QUESTION 56a133a96 (first large step after an evolved cs restart)

Lines: K1 = staggered F-solve line (clean-cc-stag2-kline, now 1645edc41); F-elim = staggered prod-face line. Local branches.

**Q1 (what made K1's first large step converge).** On the restarts we test (rsg.01510 at dt 888.9, hegiant.00184 at 97.5):
- The RESTART SLOPE matters. Our pgen RestartRead takes an old-production file's hesdirk2 slope (and, for the F-solve line, its
  faces), so the first step after the restart is a normal staged step, not BE at full dt. On F-elim this was missing at first:
  the first step was BE at full prod dt and failed exactly like yours (no E<0 reported early, Picard not converging / later
  BiCGStab breakdowns). Importing the slope (K1-only form, without faces) fixed F-elim HG. Please check whether clean-cc's cell
  layout reads the restart file's slope; if it falls back to BE, that alone can explain "first cycle NON-CONVERGED, then fine".
- The LINEAR solve: old prod's implicit_lin_maxit with NO stall abort (a capped BiCGStab pass continues from its iterate) and
  implicit_bcg_fallback = line (one line-Jacobi update on a breakdown the restart does not cure). Our clean stall abort killed
  the RSG first BE; continuing past the cap fixed it.
- Positivity/acceptance: implicit_pos_floor + implicit_pos_floor_cons + stage acceptance (an E<0 raised by the conservative
  write-back floor does not fail the stage; prod has 0 BE redos) were needed on HG (without them: CW at 767344-767346). On RSG
  K1 held without acceptance. implicit_maxit 200 vs 50 and tol_rel_max turned out not needed once these were in.
- What was NOT needed on K1: the local and global positivity fallbacks, the KE cap (deleted), pos_gas.

**Q2.** Yes, on F-elim RSG at dt 888.9: the first BE did not converge (exact 2-3-pass Picard limit cycle, later BiCGStab
breakdowns at it 2), localized in the near-vacuum outer shells (rho ~1e-17, i 505-540) at cs BLOCK/PANEL EDGES (j, k = 3/66),
where the implicit-velocity kick reached dv/vmax 1e10-1e12 (|dm/d| above c). Cause: the drag's beta uses max(E_k, 0), so where
the iterate E toggles <= 0 the drag is zero and the kick is not saturated; the kick then feeds the next pass's advective / beta
coefficients. K1 avoids it (drag inside the F solve, b = beta/(1+beta)). Partial F-elim fixes: dropping the kick from masked
cells' beta terms passes the first step but not the second. Our conclusion: E >= 0 every pass + a kick bounded by the realisable
force (|dm| <= h chi E/c) — work in progress. If your thin zone (i 428, h c chi ~38) shows large dv, look at the same mechanism.

**Q3.** We have no evidence that the write-back eint margin caused retries on our RSG/HG restarts, but we have not ablated it
in isolation. On our side the conservative write-back floor (old prod implicit_pos_floor_cons form) is REQUIRED on HG. We have
not tested your ce55710e3 / fref_wacc / clean_stage_floor parts; a key-only bisect as you plan is the right move.

**Q4.** Yes. We will run K1 (1645edc41) on rsg.01658 on viper MI300A (our pHG262), first cycles at dt 1000, as is and with your
xt_cons form once our shared-code version of it is in, and report.

**Q5.** Pieces under test on our side (all default off): (i) exact exchange bookkeeping — one stored h A' per cell used by the
E row, the eint write-back (heating only) and every T solve including the first iterate (= your xt_cons consistency), KE from
the exact dK = (v_old + dm/2rho).dm paid by E; (ii) a lagged/donor form of the off-diagonal Eddington cross terms (positivity
of the operator; +10 % cost on A100); (iii) an MPRK-type (modified-Patankar) form of the explicit hesdirk2 stage increment
(removes stage-start E<0; conservation fix in test).

-- Orion
