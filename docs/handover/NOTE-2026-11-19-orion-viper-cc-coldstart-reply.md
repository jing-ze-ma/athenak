# NOTE 2026-11-19 Orion -> viper (USER): reply to f716a219a (clean-cc cold start at dt ~1000, centred advective term)

Thanks; your diagnosis matches what we see on Orion's staggered (stag2) lines. Nothing pushed on our side either (local branches).

1. **The staggered F-solve line (K1) already upwinds the advective flux.** In rad_m1.cpp `cln_upw = cln_split`: with the staggered
   F solve, the donor-cell share of the advective flux (key set `upwind`) is on by default, everywhere (not only above Peclet 2).
   K1 + old prod's solver settings (local positivity fallback, implicit_maxit 200, implicit_lin_maxit 200 with no stall abort,
   implicit_tol_rel_max 1e-4, stage acceptance with the conservative write-back floor) holds prod's dt on our evolved restarts:
   HG hegiant.00184 767342-767361 at dt 97.5 and RSG rsg.01510 684504-684520 at dt 888.9, 0 CLEAN WARNING, 0 retries (Raven CPU,
   24 ranks). The upwinding is likely part of why.
2. **Our F-eliminated line (staggered, prod face, no F solve) does NOT upwind it** (`cln_upw` follows the F solve), and it fails RSG
   the same way: E<0 cells driven by the cell row's kick/advective terms, a Picard limit cycle at dt 888.9, BE redo not converging.
   So the same cause shows up on a different layout. We will test the donor-cell upwinding on that line and report.
3. Difference worth comparing: ours is a donor-cell share of the advective flux on every face (the diffusion part unchanged);
   yours is a local Rusanov dissipation >= |a_n| only where Peclet > 2, which keeps the centred (2nd-order) form elsewhere.
   If your linear-wave / IDORT order decks show the Peclet switch costs nothing, it may be the better form for us too.
4. Your second item (thin outer zone, i 428, h c chi ~38): on the staggered lines the remaining positivity issues we found are
   (a) the off-diagonal Eddington cross terms (not an M-matrix; we test a donor/lagged form), (b) the explicit stage increment
   of hesdirk2 making the stage's old vector negative in thin cells. Happy to compare notes when you get there.

-- Orion
