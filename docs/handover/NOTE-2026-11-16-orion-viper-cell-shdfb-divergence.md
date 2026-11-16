# NOTE 2026-11-16 Orion -> viper: cell layout Newton "divergence" on shd x0.5 with gas_feedback (finding, no fix)

This was found during the local staggered study (NOTE 438626ac0). Orion did NOT change the cell layout. The
reproducer is in handover_1116_shdfb/ (README, deck, log excerpt, debug hits).

- Deck: cc shd_FIN3_r12g2 at x0.5 with `gas_feedback = true` as the only change; clean defaults; binary d6b6ed4e2
  (= push 14). Cycle 0 BE: `pic 1.46e+295 ... solves 250` -> CLEAN WARNING + dt_retry twice. It is then accepted with
  neg 96 nreal 1722, and every later cycle prints pic ~1e287-1e291 with 180-250 solves per step.
- First bad quantity: rad_m1_clean.cpp CleanNewtonUpdate, which the cell layout calls from CleanCellNewtonUpdate
  only with gas coupling:
  `de = fabs(e - e_old)/fmax(fmax(e, e_old), 1.0e-300)`. When BOTH iterates are negative (the cell layout keeps E < 0
  for its limiter), the denominator is 1e-300 -> de ~ 1e299, so the loop never settles.
  - The ~16 cells are on the absorber rim (i 145-150, j 55-56 / 74-76).
  - Suggested one-liner: use max(|e|, |e_old|) in the denominator. Orion did this on the staggered path only
    (local commit 8e47cdc55).
- Underlying physics: the gas-coupled cell solve drives E to about -0.3 (against ~1) on the absorber rim. Without
  gas the same deck also ends 8/8 stages with neg > 0. So fixing the denominator removes the stall but not the
  negatives.
