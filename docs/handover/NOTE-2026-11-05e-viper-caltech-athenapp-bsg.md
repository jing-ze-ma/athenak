# NOTE 2026-11-05e viper -> Caltech: YES, Athena++ (Linhao's bsg_star.cpp) was run, in 1-D on Raven 10-15/16 (reply to TASK f0bd5acef)

Everything is in `docs/handover/athenapp-bsg/` (README.md first; the worker NOTEs in full; decks, patches, figures, summaries; 14 MB).
Raw run dirs: Raven `/raven/ptmp/jinma/athpp_audit_1015`, `/raven/ptmp/jinma/athpp_fix_1016`. They are kept in full: the planned rst cleanup excludes them.

## Short answer
- **What was run:** Athena++ e47fabb7 + his pgen, 1-D, 20-80 Rsun, nx1 384, nmu 4.
- **Every arm collapses at the inner wall:**
  - the wall cell cools 1e6 -> 2-3e5 K;
  - flux -100..-400 L into the bath;
  - the envelope deflates 54 -> 45 Rsun by 17 d, then blows out at ~18 d.
- **Cause, 1:** HLL tau-factor diffusion (3/(8 taucell) of L without a force).
- **Cause, 2:** the unconverged Jacobi, whose billing does not telescope. A converged line-solve fix removes that loss, but the converged discretisation still goes unstable at the wall.
- **AthenaK, same IC/BC:** 30 d stable (h2h).
- Caltech's plan (3) is therefore partly done already, for the taucell/angle_flag variants. The r^+2 vs inverse-square atmosphere was NOT varied there.
