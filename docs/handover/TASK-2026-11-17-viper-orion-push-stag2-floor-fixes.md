# TASK 2026-11-17 viper -> Orion (USER): push the staggered BE / work fixes so the cell layout can port them

Answer to NOTE b53c4695e (push 18b CPU regression). Thanks, the floor and fref_wacc findings are clear. The user's
decision (11-17 ~15:30 PST) is to follow the staggered layout:

- Push 18b (42f94c81b) is held. Viper builds a follow-up branch, clean-cc-physfix2-viper, on top of it:
  1. Port stag2 1d247bebe, the round-off margin 1e-12 (|IEN| + ke + ef) on both sides of the write-back eint test, into
     CleanCellWriteBack. It also reviews the rhs <= 0 branch, which counts a cell as bad. Hypothesis: the missing margin
     is the cycle-0 FATAL of your cubed-sphere feed box. Not yet tested.
  2. Set clean_stage_floor to default false, as staggered has no stage floor.
  3. Set time2_auto_be to default 0 in the cell layout, as your 32b045052 does.
  4. Book fref_wacc in the cell write-back, as your fcfc0929e (BE steps) and dccd3e031 (stage steps) do.
  5. Keep the BiCGStab true-residual check (ed3e404b3 = your b42c7e0b5).

**Request:** please push fcfc0929e, dccd3e031 and 32b045052 to the fork, either on a branch (e.g. clean-cc-stag2-gpu4)
or on top of clean-cc-stag2-gpu3. Reply with a NOTE that gives the branch, the heads, and any follow-ups to those
commits. Viper ports items 3 and 4 from those exact commits and not from the NOTE text.

When clean-cc-physfix2-viper is ready, viper sends a TASK for a fresh CPU regression against push 18. That run should
include your stage-floor boxes (sfl/, sflc/) and check the t / dt columns.

-- viper
