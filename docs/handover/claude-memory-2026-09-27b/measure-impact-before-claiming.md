---
name: measure-impact-before-claiming
description: "State a bug's IMPACT only after measuring it. Twice on this project I reported consequences that a cheap A/B then refuted."
metadata:
  node_type: memory
  type: feedback
---

Finding a defect by reading code is not the same as knowing what it does. Report the
defect immediately; report its **consequences** only once an A/B has shown them.

**Why.** This failure mode has now cost real work twice here:

* The whole dhj blow-up campaign ([[dhj-ck-eos-blowup]]) -- a dozen arms compared against
  each other, and a set of physical conclusions drawn, on a binary carrying a data race.
  Every comparison had to be retracted.
* [[mhd-fluxes-stale-x1-limits]] -- I found a genuinely wrong index range and told the
  user it caused a wrong flux and a live race. A build-and-compare (2 binaries, 5
  replicates, 600 cycles, ~1 h wall) showed it changes **no answer at all**: the array in
  question is zero at the affected ghost cell, and the predicted race never tripped.

**How to apply.** On this project a decisive A/B is cheap -- build two binaries in
`/viper/ptmp/jinma/claude_eos_gpu` and diff the history files bitwise -- so there is
rarely an excuse to skip it. Concretely:

1. Say what is wrong in the code (safe to state from reading).
2. Say "impact not yet measured" rather than guessing, or run the A/B first.
3. Check that the values involved are non-trivial before predicting an effect -- a range
   bug over cells where the data is identically zero is a no-op.
4. Confirm the path is even live (`flag=on` vs `off` must differ) before concluding from a
   null result; otherwise a dead branch reads as a clean bill of health.
5. If a regression test would pass with AND without the fix, do not add it -- say why it
   is missing instead. Same principle the user applies to experiments in
   [[confounded-tests-rejected]].
