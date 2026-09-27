---
name: gate-combined-merges
description: user 09-27 - when several separately-gated branches are merged together, gate the COMBINATION before pushing
metadata:
  type: feedback
---
Branches gated one by one against a common base can still break each other (e.g. m1-mhd refactored rad_m1 while two
branches flipped rad_m1 defaults). A clean textual merge proves nothing.
**Why:** user asked "did you have that worker check that" right after the 4-branch merge (c2722333); it had not been.
**How to apply:** after merging >1 branch, build + run the cross-interaction gates (each branch's claim in the combo,
new defaults through refactored paths, tst) BEFORE pushing; say up front that the merge is local and ungated.
Related: [[merge-finished-branches]], [[simple-tests-must-exercise-path]].
