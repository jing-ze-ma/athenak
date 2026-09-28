---
name: simple-tests-must-exercise-path
description: user 09-27 asked whether ke-dt was checked on simple tests; LE/radwave "identical" only because keys are no-ops without wb_arad -- not a test
metadata:
  type: feedback
---
Every new scheme needs simple tests (analytic or converged reference) that actually run the changed code path.
"Identical with keys on" on a test where the keys are no-ops proves only non-regression, not correctness.
**Why:** user asked "did you check that scheme on simple tests?" before the ke-dt-0926 default flip (09-27); the
answer was no (only He box + sph wedge exercised force_reference_work=split).
**How to apply:** before recommending a default flip or merge, list which simple tests hit the new path; if none,
design them (e.g. rad-supported static slab, radiation-acoustic wave vs linear theory, cfl sweep). See [[benchmark-on-radhydro]].
