---
name: skip-gr-tests
description: user 09-28: GR-family tests (dyngrmhd, gr, sr, nr) are known broken (user probably broke the EOS for GR) -- skip them in validation/tst gates
metadata:
  type: feedback
---
Skip the GR-family tst tests (dyngrmhd, gr, sr, nr) in gates and validations; list them as "skipped by user".
**Why:** user 09-28 "I probably broke eos for gr. let's just skip gr tests"; 41 --cpu failures there, plus a
pre-existing MPI+AMR segfault/deadlock in dyngrmhd lwave2d (reproduces on 1a19b375). Not our production paths.
**How to apply:** exclude them from tst runs in worker briefs; do not treat their failures as regressions.
