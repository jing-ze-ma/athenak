---
name: caltech-testing-only
description: User 2026-09-25 -- the Caltech cluster is for TESTING, not production; production stays elsewhere (viper) unless really needed
metadata:
  type: project
---

On 2026-09-25 the user said: "we don't need to really run production here unless really needed. for now it's for testing".

- Caltech jobs are short tests: builds, smoke tests, gates, and short timing and A/B runs.
- The GPU budget in [[carnegie-hpc-usage-policy]] is therefore not a production planning concern.
- Keep tests small, and still ask before long or multi-GPU jobs.
- Use scratch (purged after 14 d) for tests. `/resnick/groups/carnegie_poc/jingze/` is for anything worth keeping.

**Why:** production physics runs stay on viper; Caltech is mainly for porting and testing.
**How to apply:** do not propose production campaigns or GPU-hour budgets here unless the user asks.
