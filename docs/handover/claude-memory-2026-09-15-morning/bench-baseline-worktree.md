---
name: bench-baseline-worktree
description: A reusable pre-change CPU baseline binary lives in a git worktree at /viper/u2/jinma/ATHENAK/bench/base_wt; the other bench binaries are GPU builds that cannot run on the login node
metadata:
  type: reference
---

**`/viper/u2/jinma/ATHENAK/bench/base_wt`** is a git WORKTREE parked at a pre-change HEAD
with its own build (`build/src/athena`).  It reproduces every recorded baseline
digit-for-digit, so an A/B does not need a rebuild of the old side.

* The other `athena` binaries under `bench/` are **GPU builds and CANNOT run on the login
  node** -- they fail in a way that is easy to misread as a code problem.
* A fresh worktree needs `kokkos` symlinked in before it will configure.

**Why:** re-creating a baseline build by hand costs ~15 min and risks a silently
different configure, which invalidates the comparison.
**How to apply:** for any A/B against a committed baseline, run the old side out of
`base_wt` rather than checking out and rebuilding in the working tree.
Cf. [[validate-the-instrument]], [[measure-impact-before-claiming]].
