---
name: caltech-port-hold-merge
description: DONE 2026-09-25 -- caltech-port merged into rt-integration (39ca34a6) and pushed (83876228); viper must run git submodule update --init once + HIP bitwise gate
metadata:
  node_type: memory
  type: project
  originSessionId: 1f4629f7-1667-4eeb-be90-3eb8315950c0
  modified: 2026-09-26T00:19:23.415Z
---

Branch `caltech-port` (from rt-integration 265ca539; first commit 2c2368d5) holds: kokkos restored as a submodule (4.6.02 08ceff92; 3596f5b3 had made it a symlink to a viper path), nvcc fixes (bvals_cc.cpp generic lambda, captured helper lambdas in conduction.cpp and others), Caltech build/job/smoke scripts.

**MERGED + PUSHED 2026-09-25 22:07 PDT** (user: "merge it here ... verify"; "push it"): merge 39ca34a6 on top of viper 3f13cd38, handover note 83876228; gates: CPU bitwise vs 3f13cd38, GPU vs CPU t=0 identical, 2 vs 1 H200 bitwise. Branch caltech-port deleted. Earlier: **Hold the merge** (user 2026-09-25): viper may still have tasks to push to rt-integration first. Overrides the standing "merge gated branches promptly" rule ([[merge-finished-branches]]) for this branch.

**Why:** avoid conflicts with pending viper work; the kokkos change also alters how viper checkouts look (needs `git submodule update --init` once).
**How to apply:** keep gating caltech-port (diff review + CPU bitwise vs the unported binary); when the user clears it, rebase/merge onto the latest rt-integration, re-gate, merge, push, and add a handover note about the one-time submodule update on viper.
