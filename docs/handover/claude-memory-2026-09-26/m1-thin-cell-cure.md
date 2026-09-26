---
name: m1-thin-cell-cure
description: The multi-D thin-cell instability of implicit M1 chi(f) closures (and the 09-22 He-slab thin-top damage) is a lagged-closure scheme effect; cure = <rad_m1>/implicit_closure_thin_relax = 1.5 (merged b7e48012, default off); the 09-24 "lp 0.5 + offdiag none" choice does NOT cure it
metadata:
  type: project
---

tests_m1/runs_5c_thinstab/README.md. Growth |g| ~ max(chi', b/f)/tau_cell per step, dt-independent; decays for
tau_cell >= 0.25; Eddington stable, M1/Kershaw unstable; 1-D stable. m1-thin hll and trans_limit lp + offdiag none do
NOT cure it; closure_lag = pass works only to tau_cell 0.125 at 5-7x passes. implicit_closure_thin_relax = 1.5: seeds
decay over 500 steps at tau_cell 0.006-0.5 under be AND hesdirk2, 0 extra passes; He slab 09-22 input: 332 fallbacks ->
0. Limits: slow adaptation at tau_cell ~0.006 (3 %), f > 0.69 in thin cells not covered, first step after restart
unrelaxed. Fixed-tensor closures (vet_sc, vet_col, Eddington) don't need it.

**Why:** supersedes the m1-thin branch decision (keep unmerged; its thin-top cure is refuted).

**How to apply:** for M1/Kershaw closures with thin regions, name implicit_closure_thin_relax = 1.5; prefer vet_col /
vet_sc in thin tops. Default-on is a user decision. Related: [[vet-col-enough-no-sc-on-sp]].
