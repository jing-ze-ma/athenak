---
name: eos-inversion-nan-trap
description: The EOS-table NaN rows silently corrupt a (rho,e)->T root find; always validate a post-processing inversion against ck_dump_file
metadata:
  type: feedback
---

Post-processing that inverts the tabulated EOS for `T` and `p` from a binary dump must
NOT bracket the root with `(E < target).sum(axis=0) - 1`. **The table from
[[eos-table-dump]] is NaN for `log10 T < ~1.85`** (about 35 rows at `dy = 0.01`, and the
boundary is density-dependent). `NaN < t` evaluates False, so those rows are skipped in the
count and every index comes out ~35 rows low -- **T too low by exactly 0.35 dex (2.24x), and
p too low by ~3.5x.**

**Why:** the bug is silent. Nothing raises, no NaN propagates into the answer, and the
result is smooth, monotonic and physically plausible. It survived a "validation" against
the ideal-gas estimate `rho = p*mu*m_H/(kT)` because that check used the same wrong T and
mu, so it was circular.

**How to apply:** mask before counting --
`Efill = np.where(np.isfinite(Ec), Ec, -np.inf)`, then clamp the index to the first finite
row per column. And **validate against `problem/ck_dump_file`**, which writes one column's
`p[bar]` and `T[K]` straight from AthenaK. Post-fix that agrees to 4-5 digits and the
implied gravity comes out 943.7 against the input 942. **Any pressure or temperature in this
campaign that predates 2026-08-25 and was not checked against `ck_dump_file` is suspect.**

**What it did and did not change (2026-08-25, eighth session).** The x1max decision in
[[ck-limb-run]] SURVIVED and got stronger: on the retired grid 68.3 % of terminator columns
had the 1e-6 bar isobar outside the domain, not the 43 % I first reported, and both limbs
were outside. `f_stretch_theta` and the `dfloor` A/B were unaffected -- both came from cycle
counts and dt, not the EOS. **Withdrawn entirely:** a claimed factor-3.6 discrepancy between
the realized gravity and `<problem>/grav`, which was this bug; and the day/night
`H_observed/H_hydrostatic` ratios. See [[dhj-isobar-vs-shell]].
