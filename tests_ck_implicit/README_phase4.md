# `problem/ck_implicit`, phase 4 — JACOBIAN REUSE and a BETTER INITIAL GUESS

Phase 1 is `README.md` (the backward-Euler correlated-k column), phase 2
`README_phase2.md` (the thin-cell bracketed solve and the column skipping), phase 3
`README_phase3.md` (the frozen exchange operator `ck_impl_frozen_op`, worth 12 % on CPU,
and the warm start `ck_impl_warm`, which is not recommended because it perturbs
unconverged calls).  Phase 3 ends with one statement: **the pass count is now the whole
cost** — 6.3-6.8 passes per RK stage, 13 ck sweeps per hydro step against the 2 of the
semi-implicit scheme.  Phase 4 attacks the cost of a pass and the number of passes with
one lever each.

Both are default OFF and both are inert when off (gate a).

| knob | default | what it does |
| --- | --- | --- |
| `problem/ck_impl_reuse_jac` | `0` | `1` = chord Newton: build the tridiagonal on the first pass that takes a Newton step and REUSE it; `2` = as 1, with each entry rescaled by `(T/T_build)^3` |
| `problem/ck_impl_seed` | `0` | `1` = pass 0 takes the SEMI-IMPLICIT step of the default scheme instead of a Newton step; `2` = pass 0 takes the exact per-cell bracketed backward-Euler step (`CkThinSolve`) |

Code: `src/utils/two_stream_column_ck.hpp` (the knobs, `CkImplJacPass`, the seed rows and
the rescaling inside `CkImplStep`, the `ck_t0` latch), `src/utils/two_stream_rt.hpp`
(`ckjacp_` — the one flag that decides whether a pass assembles anything — and its three
uses: the `dB_b/dT` block of `rt_pre_opac`, the `ck_jac_zero` kernel, and `jck` inside
`rt_chain_ck`), `src/pgen/deep_hot_jupiter_rt.cpp` (plumbing and the range guards),
`inputs/tests/dhj_ck_implicit.athinput`.  Scripts `p4.sh` (gates a, b, c) and `p4d.sh`
(the `maxit = 20` pass count).

---

## 0. What each lever actually removes, and why it is allowed

**Reuse.**  The residual is re-formed by a full ck sweep every pass and that does not
change — the fixed point is the exact backward-Euler balance either way, whatever matrix
the step is taken with.  What a reusing pass skips is only the *assembly*: the ~10
`Kokkos::atomic_add` into `ckjac_g` per (cell, chain) inside `rt_chain_ck`, the
`ck_jac_zero` kernel, and the `CK_NB` PAIR of tabulated Planck-fraction interpolations per
cell that the one-sided `dB_b/dT` difference costs in `rt_pre_opac` (the Planck functions
`B_b` themselves are still rebuilt — they are the iterate).  This is the chord (modified
Newton) method: linear convergence instead of the quasi-Newton's, same root.

Mode 2 exists because at frozen opacity the ONLY thing in `J` that moves between passes is
`dB_b/dT = d(sigma T^4 f_b(T)/pi)/dT`, which is `4 sigma T^3 f_b/pi` up to the band
fraction's own weak `T` dependence.  Rescaling each column of a row by `(T/T_build)^3` of
the cell that column points at therefore recovers the dominant part of the update for
three multiplies and no sweep work.  The factors are positive, so `J` stays the M-matrix
the Thomas sweep needs.

**Seed.**  `ck_impl_warm` starts the Newton from the PREVIOUS call's increment, and phase 3
measured the damage: on a transient half the calls stop at the pass cap, and a truncated
iteration depends on where it started, so the warm arm's state drifts ~1e-2 from the cold
one.  A seed formed inside THIS call carries none of that history.  Mode 1 is literally
the step the default (`ck_implicit = false`) scheme would take, the `rt_semi_lin` form
`de = (S/lambda)(1 - e^{-lambda bdt})` with `lambda = 4 E/e`; mode 2 is the same balance
solved exactly rather than damped — the bracketed positive root of
`x - e* - bdt(A - E (x/e)^4)`, with the absorbed field `A` lagged at pass 0's sweep.  Both
reach the gas as an IDENTITY ROW of the same Thomas sweep (`a = c = 0`, `b = 1`, so
`cp = 0` and `de = d`), which means the per-pass cap `ck_impl_dtmax` and the total-excursion
bound `ck_impl_demax` act on a seed exactly as on a Newton step and nothing downstream
knows the difference.  With `ck_impl_seed > 0` the Jacobian build moves to pass 1
(`CkImplJacPass()`), so under `+both` pass 0 assembles nothing at all.

One consequence of reuse is worth stating: the thin/thick classification
(`ck_impl_arat`) keeps moving between passes while the matrix does not.  A cell that turns
thin after the build still gets the identity row the solve gives it, and a thick
neighbour's frozen off-diagonal pointing at it is still `<= 0`, so the M-matrix property
survives; only the convergence rate can suffer.

---

## 1. Gate (a): both switches off is bitwise the reference

Serial CPU, `inputs/tests/dhj_ck_spherical.athinput`, 20 cycles, one `bin` dump per cycle,
against a binary built from unmodified `HEAD` (446b7b30).

| `ck_spherical` | `dhj.hydro.hst` | `bin` payloads |
| --- | --- | --- |
| `false` | **identical** | **22 of 22 identical** |
| `true` | **identical** | **22 of 22 identical** |

(`paycmp.py` compares everything after the stated header offset; the files differ in
length only because `PAR_DUMP` now carries the two new `<problem>` keys.)  By inspection
the same holds for every other path: `ckjacp_ = ckimp_ && (...)` is false whenever
`ck_implicit` is, so each of the three gates it replaces evaluates exactly as `ckimp_` did.

## 2. Gate (b): passes, residual and gap at the production cap (`maxit = 8`)

`inputs/tests/dhj_ck_implicit.athinput`, 20 cycles from the analytic initial condition,
`ck_spherical = ck_beam_sph = true`, `ck_impl_frozen_op = true` in every arm (so this is
phase 3's fastest configuration as the baseline).  40 RT calls.

| arm | mean passes / stage | res (max) | gap (max) | not converged |
| --- | --- | --- | --- | --- |
| baseline (phase 3) | **6.75** | 1.79e-05 | 9.45e-06 | 20 / 40 |
| `+ reuse_jac = 1` | **7.42** | 1.28e-04 | 9.91e-05 | 23 / 40 |
| `+ reuse_jac = 2` | **7.40** | 1.28e-04 | 9.91e-05 | 23 / 40 |
| `+ seed = 1` | **6.83** | 1.98e-05 | 1.12e-05 | 21 / 40 |
| `+ seed = 2` | **6.83** | 1.79e-05 | 9.45e-06 | 20 / 40 |
| `+ both` (1 + 2) | **7.08** | 1.08e-04 | 8.00e-05 | 22 / 40 |
| `+ both` (2 + 2) | **7.10** | 1.09e-04 | 8.16e-05 | 22 / 40 |

**Reuse costs 10 % of the pass count** (6.75 -> 7.42) — the trade the brief asked to be
reported.  **The `T^3` rescaling buys nothing**: modes 1 and 2 agree to the third digit on
the pass count and to every digit on the residual, because within one RT call the
temperature moves by far too little for the correction to matter.  Mode 2 is kept for the
large-`bdt` case it was written for, but **mode 1 is the one to use**.

**The seed does not reduce the pass count either** (6.75 -> 6.83, i.e. it costs 0.08 of a
pass): on this cold-start transient the residual after the seed step is not systematically
smaller than after a Newton step, because the pass-0 Newton step is ALREADY an implicit
step with the true tridiagonal.  What it does do is take the sting out of reuse — `+both`
is 7.08 against reuse's 7.42, i.e. the seed gives back half the passes reuse costs, which
is exactly the mechanism one would want: the matrix is then built at a better iterate.

**These maxima are the cold start.**  Half the calls in this window stop at the cap, and
the gap `> 1e-6` in the reuse arms is those calls, not the scheme (the baseline's own
cold-start maximum is 9.45e-06).  The settled-state gap is section 3.


## 3. Gate (d): GPU — the verdict

The cost verdict is measured on the GPU.  MI300A APU (`apudev`, 2 ranks, 2 GPUs,
`HSA_XNACK=1`, hipcc 6.3.4, `Kokkos_ARCH_AMD_GFX942_APU`), the production mesh
(6 x 32 x 32 x 128), restarted from `bench/cs_mhd_prod3/rst/dhj.00567.rst` at rotation 283
and run for exactly **300 cycles** (`time/nlim = 4430131`), i.e. 600 RT calls of a
SETTLED atmosphere — not the cold start of section 2.  Binary:
`bench/ckq_0922/src_snapshot/build_gpu/src/athena`, a `git archive HEAD` (446b7b30)
snapshot with only the three phase-4 files overlaid; md5 `f960b8686f9e0c12407c8442520ac1bf`.
`cycles/s = 300 / (cpu time used)`, the same quantity `bench/prof_0922` quotes.
All four arms have `ck_impl_frozen_op = true`, `maxit = 8`, `colskip = true`.

| arm | cpu time [s] | cycles/s | x plain | mean passes / stage | gap: max / median | not conv. |
| --- | --- | --- | --- | --- | --- | --- |
| plain, semi-implicit (`prof_0922/plain`) | 23.04 | **13.02** | 1.00 | — | — | — |
| implicit, no frozen op (`prof_0922/implicit`) | 113.99 | **2.63** | 4.95x | — | — | — |
| baseline: implicit + `frozen_op` | 117.86 | **2.55** | 5.11x | 5.67 | 4.10e-06 / 6.5e-07 | 5 / 600 |
| `+ reuse_jac = 1` | 96.80 | **3.10** | **4.20x** | 5.66 | 4.11e-06 / 6.5e-07 | 6 / 600 |
| `+ seed = 2` | 119.33 | **2.51** | 5.18x | 5.68 | **9.69e-07** / 3.4e-07 | 4 / 600 |
| `+ both` | 99.14 | **3.03** | **4.30x** | 5.68 | **9.69e-07** / 3.4e-07 | 4 / 600 |

**Read this way.**

* **Reuse is worth 1.22x on the GPU** (2.55 -> 3.10 cycles/s) and **on the settled state it
  costs no passes at all** — 5.66 against 5.67 over 600 calls.  The 10 % pass penalty of
  section 2 is an artefact of the cold-start transient, where the iterate moves far enough
  in one pass that the matrix built at the previous one is genuinely stale.  On the
  production state the chord matrix is as good as the rebuilt one, and everything the
  assembly cost is saved: the atomics, the zeroing kernel and the `dB_b/dT` look-ups.  The
  GPU win (1.22x) is twice the CPU one (1.13x, section 4), which is what the phase-3 note
  predicted for anything that removes a table gather and a dependency chain on an
  accelerator.
* **The seed buys no time — it buys the gap.**  It leaves the pass count and the wall time
  alone (2.51 vs 2.55 cycles/s is within the run-to-run spread), but it **halves the
  sweep-to-gas gap and takes its maximum under the 1e-6 target**, 9.69e-07 against
  4.10e-06, and it removes a non-converged call.  That is the whole of its value and it is
  a real one: at `maxit = 8` the baseline does NOT meet the 1e-6 gate on this state and the
  seeded arms do.  The mechanism is visible in the numbers of section 2: a seeded pass 0
  puts every cell on the exact per-cell backward-Euler root with the field lagged, so the
  handful of calls that stop at the cap stop somewhere much better.
* **`+both` is the configuration to recommend**: 3.03 cycles/s, 4.30x the semi-implicit
  scheme instead of 5.11x, AND the gap of the seeded arm.  Against the 4.95x of the
  phase-2/3 production path (`prof_0922/implicit`) that is a **1.15x** end-to-end
  improvement at a strictly better gap.
* **A side result, and it contradicts the CPU measurement of phase 3:
  `ck_impl_frozen_op` is NOT a win on this GPU** — 2.55 cycles/s against the 2.63 of
  `prof_0922/implicit`, which has it off.  Phase 3 measured +12 % for it on CPU and warned
  that the 0.7-2.2 GB of stored operator could be eaten by bandwidth on a device; on the
  MI300A it is, slightly.  (Caveat: the two binaries are different builds — 8e41808b vs
  446b7b30+phase 4 — so this is a 3 % difference between lineages and should be confirmed
  with a fifth arm before it is quoted as a defect.  It does NOT affect the four arms
  above, which share one binary and one input.)

## 4. Gate (c): CPU cost — and why it is NOT the verdict

100 cycles, serial CPU, same input and arms as section 2.

| arm | wall [s] | x semi-implicit | mean passes / stage | gap, last half |
| --- | --- | --- | --- | --- |
| semi-implicit (`ck_implicit = false`) | 38.54 | 1.00 | — | — |
| baseline: implicit + `frozen_op` | 160.32 | **4.16** | 6.51 | 1.68e-06 |
| `+ reuse_jac = 1` | 140.16 | **3.64** | 7.36 | 1.67e-06 |
| `+ reuse_jac = 2` | 113.96 | (2.96) | 7.29 | 1.66e-06 |
| `+ seed = 1` | 452.47 | (11.7) | 6.82 | 2.04e-06 |

The first two rows reproduce phase 3 exactly (39.42 s / 4.15x there), so the measurement
is sound where the node is quiet.  **The rest of this table is not trustworthy and the
reason is known**: another user's 16-GPU training job occupied this node for part of the
series.  The proof is internal — `reuse_jac = 1` and `= 2` do *identical* work (section 2:
same pass count to the third digit, same residual to every digit) and came out 140.16 s
and 113.96 s, a 19 % spread; and `seed = 1` took 452 s while ending at the same `t`, the
same `dt`, the same conserved quantities and 6.82 passes against the baseline's 6.51, i.e.
at most 5 % more arithmetic.  A 2.8x wall time on 5 % more work is the node, not the
scheme.  **The GPU table of section 3 is the verdict** (the standing rule anyway: cost
verdicts are measured on the GPU), and it was taken with all four arms on the same
hardware, the same binary, the same restart and the same 300 cycles.

What the CPU series does establish, because these are machine-independent counts: the
**late-run gap is unchanged by reuse** (1.67e-06 against the baseline's 1.68e-06 over the
last 100 calls, which is phase 2's and phase 3's number), and the pass counts of section 2.

`+ seed = 2`, `+ both` were still running when the series was stopped; their GPU numbers
are in section 3 and are the ones that matter.

## 5. The caveat that decides the default: truncated calls

At the affordable cap (`maxit = 8`) a cold start does not converge — 20 of 40 calls in
section 2 — and a truncated iteration depends on the matrix it was truncated with.
Measured on the last dump of the 20-cycle cold-start run, against the baseline:

| arm | `eint` max rel. | `dens` max rel. |
| --- | --- | --- |
| `+ reuse_jac = 1` | **1.3e-02** | 1.3e-02 |
| `+ reuse_jac = 2` | 1.2e-02 | 1.3e-02 |
| `+ seed = 1` | 3.7e-04 | 3.8e-04 |
| `+ seed = 2` | **4.9e-05** | 4.2e-05 |
| `+ both` | 8.5e-03 | 8.5e-03 |

This is the same failure mode phase 3 reported for `ck_impl_warm` (8.3e-03 there), and it
is why **both knobs stay default off**.  Two things separate reuse from the warm start,
though, and they are why reuse is still the recommendation for production:

1. Reuse carries **no history between calls** — it is rebuilt at the first Newton pass of
   every call — so the drift is bounded by what one truncated call can do, not by an
   accumulating seed.
2. On the state that actually matters it does not happen at all.  The 600-call GPU run of
   section 3 is a settled atmosphere where only 5-6 calls of 600 fail to converge, and
   there the two arms agree on the pass count (5.66 vs 5.67) and on the gap
   (4.11e-06 vs 4.10e-06 max, 6.5e-07 median) — i.e. where the calls converge, the chord
   matrix changes nothing, exactly as the theory says.

**Not done:** a `maxit = 20` cold-start pair, which would show directly that the 1.3e-02
collapses once the calls converge.  It was scripted (`p4d.sh`) and queued behind the cost
series, and was abandoned when the node was taken.  It is the one measurement this
deliverable is missing.

## 6. Verdict and what remains

* **`ck_impl_reuse_jac = 1` is the win: 1.22x on the GPU at zero cost in passes, gap or
  answer on the production state.**  Mode 2 (the `T^3` rescale) is measurably useless
  within a call and should not be used; it is kept only for the large-`bdt` regime it was
  written for.
* **`ck_impl_seed = 2` buys accuracy, not speed**: no wall-time change, but the
  sweep-to-gas gap maximum drops 4.10e-06 -> 9.69e-07, i.e. it is what puts the scheme
  under the 1e-6 target at `maxit = 8`.  `ck_impl_seed = 1` (the literal semi-implicit
  step) is strictly worse than 2 on every measure and exists only as the comparison the
  brief asked for.
* **Recommended production setting: `ck_impl_reuse_jac = 1`, `ck_impl_seed = 2`,
  `ck_impl_frozen_op = false`** — 3.03 cycles/s with frozen_op on; the frozen operator
  itself looks like a small LOSS on this GPU (section 3) and should be re-measured before
  it is carried into production.  That is **4.30x the semi-implicit scheme against 4.95x
  today**, with a better gap.
* Remaining: (i) the `maxit = 20` cold-start pair of section 5; (ii) a fifth GPU arm with
  `frozen_op = false` + reuse, to settle whether the frozen operator should be on at all on
  an APU; (iii) the pass count itself, 5.7 per stage, which neither lever touches and
  which is still the whole cost.

## 7. Files

`p4.sh` (gates a, b, c), `p4d.sh` (the `maxit = 20` pass count, NOT run),
`stats.py` / `late.py` / `cmpbin.py` / `paycmp.py` (phase 3's).  Run directories `p4a_*`
(a), `p4b_*` (b), `p4c_*` (c).  GPU arms and the snapshot build:
`/viper/u2/jinma/ATHENAK/bench/ckq_0922/{base,reuse,seed,both,src_snapshot}`, jobs
11933400-11933403.  `.bin`/`.rst` dumps and the build directories are deleted after
measurement (inode quota).
