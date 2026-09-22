# `problem/ck_implicit`, phase 3 — the FROZEN EXCHANGE OPERATOR and the WARM START

Phase 1 is `README.md` (the backward-Euler correlated-k column), phase 2 is
`README_phase2.md` (the thin-cell bracketed solve, `ck_impl_once`, column skipping, and
the cost table that ends at **2.46x** the semi-implicit wall time and only WITH
`ck_impl_once`).  `ck_impl_once` is **off for everything in this file**: the user's rule is
that once-per-cycle splitting is not a default, so the target here is the per-stage
cadence, whose phase-2 cost is 4.7-5.7x.

What is new, all of it default OFF:

| knob | default | what it does |
| --- | --- | --- |
| `problem/ck_impl_frozen_op` | `false` | freeze the exchange operator: pass 0 stores what the sweep built out of the opacity, every later pass re-applies it to the new `B_b(T)` with no k-table look-up, no exponential and no beam |
| `problem/ck_impl_frozen_cof` | `true` | store the three half-layer coefficients as well as `kappa rho` (the memory/`expm1` trade) |
| `problem/ck_impl_warm` | `false` | start the Newton from the previous call's converged increment |

Code: `src/utils/two_stream_column_ck.hpp` (the knobs, the stored-operator arrays,
`CkWarmSeed`, the increment accumulator), `src/utils/two_stream_rt.hpp` (the store/load
inside `rt_chain_ck`, the frozen beam, the driver hook, `e^n` recorded net of the seed),
`src/pgen/deep_hot_jupiter_rt.cpp` (plumbing), `inputs/tests/dhj_ck_implicit.athinput`.

---

## 0. The idea, and why it is exact

At frozen opacity the thermal two-stream is **linear in the band Planck functions**: every
layer's source is a convex combination of the two centre `B_b`, every `step` is
`I -> (1-e0) I + cin s_in + cout s_out`, and the face mixing and the frame changes are
multiplicative factors.  So the sweep is the application of a fixed operator `M` to `B`,
and everything in `M` — the layer optical depths, `e0/cin/cout`, the `BFace` weights, the
`(1+beta)` mixing, the area ratios — depends on the **opacity and the geometry alone**,
which `ck_implicit` already freezes over the step (`ck_impl_refresh_kappa = false`).

Storing `M` itself is not an option (it is a dense per-column matrix per g-point).  What is
stored is what the sweep computes *from the opacity*:

* `kappa rho` per (cell, chain) — **the whole k-table look-up**, which the sweep pays four
  times per cell and chain (the two `ck_spherical` probe passes, the down-sweep and the
  up-sweep);
* with `ck_impl_frozen_cof`, the triple `(e0, cin, cout)` of that cell's HALF layer.  Every
  `step` in the kernel is the half of one cell taken at one chain's `mu` — the probes, the
  down-sweep and the up-sweep all cross the same half layers — so **one triple per (cell,
  chain)** describes the operator completely, and it also serves the Jacobian assembly
  (`jcof`) in double precision;
* the top boundary layer's `(1 - e^-dtau)` per (column, chain);
* the **direct beam**: at frozen opacity `Qb_g` is exactly temperature-independent, so a
  frozen pass keeps pass 0's deposit and never re-runs the pseudo-spherical ray
  integration (the O(N^2) chord walk of `ck_beam_sph`).

A frozen pass therefore re-forms the optical depths from the stored `kappa rho` **by the
same expressions** and runs the same recurrences on the new `B`, so it is not merely
consistent with the pass it replaces, it is **bitwise** that pass.  Gate (b) measures
exactly that and finds it.

**Memory (the choice asked for).**  The store is 4 Reals per (cell, chain) with
`ck_impl_frozen_cof = true`, 1 without.  On the production mesh — 6 x 32 x 32 columns, 88
chains (`ck_nquad = 1`: 11 bands x 8 g-points), ~80 correlated-k layers — that is
**1.4 GB** against **0.35 GB**; at the full `nx1 = 128` of the radial grid, 2.2 GB against
0.55 GB, plus the ghost-zone inflation the arrays carry (they are allocated over `n1`,
`n2`, `n3`).  **The default is the full store**, because the reduced one is *measured to be
slower than not freezing at all* (gate c): dropping the coefficients keeps the memory
traffic and gives back the `expm1`, which on this CPU is the cheaper of the two.  So the
answer to "too large for a GPU" is **not** `ck_impl_frozen_cof = false`; it is to store per
band block and loop the blocks, which is not built here.

---

## 1. Gate (a): `ck_implicit = false` is still bitwise the reference

Serial CPU, `inputs/tests/dhj_ck_implicit.athinput`, 20 cycles, one `bin` dump per cycle,
`ck_beam_sph = true`, this tree against a reference binary built from the unmodified tree
(`build_ckfast_ref`, snapshot of the working tree at the start of this work).

| `ck_spherical` | `dhj.hydro.hst` | `bin` payloads |
| --- | --- | --- |
| `false` | **identical** | **22 of 22 identical** |
| `true` | **identical** | **22 of 22 identical** |

(The `bin` FILES differ in length because the embedded `PAR_DUMP` now carries the three new
`<problem>` keys; `paycmp.py` compares everything after the stated header offset, which is
the data.)

## 2. Gate (b): the frozen operator against the phase-2 code path

Same input, 20 cycles with the hydro live, `ck_spherical = ck_beam_sph = true`,
`ck_implicit = true` and the phase-2 defaults (`maxit = 8`, `colskip = true`,
`arat = 2.0`), one dump per cycle.  Both paths are in the same binary and are selected by
`ck_impl_frozen_op`.

| arm | `dhj.hydro.hst` | 22 `bin` payloads | passes/call | res (max) | `ckdesum` (max) |
| --- | --- | --- | --- | --- | --- |
| phase 2 (`frozen_op = false`) | — | — | 6.75 | 1.79e-05 | 9.45e-06 |
| `frozen_op = true` | **identical** | **22/22 identical** | 6.75 | 1.79e-05 | 9.45e-06 |
| `frozen_op = true`, `frozen_cof = false` | **identical** | **22/22 identical** | 6.75 | 1.79e-05 | 9.45e-06 |

**Bitwise, not "to round-off".**  The converged temperature is the same bits, the pass
count is the same integer call by call, and the sweep-to-gas gap is the same number: the
frozen pass is the pass it replaces.  (The maxima above are the cold-start transient of
cycle 0-1; over the last half of the 100-cycle run of gate (c) the gap is **1.68e-06** on
both arms, which is phase-2's per-stage number and is what the 1e-6 target is measured
against.  It is unchanged by construction.)

That also answers gate (d) — "a 20-step run with hydro live, `T(p)` vs phase 2 within
round-off" — in the strongest possible form: the whole state, every cycle, is identical.

## 3. Gate (c): cost

100 cycles, serial CPU, `ck_spherical = ck_beam_sph = true`, one process at a time (the
node was not idle: other work of the session was running on other cores, so treat the
ratios as good to a few per cent).  `wall` is the whole run, hydro included — the same
quantity phase 1 and 2 quoted.

| arm | wall [s] | x semi-implicit | ck sweeps / hydro step | mean passes / stage |
| --- | --- | --- | --- | --- |
| semi-implicit (`ck_implicit = false`) | 39.42 | 1.00 | 2 | — |
| phase 2, per stage | 185.34 | **4.70** | 13.0 | 6.51 |
| + `ck_impl_frozen_op` | 163.79 | **4.15** | 13.0 | 6.51 |
| + `ck_impl_frozen_op`, `frozen_cof = false` | 224.81 | 5.70 | 13.0 | 6.51 |
| + `ck_impl_frozen_op` + `ck_impl_warm` | 131.76 | **3.34** | 12.7 | 6.34 |

**Read this way.**

* **The frozen operator is worth 12 %** (4.70x -> 4.15x) and costs nothing in accuracy,
  because it changes nothing at all.  Since the two arms take the same passes on the same
  states, 163.79/185.34 IS the cost of a pass, hydro included: a frozen pass is ~0.88 of a
  full one, and a little less than that once the hydro is taken out of both.  (A separate
  `maxit = 1` / `maxit = 5` pair, intended to give the marginal cost of a pass directly,
  was invalidated by two copies of its own script running at once and is not reported.)
* **The k-table look-up is NOT the bottleneck on this CPU.**  That was the premise of the
  idea and it is wrong here: the 11 x 8 k-table is small enough to sit in cache, and what a
  pass actually spends its time on is the four column recurrences themselves (the two
  `ck_spherical` probe passes, the down-sweep and the up-sweep) and their traffic over
  `Bb_g`, `Src_g`, `Fb_g` and the Jacobian atomics — none of which can be frozen, because
  all of them are the parts that depend on `B(T)`.  Storing `kappa rho` alone makes the
  run **slower than not freezing** (5.70x), which is the same statement from the other
  side: the added array traffic outweighs the look-up it removes.
* **The warm start is worth another 20 %** (4.15x -> 3.34x) — but not through the pass
  count, which it moves only from 6.51 to 6.34.  It pays through `ck_impl_colskip`: seeded
  columns cross the tolerance earlier, so the later passes sweep fewer columns.
* **The 2x target is not met.**  3.34x per stage against phase 2's 4.70x is a 1.41x
  reduction at an answer that is bitwise unchanged for the frozen part, but the remaining
  factor is the pass count (6.3 passes per stage, 12.7 sweeps per hydro step against 2),
  and no amount of freezing the opacity touches that.  The levers that would are the ones
  phase 2 already identified and that are ruled out or unbuilt here: `ck_impl_once` (a
  factor 2, forbidden as a default), a smaller `maxit` (phase 2: `maxit = 4` reaches 2.05x
  but its gap is 1.3e-05), and a line search / better-conditioned residual norm.

## 4. Gate (b'): the warm start, passes with and without

`ck_impl_maxit = 20` so that the pass count is a convergence measure and not the cap,
20 cycles from the analytic initial condition, `frozen_op = true`:

| arm | calls | mean passes | not converged |
| --- | --- | --- | --- |
| cold | 40 | 13.43 | 21 |
| warm | 40 | 14.09 | 19 |

**On this transient the warm start does not reduce the pass count** — it slightly increases
it.  The 20-cycle window is entirely the cold-start relaxation from the analytic initial
condition, where the previous call's increment is a poor predictor of the next one's (and
the two RK stages of a cycle see different `bdt`).  Over the 100-cycle run of gate (c),
where the state has begun to settle, it is worth 2.6 % of the passes and 20 % of the wall
time through the column skipping.  A verdict on a genuinely relaxed atmosphere needs the
thousands of cycles phase 1 and 2 both flagged and was not run.

**Accuracy of the warm start, and why I do not recommend switching it on yet.**  It moves
the starting iterate and not the fixed point, so where the residual test actually passes it
changes nothing.  But on this transient **half the calls do not converge** — 100 of 200 at
`maxit = 8`, and still 21 of 40 at `maxit = 20` — and a truncated iteration does depend on
where it started.  Measured, warm against cold, same binary, `frozen_op = true`:

| | worst cell, relative, `eint` |
| --- | --- |
| `maxit = 8`, after 1 cycle / after 20 | — / **8.3e-03** |
| `maxit = 20`, after 1 cycle / after 20 | **1.1e-01** / **9.7e-03** |

So the seed is not benign at the pass caps that are actually affordable: a thin top cell
whose call stopped at res ~1e-5 lands somewhere else, and the `dt` histories then diverge.
**The 20 % it buys is therefore not free**, and `ck_impl_warm` stays default off pending a
relaxed-state test where the calls converge.  The frozen operator carries none of this: it
is bitwise.

## 5. GPU layout (not measured — there is no GPU in this deliverable)

The frozen re-application inherits the parallelisation the sweep already has:
`rt_chain_ck` is launched over `(m, blk, k, j)`, so the **g-point blocks are already a
parallel dimension** (`rt_split`'s doing), 22 blocks of 4 chains here, and freezing the
operator removes no parallelism — the only serial dependency is the radial recurrence,
which is what a column solve is.  The stored arrays are laid out `(m, chain, i, k, j)`,
the same layout as `Bb_g` and `kc_g`, so that lanes of a wave — which vary in `j` — read
adjacent Reals; the layout note at the head of the kernel is the measured reason.

What I expect on an accelerator, and why it may well be better there than the 12 % measured
here: the k-table look-up is a random global-memory gather rather than a cache hit, the
`expm1` is a multi-instruction sequence on a dependency chain the kernel cannot hide (the
kernel's own note says its cost is "exp() latency on a dependency chain"), and the
`ck_beam_sph` chord walk is an O(N^2) serial loop per (face, chain) that a frozen pass
skips entirely.  Against that, the stored arrays are 1.4-2.2 GB on the production mesh and
every frozen pass streams them, so on a bandwidth-bound device the win could equally be
eaten — which is exactly what the CPU measurement of `frozen_cof = false` shows in
miniature.  **Untested; do not quote a GPU number from this file.**

## 6. What remains

1. The pass count, which is now the whole cost (12.7 sweeps per hydro step).  Nothing in
   this deliverable addresses it.
2. Per-band-block storage, which is the answer to the memory if the operator is ever to be
   frozen on a GPU-sized mesh.
3. A relaxed-state warm-start measurement (thousands of cycles), and a `maxit` large enough
   that the calls converge, so that the warm start can be judged on passes alone.
4. The conduction hand-over at `ck_pcut_bar` is still split (phase 2, section 3).

## 7. Files

`p3.sh` (gates a, b, d and the cost table), `p3b.sh` (the warm start at `maxit = 20`),
`paycmp.py` (bin payload comparison),
`cmpbin.py` (per-variable relative differences), `stats.py` / `late.py` (the verbose
line).  Run directories `p3a_*` (a), `p3b_*` (b, d), `p3c_*` (c), `p3w_*` (b').
The `.bin`/`.rst` dumps and the build directories are deleted after measurement (inode
quota).
