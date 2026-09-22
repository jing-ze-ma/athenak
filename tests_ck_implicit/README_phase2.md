# `problem/ck_implicit`, phase 2 — large `dt` and cost

Phase 1 is `README.md` in this directory: the backward-Euler correlated-k column solve,
what it is, and the two things it could not do.  This file is only about those two:

1. it **diverged at 10x and 100x** the production `dt` (`README.md` gate e);
2. it cost **5.1-5.4x the RT wall time**, 12-16 ck sweeps per hydro step against 2.

Code touched: `src/utils/two_stream_column_ck.hpp` (the thin-cell solve, the split, the
column-skip bookkeeping, `ck_impl_debug`), `src/utils/two_stream_rt.hpp` (the
classification, the skip early-returns, the emission store, the driver loop),
`src/pgen/deep_hot_jupiter_rt.cpp` (the new knobs and the once-per-step hook),
`inputs/tests/dhj_ck_implicit.athinput` (the knobs, all at their defaults).

Default is still **off**, and off is still bitwise today's scheme (gate a).

---

## 0. New parameters

| knob | default | what it does |
| --- | --- | --- |
| `problem/ck_impl_arat` | `2.0` | the two-level split: a cell is solved inside the tridiagonal only while the field it absorbs is within this factor of its own emission, `A <= arat E`.  Everything else takes a bracketed per-cell implicit step. |
| `problem/ck_impl_tau_min` | `0.0` (off) | the grey template's discriminant, the cell's own optical depth.  Kept, and available, but **measured to be the wrong one here** (section 1). Both gates must pass for a cell to be called thick. |
| `problem/ck_impl_colskip` | `true` | stop sweeping a column once its residual is under `ck_impl_tol`. |
| `problem/ck_impl_once` | `false` | run the whole radiation ONCE per hydro step with the full cycle `dt`, after the last RK stage, instead of once per RK stage. |
| `problem/ck_impl_debug` | `0` | print the worst cell of every pass with its whole row. |

---

## 1. What was ported, and what had to change

The grey mode-3 column had exactly this problem and solved it with `rt_impl_tau_min`
(deleted in `912ef43c`; the design note survives at `912ef43c~1:src/utils/two_stream_rt.hpp`
around line 763, the implementation around 5093).  **Only a cell whose own Rosseland
optical depth reaches a threshold goes into the tridiagonal; the thin ones keep a bounded
nonlinear per-cell relaxation and enter the thick cells' rows only through the field they
radiate.**  The idea is ported; two pieces of it are not, and both changes are measured,
not assumed.

### 1.1 The discriminant is `A/E`, not `dtau`

The only optical depth available *before* the sweep on the correlated-k path is the
**continuum** one — the line opacity is a function of the g index, not of the cell, and
is built inside the sweep.  A Planck-weighted continuum `kappa rho dr` under-estimates
`dtau` badly on this grid.  Measured, `ck_impl_tau_min = 1.0` on that quantity:

| configuration, first RT call, `maxit = 20`, `ck_spherical = true` | thin cells | residual |
| --- | --- | --- |
| no split | 0 | **5.70e-05** |
| `ck_impl_tau_min = 1` on the continuum `dtau` | 5184 of 8192 (63 %) | **2.58e-03** |

It calls 63 % of the production column thin, the iteration becomes per-cell Picard almost
everywhere, and the residual after 20 passes is **45x worse** than with no split at all.
So that gate is kept but defaults to off.

What the divergence is actually governed by is in the apply block's own note: the
linearised step asks for `de = e (A - E)/(4 E)` where the true backward-Euler root is
`e ((A/E)^(1/4) - 1)`, so the overshoot is a function of **`A/E` alone** — 1.3x at
`A/E = 2`, 3.8x at 16, unbounded as `E -> 0`.  `dtau` is only ever a proxy for it (a
thermalised cell has `A/E = O(1)`).  `A` and `E` are both stored by the apply block
already (`cksrc_g` and the new `ckem_g`), so the split can test the thing itself.  A
*cooling* cell (`A < E`) is never split out: there the linearisation under-steps.

The classification carries a **one-pass lag** by construction: the apply of pass `p`
writes it, `rt_pre_opac` of pass `p+1` freezes it into a second array and acts on it, and
the matrix assembly and the solve of that pass both read the frozen copy.  That is what
guarantees the decoupling and the rows can never disagree.  Zero-initialised = every cell
thin, which is the safe state the first pass of a run starts from.

### 1.2 The thin cell gets a bracketed implicit solve, not a relaxation

The grey template hands a thin cell the semi-implicit `(1-e^-x)/x` relaxation.  That is
bounded, but its fixed point is *not* the backward-Euler balance, so the residual test
could never converge on it.  Instead `CkThinSolve` solves, by bisection on a bracket that
is guaranteed by construction,

```
x - e* - bdt (A - Em (x/e)^4) = 0,    A = S + Em frozen at the current iterate
```

`f` is strictly increasing, `f(0) = -(e* + bdt A) < 0`, `f(max(e* + bdt A, e)) > 0`: the
root is unique, bracketed and **positive**, so the cell can never be driven through zero —
which is the whole failure mode at large `bdt`.  Small `bdt` gives `e* + bdt S` to
`O(bdt^2)`; large `bdt` gives radiative equilibrium `E = A`.  Both limits are exact, and
the fixed point over the Newton passes IS the backward-Euler balance, so thin cells are
still measured by the same residual test.  The same solve replaces the old bare-explicit
column fallback.

---

## 2'. Gate (b/e): behaviour at 10x and 100x the production dt

Serial CPU, `ck_spherical = true`, `ck_beam_sph = true`, `ck_impl_maxit = 20`,
`ck_impl_demax = 0` (the total-excursion bound RELEASED, so nothing is hiding behind it),
first RT call of the run from the analytic initial condition -- the same protocol as
phase-1 gate (e), and the hardest state there is.  `bdt` = 28.98 s / 289.8 s / 2898 s.

| dt | split | per-pass cap | residual after 20 | max step | capped | fallbacks |
| --- | --- | --- | --- | --- | --- | --- |
| 1x | off | 0.5 | 6.4e-05 | 0.037 | 0 | 0 |
| 1x | on | 0.5 | 1.95e-04 | 0.101 | 0 | 0 |
| 10x | off | 0.5 | 1.11e+00 | 0.5 (capped) | 24 | 0 |
| 10x | on | 0.5 | 1.12e+00 | 0.5 (capped) | 28 | 0 |
| 100x | off | 0.5 | 1.15e+02 | 0.5 (capped) | 64 | 0 |
| 100x | on | 0.5 | 9.85e+01 | 0.5 (capped) | 64 | 0 |
| 10x | off | **released** (1e30) | **2.02e+04** | **5.6e+04** | 0 | 0 |
| 10x | on | **released** (1e30) | 2.92e+04 | **1.44** | 0 | 0 |

**What changed and what did not.**

* **The step is bounded.** With the caps released -- the configuration in which phase 1
  reported the iteration *diverging* to 3.3e+05 -> 2.6e+07 at 10x and 3.7e+12 at 100x --
  the largest Newton step the solver now asks for at 10x is **1.44 e**, against
  **5.6e+04 e** without the split: four and a half decades.  No cell is ever driven
  through zero, the temperature stays positive everywhere, and there are **no fallbacks
  and no non-finite steps in any arm**.  That is what `CkThinSolve` buys, and it is the
  difference between a solver that can be run at large `dt` and one that cannot.
* **It still does not converge at 10x and 100x from the cold start, and the caps are
  still needed.**  The residual is bounded and oscillates at O(1)-O(100) instead of
  exploding, but it does not fall to `ck_impl_tol`.  Target not met.

**Why, measured with `ck_impl_debug = 1`.**  The worst cell at 100x is not a thin cell at
all: `i = 96`, `e = 4.75e+02`, `e* = 2.57e+02`, `S = -20.4`, `E = 294`, `A/E = 0.93`,
i.e. a THICK, COOLING cell one per cent from its own radiative equilibrium.  Its residual
is 1.0e+02 only because `bdt |S| = 5.9e+04` is 124x its own energy: at `bdt` = 2898 s the
quantity `bdt (dE/de)/e` is **7.2e+03**, so `|R|/e = 7.2e+03 (de/e)` and the stated
tolerance of 1e-8 in that norm demands the energy to 1.4e-12 relative -- below the
round-off of the sweep that evaluates it.  **At 100x the production `dt` the residual
norm `|R|/(e + eps e_max)` is not a reachable measure**, independently of the solver.  On
top of that the backward-Euler solution at 100x is a *different atmosphere* (most of the
column relaxes to radiative equilibrium in one step) and 20 damped passes cannot walk
there from the analytic start.

**What is therefore still open, and what would fix it** (not done, measured enough to
name): (i) test the residual in the **diagonally scaled** norm `|R|/(b_ii (e + eps
e_max))` with `b_ii = 1 + 4 bdt E/e`, which is the Newton step in energy units and is
`dt`-robust where `|R|/e` is not; (ii) a **line search** on the ck residual (each trial
costs a sweep, so at most one backtrack per pass); (iii) start the large-`dt` test from a
*relaxed* state rather than the cold start.

**The relaxed-start arm was attempted and did not measure what it was meant to** (scripts
`bigdt.sh`, directories `k_*`): restarting a 200-cycle state with `time/cfl_number` raised
to 3 and 30 produced **byte-identical** `passes/res/dstep` for the two, because the
timestep after a restart is re-limited by the code's own dt controls and never reached
10x, let alone 100x.  Those six runs are therefore not evidence either way and are not
tabulated.

---

## 4. Gate (f): cost

100 cycles, `ck_spherical = true`, `ck_beam_sph = true`, serial CPU, one process at a
time on an otherwise idle node.  `wall` is the whole run, which is the same quantity
phase 1's 5.1-5.4x was; the RT-only ratio is necessarily a little *larger*, since the
hydro is in both numbers.  `gap` is the largest `|ckdesum|` over the last 60 RT calls --
the sweep-to-gas gap, i.e. what the deliverable's 1e-6 target is on.

| arm | wall [s] | x off | ck sweeps / hydro step | gap |
| --- | --- | --- | --- | --- |
| semi-implicit (`ck_implicit = false`) | 33.15 | 1.00 | 2 | — |
| **phase 1** (no split, no colskip, per stage) | 187.43 | **5.65** | 12.9 | 1.7e-06 |
| + the `A/E` split | 189.61 | 5.72 | 13.0 | 1.7e-06 |
| + `ck_impl_colskip` | 159.03 | 4.80 | 13.0 | 1.7e-06 |
| + `ck_impl_once` (`maxit = 8`) | 103.10 | 3.11 | 8 | 2.2e-07 |
| + `ck_impl_maxit = 6` | 89.54 | **2.70** | 6 | 1.8e-07 |
| + `ck_impl_maxit = 5` | 81.48 | **2.46** | 5 | 5.6e-07 |
| + `ck_impl_maxit = 4` | 67.86 | **2.05** | 4 | 1.3e-05 |
| `once`, `maxit = 6`, colskip OFF | 98.21 | 2.96 | 6 | 1.8e-07 |

**Read this way.**

* **`ck_impl_once` is the whole win**: 5.65x -> 3.11x, exactly the factor 2 in sweeps that
  leaving the RK stage buys, and it makes the gap *better* (2.2e-07 against 1.7e-06),
  because one solve over the full `dt` is a better-conditioned problem than two over half.
* **`ck_impl_colskip` is worth 8-16 %** (189.6 -> 159.0 per stage; 98.2 -> 89.5 once), and
  costs nothing in accuracy: a retired column's residual is already below tolerance.
* **The split costs 1 % and buys nothing at the production dt** -- at 1x almost no cell has
  `A > 2E`, so almost nothing is split out.  It is insurance for large `dt`, not a speed-up.
* **Loosening the tolerance is not available.**  Measured, the sweep-to-gas gap is about
  **1e2 times `ck_impl_tol`**: at `tol = 1e-6` the relaxation run's `ckdesum` runs to
  **-1.3e-04**, a hundred times over the deliverable's own target.  The gap is a
  `dx`-weighted volume integral while the residual is a max over cells normalised by
  `e + 1e-3 e_max`, and the two differ by that factor on this grid.  So `ck_impl_tol`
  stays at 1e-8 and the pass count is controlled with `ck_impl_maxit` instead.

**Against the 2x target.**  The best configuration that keeps the gap at or under 1e-6 is
`ck_impl_once = true`, `ck_impl_maxit = 5`: **2.46x**, gap 5.6e-07.  `maxit = 4` reaches
**2.05x** but its gap is 1.3e-05, thirteen times the target.  **The target is therefore
not met: 2.46x at the required accuracy, against 5.65x before -- a 2.3x reduction, but
not to 2x.**  The untried lever is a **warm start** (seed pass 0 with the previous step's
accepted `de`); it is not implemented here and would plausibly remove one pass, which is
the remaining factor.

## 5. Gate (d): the two schemes' relaxed states

200 cycles from the analytic initial condition with `ck_spherical = true`, one arm
semi-implicit and one `ck_implicit` (`tol = 1e-6` for the implicit arm), then the
`ck_dump_file` column out of each final restart, for three columns.  `max |T_impl -
T_semi|` by pressure range, out of 2000-6400 K:

| column | `p > 10 bar` | `1-10 bar` | `0.01-1 bar` | `1e-4 .. 1e-2 bar` | `p < 1e-4 bar` |
| --- | --- | --- | --- | --- | --- |
| night, `mu0 = -0.92` | 0.03 K | 0.13 K | 12.8 K | 49.7 K | 38.8 K |
| day, `mu0 = +0.92` | 0.03 K | 0.07 K | 4.0 K | 9.2 K | 14.4 K |
| day, `mu0 = +0.38` | 0.03 K | 0.11 K | 12.7 K | 13.5 K | 8.7 K |

Same picture as phase 1's 600-cycle table, at a third the amplitude because the state is
less far along: below 1 bar the two schemes agree to 0.13 K out of 6400 K (2e-5 relative),
and they separate only in the thin, stiff region above 0.01 bar where the semi-implicit
apply's damping and `rt_de_max` are active.  **Same caveat as phase 1: 200 cycles is not
relaxed**, the two arms' `dt` histories differ slightly, and part of the upper-atmosphere
difference is a time offset this table cannot separate from a scheme difference.

## 6. Gate (f'): once-per-step against per-stage, on accuracy

Both arms implicit, both relaxed 200 cycles, `maxit = 6` for the once arm.  `max |T_once -
T_per-stage|`:

| column | `p > 10 bar` | `1-10 bar` | `0.01-1 bar` | `1e-4 .. 1e-2 bar` | `p < 1e-4 bar` |
| --- | --- | --- | --- | --- | --- |
| night | 0.77 K | 0.22 K | 0.25 K | 5.9 K | 4.9 K |
| day, `mu0 = +0.92` | 0.76 K | 0.23 K | 0.16 K | 1.0 K | 8.6 K |
| day, `mu0 = +0.38` | 0.77 K | 0.22 K | 0.26 K | 5.2 K | 4.9 K |

**Once-per-step is a smaller perturbation than the choice of scheme.**  Everywhere below
0.01 bar the two cadences agree to 0.8 K, and in the thin region they differ by 1-9 K
against the 4-50 K that separates implicit from semi-implicit in the same cells.  Its
sweep-to-gas gap is also the better of the two (2.2e-07 against 1.7e-06).  The box's
`rt_col3_once` was rejected on accuracy grounds for the He box; on this column, measured,
it is not the accuracy-limiting choice.  (The same caveat applies: 200 cycles, and part of
the difference is a `dt`-history offset.)

## 7. Recommended production configuration

```
problem/ck_implicit     = true
problem/ck_impl_once    = true
problem/ck_impl_maxit   = 5        # 6 if the gap must be <= 2e-07
problem/ck_impl_tol     = 1e-8     # do NOT loosen: the gap is ~1e2 x tol
problem/ck_impl_colskip = true
problem/ck_impl_arat    = 2.0
problem/ck_impl_demax   = 0.5      # keep the belt; never reached at the production dt
```
2.46x the semi-implicit wall time, sweep-to-gas gap 5.6e-07.

---

## 2. Gate (a): default off is still inert

Serial CPU, `inputs/tests/dhj_ck_implicit.athinput`, 20 cycles, `ck_beam_sph = true`,
this tree against a pristine copy of `ac42334d`, both reading the same input:

| `ck_spherical` | `dhj.hydro.hst` | `bin` payloads |
| --- | --- | --- |
| `true` | **identical** | **identical** (2 dumps) |
| `false` | **identical** | **identical** (2 dumps) |

The `bin`/`rst` FILES differ by 393 bytes because the embedded `PAR_DUMP` now carries the
five new `<problem>` keys; the binary payload after `header offset` is byte-for-byte the
same, and the history file is identical byte-for-byte.  Nothing else in the tree changed:
the phase-1 section-0 apply correction is untouched and is already in `ac42334d`.

## 3. The conduction hand-over at `ck_pcut_bar` — still split

Unchanged from phase 1, and deliberately: the radial radiative conduction
(`<hydro>/rad_implicit_x1`, the `rad_tau_lo/hi` blend) stays a separate operator, so the
blend's handover flux enters `S_i` explicitly, lagged by one operator, exactly as with
`ck_implicit` off.

**What folding it in would take.**  `Conduction::ImplicitRadialUpdate` already owns a
radial tridiagonal in the same unknown (the cell internal energy) over the same column,
and the two matrices would simply add — the blend weight `w` already partitions the flux
between them, `w` to the conduction and `1-w` to the two-stream, so the merged row is
`I - bdt (J_ck + J_cond)` with no double counting.  Three concrete pieces of work: (i)
`ImplicitRadialUpdate` builds its rows from a grey Rosseland conductivity and would have
to accept externally supplied off-diagonals (`ck_jac`) and an externally supplied residual
instead of forming its own; (ii) the two operators currently run in different tasks at
different points of the stage, so the merged solve has to be hoisted into one of them and
the other's task made a no-op (the grey mode-3 path does exactly this, and `rt_col_active`
is the flag that already exists for it); (iii) the residual of the merged system must be
re-evaluated by re-running *both* operators per Newton pass, which makes a pass cost a ck
sweep plus a conduction flux build rather than a ck sweep — about 15 % more per pass on
this grid.  It is a self-contained second deliverable, not a change to what is here.
