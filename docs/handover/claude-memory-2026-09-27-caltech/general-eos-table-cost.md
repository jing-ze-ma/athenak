---
name: general-eos-table-cost
description: "SUPERSEDED NUMBERS at the top -- as of 2026-08-17 the general EOS costs 3.07x per cycle / 2.32x per simulated second, not 4.4x/2.9-3.4x. The profile below still holds; the cost fell because of warm-start fixes, not the table."
metadata:
  type: project
---

## 2026-08-17: THE HEADLINE NUMBERS BELOW ARE OUT OF DATE

Re-measured after `2a154e7d` + `d1c289dd` (warm-starting the temperature inversions that
the pgen, srcterms and the well-balanced background were doing cold). Identical inputs,
identical hardware, `ohmic_resistivity = eos` in BOTH, only the EOS switched. dhj
64x56x128, bbot = 10 G, 8 ranks x 4 threads on a DEDICATED node:

| | s/cycle | dt | wall s / simulated s |
|---|---|---|---|
| `eos = ideal` + EOS x_e | 0.0448 | 11.52 | 0.00389 |
| `eos = general` + EOS x_e | 0.1376 | 15.24 | 0.00903 |

**3.07x per cycle, 2.32x per simulated second** (was 4.4x / 2.9-3.4x). The general EOS earns
a 32% LARGER dt because energy goes into dissociation instead of thermal motion, lowering
the sound speed.

Cumulative effect of the fixes on the general+eos-resistivity run: 0.1755 -> 0.1377
s/cycle, **21.5%**.

**RESTRICT ANY SUCH COMPARISON TO THE HYDRO-LIMITED WINDOW.** Averaged over all cycles it
looks like 3.22x/1.32x, and the 1.32x is an artefact: the IDEAL run's dt collapses to
exactly 1.1731 s partway through -- that is dt_diff at max_eta = 1e14, i.e. it becomes
RESISTIVITY-limited -- while the general run stays hydro-limited at 15.23 throughout. The
table above is t = 500..1500 where both are genuinely hydro-limited.

Why the ideal run falls in and the general one does not is the mu/T mismatch again: the
ideal run's log shows dfloor firing ~1e6 times while the general run's event log is
**completely EMPTY**. Fixed mu = 1.810 drives it into cold floored states the general EOS
never reaches. So a well-tuned ideal run needs `max_eta = 1e13` (see
`docs/ideal_gas_resistive.md`) to keep its advantage.

---

## Original measurement, 2026-08-16 (still correct as an attribution, wrong as a total)

Measured 2026-08-16 on dhj 64x64x128, 16x7, 500 cycles, bbot=3, by running the
intermediate configurations rather than guessing.

| configuration | ms/cycle | x ideal |
|---|---|---|
| ideal, no resistivity | 18.82 | 1.00 |
| ideal + perna | 20.31 | 1.08 |
| general(**gamma**) + perna | 23.13 | 1.23 |
| general(**table**), no resistivity | 80.75 | 4.29 |
| general(table) + perna | 88.17 | 4.68 |
| general(table) + **eos** resistivity | 88.15 | 4.68 |

Attribution per cycle: ideal baseline 18.82; +1.48 perna; **+2.83 the general-EOS machinery
itself** (reconstructing p and Gamma_1 through `wder`, the `wtemp` cache, solvers reading
reconstructed values); **+65.04 the TABLE**; **-0.02 for the x_e surface instead of the
perna fit** — i.e. `ohmic_resistivity = eos` is FREE, a 4th Hermite patch is cheaper than
perna's transcendentals, confirming the benchmark in [[resistivity-perna-uhj]].

Per SIMULATED second it is only ~3.5x, because the table gives a larger dt (13.41 vs 10.14 s
at bbot=3: different mu and Gamma_1, lower sound speed).

## ~3x of the table cost is redundant, and removing it is BITWISE safe

Per cell, `SingleC2P_General{Hyd,MHD}` currently costs `3N + 9` bicubic Hermite patch
evaluations, where N = root-find iterations (2-3 warm-started):

1. `Temperature()` -> `SolveTemperature` -> N x `EvalResidual` -> N x `Eval`, and `Eval`
   always interpolates ITE **and ITP and ITMU**. Mode 0 uses only `s.e` and `s.dlne_dlnt`,
   both of which come from the ITE patch alone (`Interpolate` already returns f, fx, fy).
   **ITP and ITMU are computed and discarded N times.**
2. `BelowPressureFloor(d,e,temp)` -> `Pressure()` -> a full `Eval` (3 patches).
3. `Pressure(d,e,temp)` -> a full `Eval` at the SAME (d,T) -- an exact duplicate of 2.
4. `Gamma1(d,e,temp)` -> a full `Eval` at the SAME (d,T) again -- a third.

Optimal is `N + 3`: an energy-only residual, then ONE shared `Eval` feeding the floor test,
p and Gamma_1. With N=3 that is 18 patches -> 6.

**Both changes are bitwise identical by construction** -- skipping ITP/ITMU does not alter
the residual, and reusing one `Eval` returns literally the same numbers -- so the fix can be
verified exactly against a reference run. If patch evaluation dominates the +65 ms this
should take the table cost to ~22-30 ms/cycle, i.e. **4.7x -> ~2.5x per cycle and ~3.5x ->
~1.9x per simulated second.**

## DONE 2026-08-16, commit `2dedcbdb` — and the prediction above was WRONG

Implemented as `EvalEOnly`/`EvalPOnly` in `eos_table.hpp` and `PressureAndGamma1` in
`eos.hpp`, with the two `general_c2p_*.hpp` evaluating once before the floor test and a
`stale` flag for the branches where a floor moves the state.

**Measured: 88.0 -> 79.0 ms/cycle, 10.3% overall, 12% off the table's share (4.73x -> 4.43x
ideal). NOT the ~2.5x predicted.** So bicubic patch evaluation is only ~a quarter of the
table cost. Most likely GCC was already CSE-ing the three identical inlined `Eval` calls at
the same (d,T), so only the root-find saving was real.

**Verification method worth reusing:** at the default `-O3` the results were NOT bitwise --
`dens`, `eint`, `bcc1`, `bcc2` identical, but the near-cancelling residuals (velocities,
`bcc3` which is 3e-6 of `bcc1`) moved by ~5e-8 of scale. Rebuilding BOTH versions with
`-ffp-contract=off` made all 42 dumps byte-identical, proving the change is algebraically
exact and the difference was GCC choosing different FMA fusions in the restructured code.
**Do this before concluding a refactor changed the answer.** All six general-EOS regression
tests pass.

## STILL ON THE TABLE: the transcendentals, which are NOT bitwise-safe to remove

The root find works in `z = log10 T`, then `EvalResidual` does `Pow10(z)` to get T, and
`Interpolate` immediately does `log10(t)` to get back to `z`. It also recomputes
`log10(rho)` every iteration although rho is FIXED throughout the solve. And
`g = log10(e) - ltarget` with `e = rho*Pow10(ev)` is analytically `log10(rho) + ev`.

So the mode-0 residual could in principle be one patch evaluation plus additions, with
`log10(rho)` hoisted out of the loop — removing three or four transcendentals per
iteration. That needs an `Interpolate` entry point taking `(x,y) = (log10 rho, log10 T)`
directly. **It is mathematically exact but NOT bitwise**, so it cannot be verified the same
way. Unprofiled: nobody has yet confirmed the transcendentals are actually where the time
goes. Profile first (`perf`), then decide.


## PROFILED 2026-08-16 (perf, single thread, startup excluded with -D 7000)

Build `-O3 -g -fno-omit-frame-pointer`; dhj 64x16x64, 80 cycles. Startup MUST be excluded or
`EOSCompositionModel::EvaluateAt` (table build, one-off) shows up at 10%.

| category | share of the time loop |
|---|---|
| **libm — log / exp / log10 / fmax** | **48.7%** |
| EOS table structural code | 32.4% |
| **actual fluid solver + RT** | **14.0%** |
| everything else | 4.4% |

Top symbols: `__ieee754_log` 16.4, `__ieee754_exp` 11.2, `EOS_Data::Temperature` 11.1,
`__log10_finite` 11.1, `EOSTable::HermitePatch` **6.0**, `EOS_Data::Pressure(double,double)`
4.0, `GeneralMHD::ConsToPrim` 3.9, `mhd::HLLD` 3.8, `exp` 3.2, `__fmax` 3.1.

**HermitePatch is only 6%** — which is exactly why `2dedcbdb` bought 10% and not 2.5x. The
patch-count model was the wrong model.

### Where the transcendentals come from, and why they are removable

The root find works in `z = log10 T`. Per residual evaluation it does:
`Pow10(z)` (an **exp**) to make T, then `Interpolate` immediately does **log10(t)** to get
`z` back; **log10(rho)** recomputed every iteration although rho is FIXED for the whole
solve; and `g = log10(e) - ltarget` where `e = rho*Pow10(ev)`, so analytically
`log10(e) = log10(rho) + ev` — a fourth transcendental that cancels.

So a mode-0 residual could be one Hermite patch plus additions, with `log10(rho)` hoisted
out of the loop. Needs an `Interpolate` entry taking `(x,y) = (log10 rho, log10 T)`
directly. **Mathematically exact but NOT bitwise**, so it cannot be verified the way
`2dedcbdb` was (see that entry for the `-ffp-contract=off` trick).

`__fmax` at 3.1% is a real libm CALL, not an inlined instruction — the `fmin`/`fmax` clamps
in `Interpolate`/`SolveLog` could be ternaries.

### Secondary: a cold-start root find on the boundary, ~4%

`EOS_Data::Pressure(double,double)` — the TWO-argument form, documented "setup-time use
only" because it solves for T itself with no warm start — is 4% of the time loop. It comes
from `wb_background.hpp:201`, inside `WBAdvance`, which the dhj outer-x1 boundary calls for
every ghost cell every stage. (`SetWbBackgroundPressure` is NOT the culprit: it early-returns
unless `use_wellbalance_static`, and runs once at init.) Passing a temperature guess there
would recover most of it.

**ESTIMATE, and I got the last one wrong:** removing the root-find transcendentals plausibly
takes 25-35% off the total. Verify, do not trust.

## THE ANSWER: tabulate on (log rho, log e), not (log rho, log T)

Measured 2026-08-16. Everything short of this is single digits; this is the only change
worth ~2x.

**The ceiling, measured not estimated.** Forcing the root find to exit after ONE Newton step
(`logtol = 1e30`, physics wrong, cost representative) is a lower bound on what a
no-inversion EOS would cost, and it is generous because that proxy still pays one residual
plus the final Eval:

| | ms/cycle | x ideal |
|---|---|---|
| ideal | 17.68 | 1.00 |
| general, ONE-iteration proxy | 45.97 | **2.60** |
| general, current | 79.95 | 4.52 |

So killing the inversion is worth **43%**, and a real (rho,e) table should beat the proxy --
it removes the residual and the Pow10/log10 round trip entirely -- landing near **2.0-2.3x**.

**Why it works:** the code stores (d,e) but tabulates on (rho,T), so every cell every stage
inverts e(rho,T) for T. Tabulate T, p, Gamma_1, mu, x_e on (log rho, log e) instead and
ConsToPrim becomes a direct lookup: 2 log10 plus 3-4 Hermite patches, no iteration.

Notes for whoever does it:
- e is monotonic in T at fixed rho (already relied on for the bracket), so the inverse table
  is well defined.
- Accuracy should IMPROVE across ionization zones. The current T grid needs `dlogt = 0.2 x
  dlogd` precisely because e changes fast at nearly constant T there; on an e grid that
  region is naturally resolved.
- The (p,T) and (rho,T) directions are still needed by pgens, `DensityFromPressureTemperature`
  and the WB background -- keep the analytic model or a second table for setup, where cost
  does not matter.
- The pressure floor needs e(rho,p), still an inversion, but only in floored cells.
- Gamma_1/chi_rho/chi_T are defined at constant T; tabulate them directly as their own
  surfaces rather than deriving them.

**Everything else, for completeness (measured unless marked):**

| change | effect |
|---|---|
| remove redundant Hermite patches (`2dedcbdb`, DONE) | 10% -> 4.43x |
| `logtol` 1e-13 -> 1e-8 | 7.6%, and perturbs the solution ~1e-4 over 500 cycles |
| remove the residual's transcendentals | ESTIMATED 25-35%, not bitwise |
| `fmax` -> ternary (it is a real libm call) | ~3% |
| warm-start `WBAdvance`'s `Pressure(d,e)` on the boundary | ~4% |

**Floor:** with a completely free EOS the general path still costs ~1.15x ideal (the +2.8
ms/cycle of reconstructing p and Gamma_1 through `wder`). So the useful range is
4.52x today -> ~2.1x achievable -> 1.15x unreachable.
