---
name: cs-shocks-through-seams
description: Shocks on the cubed sphere -- a shock CROSSES a seam fine, but strong structure sitting on one turned the whole domain to NaN silently until the along-seam resample was limited (b4f7aa25); records the new iprob=12 blast test and the three limiters that failed; FOFC is DECIDED not needed
metadata:
  type: project
---

Done 2026-08-31: `da426248` (the test), `b4f7aa25` (the fix), `f493c209` (a GPU build bug
in the test). Asked for because the user's applications put shocks through seams/vertices.

## THE TEST: iprob = 12, `inputs/tests/cubed_sphere_blast.athinput`

A geodesic cap of overpressure, uniform in RADIUS so the shock runs TANGENTIALLY and
therefore through the seam. **The gate needs no exact solution**: on a sphere the problem
is rotation-invariant, so the same blast centred in a panel INTERIOR, on a SEAM and on a
cube VERTEX must agree. The interior run is a control differing only in where the grid
sits under identical physics. Every other cubed-sphere test in the repo is SMOOTH, several
are exact equilibria, so none of them said anything about shock capturing.

## WHAT IT FOUND

    blast survival, pressure contrast vs p0=1    before        after
      panel INTERIOR                             fails 300     survives 1000
      SEAM                                       fails  70     survives 1000
      cube VERTEX                                fails  70     survives 1000

1. **A SHOCK CROSSES A SEAM PERFECTLY WELL**, even unlimited -- the interior blast drives
   its front past the seam (1.163 rad against a seam at 0.785, p ~ 28 there) and finishes
   clean. Seams are not a barrier to a PASSING shock. I predicted the opposite.
2. **Every failure came from ONE thing**: the unlimited along-seam resample on the
   CELL-CENTRED variables. The interior case failed only once its own shock reached a
   seam. There is **no separate "scheme limit"** -- I claimed one and it was wrong.
3. **THE FAILURE IS SILENT**: nothing aborts, no floor is reported, the run prints
   "Terminating on time limit" and every cell is NaN. Hence the NON-FINITE counter.

## THE FIX, and why the two paths differ

* **CELL-CENTRED: clamp UNCONDITIONALLY.** Both guards used on the fc side let the VERTEX
  failure through. The one that mattered is the monotonicity test: **a blast cap has a
  FLAT TOP, so a stencil astride its edge is NOT monotone**, and a guard written to protect
  smooth extrema skipped exactly the cells needing it. **NON-MONOTONE DOES NOT IMPLY
  SMOOTH.** Costs ~6% on rigidrot v_xi.
* **FACE-CENTRED: clamp only when MONOTONE and INTERPOLATING.** Unconditional costs 3.5x
  there and the seam case is already fixed by the guarded form.
* The asymmetry is load-bearing: cc carries density/pressure, where an inconsistent ghost
  goes straight to a negative state; fc carries B, with no such positivity constraint.

## THREE LIMITERS MEASURED AND REJECTED -- do not retry

| | why it failed |
|---|---|
| unconditional clamp on FC | 3.5x on the smooth halo; the resample sometimes EXTRAPOLATES, where leaving the node range is CORRECT |
| 2nd-difference test x relative span `(hi-lo) > 0.1*(|hi|+|lo|)` | fires 10^7 times on a SMOOTH run -- the span test breaks wherever the stencil straddles ZERO |
| van Leer, exactly as `PLM()` uses it | 4.3x worse -- a limited LINEAR slope discards the quadratic term the resample's accuracy depends on. The right limiter for a face value from cell averages is the WRONG shape for interpolating at a fractional offset |

Also rejected for the CORNER extrapolation: a monotone clamp (didn't fix the vertex, 15x
on the smooth halo) and a blow-up fallback (didn't fix it, 11x; span->0 makes it fire on
near-constant data). And the geometric cube-vertex fill did NOT fix the vertex either --
which is what finally pointed at the resample rather than the corner.

## STILL MISSING for shock work

**FOFC is a startup FATAL on the cubed sphere** ({hydro,mhd}_fofc.cpp have no gnomonic
form), so there is no first-order fallback to switch on.

**DECIDED 2026-09-02: FOFC IS NOT NEEDED. Do not re-propose it.** The user's call, after
the seam and conservation work landed. It had been written up here as "the highest-value
remaining piece"; that judgement is withdrawn. The blast test (iprob=12) is stable and
positive without it once the along-seam resample is clamped, which is what the fallback
would have been there for. The startup fatal stays -- refusing an unimplemented option is
the right behaviour, and nothing needs to change in the code. If a future problem does put
a strong shock on a cubed-sphere grid and the floors start firing, THAT is the trigger to
revisit, not the absence of the feature.

## TWO PROCESS TRAPS THIS COST

* **`make | grep | head -3` kills make with SIGPIPE** and sbatch runs anyway -- a GPU job
  reported PRE-fix numbers from a 4-hour-stale binary and nearly passed as a clean GPU
  verification. Never pipe a build through `head`; check the binary timestamp.
* A **CPU build accepts a device lambda reading file-scope HOST globals**; hipcc rejects
  it. iprob=12 shipped in da426248 not building for GPU at all. Copy into locals before
  the `par_for`, as CLAUDE.md says.
