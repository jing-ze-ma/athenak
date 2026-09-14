---
name: red-giant-conduction-dt-collapse
description: "OPEN as of 2026-09-08: the red giant run loses its timestep to RADIATIVE CONDUCTION, twice; the RT source fix moved it from t=2.9e4 to 1.86e5 but did not remove it. Deepening the tau blend to 10/100 makes it NaN. RESUME HERE."
metadata:
  type: project
---

Continues [[red-giant-envelope-project]] and [[rt-source-semi-implicit]]. Runs live in
`/orion/ptmp/jinma/Athenak/red_giant/`.

## Where it stands

The `red_giant` cubed-sphere run (`inputs/hydro/red_giant_cs.athinput`) loses its
timestep to the **radiative conduction** dt, not hydro. Two events so far:

| | time | dt before -> after | RT clipping? |
| --- | --- | --- | --- |
| before the RT fix | 2.938e4 s | 32 -> 2e-3 s | yes, 144 cells |
| after 048dff30 | 1.859e5 s | 29 -> 0.70 s | **no** |

The semi-implicit RT source (048dff30) bought a factor 6 in time and turned a death
(2e-3 s) into a crawl (0.70 s), and the RT source limiter no longer fires at all. So the
stiff radiative source was real and is fixed, but it was not the whole story.

## The user's correction, which is right and which I had blurred

**The optically thin top of the atmosphere contributes NO conduction constraint.** The
kernel in `Conduction::NewTimeStep` returns early wherever the tau blend gives the cell's
faces zero weight (`if (wmax == 0.0 && blend_r) return;`). With `rad_tau_lo = 1` that is
everything above tau = 1. So the cell that sets the collapsing dt is at **tau > 1**, at or
below the photosphere -- it CANNOT be one of the cold/hot oscillating cells at i = 183-191
that I pointed at. My "top of the atmosphere" wording conflated the CAUSE (an RT
instability up there) with the LOCATION of the limiting cell, and the causal link between
them was never verified.

Commit 121acbdc exists so this stops being guesswork: the conduction dt now reduces with
`Kokkos::MinLoc` and `Mesh::NewTimeStep` prints the winning `(m,k,j,i)` on the line after
a collapse it attributes to conduction.

## The collapse is PROGRESSIVE, not a single event

Only ONE `dt COLLAPSE` line fires (the >4x-drop trigger), at cycle 6531, but dt keeps
degrading after it without ever dropping 4x in one cycle again:

    cycle  6531   t = 1.862e5   dt = 5.96    <- the reported collapse
    cycle  7300   t = 1.880e5   dt = 0.70
    cycle 49200   t = 1.917e5   dt = 0.072

Two decades over 43000 cycles. So whatever is happening at the photosphere is running
away slowly, and a fix that only stops the first drop will not be enough -- gate any
candidate fix on dt STILL being healthy at t = 2e5, not just on surviving 1.86e5.

## Why the tau=10/100 experiment failed, in hindsight

Moving the handover to tau 10-100 takes the conduction constraint OFF the i = 126-130
cells and hands that whole band to the two-stream alone. dt then jumps to the hydro limit
(~800 s), and tau ~ 1-10 is exactly where the explicit radiative source is stiffest. So
the two changes fight: the blend window is where conduction is expensive AND where the RT
needs the small timestep. Any fix has to address the cost, not move the boundary.

Directions not yet tried, in order of promise:
1. Check whether the FLUX LIMITER is actually biting at i ~ 127. At tau ~ 1 the flux is
   near free-streaming, so s = kappa|grad T|/(2 F_free) should be O(1) and the limiter
   should already be cutting kappa hard. If it is not, that is the bug.
2. Super-time-stepping for conduction. The code already has RKG STS for resistivity
   (`src/diffusion/resistivity.cpp`, the `*_rkg_*` task lists); conduction has none.
   This is the standard cure for a parabolic dt and would be reusable well beyond this run.
3. Accept dt ~ 6 s. That is 3.3e5 cycles for one convective crossing -- too slow.

## REFUTED: deepening the tau blend

`rad_tau_lo = 10, rad_tau_hi = 100` (run `relax_tau10`, job 193599) **NaNs by cycle 100**,
at t = 8.4e4 with 470400 non-finite cells. Removing the conduction constraint let dt jump
to the hydro limit (~800 s average over those cycles) and the scheme blew apart. So the
conduction dt is doing real work and cannot simply be pushed out of the way. Do not retry
this without also capping dt.

## Numbers so far

Flux gate (`tools/grid/flux_gate.py`), WALL run with the RT fix, net longwave at the top
over L/(4 pi r^2):

    t = 0        0.877   (implied Teff 3713 K)
    t = 2.4e4 s  0.702   (implied Teff 3512 K)

Falling because the atmosphere is still relaxing off its grey-Eddington initial condition
onto the correlated-k structure. **2.4e4 s is far too early to read** -- a convective
crossing of this domain is ~1.8e6 s. No open-boundary gate number exists yet.

## What is running / what to do first

- `loc_diag` (job 193600, 6 ranks, `build_rg_diag`, tlim 1.95e5): the instrumented rerun.
  **Read its `### dt COLLAPSE` line and the `conduction dt is set by cell` line under it.
  That names the cell. START THERE.** It reaches t = 1.86e5 about 400 s into the run.

  **ANSWERED 2026-09-08.** Both collapses are set by a cell INSIDE THE TAU BLEND WINDOW,
  just below the photosphere -- not deep, and not in the thin top either:

  | event | cell (m,k,j,i) | r | tau at t=0 | w at t=0 |
  | --- | --- | --- | --- | --- |
  | t = 362 s (the normal first drop, 363 -> 17.3 s) | (0,2,6,127) | 3.4176e12 | 3.31 | 0.53 |
  | t = 1.86e5 s (the collapse, 30.6 -> 5.96 s) | (0,16,3,130) | 3.4263e12 | 0.76 | 0.00 |

  tau = 2/3 sits at r = 3.4269e12, so i = 127 is one or two cells UNDER the photosphere
  and i = 130 is essentially AT it (its w was 0 at t = 0 and had risen above 0 by
  t = 1.86e5, i.e. the photosphere moved inward onto it).

  **TRAP: do not read depth off the radial index.** The fitted polynomial stretch packs
  cells at the photosphere, so i = 127 of 192 active cells is 66 % through in index but
  **95 % of the way out in radius**. I called it "deep, far below the photosphere" on the
  index alone and that was wrong. Always convert with the column dump's own r column.

