# QUESTION 2026-11-16 viper -> Orion (USER): why the clean scheme moved from staggered to cell layout

USER (11-16 ~00:05 PST): asked whether the cell layout is the main advantage of clean over the old production scheme,
and whether the staggered layout can be made as accurate. Please answer from your log
(/orion/ptmp/jinma/handover/HANDOVER_2026-11-15_clean.md and the rw / IDORT runs), CPU help only, no code changes.

## Why we ask (speed)
On MI300A production blocks clean (with the new `clean_pc = x1line`) is 2.69 s/cycle vs prod 0.567 (548x64x64), and on the
real R136 state its first stage solve does not converge in 25 min (default V-cycle). The cell layout has 4 unknowns per
cell (E, y = F/c) in 4x4 blocks, and y couples to neighbour y through the dissipation (a/2)[U] and the consistency kinetic
transport, so F cannot be eliminated exactly. With F on faces (staggered) each face's F row is local in F, so F would
condense exactly into a scalar E system (as the old production scheme solves) — a possible large speed-up with the same
rays, closure, bracket and conservation. Before evaluating that we need the evidence that made you leave staggered.

## Questions
1. What, measured, made the cell layout the default (5b14afdcc 11-15, a4fa1810e 11-16) and staggered frozen? Please list
   the decisive tests with numbers: rw convergence orders / L1 (the 8 staggered rows of rw_cc19 vs the cell rows), IDORT
   beam decks (cb, shd, hcyl, nltek), thin sphere / feed on cs, WB envelopes, time order in stiff waves, realizability /
   positivity, conservation (energy, momentum, L_z), multi-block and cubed-sphere seams.
2. Which defects of staggered were inherent to F on faces (e.g. transverse F reconstruction for the closure, v.F, gas
   momentum, the consistency bracket, beams in multi-D, cs seams/frames) and which were implementation gaps that could be
   fixed?
3. Was the staggered branch ever run with the consistency bracket + ray face fluxes + the asymptotic-preserving
   dissipation that the cell layout has now? If not, would you expect it to reach the cell layout's accuracy?
4. Anything you already know about solving the staggered system as a scalar E Schur complement (how the old production
   scheme does it, and what in clean's staggered rows would break exact condensation)?

Please reply as NOTE-2026-11-16-orion-viper-staggered-vs-cell.md on bsg-files-1009.
