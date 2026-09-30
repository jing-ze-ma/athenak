---
name: two-phase-relax-then-accurate
description: USER POLICY 09-26 - cheap schemes allowed during relaxation/spin-up; the steady-state (measured) phase uses accurate settings (round-off-neutral speed-ups only, cfl 0.9 kept)
metadata:
  type: feedback
---
User 09-26: "keep cfl 0.9, we can use lame schemes during relaxation but should use accurate schemes at the steady state."

- **Relaxation / spin-up phase.** Any cheap setting that is stable and does not bias the state the run relaxes to is allowed: multi-rate (implicit_mr_every), ck e8 + x16, vet_col_every N, looser Picard tol, cruder top treatment.
  - The check is not round-off. It is: switch to accurate settings at the end, and confirm that the diagnostics do not drift beyond a few thermal times of the top layers. If the deep T or flux drifts, the cheap setting biased the deep state.
- **Steady-state / measurement phase.** Only round-off-neutral speed-ups (see [[no-accuracy-sacrifice]]), with ONE exception: the step size. hesdirk2 at cfl 0.9 is kept, because its time error matches the previous production scheme (be at 0.3).
- The accuracy region rules still hold in both phases: [[accuracy-region-interior-only]] (massive stars) and [[resolution-only-below-1e-6-bar]] (dhj).

**How to apply:** production plans have two phases with explicit key sets and a switch-over check. Judge a speed-up test by "phase-1 eligible" (stable, relaxes to the same state) vs "phase-2 eligible" (round-off-neutral).
