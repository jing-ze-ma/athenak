# NOTE viper -> caltech (USER): CORRECTION to de4be1b63 (vet_gd theta wall): the accretor walls are CLOSED hse walls; mirror is the right model, vacuum is wrong

This corrects two statements in NOTE-2026-10-27-viper-caltech-vetgd-theta-wall.md. Everything else there stands: the bug, its size, who is affected, and the read-only check. Still nothing to change in the running jobs.

## What was wrong
1. **"The theta walls at ±20 deg are FLUX/outflow walls."** Wrong. In the run12 deck the gas walls are closed hydrostatic walls (`ix2_bc = ox2_bc = user`, pgen `theta_bc = hse`). The M1 x2 faces are `implicit_bc_x2min/max = flux`, which imposes the star's equilibrium lateral flux. Beyond the wall there is more star, not empty space.
2. **"Do NOT use `mirror` for the accretor; `vacuum` in progress."** Reversed. Mirror is the consistent model for these walls, and it is also what the ray reference raynz/raycheck assumes: it folds theta about the wall (`vetgd_acc_1021/jobs/raycheck.py`, `fold()`).

## Evidence (viper, X6 restart, cycles 682 -> 882, production keys; `fork/wallbug-1027` = 917a084cf)
- **`vet_gd_theta_wall = vacuum`** (new, opt-in) matches an exact long-characteristics reference to 0.5 % for a wall open to vacuum. On the accretor it gives E in the wall rows, outer thin region, **+959 %** against periodic: the dark sky disagrees with the imposed wall flux, and radiation piles up. Not usable here.
- **`extrap`**: -95 % in the same rows, 5 stage fallbacks, +14-34 % cost. Not usable here.
- **`mirror`**: the same rows change by about 5 % (median) against periodic.
- Why the old periodic wrap is close but not right:
  - The wedge is equatorially symmetric (x2min + x2max = pi), so the periodic wrap reads the mirrored position.
  - But it keeps the ray direction instead of reflecting it (z -> -z), and the stream breaks the symmetry.

## Open (USER decision, not yet made)
- With the production GL10x12 staggered direction set, mirror is approximate: the error near the phi walls grows with resolution (p90 0.10 -> 0.43 on the CPU test). It is exact only with an odd `vet_gd_gl_nmu` (e.g. 11) or `vet_gd_gl_stagger = false`.
- Whether run12 switches to mirror, and with which direction set, comes in a later TASK. No binary or deck change until then.
