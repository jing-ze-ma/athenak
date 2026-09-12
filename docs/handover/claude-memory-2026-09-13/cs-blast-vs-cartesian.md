---
name: cs-blast-vs-cartesian
description: MEASURED -- a cubed-sphere panel seam imprints no more on a blast than a Cartesian grid's own anisotropy (4.1% vs 4.5%); only the cube VERTEX is genuinely worse, ~2x
metadata:
  type: project
---

**The question "is the cubed sphere good enough for shocks" now has a NUMBER, and the
reference is a Cartesian grid running the same blast.**

At matched shock strength (NOT matched time -- front width depends on strength), with
reconstruction, solver, gamma, CFL, pressure contrast, field strength and
cells-per-blast-radius all held equal:

```
run                            peak/ambient   front width   departure from one curve
Cartesian 3-D, 100^3               4.97        1.83 cells          4.6%
Cartesian 2-D, 200^2               5.28        1.81 cells          3.3%
cubed sphere, panel interior       5.29        1.60 cells          4.1%
cubed sphere, cube vertex          5.22        1.35 cells          9.7%
```

The last column is the rms departure from a single curve through the front as a fraction
of the shock jump. The blast is spherically symmetric but no grid is, so **a Cartesian
grid has its own version of this defect** -- that is the honest reference, and against it
a panel **seam costs nothing measurable**. Only the **cube vertex (~2x) is genuinely
worse**.

**The Cartesian 4.6% is stable, not one dump's fluke: 4.5% +/- 0.2% (range 4.2-4.8%)
across the 17 dumps that still hold a shock**, drifting down slowly as it weakens.
The 2-D value is LOWER than the 3-D (3.3% vs 4.6%) because a 2-D grid has fewer
directions in which to be anisotropic -- which is why the 3-D run is the right reference.

**Do NOT over-read the front-width column:** the Cartesian width wanders between 1.5 and
2.3 cells depending on where the front sits relative to the grid, so 1.35-1.83 means
COMPARABLE, not "the cubed sphere is sharper". And the geometries genuinely differ (3-D
box vs spherical shell), so this compares front sharpness and scatter, NOT the profiles
themselves -- a 3-D blast decays faster than a quasi-2-D one on a shell.

Artifact with the plots (3 panels, incl. all three fronts aligned and normalised):
https://claude.ai/code/artifact/fbe368bd-5960-40c7-a951-0cb8c27fe154
Cartesian runs used `PROBLEM=blast`; `athena_read` has NO `vtk` function, so the VTK
dumps need a hand-written parser. See [[cs-shocks-through-seams]], [[cs-mhd-blast]].
