---
name: reconstruction-what-actually-runs
description: WHAT RECONSTRUCTION ACTUALLY RUNS (checked 2026-09-06 in mhd_fluxes.cpp / hydro_fluxes.cpp) -- on spherical polar ALL THREE directions are hard-wired to the stretched-grid PLM (GridPiecewiseLinearX1/X2/X3) and `<mhd>/reconstruct` is IGNORED; on the cubed sphere x1/x2/x3 honour `reconstruct`, but x1 treats the STRETCHED radial grid as uniform
metadata:
  type: project
---

    grid   x1 (radial)                          x2, x3
    sp     GridPiecewiseLinearX1 (non-uniform,   GridPiecewiseLinearX2/X3 (non-uniform PLM)
           x1v/x1f-aware, wb options) ALWAYS     ALWAYS -- `reconstruct` ignored
    cs     the selected method (plm/ppm4/ppmx/   the selected method
           wenoz) with UNIFORM spacing -- the
           radial stretch is NOT accounted for
           (mhd_fluxes.cpp:160-193: only sp
           takes the Grid* branch)

Consequences: (1) every sp claim about "reconstruct = X" was PLM; (2) the sp polar-row fix
had to be put INSIDE GridPiecewiseLinearX2 ([[sp-polar-row-reconstruction]]); (3) on cs with
use_grid_stretch_r_poly the radial reconstruction carries an O((dx_{i+1}-dx_i)/dx) error
per cell (~0.7 % stretch per cell in production) -- pre-existing, every cs production run has
it; the fix would be to route cs x1 through GridPiecewiseLinearX1 (caps radial at PLM) or a
non-uniform PPM. Flagged by the user 2026-09-06, NOT yet fixed. cs bcc DOES use the
stretch-aware interpolation (CellCenteredRadialFld) -- only the face-state reconstruction
does not.
