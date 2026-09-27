---
name: cs-vertex-limiter-clipping
description: THE cube-vertex residual is LIMITER CLIPPING (2026-09-06) -- same root cause as the sp pole. WENO-Z cuts the vertex worst-cell force 6.5x (3.22e-2 -> 4.96e-3) and the 2x2 mean 3x; the Cartesian momentum update (d166238b, cs_cart_momentum, default OFF) does NOTHING at the vertex (6 %). On cs `reconstruct` IS honoured (unlike sp)
metadata:
  type: project
---

Chain of exclusions at the cube vertex (cs_test iprob=8, one step, |F| in (B^2/2)/r, 32/panel,
GS07 + WB, raw HLLD): exact ghosts (oracle) <1 %  [[cs-vertex-oracle-halo-innocent]]; the
face-sum WB source: vertex unchanged; the CARTESIAN momentum update (no geometric source at
all, exact closure -- null test round-off, hydro rigid rotation unchanged, interior -15 %):
vertex 2x2 7.06e-3 -> 6.66e-3, worst 3.22e-2 -> 3.17e-2. Then the reconstruction switch:

    recon    interior   edge     vertex 2x2   vertex max
    dc       1.42e-3    4.83e-3  1.17e-2      3.15e-2
    plm      2.09e-3    3.63e-3  7.06e-3      3.22e-2
    ppm4     1.66e-3    2.68e-3  4.26e-3      3.28e-2   (PPM still clips the extremum)
    wenoz    8.98e-4    1.17e-3  2.36e-3      4.96e-3   (6.5x on the worst cell)

So the vertex cell is an EXTREMUM of the component profiles in the extended chart and every
TVD limiter flattens it -- the twin of [[sp-polar-row-reconstruction]] (there the mirror
ghost makes every even profile an extremum; here the chart continuation does). Exact
ghosts cannot help because the extremum is real. The cure on cs is a non-clipping
reconstruction: WENO-Z (needs nghost = 3), or an extremum-preserving PPM. The vertex is
then ~2x the 2x2 mean and ~5x the interior mean, i.e. ordinary truncation.

Gates in flight (scratchpad sp_src / wb_speed, task bn9zvmpt2): iprob=9 convergence with
wenoz, the vertex MHD blast, the strat cost (nghost 3, WENO), and low-beta stability with
raw HLLD. If they pass, the production cs choice becomes reconstruct = wenoz.
The `cs_cart_momentum` kernel stays as a measured negative (keep off).

## ppmx IS the extremum-preserving PPM (Colella-Sekora `extrema` flag) -- gates 2026-09-06
One-step uniform field (iprob=8), GS07+WB, raw HLLD:
    n/panel  recon  interior   edge      vertex 2x2  vertex max
    32       plm    2.09e-3    3.63e-3   7.06e-3     3.22e-2
    32       ppmx   8.99e-4    1.17e-3   2.35e-3     4.94e-3
    64       plm    5.70e-4    1.23e-3   3.25e-3     1.47e-2
    64       ppmx   2.27e-4    3.98e-4   1.04e-3     2.36e-3
    EXACT ghosts with ppmx: identical to 1 % -> halo depth (ng=4) is NOT limiting.
iprob=9 convergence: ppmx ABSOLUTE errors 1.8x (beta 2) to 5.7x (beta 0.1, vertex) below
PLM at 64/panel in every region, but the measured ORDER drops (per step 2.6 -> 2.2-2.3,
vertex 0.9 at 32->64; finite-time 1.4 -> 0.9) -- same for wenoz. Reading: the higher-order
face states expose the floor set by the radial BC / RK2 error ([[cubed-sphere-mhd-convergence]]
"phase lag, one law in dt"), not a ppmx defect (uniform-field test converges 4x per doubling,
exact ghosts change nothing). Vertex blast with ppmx: clean (285 cycles, 0 non-finite).
Cost and low-beta stability: see the sweep (task bhgq64u5k) -- ppmx needs nghost = 4.

## 2026-09-06: the production ppmx arm DIED, then resubmitted angular-only
bench/cs_deep_ablate/ppmx on the pre-fix binary froze at dt = 0 at 0.105 rot with the
energy/temperature floors hit 26x more than the control: PPM ran in the RADIAL direction
on the stretched grid with the uniform stencil. Since 979edada ([[cs-radial-unification]])
x1 on cs is ALWAYS the position-aware Grid PLM and `reconstruct` governs the angular
sweeps only, so the resubmitted arm (job 11442987) is ppmx-angular / PLM-radial, the
combination the vertex fix actually needs. Control = rcmfix (same binary, plm).
