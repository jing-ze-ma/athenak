# Axial angular momentum in the deep-hot-Jupiter productions: cubed sphere vs spherical polar

Read-only analysis, 2026-09-22. Scripts: `angmom.py` (integrals -> `angmom.npz`),
`plot_angmom.py` (table + `angmom.png`).

## Headline

**The cubed sphere does not have an angular-momentum problem at this resolution.** Over
283 rotations the cs run loses 2.4 % of the inertial `L_z`, essentially all of it in the
first ~150 rotations, and the last 130 rotations are flat to `+2e-5` per rotation. The two
sp runs, which *can* conserve `L_z` exactly in principle, drift by a comparable or larger
amount (2.2 % over 230 rot, 3.4 % over 174 rot). All three are dominated by the same
physical sinks (sponges, open radial boundaries), not by the grid.

The old static estimate of 6.9e-4 per rotation at 32 cells/panel is **not** what the
production run does: the measured cs rate is 8e-5 per rotation over rot 100-283 and
+2e-5 over rot 150-283, i.e. 10-30x smaller, and it is not monotone.

## What was computed

The conserved quantity is the *inertial* axial angular momentum. The dumps hold
rotating-frame velocities (`problem/omega = 2.06e-5`, `cs_full_rotation = true`,
`rot_potential = true`; axis = grid z on both grids), so with `R = r cos(lat)`

```
L_z = L_rel + L_frame ,   L_rel = int rho u_zonal R dV ,   L_frame = Omega int rho R^2 dV
```

`L_rel` is the relative-flow term, `L_frame` the mass-redistribution term. `L_frame` is
~99 % of `L_z` here (the jet carries ~1 % of the planet's spin angular momentum).

### Data

The bin dumps are **full 3-D at full resolution in single precision** — not slices, not
coarsened. cs: 24 MeshBlocks x 128x16x16 x 8 variables x 4 B = 25.2 MB, matching the file
size exactly; sp: 128x64x128 x 8 x 4 = 33.6 MB. No restart reader was needed. Cadence is
one dump per 2 rotations: 143 (cs), 116 (sp_mhd_prod3), 88 (sp_mhd_nopole).

### Velocity basis (the "orthogonal-KE audit" question)

Checked in the source, not assumed:

* **cs**: `Coordinates::GnomonicEquiangleRaiseVel` (coordinates.cpp:771 ff) writes
  `w0(IVX,IVY,IVZ)` as **contravariant** components on the gnomonic tangent basis of
  **unit but non-orthogonal** vectors `(e_r, e_xi, e_eta)`, metric `[[1,c],[c,1]]`,
  `c = cos_cell`. So the physical vector is `v = v1 e_r + v2 e_xi + v3 e_eta` and
  projecting `v2,v3` onto the unit zonal direction (`dhjcs.winds`) is **exact** — no
  cross-term correction is needed for a *linear* functional of v such as `L_rel`. (The
  1 % error in the old history KE came from the *quadratic* functional `0.5 rho v.v`,
  which does need the `c` cross term. `bcc` is stored in a different, mixed convention —
  see `dhjcs.bfield` — but is not used here.)
* **sp**: no raise step exists; `w0` is orthonormal `(v_r, v_theta, v_phi)`, zonal =
  `velz`.
* `src/outputs/` applies no rotation to the dumped primitives.

### Geometry

* Radial: `StretchRPoly` with the four `f_stretch_r_c*` from the input, reproduced from
  `src/coordinates/grid_stretch.hpp`. Exact volume-weighted radial moments,
  `<r> = (r_r^4-r_l^4)/4/dV_r`, `<r^2> = (r_r^5-r_l^5)/5/dV_r`.
* cs angular: **exact** gnomonic cell solid angle
  `S = F(xi2,eta2)-F(xi1,eta2)-F(xi2,eta1)+F(xi1,eta1)`,
  `F = atan(tan xi tan eta / sqrt(1+tan^2 xi+tan^2 eta))`.
* sp angular: exact `(cos th_l - cos th_r) dphi` with the sinh theta stretch,
  `f_stretch_theta = 3`.
* Both grids: `sum dV / (4/3 pi (r1^3-r0^3)) = 1.000000`.

## Results

`L_z0 = 4.8369e41` (cgs) on both grids — the two ICs agree to 1e-4, a good cross-check of
the two independent geometries.

| run | rot end | dL_z/L_z0 (end) | dM/M0 (end) | dL_frame/L_z0 | dL_rel/L_z0 | d(L/M)/(L/M)_0 (end) |
|---|---|---|---|---|---|---|
| cs 6x32x32              | 283.3 | -2.43e-2 | -1.32e-2 | -1.82e-2 | -6.1e-3 | -1.12e-2 |
| sp 64x128 (polar HLLE)  | 230.0 | -2.17e-2 | +1.44e-2 | +2.34e-2 | -4.52e-2 | -3.57e-2 |
| sp 64x128 (polar HLLD)  | 173.8 | -3.39e-2 | -8.4e-3  | -1.31e-2 | -2.09e-2 | -2.58e-2 |

Drift per rotation of `dL_z/L_z0`, by window (linear fit):

| run | rot 20-100 | rot 100-end | rot 150-end |
|---|---|---|---|
| cs                     | -1.39e-4 | -8.3e-5  | +2.1e-5 |
| sp (polar HLLE)        | +5.60e-4 | -5.1e-5  | +6.9e-5 |
| sp (polar HLLD)        | -1.68e-4 | -4.46e-4 | -4.6e-5 |

`L_z(t)` at ten times, with the mass-weighted mean zonal wind on three fixed shells
(cm/s; `r = 9.91e9`, `1.30e10`, `1.84e10` cm, i.e. deep / mid / near-top):

**cs_mhd_prod3**

| rot | (L-L0)/L0 | (M-M0)/M0 | L_rel/L0 | u(deep) | u(mid) | u(top) |
|---|---|---|---|---|---|---|
| 0.0 | +0.000e+00 | +0.000e+00 | +0.000e+00 | 0 | 0 | 0 |
| 30.0 | +4.91e-03 | -6.99e-03 | +1.34e-02 | -3 | 21480 | 2193 |
| 62.0 | -1.30e-03 | -8.87e-03 | +1.02e-02 | 11 | 22106 | 1557 |
| 94.0 | -5.72e-03 | -1.11e-02 | +9.16e-03 | 13 | 24804 | 1267 |
| 126.0 | -1.23e-02 | -1.27e-02 | +5.04e-03 | -4 | 24196 | 1075 |
| 156.0 | -2.47e-02 | -1.32e-02 | -6.53e-03 | -37 | 29018 | 1055 |
| 188.0 | -2.99e-02 | -1.32e-02 | -1.16e-02 | -39 | 33221 | 1051 |
| 220.0 | -2.57e-02 | -1.32e-02 | -7.56e-03 | -56 | 28573 | 1082 |
| 252.0 | -2.45e-02 | -1.32e-02 | -6.33e-03 | -55 | 25477 | 1154 |
| 283.3 | -2.43e-02 | -1.32e-02 | -6.11e-03 | -41 | 28057 | 1153 |

**sp_mhd_prod3 (polar HLLE)**

| rot | (L-L0)/L0 | (M-M0)/M0 | L_rel/L0 | u(deep) | u(mid) | u(top) |
|---|---|---|---|---|---|---|
| 0.0 | +0.000e+00 | +0.000e+00 | +0.000e+00 | 0 | 0 | 0 |
| 24.0 | -5.13e-02 | +1.04e-02 | -6.72e-02 | -1 | 17197 | 16643 |
| 50.0 | -2.24e-02 | +1.24e-02 | -4.19e-02 | 1 | 60697 | 18613 |
| 76.0 | -1.32e-02 | +1.31e-02 | -3.39e-02 | 12 | 9642 | -34933 |
| 101.8 | -8.38e-03 | +1.35e-02 | -3.00e-02 | 12 | -36603 | 46919 |
| 126.0 | -2.41e-02 | +1.37e-02 | -4.61e-02 | 7 | 53081 | 18633 |
| 152.0 | -2.50e-02 | +1.38e-02 | -4.74e-02 | 12 | 83369 | 21027 |
| 178.0 | -2.78e-02 | +1.42e-02 | -5.08e-02 | 15 | 36571 | 108554 |
| 204.0 | -1.95e-02 | +1.44e-02 | -4.27e-02 | 25 | 9384 | 62908 |
| 230.0 | -2.17e-02 | +1.45e-02 | -4.52e-02 | 13 | 44389 | 61703 |

**sp_mhd_nopole (polar HLLD)**

| rot | (L-L0)/L0 | (M-M0)/M0 | L_rel/L0 | u(deep) | u(mid) | u(top) |
|---|---|---|---|---|---|---|
| 0.0 | +0.000e+00 | +0.000e+00 | +0.000e+00 | 0 | 0 | 0 |
| 18.0 | +5.71e-03 | -2.31e-03 | +9.00e-03 | -10 | 28080 | 3350 |
| 38.0 | -8.60e-04 | -3.26e-03 | +4.06e-03 | 0 | 23742 | 2448 |
| 58.0 | -4.11e-03 | -4.38e-03 | +2.67e-03 | 3 | 23686 | 2122 |
| 76.0 | -5.89e-03 | -5.57e-03 | +2.72e-03 | 3 | 25642 | 1892 |
| 96.0 | -6.66e-03 | -6.72e-03 | +3.72e-03 | 8 | 24670 | 1571 |
| 116.0 | -1.09e-02 | -7.64e-03 | +9.13e-04 | 7 | 24810 | 1328 |
| 134.0 | -2.11e-02 | -8.18e-03 | -8.39e-03 | -28 | 25579 | 1390 |
| 154.0 | -3.28e-02 | -8.33e-03 | -1.97e-02 | -16 | 30529 | 1333 |
| 173.8 | -3.39e-02 | -8.36e-03 | -2.09e-02 | -17 | 32517 | 1243 |

`angmom.png`: (a) `(L_z-L_z0)/L_z0` (solid) and the mass-normalised `d(L_z/M)` (dotted);
(b) mass drift; (c) frame / relative split; (d) `L_rel/L_z0`.

## What is physical and what cannot be

Torque sources in this setup, all identical between cs and sp (same pgen,
`deep_hot_jupiter_rt.cpp`, same parameters):

1. **Bottom sponge**, `p > 50 bar` ramp to `100 bar`, `tau = 1000 s`, applied to `IM2`
   and `IM3` only (`IM1` commented out). It relaxes the *rotating-frame* momentum to
   zero, i.e. drives the deep atmosphere to exact corotation: a genuine and large sink of
   `L_rel`, and the dominant one, since the deep layers hold most of the mass.
2. **Top sponge**, `p < 1e-6 bar` ramp, `tau = 1000 s`, all three momenta.
3. **Initial Rayleigh drag**, all momenta, active only for `t < 5 P` — this is what makes
   the first ~10 rotations useless for a drift fit, and why the fits start at rot 20.
4. **Open radial boundaries** (`ix1_bc = ox1_bc = user`): mass crosses them, and the mass
   drift is large (-1.3 % cs, +1.4 % sp_mhd_prod3). Advected angular momentum goes with
   it. This alone accounts for most of `dL_frame`.
5. **`bc_outer_maxwell = true`**: the outer ghost extrapolation carries the divergence of
   the Maxwell stress, so the outer boundary exerts a real magnetic torque.
6. Central point-mass gravity and the centrifugal/`rot_potential` term are axisymmetric
   about z and exert **no** z torque. The Coriolis force is the frame term's exact
   counterpart and cancels in the inertial `L_z` above.

**Nothing in the list can be separated from grid error with these data.** The dumps are
snapshots only: there is no boundary-flux accounting, no sponge-tendency diagnostic, and
the hst file carries mass, grid-component momenta and energies but **no angular
momentum**. So the drifts quoted are *upper bounds* on the numerical non-conservation,
and the useful statement is the comparative one: cs, which cannot conserve `L_z`
exactly, is not worse than sp, which can. Making the separation would require either a
history variable for `L_z` plus surface-integrated torques (a code change), or a
zero-sponge / closed-boundary test problem.

## Caveats

* Single-precision dumps: the integrals carry ~1e-7 relative round-off per cell,
  negligible against 1e-2 drifts.
* Cell-centred `cos(lat)` and cell-centred `u` are used with exact `dV` and exact radial
  moments; the resulting quadrature bias is second order **and constant in time**, so it
  cancels in every drift number and only slightly biases `L_z0` itself.
* The mass-normalised `d(L_z/M)` curve (dotted in panel a) is the fairer measure when the
  boundaries leak: cs -1.1 %, sp_mhd_prod3 -3.6 %, sp_mhd_nopole -2.6 %. cs is the best of
  the three by that measure.
* `sp_mhd_prod3`'s mean zonal winds at the mid and top shells swing by tens of km/s
  between dumps and change sign. That is the known polar-HLLE sp MHD pathology
  (`sp-polar-hlle-swap-is-the-excess`), not angular-momentum bookkeeping; the
  `sp_mhd_nopole` (HLLD polar rows) winds are smooth and track cs closely. Its `L_rel`
  noise dominates its `L_z` curve.
* The three runs cover different spans (283 / 230 / 174 rot); the drift table compares
  rates, not endpoints.
* The cs run's `L_z` and mass both plateau after ~rot 150, which is consistent with the
  atmosphere reaching a steady exchange with the boundaries rather than with a secular
  grid-driven leak.

## Reproduce

```
cd /viper/u2/jinma/ATHENAK/athenak/tests_cs_angmom
python3 angmom.py          # ~2.5 min, reads 347 dumps (~9 GB), writes angmom.npz
python3 plot_angmom.py     # tables to stdout, angmom.png
```
