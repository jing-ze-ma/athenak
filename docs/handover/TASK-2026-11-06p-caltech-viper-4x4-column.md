# TASK 2026-11-06 Caltech -> viper (USER): `his1d_g` on exactly our column setup, 4 x 4 lateral cells

The USER asks you to run real Athena++ on exactly the setup of our AthenaK column `Lr` (`vstarfix_1106/col/Lr`), so that the geometry and angular setup match too: same mesh shape, not a true 1-D line.

## Why
With grid and dt matched, comoving L/L_* agrees to ~0.01-0.05 between `Lr` (AthenaK) and `his1d_g` at 0.5-1 d, at 8 radii from 20.16 to 26.41 Rsun. But KE and v_r differ: at 1 d, v_r at 22.66 Rsun is -42 m/s in ours and -122 m/s in yours, and at 3 d your KE is 3.7e40. `Lr` may or may not reach it; ~4 d results in about 30 min. Your line has 1-D angular geometry (nzeta 5, no npsi) and ours is a 3-D mesh. We want to remove that difference.

## Setup, as `Lr`
- **Domain:** spherical-polar wedge theta pi/3..2pi/3, phi 0..pi/3, **nx2 = 4, nx3 = 4**, periodic x2 and x3 as his deck, one uniform lateral grid at every radius. No refinement laterally.
- **Radial faces:** exactly his 256: 48 x 0.3125 Rsun (20-35), 128 x 0.15625 (35-55), 80 x 0.3125 (55-80). SMR would refine laterally too (8 x 8 at 35-55), so please use a user mesh generator or explicit faces for non-uniform x1 if it's cheap. If it isn't, use SMR with 8 x 8 at 35-55 and tell us; it's then not exact there, but the wall is the same.
- **Radiation:** his 3-D angular set: nmu 4, angle_flag 1, nzeta 5, npsi 6, taucell 4, error_limit 1e-6, nlimit default.
- **Everything else as `his1d_g`:** his units, original tables, GM 260.795, cfl_number 0.3, cfl_rad 0.025, floors, BC, **no seed** (laterally uniform IC, `problem/seed_amp` off). Report if lateral motion appears above round-off.
- Ours has radiation signal speed off in the CFL, so dt should be the same as yours (~10.6 s at the start). Note yours.

## Please report, to at least 8 d; to 25 d if cheap
Same columns as `his1d_g` so that we can compare directly:
- `*_ke.txt`: KE at hst cadence, plus KE_r and KE_lat separately.
- `*_wall_profiles.txt`: shell-mean first 24 cells, at 0.5 d steps to 5 d, then each 1 d to 10 d: r, rho, Pgas, Er(lab), T, vr, Fr1(lab), comoving Fr01 if you have it, N^2.
- dt at those times; Jacobi stats.
- Does it agree with `his1d_g`? In particular KE at 2-4 d and v_r at 1 d: this tells us whether the 1-D line vs 3-D angular mesh matters at all.

Short NOTE on this branch, first results as soon as they exist. Don't delay Orion.
