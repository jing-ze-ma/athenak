# NOTE 2026-10-17 (Raven): fin51-r-fix moved b16b66db -> 95fb04db

For: Caltech, Viper, DeltaAI (all sites building production from `fin51-r-fix`).

## What changed
`fin51-r-fix` was fast-forwarded (user decision, 2026-10-17) to **95fb04db** = b16b66db + one commit
(branch `vgd-latuni-1017`, 3 files, +60 lines: `src/rad_m1/rad_m1.cpp`, `rad_m1.hpp`, `rad_m1_vetlat.cpp`).

New startup FATAL in `VetLatInit` when **vet_gd** or **vet_col_lat** (lateral SC) runs on a mesh that is
not uniform in theta (x2) or phi (x3). Reason: the lateral SC upwind interpolation is bilinear in (theta, phi)
INDEX space and assumes uniform lateral spacing. Radial stretch (e.g. `use_grid_stretch_r_poly`) is fine.

Refused:
1. `mesh/use_grid_stretch_theta = true` (already refused for all M1 runs in `rad_m1.cpp`; repeated here).
2. `mesh/x2rat` or `mesh/x3rat` present in the deck and != 1. NOTE: AthenaK has no such keys; before this
   change they were silently ignored in a deck file (grid stayed uniform). Now they FATAL with vet_gd / vet_col_lat.
3. Built `xx2f`/`xx3f` active cell widths varying by > 1e-9 relative (1e-3 single precision), checked over all
   MeshBlocks and ranks (catches any future lateral map).

## Is it bitwise?
Yes, for every valid deck. Gates (Raven, /raven/ptmp/jinma/vgdlatuni_1017):
- A100, production decks bsg / e1 / hg (vet_gd + r-poly stretch), 4 GPUs, 10 cycles: every hst/bin/rst md5
  identical to b16b66db, Picard means identical (job 31245454).
- CPU, IDORT 20-dir subset (10 Cartesian, 10 sp vet_gd), 10 cycles: final restart bitwise (job 31245453).
- FATAL checks: x2rat, x3rat, both; vet_gd and vet_col_lat-only; CPU and A100 (jobs 31246202, 31246203).

## What to do
- Running production jobs: nothing. Do not rebuild or restart them for this change.
- Next build from `fin51-r-fix`: expect 95fb04db. If a deck now FATALs at startup naming x2rat / x3rat /
  use_grid_stretch_theta, remove that key (it never had an effect) or set it to 1.
- On CPU the FATAL prints after the `cycle=0` line (VetLatInit runs lazily in the first step), before any
  step completes or output is written.
