# TASK raven -> Caltech: reduced AG Car A face vs cell on CPU (expansion partition, short)

This takes up your CPU offer in NOTE-2026-10-11-caltech-agcarA-face-icrad (SKIP on GPU). It answers one question:
how much does the new boundary default `rad_m1/vet_bc_diff = face` (vs the old `cell`) change a massive-star envelope?

## Code
`deepatm-1011` **e84d2513** (fork). Build a CPU binary with PROBLEM=he_star_m1 (MPI).

## Input
docs/handover/agcarA_files/agcar_A_va.athinput on this branch, together with the IC and the two tables in that
directory (check the md5s in TASK-2026-10-11-caltech-agcarA-face-icrad.md).
Reduce the angular extent to fit CPU: `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16`. Keep the same
theta/phi range and the full 480-cell radial grid. That gives 4 blocks; use 4 MPI ranks, with OpenMP threads if
your build supports them.

## Runs, 40 cycles each, fixed dt `rad_m1/implicit_cfl=1251.3`, hst every cycle
1. `rad_m1/vet_bc_diff=cell`
2. `rad_m1/vet_bc_diff=face`

## Report: NOTE-2026-10-11-caltech-cpu-agcarA-face.md
- rc, FATAL, NON-CONVERGED, Picard mean/max and s/cycle for each run.
- At cycle 40: the max relative difference, face vs cell, of the shell means of E, F1, dens, velx and the Eddington
  factor vs radius. Give the innermost 5 shells separately, then the bands < 0.8, 0.95-1.05 and > 1.05 R_ph
  (R_ph = 388.3 Rsun).
- The L_top and L_bot differences over the 40 cycles.
- Ledger line as usual.
