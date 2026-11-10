# TASK 2026-11-10 viper -> ORION: vgdh-geo -- more vet_gd directions (geodesic level 3, 4) on the cs / sp order ladders

Addendum to TASK-2026-11-10-viper-orion-vgdh.md (e9b71e840): same harness (raycheap18 st17), same binary, same CAND,
rotation OFF. Run AFTER the vgdh arms if those are still running. Answer as NOTE vgdh-geo on bsg-files-1009.

## Why (user 13:45 PST)
What worked on Cartesian (vet_sc, 512 rays, ~2nd order with the face-flux fixes) should carry over; the face fixes are
in CAND already. What differs on cs / sp is the closure: vet_gd with the 42-direction geodesic set (vet_gd_geo = 2). If
the space order goes toward 2 with more directions, the cap is the angular resolution and vet_gd_hr_latfit is a stopgap.

## Code
fork/vgdh-1110 92d4621fd (as in the vgdh TASK); `vet_gd_geo = n` is an existing key (geo-1031: geodesic level n,
10 n^2 + 2 directions: 2 -> 42, 3 -> 92, 4 -> 162; needs vet_gd_gl_nmu = 0, which CAND has). CPU MPI build line as in
NOTE-2026-11-08-orion-viper-consfix-ledger.md (gcc/13 openmpi/4.1 cmake/3.30, Kokkos SERIAL):
`cmake -B build -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build -j32`
(no PROBLEM: the pulse decks use the built-in generators, as for hr2cand). viper CPU smoke: CAND + geo4 pulse_cs C_32,
6 ranks, 3 cycles, rc 0, no FATAL, "geodesic level 4 (162 directions)" in the log.

## Arms (keys ADDED to CAND; no latfit)
- **CAND-geo3**: `vet_gd_geo = 3`
- **CAND-geo4**: `vet_gd_geo = 4`
- (optional, if 1-2 finish early) **CAND-geo4-lf2**: `vet_gd_geo = 4`, `vet_gd_hr_latfit = 2`

## Runs, in this order
1. pulse_cs C R3 n = 32 64 128 256: CAND-geo3, CAND-geo4 (viper runs 32..128 locally; the 256 level is what we need).
2. pulse_u full refinement C 64 128 256 512 (hr2prod spfull decks = CAND on sp): CAND-geo3, CAND-geo4.
3. Cost: wall per step (or rank-seconds per cycle) of CAND, CAND-geo3, CAND-geo4 at pulse_cs C_128 and pulse_u C_256,
   same node, same ranks.

## Expected table (fill in)
| arm | pulse_cs E p (32-64-128, 64-128-256) | F1 p | lateral spread per level | pulse_u E p (64..512) | F1 p | s/step |
|---|---|---|---|---|---|---|
| CAND (42) | 0.91 1.72 (Raven A100) | 1.11 1.56 | | 1.58 1.03 | 1.22 0.56 | |
| CAND-geo3 (92) | | | | | | |
| CAND-geo4 (162) | | | | | | |
Report rc / FATAL / NON-CONVERGED counts as in hr2cand. NOTE back (NOTE-2026-11-10-orion-vgdh-geo.md).
