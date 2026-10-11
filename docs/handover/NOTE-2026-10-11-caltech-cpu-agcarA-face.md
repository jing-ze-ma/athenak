# NOTE Caltech -> raven: reduced AG Car A face vs cell on CPU (TASK-2026-10-11-caltech-cpu-agcarA-face)

Picked up by Caltech 2026-10-10 19:10 PDT; done 19:50 PDT.

## Setup
- Binary: CPU, PROBLEM=he_star_m1, MPI only (gcc 13.2 + openmpi 5.0.1, no OpenMP), fork deepatm-1011 **e84d2513**,
  `athena_cpu_he_star_m1_e84d2513` md5 **a51380cf40132c5e343655f707ac2275** (incremental build, expansion partition).
- Inputs: all four md5s match (deck 255f7a17..., IC 4ab9e3f3..., Rosseland b28e97c5..., Planck 493f9d01...).
- `vet_bc_diff` cannot be added from the command line (AthenaK only overrides keys that exist: FATAL "not found"), so each
  run uses a copy of the deck with `vet_bc_diff = cell|face` added to <rad_m1>; nothing else in the deck changed.
- Command line: `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16 time/nlim=40 rad_m1/implicit_cfl=1251.3`,
  table/IC paths, bin outputs at dcycle 40 (output6 switched to m1_vet for the Eddington factor m1_d11). 4 blocks of 480x16x16,
  4 MPI ranks (srun --mpi=pmix), 1 node of expansion (both runs on hpc-23-18 at the same time). dt = 1239.996 s, t_end = 49599.8 s.

## Runs
| run | job | rc | FATAL | NON-CONVERGED | Picard mean/max | s/cycle (main-loop cpu time/40) | wall |
|---|---|---|---|---|---|---|---|
| cell | 4333486 | 0 | 0 | 0 | 4.0 / 9 | 15.3 | 637 s |
| face | 4333487 | 0 | 0 | 0 | 4.0 / 9 | 15.3 | 635 s |
(3-cycle timing check first: 23 s/cycle incl. startup cycles; projection 16 min.)

## Face vs cell at cycle 40: max relative difference of the shell means (sin theta-weighted means over theta, phi)
rel = |face - cell| / |cell| per shell, max over the band (nrm = max|diff| / max|cell| in the band, for quantities near zero).
Bands in R_ph = 388.3 Rsun; the grid starts at 0.338 R_ph, so "< 0.8" has 139 shells, 0.95-1.05 has 86, > 1.05 has 170.

| quantity | innermost 5 shells (0.338-0.347 R_ph) | < 0.8 R_ph | 0.95-1.05 R_ph | > 1.05 R_ph |
|---|---|---|---|---|
| E | 0 (bitwise) | 0 (bitwise) | 8.9e-4 at 0.977 (nrm 1.3e-4) | 3.4e-7 |
| F1 | 0 (bitwise) | 3.6e-11 | 8.0e-7 | 3.4e-7 |
| dens | 0 (bitwise) | 0 (bitwise) | 2.3e-4 at 0.977 (nrm 6.9e-5) | 3.0e-9 |
| velx | <= 1.1e-10 | 3.3e-10 | 1.8e-2 at 0.976 (velx near 0 there; nrm 1.7e-4) | 5.5e-7 |
| f_rr (m1_d11) | 0 (bitwise; = 1/3) | 0 (bitwise) | 9.1e-4 at 0.977 (nrm 5.9e-4) | 2.0e-9 |

The same pattern holds at cycle 3: the difference is present from the first steps, all of it at ~0.977 R_ph, and does not grow
(E 8.90e-4 at cycle 3 vs 8.90e-4 at 40; dens grows from 1.4e-5 to 2.3e-4; velx 8.6e-2 -> 1.8e-2).

## Luminosities over the 40 cycles (hst, every cycle)
- L_bot: face = cell to round-off (max rel 6.4e-15). L_in: bitwise.
- L_top: max rel 5.0e-6 (cycle 3), 3.4e-7 at cycle 40. L_top cell goes 2.843095e38 -> 2.958946e38 (cycle 5 transient) -> 2.844125e38;
  face 2.844124e38 at cycle 40.

## Reading
On this reduced (32x32) AG Car A wedge over 40 fixed-dt cycles, face vs cell changes nothing at the inner boundary or in the
interior (bitwise or round-off below 0.8 R_ph, incl. the innermost 5 shells). The only change is a localized one near
0.977 R_ph (photosphere band): E and f_rr ~9e-4, dens 2e-4, and L_top 3e-7 to 5e-6. Why the switch shows up there and not at the wall was not investigated (f_rr = 1/3 exactly in the inner shells). Cost is identical (15.3 s/cycle both).

## Files (Caltech)
/resnick/groups/carnegie_poc/jingze/agcarA_face_1011: run.sh, files/agcar_A_va_{cell,face}.athinput, run_{cell,face}_n40/,
ana/cmp.py, logs/ana.4333595.out. CPU use: 3.4 CPU-h total (build 1.7, runs 1.5, tests and failed launches 0.2), 0 GPU-h.
