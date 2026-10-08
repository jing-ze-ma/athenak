# NOTE 2026-10-08 (orion): AG Car on Orion CPU, smoke B passes, throughput ~1/9 of a viper node

Answer to TASK-2026-10-08-orion-agcar.md. Report only; no production run started.

## 1. Binary and files
- Branch **he-ic-eint-from-t** at b40c09f1 (code = 2fbc6aa1, the viper GPU smoke binary's code; user asked
  for this branch rather than rt-integration). CPU build `build_agcar`: PROBLEM=he_star_m1, MPI on, Kokkos
  OpenMP on, Release, gcc 13.1 + openmpi 4.1.8.
- Binary `/orion/ptmp/jinma/agcar_1008/bin/athena_he_cpu_2fbc6aa1`, md5 **d2fefd7ef054eec01fa5346630a28f7f**.
- Files in `/orion/ptmp/jinma/agcar_1008/files`, `SETUP.sh` md5 check OK. Inputs unchanged.
- Job script `/orion/ptmp/jinma/agcar_1008/shake/shake_orion.sh` (`NB` = meshblock nx2 = nx3).

## 2. Layout: only the original 4 blocks (480 x 64 x 64) run
vet_gd needs a lateral band of **57 cells** (near-tangent rays) and fills its band only from the immediate
neighbour block, so meshblock nx2, nx3 must be >= 57:

| job | meshblock | ranks x threads | result |
|---|---|---|---|
| 213432 | 480 x 16 x 16 (64 blocks) | 16 x 7 | FATAL "clamped lateral reads (MeshBlocks too small for the reach)", then segfault, rc 139 |
| 213433 | 480 x 32 x 32 (16 blocks) | 16 x 7 | same FATAL, rc 139 |
| 213434 | 480 x 64 x 64 (4 blocks) | 4 x 28 | **rc 0** |

So the TASK's "split theta/phi on CPU" does not hold for vet_gd: at most 4 ranks, so in practice one node
(4 nodes x 1 rank x 112 threads is the only way to use more, not tested; pure OpenMP across both sockets
measured ~15 % slower than 16 x 7 on another problem). `vet_gd_allow_clamp = true` would run smaller blocks
but makes the near-tangent rays wrong: not used.

## 3. Smoke B (job 213434, 1 node p.exclusive, 2 x Xeon 8480+, 112 cores) vs viper 12130898

| | orion CPU | viper GPU | rel. diff |
|---|---|---|---|
| t after 10 cycles | 2.4186988199011507e+03 | 2.4186988199012844e+03 | 6e-14 |
| dt at cycle 10 | 2.4186976857005658e+02 | 2.4186976857007841e+02 | 9e-14 |
| mass (col 3) | 1.7891809043880330e+31 | 1.7891809043883109e+31 | 2e-13 |
| tot-E (col 7) | 1.0633773260125708e+46 | 1.0633773260127629e+46 | 2e-13 |
| IC column T(rho,eint)/T_col | 5.3e-14 | 6.4e-14 | |
| he_ic_balance cells max | 1.07898e-07 | 1.07898e-07 | |
| Picard mean / max | 6.9 / 24 | 6.9 / 24 | |
| rc / FATAL / NON-CONVERGED | 0 / 0 / 0 | 0 / 0 / 0 | |
| zone-cycles/cpu_second | 3.66e5 | 3.18e6 | |
| s/cycle (cpu time / 10, incl. setup) | ~21.5 | ~2.5 | |

**PASS.** Per core: 3.3e3 zone-cycles/s. Smoke A not run (same layout and per-cycle cost class; the
throughput answer does not change).

## 4. Projection (lower limits, dt falls once convection starts)
- B: >= 1,900 cycles x 21.5 s = **>= 11 h on 1 node** (~11 node-h).
- A: >= 6,500 cycles x 21.5 s = **>= 39 h on 1 node**, likely more (Picard mean 13.5 vs 6.9 for B on viper).
- Comparison: viper node 2.5 s/cycle, Raven node 1.4-1.55 s/cycle. One orion node ~ 1/9 of a viper node,
  consistent with memory bandwidth (~0.6 TB/s vs ~10 TB/s for 2 MI300A); CPU-specific tuning could buy
  maybe 10-30 %, not a factor of 9.
- Nodes: p.exclusive had idle nodes today (96/98 allocated), but the 4-block limit means extra nodes do not
  help this run.

## 5. Code change for CPU? (not done)
Smaller blocks need a multi-hop vet_gd lateral halo (band filled from blocks 2-4 away): a substantial rewrite
of rad_m1_vetgd.cpp, and more blocks add block-Jacobi lag at block boundaries (different convergence, not
comparable with the viper/Raven copies). Even with perfect scaling ~9 orion nodes = one viper node.

## 6. Recommendation
Orion is usable only as a 1-node backup for **B** (>= 11 h); keep A on viper/Raven. The user decides.
