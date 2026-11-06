# NOTE 2026-11-06 viper -> Caltech: cfl_rad is used, the wall film goes away at his dt; 3-D plan and cost; seed correction (TASK-2026-11-06c, NOTE-2026-11-06d)

**1. `cfl_rad` IS used, against the assumption in NOTE-2026-11-06d, and it changes the 1-D answer.** `src/hydro/new_blockdt.cpp:200` multiplies the implicit-radiation dt by `cfl_rad`, default 1. It is a plain multiplier on the hydro dt, not a light-speed limit. His deck has 0.025, so dt is ~5.3 s (effective hydro CFL 0.0075). **All our earlier 1-D and audit arms ran the default, dt ~200 s.** That includes A_af1, A_tc4, the ledger runs and NOTEs 11-06, 11-06b, 11-06c.

**2. So the one-cell wall film in NOTE-2026-11-06b comes from cfl_rad = 1, not from his setup.** With his cfl_rad it has not formed by 2-3 d, the time it is established at cfl_rad = 1. Interim; the cfl_rad = 0.025 arms are still running to 10 d (15 d if useful) and could still develop it later.

## 1-D arms (384 uniform, seeded; serial binary with his pgen)
`his1d` = his units and original tables, nzeta 5 (npsi is impossible in 1-D), cfl_rad 0.025. `his1d_cfl1` = the same with cfl_rad 1. `af1_cflrad` = earlier A_af1 deck (padded tables) with cfl_rad 0.025.
Faces are as in 11-06b; N^2 in 1e-10 s^-2, full-sphere dE, KE in the pi/3 sr wedge, Lc0 = cell-centred 4 pi r^2 F_lab in cell 0 (mostly advective).

| arm | t [d] | T/IC cells 0, 1 | rho/IC cell 1 | N^2 faces 0-1, 1-2, 2-3 | rms N^2 cells 0-23 | KE [erg] | dE [1e45] | Lc0/L_* |
|---|---|---|---|---|---|---|---|---|
| his1d | 0.5 | 1.020, 0.993 | 1.01 | -5, +4, +8 | 12 | 1.2e40 | -0.03 | -6 |
| his1d | 1 | 1.022, 0.991 | 1.02 | -10, +7, +5 | 11 | 5.2e39 | -0.10 | -3 |
| his1d | 2 | 1.024, 0.990 | 1.02 | -19, +10, -1 | 12 | 1.2e40 | -0.36 | -2 |
| his1d_cfl1 | 0.5 | 1.033, 0.994 | 1.01 | -35, -11, +1 | 17 | 9.5e39 | -0.17 | -50 |
| his1d_cfl1 | 1 | 1.039, 0.993 | 1.01 | -55, -25, -3 | 22 | 1.3e40 | -0.54 | -56 |
| his1d_cfl1 | 2 | 1.059, 0.987 | 1.04 | -126, -129, -22 | 70 | 1.1e41 | -2.65 | -124 |
| his1d_cfl1 | 3 | 1.107, 0.645 | 1.57 | -745, +140, +306 | 221 | 7.2e41 | -9.75 | -213 |
| his1d_cfl1 | 5 | 1.171, 0.385 | 2.52 | -1484, +1060, +376 | 554 | 7.0e41 | -25.67 | -75 |
| his1d_cfl1 | 10 | 1.271, 0.287 | 3.27 | -1938, -158, +400 | 528 | 9.7e41 | -64.20 | -218 |
| af1_cflrad | 1 | 1.045, 0.985 | 1.04 | -80, +34, -14 | 22 | 5.4e40 | -0.30 | -13 |
| af1_cflrad | 3 | 1.050, 0.985 | 1.04 | -92, +38, -42 | 27 | 2.9e40 | -0.76 | -12 |
| A_af1 (11-06b) | 3 | 1.156, 0.420 | 2.28 | -1319, +676, +522 | 405 | 1.3e42 | -15.14 | -142 |

Reading:
- Jacobi: his1d_cfl1 hits the 101-iteration cap in 7304 of 9062 reported solves, his1d in 0 of ~70700 so far.
- At 2 d his1d has about 1/7 of the cfl_rad = 1 energy loss, a wall flux ~60x smaller, and first-face N^2 of -19e-10. That last value is close to your layer's -16 to -20, but on one face of half your dr and with the next face positive. rms N^2 stays at the IC seed level.
- dt alone is not everything: af1_cflrad (padded tables, nzeta 4, floors differ) has N^2 ~4x larger than his1d at 1 d, though both are quiet.
- **Not settled:** whether his1d reaches a sudden inflow later, and 3-D.

## What this changes in earlier NOTEs
- **11-06b, 11-06c:** the film, the billing magnitudes and "the wall instability is a property of his discretisation" hold for cfl_rad = 1 only. Please do not read them as his run. The audit's C_Bd1/F_Bd1 "same pattern with zero billing" were also at dt ~200 s with a closed wall; convergence alone may not be the point, and dt may be.
- **11-06 (early KE):** "no kick" still needs checking at his dt (hst at 100 s not yet read). The KE and wall-drain numbers there are cfl_rad = 1. his1d KE at 0.5 d is 1.2e40, in line with the 1e40 level you quote from his Fig. 2.
- **11-05e and your toy:** if the toy or v* mimic targets Athena++ at dt ~200 s, it targets something he did not run. Your v* key should probably be compared at dt ~5 s or against his1d.

## Seed correction
**NOTE-2026-11-06 said "no density seed (1-D)". That is wrong: the pgen seeds in 1-D too.** At t = 0 rho scatters by +-2.6 % (rms 1.37 %) and P by 3.7e-4, matching amp 0.05 uniform. That also explains the "noisy IC" in 11-06b and the IC N^2 rms 22 here.
Unseeded 1-D (added `problem/seed_amp`, read only when named; absent = byte-identical dumps) at cfl_rad = 1 gives IC N^2 +9, +9, +9, matching your table, and the same film at about the same time and size (A_af1 first-face N^2 at 3 d -982 vs -1319; 5 d -1856 vs -1619). So at cfl_rad = 1 the seed does not cause the film. Not tested at cfl_rad = 0.025.

## 3-D plan (Raven CPU; viper is GPU only)
Directory `/raven/ptmp/jinma/athpp_repro3d`, run `prod_seeded`, SEEDED arm first as you asked.
- **Binary:** public Athena++ e47fabb7 (his exact commit is not available) + his `bsg_star.cpp`. One edit: 4 lines of `phydro->rad_source_term` -> `pnrrad->pradintegrator->rad_source`, because the former is not public. Diagnostic user output only (`pgen_changes.diff`). gcc/13, openmpi/5.0, hdf5-mpi, -O3, hllc, spherical_polar.
- **Inputs:** his deck `athinput.bsg_star.rusty` byte-identical, including the 0.2 code dump cadence (USER). His original tables from `BSG_Linhao.zip`. `MESA_profile_combined_cgs.txt` converted to his code units with the audit's formula (r=20: T 0.998, rho 0.993, consistent with "values at the bottom"). Effective settings: nmu 4, angle_flag 1, nzeta 5, npsi 6, **x2 periodic**, nlimit default 100, seed as coded.
- **Resources:** 32 nodes, 2304 ranks for 2560 16^3 blocks; two 24 h links queued (31975424, 31975425), restart chain. The scheduler's estimate for link 1 is ~25 min after this note; estimates have moved by many hours today.
- **Cost, scaled from a 1/16-lateral probe of this deck (72 ranks, 16 cycles, early times only):** 2.4e5 zone-cycles/s, i.e. 3350 per core-s; 0.44 GB/block. Full grid not tested before launch, it cannot fit on a dev node. ~78000 cycles per 5 d, so ~30 h per 5 d (68000 core-h) and ~150 h for 25 d if the rate holds. It may drop when iterations rise or at 1 block per rank.
- **Disk:** ~2.4 GB per dump, so ~0.9 TB per 5 d and ~4.5 TB for 25 d; ~20 GB per restart.
- **Not queued:** the NO SEED arm (optional per 11-06d). It would cost the same; say if you want it.

Results for the first 5 d will come as they arrive, in the TASK-2026-11-06c columns. 1-D summaries: `docs/handover/athenapp-bsg/cflrad_1106/`.
