# NOTE 2026-10-21 viper + raven -> caltech: f55be371 (rank-ordered sums) reproducibility check: PASS on MI300A and A100

Re TASK-2026-10-21-caltech-viper-raven-detsum-check (11e4e0f3 + addendum 59d41b7a).

## Builds (each site's production recipe, incremental from the 4de8168d build dirs)
- Viper MI300A: gcc/16 rocm/7.1 openmpi_gpu/5.0, hipcc roc-7.1.0, Kokkos_ARCH_AMD_GFX942_APU, MPI, Release, PROBLEM=he_star_m1
  (same recipe as he_apu 4de8168d_r71; no extra CXX flags) -> athena_he_apu_f55be371_r71 md5 fa0efe0162bbf60b1c774fe3bf1168db.
- Raven A100: gcc/13 cuda/12.6 openmpi_gpu/5.0, nvcc_wrapper, AMPERE80, MPI, Release, PROBLEM=he_star_m1
  -> athena_he_a100_f55be371 md5 1c2f2a0990014e3c294031da27daded2.
- Current production binaries: viper 4de8168d_r71 (0a0d8746), raven 4de8168d (9af6efbb).

## Decks (deck A analogue, from each site's own production deck)
RSG E: production rsg_E_start1020(r)_clip_rst5d deck; HG: production hgglob_fresh_start1020(r)_clip deck. Changes only:
cs N=32 (mesh nx2=nx3=32), 24 blocks NX1x16x16, nlim 200, outputs by cycle (hst every cycle, rst + log every 100,
hydro_w/m1/m1_vet bins at 200). C, C2 fresh from t=0; R = restart from C's cycle-100 rst.

## Env (production env, collectives at the MPI default)
- Viper: HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1 OMP_NUM_THREADS=1 UCX_RNDV_THRESH=inf, srun --cpu-bind=cores; 2 nodes x 2 MI300A (4 ranks).
- Raven: KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1, srun --cpu-bind=cores; 2 nodes x 4 A100 (8 ranks).
- No HCOLL at either site: Open MPI 5 coll components are self libnbc inter [cuda] han ftagree adapt sync tuned monitoring basic sm.

## Results (all hst rows of hydro.hst + user.hst, and the final hydro_w/m1/m1_vet bins)
| site | deck | binary | C2==C | R==C | s/cycle |
|---|---|---|---|---|---|
| Viper | RSG E | f55be371 | IDENT (201 rows + bins) | IDENT (100 rows + bins) | 0.72 |
| Viper | RSG E | 4de8168d | IDENT (201 rows + bins) | - | 0.78 |
| Viper | HG    | f55be371 | IDENT (201 rows + bins) | IDENT (100 rows; bin data IDENT, header param padding only) | 0.92 |
| Raven | RSG E | f55be371 | IDENT (201 rows + bins) | IDENT (100 rows + bins) | 0.55 |
| Raven | RSG E | 4de8168d | IDENT (201 rows + bins) | - | 0.64 |
| Raven | HG    | f55be371 | IDENT (201 rows + bins) | IDENT (100 rows; bin data IDENT, header param padding only) | 0.72 |

- Viper and Raven were NOT affected with the current binaries: neither MPI stack has HCOLL, and Open MPI's own collectives
  (tuned/han) gave the same summation order run to run.
- f55be371 vs 4de8168d: round-off change as expected (first hst difference at cycle 2-4); f55be371 (including 63a7d28b's faster
  vet_gd exchange) is about 8 % faster per cycle on MI300A and 14 % on A100 for RSG E N=32 (separate jobs, same partition).
- The HG global deck (vet_col + vet_gd blend) is reproducible run to run and across a restart with f55be371 on both sites,
  so the unconverted MPI_SUM list in the addendum does not break reproducibility at either site with these collectives.
- Not done: production chains were not swapped (USER cancelled the global chains on viper and Raven; the accretor has priority
  there, globals run on Caltech).
Jobs: viper 12315481-88, raven 31378662-69. Dirs: /viper/ptmp2/jinma/detsum_1021, /raven/ptmp/jinma/detsum_1021.
