# TASK 2026-11-09 viper -> ORION: raycheap2-1109 CPU gates + one convergence check (CPU only)

**Purpose:** correctness gates on CPU for the test-only key added on branch `raycheap2-1109`, plus one
self-convergence check of the physical-clock EMA closure at the production-like rotation cadence.
Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap.md` on this branch (xx = the day you finish).
Everything needed is on fork branches; no files from viper/Raven.

## 1. Code
- NEW = `fork/raycheap2-1109` head `a39c9ec49` (= `fork/prod-1112-cand` 3330bc0dd + 4 commits):
  - `367d2e8b3`/`f5eac602b` rad_m1/vet_rotate_avg_ema_dema (double EMA; opt-in) + rad_m1_beam sph_shell
    problem/flux_mod_amp|period (old test hook, not used here);
  - `063d6d96e` rad_m1/implicit_flux_x1min_mod_amp|period (test only, read only when named): the inner face
    flux F0 (1 + A sin(2 pi t/P)) set before each hesdirk2 stage solve at t_n + dt (Time2FormStage);
  - `a39c9ec49` tests_raycheap/ (decks + scripts used below; inputs only).
- BASE = `fork/prod-1112-cand` head `3330bc0dd`.
- Build both, CPU MPI, PROBLEM none, from git archives (any CPU flags; the same for both):
  ```
  git archive <sha> | tar -x -C src_<sha>; (kokkos submodule at the sha's gitlink)
  cmake -S src_<sha> -B b_<sha> -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off
  make -C b_<sha> -j
  ```
  (viper/Raven used gcc 13 + OpenMPI 5.0; CC/CXX must be gcc/g++, not the system gcc 7.)

## 2. Gates (all decks and scripts in `tests_raycheap/` of the NEW tree)
`bash tests_raycheap/run_gates.sh <NEW athena> <BASE athena> <outdir>` (uses `mpirun -np 1|4`; adapt the launcher
line in `run()` if Orion needs srun). Decks: cubed-sphere test star (cs lo, 32x48x48 per panel, 6 MeshBlocks), 40 cycles.

| gate | runs | PASS |
|---|---|---|
| G1 keys absent | NEW vs BASE on G_cs_lo_shd_geo2_k8e2 (production K8), ..._ravg_ema16 (EMA16), G_cs_lo_shd_ema16_h2 (hesdirk2) | `m1 bitwise 4/4`, `m1_lat bitwise 12/12` on all three |
| G2 path exercised | NEW: G_cs_lo_shd_ema16_h2_mod vs G_cs_lo_shd_ema16_h2 | must DIFFER (bitwise < all vars) |
| G3 restart | NEW mod deck: 20 + restart 20 vs 40 | bitwise 4/4 and 12/12 |
| G4 ranks | NEW mod deck: 4 ranks vs 1 rank | bitwise 4/4 and 12/12 |

Every run must end rc=0 with 0 FATAL (the script prints rc and FATAL counts).

## 3. Convergence check (CPU, ~1-2 h per deck at 24 ranks on viper/Raven CPUs)
Decks `tests_raycheap/Z_Q128T8_c{0.15,0.075,0.0375,0.01875}.athinput`: test star "hot" (24 MeshBlocks: run with
24 ranks), hesdirk2 + extrap2, inner flux F0 (1 + 0.1 sin(2 pi t/400)), t_end 400, EMA128 on the physical clock
(vet_gd_rotate_dt 21.44 = 8 coarse steps, align on), dt 1.34 / 0.67 / 0.335 / 0.167.
Run each with the NEW binary (`mpirun -np 24 athena -i <deck> -d <dir>`), then
`python3 tests_raycheap/conv.py <dir_c0.15> <dir_c0.075> <dir_c0.0375> <dir_c0.01875>`.
**Expected (viper/Raven CPU, same scheme at T_rot 10.72, EMA16/64):** d ~ 1e-4, 3e-5, 8e-6 and p -> 2.0-2.2 at the
finest pair. PASS = finest p >= 1.8 and all runs rc=0, 0 FATAL. Report d, p and t_end.

## 4. NOTE back
`docs/handover/NOTE-2026-11-xx-orion-raycheap.md`: build shas + compiler, the run_gates.sh summary lines (G1-G4),
the conv.py output, wall times, anything that did not run. Do not change code; if a gate fails, report the
first differing variable and stop.
