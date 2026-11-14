# NOTE 2026-11-14 Orion -> viper: unify-1114 CPU batch (TASK-2026-11-14-orion-unify.md)

All times are PST, taken from sacct/date (the Orion clock is CET = PST+9).

## Summary
- **Part 1 (Cartesian ladders):** all 69 runs finished with rc 0, 0 FATAL and NON-CONVERGED=0.
  - new and mlsf converge at second order on every case: finest-pair E orders thick 2.01/2.01, pulse_edd 1.99/2.00, fsm 2.00/2.00. Time-only (thick_T) gives new = mlsf 2.05/2.04/2.03, close to your ~2.02.
  - old is below second order on thick (finest E 1.76) and fsm (1.56). pulse_edd old = new to every printed digit (Eddington closure, so the vet_col defaults do not act).
  - **My old numbers do not reproduce yours** (old 1.84/1.41/1.16 space+time, 1.68/1.55/1.45 time-only). Orion thick_T old gives 1.30/1.18/1.10. Details are in Part 1.
- **Part 2 (Milne n=256, steady):** old = new L1 2.044e-5, msf L1 1.023e-5 (your n=128: 6.49e-5 / 2.04e-5). msf halves L1, as it does at n=128.
- **Part 3 (tst suite):** the failure set is **identical on the branch and on prod-1118**, so the branch adds no failures.
  - cpu: 52 failed / 216 passed / 18 skipped on both.
  - mpicpu, run file by file: the same 16 failures and 4 timed-out files on both.
  - The rad_m1 tests (slab_cpu, opcheck_mpicpu, restart_mpicpu) all SKIP: the He-slab data set is not on Orion, so there is **no slab_cpu value**.

## Builds
| arm | commit | md5 of build_u/src/athena |
|---|---|---|
| branch (new, mlsf, msf) | origin/unify-1114 = 7a3f4e23d89533ad338fa19c39d50f311f9746af (as the TASK expects) | 9a97e4c038b434518155dcc632c03df8 |
| base (old) | origin/prod-1118 = b2d852bb6f5cffe45d9454bde042e86560c3b9a8 | 03a55c10aaecde8eab745396b7b24b81 |

- Sources: `git archive <commit>` into /orion/ptmp/jinma/unify_1114/{u,p}/src. Kokkos was copied in from src_vfix (git archive omits the submodule). No worktree was used and the main checkout was not touched.
- CMake: `-D Athena_ENABLE_MPI=ON -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off`, then `make -j32`.
- Modules: gcc/13 (13.1.0), openmpi/4.1 (4.1.8), cmake/3.30.
- Both binaries are MPI builds run as singletons (direct exec, 1 rank), which is the established Orion practice.
- The tst suite's own MPI-off build (`--cpu`) did compile on this branch.
- `-D PROBLEM=` (empty, as written in the TASK) fails CMake configure on Orion (src/CMakeLists.txt:282, "if given arguments: STREQUAL deep_hot_jupiter_rt"). I left it out, which gives the default PROBLEM=built_in_pgens.

## Part 1: Cartesian ladders
Script edits (copies in /orion/ptmp/jinma/unify_1114/scripts; the originals are from bsg-files-1009 3e0e3c16a):
- solib.py: `HERE = /orion/ptmp/jinma/unify_1114/u/src/tests_m1/runs_5q_sporder2b`, `W = /orion/ptmp/jinma/unify_1114/cart`.
- cst.py: `BIN['u0'] = .../p/build_u/src/athena` (prod-1118), `BIN['new'] = .../u/build_u/src/athena` (7a3f4e23d). Your cst.py had new = a49481fa_u2b.

There are no other edits. Levels were 32..1024 as in the TASK commands, and thick_T used the default levels 1, 2, 4, 8, 16 (NT = 128). All 9 ladders and 3 thick_T ladders ran at once on one exclusive node (orion060), each level a singleton.

Run status:
- All 69 runs: rc 0, FATAL 0.
- The "NONCONV 1" that cst.py prints for every run is the summary line `NON-CONVERGED=0.000000e+00`; all 69 logs carry the value 0.
- Wall time per ladder: thick_T 18-20 s, pulse_edd 185-204 s, thick 216-238 s, fsm 295-337 s.

Eval output, complete:
```
thick      old      E   t=1       L1 6.05e-03 1.50e-03 3.88e-04 1.06e-04 3.13e-05 | orders  2.01  1.95  1.87  1.76
thick      old      F1  t=1       L1 7.61e-03 1.95e-03 5.02e-04 1.48e-04 4.94e-05 | orders  1.96  1.96  1.76  1.58
thick      new      E   t=1       L1 6.31e-03 1.50e-03 3.60e-04 8.78e-05 2.17e-05 | orders  2.07  2.06  2.04  2.01
thick      new      F1  t=1       L1 7.74e-03 1.87e-03 4.50e-04 1.10e-04 2.73e-05 | orders  2.05  2.05  2.03  2.01
thick      mlsf     E   t=1       L1 6.31e-03 1.50e-03 3.60e-04 8.78e-05 2.18e-05 | orders  2.07  2.06  2.04  2.01
thick      mlsf     F1  t=1       L1 7.72e-03 1.87e-03 4.50e-04 1.10e-04 2.72e-05 | orders  2.05  2.05  2.03  2.01
pulse_edd  old      E   t=1       L1 4.79e-02 1.20e-02 2.98e-03 7.46e-04 1.88e-04 | orders  1.99  2.01  2.00  1.99
pulse_edd  old      F1  t=1       L1 3.69e-02 9.95e-03 2.54e-03 6.42e-04 1.64e-04 | orders  1.89  1.97  1.98  1.97
pulse_edd  new      E   t=1       L1 4.79e-02 1.20e-02 2.98e-03 7.46e-04 1.88e-04 | orders  1.99  2.01  2.00  1.99
pulse_edd  new      F1  t=1       L1 3.69e-02 9.95e-03 2.54e-03 6.42e-04 1.64e-04 | orders  1.89  1.97  1.98  1.97
pulse_edd  mlsf     E   t=1       L1 4.90e-02 1.22e-02 3.04e-03 7.61e-04 1.90e-04 | orders  2.00  2.01  2.00  2.00
pulse_edd  mlsf     F1  t=1       L1 3.61e-02 9.77e-03 2.49e-03 6.24e-04 1.56e-04 | orders  1.88  1.97  2.00  2.00
fsm        old      E   t=1.5     L1 6.34e-02 1.35e-02 3.66e-03 1.09e-03 3.71e-04 | orders  2.23  1.88  1.75  1.56
fsm        old      F1  t=1.5     L1 9.63e-02 1.20e-02 2.87e-03 1.05e-03 4.45e-04 | orders  3.01  2.06  1.46  1.24
fsm        new      E   t=1.5     L1 3.82e-02 9.85e-03 2.46e-03 6.14e-04 1.54e-04 | orders  1.96  2.00  2.00  2.00
fsm        new      F1  t=1.5     L1 3.45e-02 9.00e-03 2.27e-03 5.73e-04 1.44e-04 | orders  1.94  1.99  1.99  1.99
fsm        mlsf     E   t=1.5     L1 3.86e-02 9.91e-03 2.48e-03 6.19e-04 1.55e-04 | orders  1.96  2.00  2.00  2.00
fsm        mlsf     F1  t=1.5     L1 3.41e-02 8.94e-03 2.26e-03 5.68e-04 1.42e-04 | orders  1.93  1.98  1.99  2.00
thick_T    old      E   t=1       L1 3.45e-04 1.41e-04 6.23e-05 2.91e-05 | orders  1.30  1.18  1.10
thick_T    old      F1  t=1       L1 7.61e-04 3.08e-04 1.37e-04 6.40e-05 | orders  1.30  1.17  1.09
thick_T    new      E   t=1       L1 1.04e-04 2.52e-05 6.13e-06 1.51e-06 | orders  2.05  2.04  2.03
thick_T    new      F1  t=1       L1 3.51e-04 8.45e-05 2.06e-05 5.07e-06 | orders  2.05  2.04  2.02
thick_T    mlsf     E   t=1       L1 1.04e-04 2.52e-05 6.13e-06 1.51e-06 | orders  2.05  2.04  2.03
thick_T    mlsf     F1  t=1       L1 3.51e-04 8.45e-05 2.06e-05 5.07e-06 | orders  2.05  2.04  2.02
```

Against your reference:

| quantity | viper old -> new | Orion old / new / mlsf |
|---|---|---|
| space+time E orders, finest pair (512->1024) | 1.84/1.41/1.16 -> 1.97/1.97/1.95 | thick 1.76/2.01/2.01, pulse_edd 1.99/1.99/2.00, fsm 1.56/2.00/2.00 |
| space+time E orders, 256->512 pair | (same reference) | thick 1.87/2.04/2.04, pulse_edd 2.00/2.00/2.00, fsm 1.75/2.00/2.00 |
| time-only (thick_T) E orders | 1.68/1.55/1.45 -> ~2.02 | old 1.30/1.18/1.10, new 2.05/2.04/2.03, mlsf 2.05/2.04/2.03 |

- **new and mlsf:** they reproduce your "new" picture, second order everywhere.
- **old, space+time:** I cannot map your old triple 1.84/1.41/1.16 onto these ladders under any of these readings: one value per case (thick/pulse_edd/fsm), or a sequence of orders for one case.
  - pulse_edd old is identical to new on Orion, which rules out pulse_edd being one of three old values below 1.9.
  - Was the triple perhaps the `pulse` (vet_col) case, or a different level set?
- **old, time-only:** Orion gives 1.30/1.18/1.10 against your 1.68/1.55/1.45 from the same b2d852bb6 source. The differences from viper are the CPU build flags (-ffp-contract=off, MPI build run as a singleton) and the compiler.
- **Unresolved:** please send the ladder and command behind the old reference.

## Part 2: Milne slab n=256
- Decks are copies of tests_m1/runs_5e_vetcol2/inp/milne_vc.athinput.
- Keys absent from the deck were inserted with splitkeys.py:
  - every arm: `time_scheme = be`
  - msf only: `implicit_marshak_face = linear`, `vet_col_surface_face = true`
- `vet_col_surface_q` is already in the deck (false), so it went on the command line as `rad_m1/vet_col_surface_q=true`, together with `mesh/nx1=256 meshblock/nx1=256 time/nlim=12800`.
- Binaries: old = prod-1118, new and msf = branch.

Status:
- All three arms: rc 0, 0 FATAL, NON-CONVERGED=0, cycle 12800 reached.
- Wall time: old 461 s, new 519 s, msf 504 s.
- Steady state: the m1 tab at cycle 12800 equals the one at cycle 12600 exactly (max relative change 0) in all three arms.
- old and new profiles differ by at most 1.2e-9 relative in m1_e. Their final times differ (84.03 vs 83.57) because the dt sequences differ.

`milne_v.py /orion/ptmp/jinma/unify_1114/milne old new msf`:
```
# Hopf q(0) 0.577350 (1/sqrt3 0.577350), q(0.01) 0.588235 (0.588236), q(0.1) 0.627919 (0.627919), q(1) 0.698539 (0.697729), q(inf) 0.710446 (0.710446)
old      n=256 t=84.03  L1 2.044e-05 Linf 4.612e-04  top E/E_ref-1 +1.433e-03
new      n=256 t=83.57  L1 2.044e-05 Linf 4.612e-04  top E/E_ref-1 +1.433e-03
msf      n=256 t=83.56  L1 1.023e-05 Linf 6.283e-04  top E/E_ref-1 +2.515e-03
```
Against your n=128 (old = new 6.49e-5, msf 2.04e-5):
- old = new again, L1 down 3.2x from n=128 to n=256.
- msf L1 is half of old/new, against 3.2x lower at n=128.
- msf's Linf and its top-cell error are larger than old/new (6.28e-4 vs 4.61e-4, and +2.5e-3 vs +1.4e-3).

## Part 3: tst suite
Setup and how it ran:
- Run in the branch copy (u/src/tst), and for comparison also in the prod-1118 copy (p/src/tst).
- `ATHENAK_M1_DATA` was not set. A search of /orion/u/jinma and /orion/ptmp/jinma finds rosseland_he_x0.0_z0.02.txt and planck_he_x0.0_z0.02_ferg+tops.txt (e.g. /orion/ptmp/jinma/consfix_1109/inputs). It does **not** find ic_m1_V3edd_pgen.txt or arad_V3edd.txt, so m1_common.HAVE_DATA is False and all three rad_m1 tests SKIP. **slab_cpu: SKIPPED on both builds, no value.**
- With python-waterboa loaded, CMake picks the system GNU 7.5, which Kokkos rejects ("Compiler not supported by Kokkos for C++17"). I exported `CC`/`CXX` to the gcc/13 binaries.
- `--cpu`: `python run_test_suite.py --cpu` as written.
- `--mpicpu`: as written, the run hung in its first file (dyngrmhd lwave2d_amr) for more than 3 h, and job 229584 hit its time limit. So I built exactly as `--mpicpu` does (`testutils.clean_make(-D Athena_ENABLE_MPI=ON)`) and ran pytest on each `*_mpicpu.py` file under `timeout 1200`.

### --cpu
- Branch: 52 failed, 216 passed, 18 skipped (46 min).
- prod-1118: 52 failed, 216 passed, 18 skipped (46 min).
- **The two FAILED lists are identical.** The failures, with their messages:

| test | count | message |
|---|---|---|
| dyngrmhd nqt_shocktube | 15 | test_run: athena FATAL dyn_grmhd.cpp:92 "<mhd> dyn_eos = 'compose' not implemented for GR dynamics". test_convergence: KeyError ('mhd','hlle','wenoz',256) |
| dyngrmhd tab_shocktube | 15 | same as nqt_shocktube (compose EOS; the table inputs/tables/SFHo_reduced.athtab is also absent) |
| dyngrmhd lwave1d [mhd-ppmx], [mhd-wenoz] | 2 | 0-wave error 5.05082e-10 > 5e-10; 4.92708e-10 > 4.9e-10 |
| gr lwave1d [hydro-ppmx], [mhd-ppmx], [mhd-wenoz] | 3 | 4-wave 1.96998e-11 > 1.1e-11; 0-wave 5.05085e-10 > 5e-10; 4.92708e-10 > 4.9e-10 |
| sr lwave1d [hydro-ppmx], [mhd-plm], [mhd-ppmx], [mhd-wenoz] | 4 | 2.00643e-11 > 1.2e-11; 6.00064e-08 > 6e-08; 5.19201e-10 > 5.1e-10; 4.40788e-10 > 4.4e-10 |
| nr lwave1d [mhd-ppm4-rk2], [mhd-ppm4-rk3], [mhd-wenoz-rk3] | 3 | 1.77437e-08 > 1.5e-08; 1.13157e-08 > 7.4e-09; error ratio 0.138994 > 0.1 |
| mhd fofc_blast no_flagged_cell_is_bitwise_identical (cart_x1, cart_x3_3d, spherical_polar) x (llf, hlle) | 6 | "fofc=true changed the mass budget with no cell flagged" (e.g. 1.0000000000000016 vs 0.9999999999999933) |
| mhd fofc_cs no_flagged_cell_is_bitwise_identical [llf], [hlle] | 2 | same (13.613568165555678 vs 13.61356816555579) |
| mhd ohm_drag | 1 | RKG torque ratio 0.96616 vs explicit 0.97312 (diff 0.006958 > 0.006) |
| rad cs_implicit_ang | 1 | "sts: the rate error did not improve with resolution" (0.0028431445 vs 0.0028426528) |

### --mpicpu (file by file, 1200 s timeout)
The results are **identical on the branch and on prod-1118**.

| file | result |
|---|---|
| gr_bondi, hydro_fofc_blast (2), in_cshock2d (2), mhd_fofc_blast (2), mg_binary_gravity (4), mg_poisson3d (6) | pass |
| mg_jeans3d | 1 pass, 1 FAIL: test_jeans_amr_mpicpu "Jeans AMR MPI run failed" (rc 139, segfault) |
| nr lwave2d_amr | 8 FAIL: athena segfault in MeshBlockTree::AddNodeWithoutRefinement <- MeshRefinement::UpdateMeshBlockTree <- AdaptiveMeshRefinement (mpirun -np 16) |
| nr lwave3d_amr | 4 FAIL: the same AMR segfault |
| rad lwave2d_amr | 1 FAIL: segfault (np 8, AMR) |
| sbox hydroshwave, sbox mhdshwave | FAIL: the run reaches tlim, then "free(): invalid pointer", abort (signal 6) at exit |
| rad_m1 opcheck, rad_m1 restart, rad dhj_ck | SKIP (data absent) |
| dyngrmhd lwave2d_amr, gr lwave2d_amr, sr lwave2d_amr, z4c lwave2d_amr | TIMEOUT 1200 s with no output (testutils buffers until exit); probably the same AMR-with-MPI crash leaving mpirun hung, not checked |

The whole tst failure set is the same at b2d852bb6 and 7a3f4e23d. These failures are common to the fork on Orion: AMR with MPI segfaults, the compose EOS is absent, and lwave thresholds are exceeded by at most 2x. None of them is a rad_m1 failure.

## Jobs and wall times (PST)
| job | what | partition / node | start - end (PST) | state |
|---|---|---|---|---|
| 229567, 229570 | first submissions | - | - | cancelled by me before start (a wider job replaced the first; the second had the `-D PROBLEM=` error) |
| 229572 | both builds (203 / 212 s) + Part 1 + Part 2 + evals; its tst attempt failed on the GNU 7.5 compiler | p.exclusive, orion060 | 11-14 03:32:55 - 03:45:16 | COMPLETED (12 min) |
| 229572.0 | tst rerun as an overlap step | orion060 | 03:37:24 - 03:45:17 | killed when 229572 ended |
| 229584 | tst --cpu (branch, 3019 s, rc 1) + --mpicpu as written (hung in first file) | p.exclusive, orion060 | 03:45:39 - 07:45:50 | TIMEOUT 4 h |
| 229604 | branch --mpicpu file by file | p.exclusive, orion036 | 07:47:34 - 09:14:23 | COMPLETED |
| 229605 | prod-1118 --cpu + --mpicpu file by file | p.exclusive, orion060 | 07:47:34 - 10:04:05 | COMPLETED |

Queue wait for 229572 was about 3.8 h (submitted 11-13 23:45:49 PST).

## Deviations from the TASK
1. Sources come from `git archive` copies plus kokkos, not `git worktree add`. The builds are MPI-on with `-ffp-contract=off`, run as singletons.
2. `-D PROBLEM=` was dropped: it breaks CMake configure here.
3. The cst.py and solib.py edits are listed above (paths and BIN only).
4. The Milne keys are in the decks (splitkeys.py); only `vet_col_surface_q` went on the command line.
5. tst: `CC`/`CXX` were exported to gcc/13. `--mpicpu` was run file by file under `timeout 1200`, because the plain suite hangs on the AMR files. The suite was also run on prod-1118 for comparison. `ATHENAK_M1_DATA` was unset because the data is absent, so slab_cpu has no value.

Files are in /orion/ptmp/jinma/unify_1114:
- eval_ladders.txt, eval_thick_T.txt, run_status.txt, milne_eval.txt
- tst2_cpu.log (branch cpu), tst3_p_cpu.log, tst3_{u,p}_mpicpu_*.log
- job_all_229572.out, tst_full_{u,p}_*.out
