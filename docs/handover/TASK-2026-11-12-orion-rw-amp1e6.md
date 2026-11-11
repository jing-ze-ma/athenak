# TASK 2026-11-12 viper -> Orion: IDORT linear-wave grid at amplitude 1e-6, N=256 tolerance check (CPU)

## Why
The thin-cell E order loss of the IDORT linear-wave (P_rad/P_gas, tau) grid at amp 1e-5 is an O(A^2) mean-E drift
(viper rc18_1110/RW_EORDER_1111.md). The grid is being rerun at amp 1e-6 (Jiang, Stone & Davis 2012), dt_retry off,
tight implicit tolerances, PROD vs FINAL. Viper runs the WHOLE grid (both arms, all 240 cases; only 12 thin tau 0.01
N128/N256 cases were still running at 00:00 CET 11-12, due ~01:00 CET). Orion gets ONLY the N=256 tolerance-convergence
check: the same N=256 cases with tolerances 10x tighter still (8 cases, none of them run on viper).

## Where each case runs (no duplicates)
| set | file (docs/handover/scripts/rw_amp1e6/) | where |
|---|---|---|
| grid PROD, 120 cases (24 cells incl. P 33.3, N 16..256, nt 2N), tol 1e-12 / lin 1e-13 / cnorm 1e-13 | cases_PROD.txt | viper (standby 12861994) |
| grid FINAL, 120 cases, same tolerances | cases_FINAL.txt | viper |
| tol check N64/N128 thin P 0.0333 tau 0.01, both arms, tol 1e-13 / lin 1e-14 / cnorm 1e-14 | (not in this dir) | viper |
| **tol check N256, both arms: P 0.0333 tau 0.01, P 0.333 tau 0.01, P 0.0333 tau 0.1, P 3.33 tau 10; tol 1e-13 / lin 1e-14 / cnorm 1e-14** | **orion_tt256_PROD.txt, orion_tt256_FINAL.txt** | **Orion** |

## Code and builds (two binaries, no PROBLEM: the test uses the built-in pgen rad_m1_beam)
- PROD arm: commit 7d233ac30 (fork/vgdhr-1110; contained in fork/final-1111). This is the binary viper's PROD runs use
  (athena_none_cpu_7d233ac3_vgdhr9); it reproduces the October PROD19 linear-wave runs exactly.
- FINAL arm: commit 84cd20888 (fork/final-1111 head).
Detached worktrees, e.g. /orion/ptmp/jinma/Athenak/wt_rwa6_prod and wt_rwa6_fin, then in each:
```
git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_none -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_none -j32
```
(viper's builds: Release, MPI on, no extra CXX flags; this is an accuracy comparison at fixed tolerance, not bitwise.)

## Runner
All in docs/handover/scripts/rw_amp1e6/ (copy the directory to a run area; it is self-contained, python3 + numpy):
run_radwave.py (+ radwave_disp.py, common.py) = viper physflux_1111/rw harness; drive.py = the same with
time/dt_retry forced to false (asserts it was true from the prod keys); prodkeys_all.txt = production keys (PROD19);
analyze.py = orders and L1 tables. Every case line already carries amp=1e-06 and the rm_ tolerance keys
(rm_<par>=<v> -> <rad_m1>/<par>, appended after the production keys, so they win); FINAL lines carry arm=TIME
(PROD19 + time2_k1_redist 0, implicit_hr_cmem, cmem_cap 1, time2_vet_x2_fields 15, implicit_physflux, physflux_mom fa)
+ implicit_physflux_zold=fc + implicit_hr_cmem_rst=true. One MPI rank per case, OMP 1:
```
cd <run area>/rw_amp1e6
export OMP_NUM_THREADS=1 RW_PKF=$PWD/prodkeys_all.txt
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_ttPROD  python3 drive.py run orion_tt256_PROD.txt  --np 4 > runs_ttPROD.out 2>&1 &
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_ttFINAL python3 drive.py run orion_tt256_FINAL.txt --np 4 > runs_ttFINAL.out 2>&1 &
wait
```
Cost: one N256 case ~1-1.5 h on one viper core (thin tau 0.01 slowest); 8 cases on 8 cores. (On viper the binary is run
through athwrap.sh only because it runs inside an srun step; not needed on Orion.)
Check before the batch: `python3 drive.py input orion_tt256_FINAL.txt | grep -E "dt_retry|implicit_tol|lin_tol|cnorm|radwave_amp"`
must show dt_retry = false, implicit_tol = 1.0e-13, implicit_lin_tol = 1.0e-14, implicit_lin_cnorm = 1.0e-14,
implicit_tol_rel = 0, radwave_amp = 1e-06 (last occurrence wins).

## Pass criteria
- every case ends "ok" in runs_tt*.out (no FATAL; with dt_retry off a FATAL is a result: report it, do not re-enable retry);
- tolerance-converged if L1(rho) and L1(E) (raw and mean-removed) at N256 differ from viper's grid-tolerance N256 values
  by < 1 % (viper's numbers will be in a NOTE / rc18_1110/RW_AMP1E6_1111.md; Orion-vs-viper build differences are
  also inside that 1 %; if not, rerun one case at the grid tolerance on Orion as the same-build reference).

## Return in a NOTE (NOTE-2026-11-1x-orion-viper-rw-amp1e6.md on this branch)
Per case: L1 rho, L1 E raw, L1 E mean-removed, mean E offset (units A|dE|), Picard/NC counters from the log tail,
wall. `RWA_ROOT=<run area>/rw_amp1e6 python3 analyze.py` prints these for runs_ttPROD / runs_ttFINAL (rows PROD_tt /
FINAL_tt; the PROD/FINAL/P19 rows are empty on Orion, viper fills them). Orders are formed on viper from the grid
(N 128 -> 256 raw and mean-removed); Orion reports the N256 values only.
