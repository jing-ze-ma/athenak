# clean-cc-pc progress (stronger preconditioner for the cell layout), 11-15 PST

Worktree /viper/ptmp2/jinma/clean_1114/wt_pc, branch clean-cc-pc off clean-cc-bstress 2a68179c7 (local, not pushed).
Scripts in /viper/ptmp2/jinma/clean_1114/pc: run_pc.sh (= suite/run_ladder.sh with a cycle limit), launch.sh ARM BIN
"keys" (deck set), stats.py ARM (iterations/stage, wall), cmp.py ARM REF (rel. L1 of the last m1 output),
gpu/time.sh (apudev timing, arms interleaved, repeats).  Runs: suite/<deck>/<ARM>_x<f>, pc/rw/<ARM>.
Common keys "scheme = clean;clean_layout = cell;clean_closure = extrap;clean_maxpic = 4".
Deck set: shd x0.25 (2 cycles), shd x0.5 (1), ws66 x0.5 (3), cbg x0.5 (3), nltek1 x1 (60), rw P3.3 tau10 N128 (1 period).

## Code (final 1f565b834; binaries athena_none_cpu_680023ee_pcf, athena_none_gpu72_680023ee_pcfg; 680023ee vs
## 1f565b83 = one comment line)
The cell-layout preconditioner is a multigrid V-cycle on the first-order 4x4 block operator:
- levels: halve every direction with an even number >= 4 of cells (single MeshBlock), down to where none halves;
- coarse blocks by aggregation (Galerkin with piecewise-constant prolongation, summing restriction): sum of the fine
  blocks, intra-aggregate couplings folded into the diagonal block (CleanCellCoarseSetup);
- smoother: the existing block line relaxation x1 -> x2 -> x3 (clean_pc_sweeps sweeps, default 1, pre and post), the
  mesh-level residual with the exact operator, coarse levels with their blocks; coarsest level 4x the sweeps;
- operator, BiCGStab and tolerances unchanged.  CleanCellLineSolve generalised to a CleanLevel (rad_m1_clean.hpp).
- CLEAN line now ends with "solves N" (linear solves of the stage).
No new key (clean_pc line|mg existed during the tests and was removed: mg won everywhere).
rad_m1_clean_cell.hpp untouched.

Rejected variants (ws66 x0.5 its/stage, PM1 = final 526): coarse-correction factor 1.5 -> 637, 2.0 -> 4791 (unconverged);
GMRES(30) + mg 873 GMRES its (= fewer PC applications than BiCGStab's 2x526, same wall); GMRES(30) + line PC 5032;
W-cycle 410 (ws66) / 406 (cbg, V 456) but shd x0.25 436 (V 413) and more work per iteration; lagged PLM jumps in the
blocks 502 on ws66 but 1242 on shd x0.25 and up to 178 per stage on nltek1 (V 22) -> rejected.
With clean_recon = none (operator = PC operator) mg gives ~25 its per solve, as with PLM (~24): the floor is the
aggregation V-cycle itself, not the reconstruction mismatch.

## Iterations (CPU, BiCGStab its per stage summed over the solves of the stage; PC0 = baseline 2a68179c_cc13x,
## clean_pc_sweeps 2; PF = final)
| deck | PC0 its/stage | PF its/stage | ratio |
|---|---|---|---|
| shd x0.25 | 2070 | 413 | 5.0 |
| shd x0.5 | 6655 (1 stage) | 566 | 11.8 |
| ws66 x0.5 | 1988 (22 solves/stage: ~94 per solve) | 526 (~24 per solve) | 3.8 |
| cbg x0.5 | 1937 | 456 | 4.2 |
| nltek1 x1 | 3.0 (max 49) | 1.5 (max 22) | 2 |
| rw | 3.0 | 3.0 | 1 |

## GPU wall (apudev, binary athena_none_gpu72_5a12bcf5_pc4g with both arms, interleaved, 2 repeats; jobs
## 13166281-4, 13166509; final binary check 13166520: ws66 41 s = GV)
| deck | line (GL) s | mg (GV) s | speed-up |
|---|---|---|---|
| shd x0.25 | 76, 75 | 30, 30 | 2.5 |
| shd x0.5 | BiCGStab breakdown at it 177 in the first stage, both repeats (stage unconverged, step accepted) | 30, 30 | - |
| ws66 x0.5 | 79, 78 | 41, 41 | 1.9 |
| cbg x0.5 | 75, 75 | 36, 36 | 2.1 |
| nltek1 x1 | 209, 209 | 208, 207 | 1.0 |
(rw not timed on GPU.)

## Agreement with the baseline
clean_tol 1e-10 (PF/PM1 vs PC0): max rel L1 over m1 variables shd x0.25 6.4e-9, ws66 1.7e-9, cbg 3.5e-9 (worst F2);
nltek1 E 5.0e-7, F 2.7e-6; rw E 8e-16 (1.1e-8 of the wave amplitude; rw runs at clean_tol 1e-14).
clean_tol 1e-13 (PT1 vs PT0, same binary): shd x0.25 1.3e-12, ws66 1.3e-12, cbg 3.1e-12; nltek1 E 5.2e-10, F 2.0e-9.
The difference scales with the solver tolerance (operator unchanged); same marks/negatives/cycles/dt in every run.
nltek1 vs exact transport: PC0 9.4796e-5, PF 9.4806e-5; at clean_tol 1e-13 both 9.4110e-5.
Final binary PF vs the test binary's mg arm (PM1): bitwise identical on all six runs.
