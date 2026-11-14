# TASK 2026-11-14 viper -> Orion: grid-unification CPU batch (branch fork/unify-1114 = 7a3f4e23d) (USER GO 11-13 ~23:30 PST)

Context: unify-1114 makes the rad_m1 method the same on Cartesian / sp / cs (inventory: 27 historical differences). Steps 1-3 are done
(free fixes; sp/cs column-closure defaults on Cartesian: vet_col_order2, fk_min 0, time2_vet extrap2, surface_face; one face-limiter default
on all grids). Production decks are bitwise on viper GPU. This batch = CPU correctness only. Record hashes and md5s.

## Build (new worktrees; never touch main; no destructive git in existing worktrees)
Branch: `git worktree add --detach wt_u fork/unify-1114` (expect 7a3f4e23d). Base for the 'old' arms: `wt_p18b` at fork/prod-1118 = b2d852bb6.
Each: `cmake -B build_u -D PROBLEM= -D CMAKE_BUILD_TYPE=Release` (+ your usual MPI/compiler options) `&& make -C build_u -j` -> build_u/src/athena.

## Scripts: docs/handover/task-unify-1114/{cst.py,solib.py} (this commit)
solib.py imports from <repo>/tests_m1/runs_5q_sporder2b and vis/python: set HERE in solib.py to <repo>/tests_m1/runs_5q_sporder2b and W to your work dir;
in cst.py set the BIN paths: 'u0' = base binary, 'new' = branch binary.

## 1) Cartesian ladders (each runs its levels in parallel, 1 rank each)
```
for c in thick pulse_edd fsm; do for a in old new mlsf; do python3 cst.py run  $c $a 32 64 128 256 512 1024; done; done
for c in thick pulse_edd fsm; do for a in old new mlsf; do python3 cst.py eval $c $a 32 64 128 256 512 1024; done; done
for a in old new mlsf; do python3 cst.py run thick_T $a; python3 cst.py eval thick_T $a; done      # time-only ladder
```
Viper reference (to n=512/256): space+time E orders old 1.84/1.41/1.16 -> new 1.97/1.97/1.95; time-only old 1.68/1.55/1.45 -> new ~2.02.

## 2) Milne slab n=256 to steady state (deck tests_m1/runs_5e_vetcol2/inp/milne_vc.athinput; put the keys INTO the deck: fresh runs reject absent cmdline keys)
- old = base binary + `rad_m1/vet_col_surface_q=true rad_m1/time_scheme=be`
- new = branch binary + the same
- msf = branch binary + the same + `rad_m1/implicit_marshak_face=linear rad_m1/vet_col_surface_face=true`
Run: `athena -i deck_<arm>.athinput -d <dir> mesh/nx1=256 meshblock/nx1=256 time/nlim=12800`;
score: `python3 tests_m1/runs_5e_vetcol2/milne_v.py <root> <run dirs>`. Viper n=128: old = new L1 6.49e-5, msf 2.04e-5.

## 3) tst suite at the branch HEAD
`cd tst && python run_test_suite.py --cpu && python run_test_suite.py --mpicpu` (rad_m1 tests need ATHENAK_M1_DATA).
slab_cpu is EXPECTED to fail on both prod-1118 and the branch (stale REF table; leave it) - report its value. Report every other failure.

## Return: NOTE-2026-11-14-orion-viper-unify.md on bsg-files-1009, part by part. Times in PST. Lower priority than the open thincomb TASK.
