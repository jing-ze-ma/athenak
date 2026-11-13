# TASK 2026-11-13 viper -> Orion: FULL IDORT suite (CPU) for D6 on the d6int-1113 build

The gated object is now **d6int** = the D6 keys on the build of **the `fork/d6int-1113` tip at the time of the paste; at least
62abe7ec1 (contains 35af2e3bc, the MPI-off build fix your CPU builds need)** (= xffuse-1113 + the opt-in radlz L_z key (default
off, the D6 keys do not set it) + the dhj speed-ups + the merge of the dhj/production line). Record `git rev-parse fork/d6int-1113`
(the commit you built) in the NOTE. The user (11-12) wants d6int to pass the FULL IDORT suite on Orion
(CPU, correctness only; CPU wall times are not results). It replaces the d6fast build of TASK-2026-11-12-orion-idortd6-UPDATE-d6fast
for the verdict; if d6fast runs are in flight, let them finish and report them as "d6fast", but the verdict that counts is d6int's.

Three parts: (a) linear-wave grid edd + vetf, (b) the 5 sp/cs curv cells, (c) the IDORT set that ran on Raven for d6fast
(153 rows: 41 Cartesian main rows + 12 rsx shock rows + 100 spx/xtra sp rows, + 7 OpenMP rows). For (a)/(b) the references are
your own FIN runs (wt_fin, 421a3d024). For (c) you have no FIN runs: (c) runs a FIN arm next to the D6 arm on the same node.

## Build (new worktree; never touch main; no destructive git in any existing worktree)
```
cd /orion/ptmp/jinma/Athenak && git fetch fork
git rev-parse fork/d6int-1113                            # record it; must contain 62abe7ec1: git merge-base --is-ancestor 62abe7ec1 fork/d6int-1113
git worktree add --detach wt_d6i fork/d6int-1113
git -C wt_d6i submodule update --init --recursive        # kokkos, in the new worktree only (as for wt_d6f)
```
1. `wt_d6i/build_none`: exactly the CMake options of `wt_fin/build_none` (copy from its CMakeCache.txt: CMAKE_CXX_COMPILER (gcc),
   Athena_ENABLE_MPI, Athena_ENABLE_OPENMP, CMAKE_BUILD_TYPE, CMAKE_CXX_FLAGS including `-ffp-contract=off`), no PROBLEM
   (every deck of this TASK uses the built-in `pgen_name = rad_m1_beam`; no he_star_m1 build is needed). `make -j`, record md5.
   EXE=/orion/ptmp/jinma/Athenak/wt_d6i/build_none/src/athena
2. OpenMP builds for the 7 OMP rows of part (c), only if `wt_fin/build_none` has Athena_ENABLE_OPENMP=OFF:
   `wt_d6i/build_omp` = the same options + `-D Athena_ENABLE_OPENMP=ON`, and `wt_fin/build_omp` the same on the existing wt_fin
   tree (a new build dir only, no git command in wt_fin). If build_none already has OpenMP ON, use the build_none binaries for
   EXE_OMP_*. Record both md5 and `athena -c | grep -i openmp`.
3. FIN binary for the part (c) FIN arm: the existing `/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena` (421a3d024,
   md5 e151dc3b...), the binary of your idortfin / idortd6 FIN references.

Fresh checkout of this branch for the decks: `F=<checkout>/docs/handover/scripts`.

## (a) Linear-wave grid, arm D6 on d6int (180 cases: edd + vetf)
Exactly part (1) of TASK-2026-11-12-orion-idortd6.md (1d5716074): same `fin_D6_<c>.txt` / `fins_D6_<c>.txt` lists (make them as
there if they do not exist; the `grep -c` check must give 0 0) and the same pre-check of one fins_D6_vetf deck, with
`RW_EXE=$EXE` and NEW run dirs tagged `d6i`: `runs_finD6i_<c>`, `fins_D6i_<c>.out`, `fin_D6i_<c>.out`. Smoke first (N16 lists,
`--np 18`), then the full lists (`--np 90`, vetf first). A FATAL is a result: report it, do not add keys or re-enable dt_retry.

## (b) The 5 sp/cs curv cells, arm D6 on d6int
Exactly part (2) of 1d5716074: decks `$F/idortfin_1112/curv/<cell>_D6/in.athinput`, cells T2cs_shdC, T2sp_shd, T1sp_star,
T3cs_hot5V, T3sp_hot5, ranks cs 6 / sp 1, `EXE=$EXE`, run dirs `$O/curv/<cell>_D6i`, output `ev_all_d6i.txt`
(`ATHENAK_VIS=/orion/ptmp/jinma/Athenak/wt_d6i/vis/python python3 $F/idortfin_1112/curv/ev.py */ > ev_all_d6i.txt`).
Smoke: one cs and one sp cell with `time/nlim=3` in a separate dir before the 5 runs.

## (c) The Raven IDORT set (D6 arm p on d6int + FIN arm f on wt_fin, same node, same decks except the method keys)
Everything is in `$F/idortd6int_1113/cart/` (C below):
- `decks_d6/<row>/in.athinput` and `decks_fin/<row>/in.athinput`: the PRE-RESOLVED p-arm decks (source deck + args.txt + PROD keys
  `tool/resolve/prodkeys_all.txt` + TIME keys `tool/resolve/timekeys.sh`, with CAND=d6 resp. CAND=fin), each with its
  `timekeys.txt` (keys applied) and `prod_drop.txt`. No args.txt / two.txt / command-line keys are needed (no deck has a two.txt).
  Verified on viper: all 41 d6 decks Raven ran (and the 2 OMP decks checked) are identical to Raven's d6fast decks, and 152 of 153
  fin decks are identical to Raven's FIN decks (idort_x1_1112), up to the shock-reference path. The one exception, Raven's
  cart/slabc/FIN3_n32 FIN row, was a stale X1 deck on Raven; the deck here is the correct FIN deck.
  The resolve tools are in `tool/resolve/` for provenance only (their paths point to the viper scratch).
- `refs/rsx/ref_m{1.05,1.2,2,3,5,27,50}.txt.gz` (rsx rows) and `refs/caltech/ref_m{1.2,2,3}.txt.gz` (rs_m* rows), md5 of the
  uncompressed files in `refs/md5_uncompressed.txt`. The decks read them at `/orion/ptmp/jinma/idortd6int_1113/refs/{rsx,caltech}/`.
- `lists/list_mainrk.txt` (41), `list_rsx.txt` (12), `list_spxx.txt` (100), `list_omp.txt` (7), `list_smoke.txt` (5):
  `<row> <cores> <est s>` (est = Raven seconds; sched.py starts the longest first). The Raven `_rev` lists had the same rows.
- `tool/sched.py` (Raven scheduler, unchanged), `tool/one_o.sh` (one row, arms p and f concurrently, serial binary, no srun),
  `tool/omp_o.sh` (the 7 OMP rows, 10 threads, both arms), `tool/cmp_o.sh` (rows + FLAGs + cmpout), `tool/cmpout.py`.
- `raven_ref/`: Raven summaries for comparison: `sum_runs_idort_x1_1112.txt` (Raven FIN, 153 rows),
  `sum_runs_idort_d6_1113.txt` (Raven d6fast, the 41 main rows; rsx/spxx did not run for d6fast on Raven),
  `sum_omp_*.txt` (7 OMP rows each), `cmp_d6fast_vs_fin.txt` (Raven jobs/cmp_d6.sh output; its section (3) failed on Raven).

```
W=/orion/ptmp/jinma/idortd6int_1113; C=$F/idortd6int_1113/cart; mkdir -p $W/refs
cp -r $C/refs/rsx $C/refs/caltech $W/refs/ && gunzip $W/refs/*/*.gz && (cd $W/refs && md5sum -c $C/refs/md5_uncompressed.txt)
export C W OMP_PROC_BIND=false EXE_P=$EXE EXE_F=/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena ARMS="p f"
export EXE_OMP_P=<wt_d6i omp or build_none binary> EXE_OMP_F=<wt_fin omp or build_none binary>
export BINCONV=/orion/ptmp/jinma/Athenak/wt_d6i/vis/python
# smoke (nlim 3, own dir): 5 rows x 2 arms + the 7 OMP rows x 2 arms; every row rc=0 fatal=0, timekeys.txt of every p deck has implicit_blend_dtc=1e-6
W=$W/smoke NLIM=3 python3 $C/tool/sched.py $C/lists/list_smoke.txt 40 $C/tool/one_o.sh > $W/smoke_main.out 2>&1
W=$W/smoke NLIM=3 bash $C/tool/omp_o.sh > $W/smoke_omp.out 2>&1
grep -c "p rc=0 fatal=0 .* f rc=0 fatal=0" $W/smoke_main.out; grep -c " rc=0 fatal=0" $W/smoke_omp.out   # must be 5 and 14
grep -L "implicit_blend_dtc=1e-6" $W/smoke/runs*/*/*/*/p/timekeys.txt $W/smoke/runs*/*/*/p/timekeys.txt    # must print nothing
# full set: one node per list (NC = cores of the node; each row takes 2 x <cores>), OMP rows on a node of their own
python3 $C/tool/sched.py $C/lists/list_mainrk.txt <NC> $C/tool/one_o.sh > $W/main.out 2>&1 &
python3 $C/tool/sched.py $C/lists/list_rsx.txt    <NC> $C/tool/one_o.sh > $W/rsx.out  2>&1 &
python3 $C/tool/sched.py $C/lists/list_spxx.txt   <NC> $C/tool/one_o.sh > $W/spxx.out 2>&1 &
bash $C/tool/omp_o.sh > $W/omp.out 2>&1 &
wait; bash $C/tool/cmp_o.sh $W > $W/cmp_d6i_vs_fin.txt
```
The longest Raven rows (pipe/FIN3_c220 ~ 8 h FIN, ddf/ddiff/ddn ~ 2 h serial, rsx n2048 / sp rsx ~ 1.5 h) set the wall time;
a re-run of a list skips rows with rc=0 (stale `.lock_<arm>` dirs of a killed run must be removed by hand).
Metrics: run your IDORT row analyzers (the ones of the wavefix-1013 / ridort-1014 / idortrs-1017 / sp dtree-1019 returns:
eq/eqts, shockcmp vs the refs, hcyl/shd/slab/cb/pipe/marsh, the sp evalav set, ladder od.py) on `runs/` and `runs_omp/` for both
arms (arm dirs `p` and `f`). If you have no analyzer for a row family, say so and use cmp_o.sh section (2) for it.

## Pass rules (each part)
PASS per row / case / cell: D6 rc=0, 0 FATAL, Picard NON-CONVERGED=0, and every metric <= 1.05 x the FIN value of the same row
(part (b): <= 1.05 x X1 of the same cell, and give D6/FIN per metric as in 1d5716074). FLAG anything else, including FLAG_PICARD
(D6 Picard mean > 1.25 x FIN) and any row where FIN is clean and D6 is not. A FATAL is a result: no added keys, no dt_retry.
Part (a) also keeps the 1d5716074 flags (lower order 128 -> 256 or N256 L1 rho/E > 1 % above FIN; P3.33 tau10 fix and tau0.01
flag; face_mmat limited-face count).

## Return: NOTE-2026-11-1x-orion-viper-idortd6int.md on bsg-files-1009 (NOTE only), post each part as it finishes
Header: wt_d6i HEAD = `git rev-parse fork/d6int-1113` at build time (must contain 62abe7ec1 and 35af2e3bc), md5 of every binary used (d6int build_none, omp builds, wt_fin), CMake options used.
(a) as part (1) of 1d5716074 with a d6int column (and d6fast beside it if you have it); PASS / FAIL.
(b) ev_all_d6i.txt whole + rc / FATAL / NON-CONVERGED of the 5 runs, FIN and X1 beside each D6 metric; PASS / FAIL.
(c) smoke counts; `cmp_d6i_vs_fin.txt` whole; the analyzer tables D6 vs FIN with every metric ratio > 1.05 FLAGged; the FLAG
list; your FIN arm vs `raven_ref/sum_runs_idort_x1_1112.txt` (rc / Picard mean per row, the expected agreement is round-off
level, gcc vs Raven gcc 13) and D6 vs `raven_ref/sum_runs_idort_d6_1113.txt` for the 41 main rows; PASS / FAIL.
Times in the NOTE in PST.
