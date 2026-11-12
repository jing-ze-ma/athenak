# TASK 2026-11-12 viper -> Orion: IDORT linear-wave grid at amp 1e-6, FINAL METHOD (FIN) vs FINAL, ALL CELLS AT ONCE (CPU)

**Supersedes TASK-2026-11-12-orion-idortfc-rw (ce9aa8076 + update 495210d48): do not run cmem (7bc19d109 / aea1cae35) or
c6e42d5f1.** The final method is decided (RSG verdict, 11-12). USER: make the final gates fast -> run every cell in parallel,
as many nodes as you can get (360 singleton cases; the slowest single case is vetf tau 0.01 N256, ~5000 s, so with >= 360 cores
the whole grid takes ~1.5 h). Same lists, drivers, analyzers and return format as 05e5a07c9; only binary, arms and parallelism change.

## Binary: combo-all-1112 421a3d024, MPI ON, run as singletons (as your NOTE 3ac548804)
```
git -C /orion/ptmp/jinma/Athenak fetch fork
git -C /orion/ptmp/jinma/Athenak worktree add --detach /orion/ptmp/jinma/Athenak/wt_fin 421a3d024
cd /orion/ptmp/jinma/Athenak/wt_fin && git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_none -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_none -j32
```
Report md5 and `git rev-parse HEAD` (421a3d024bb3...). Singletons: OMPI_MCA_hwloc_base_binding_policy=none, OMP_NUM_THREADS=1.

## Arms (both on the ONE 421a3d024 binary)
- FINAL = all_FINAL.txt as posted (TIME + implicit_physflux_zold=fc + implicit_hr_cmem_rst=true)
- FIN (final method) = TIME + time2_stage_local_be + time2_face_fct + time2_k1_redist=0 + implicit_face_mmat + implicit_cs2_wbsolve +
  implicit_am_cs_fold=applied, i.e. FINAL WITHOUT zold/cmem_rst. NOT in FIN: implicit_hr_cmem_fainit, implicit_hr_cmem_implicit.
  (sp/cs-only keys vet_col_lat_posfix=lag, implicit_cs_e0_metric are not added: Cartesian 1-D; vet_gd_rotate_xfade: no vet_gd here.)
```
mkdir -p /orion/ptmp/jinma/rwfin_1112 && cd /orion/ptmp/jinma/rwfin_1112
cp -r <bsg-files-1009 checkout>/docs/handover/scripts/{rw_amp1e6,rw_edd,rw_wbs} .
cd rw_wbs
X="rm_time2_stage_local_be=true rm_time2_face_fct=true rm_time2_k1_redist=0 rm_implicit_face_mmat=true rm_implicit_cs2_wbsolve=true rm_implicit_am_cs_fold=applied"
for c in edd vetf; do grep " $c " all_FINAL.txt > fin_FINAL_$c.txt
  sed -e "s/ rm_implicit_physflux_zold=fc//" -e "s/ rm_implicit_hr_cmem_rst=true//" -e "s/\$/ $X/" fin_FINAL_$c.txt > fin_FIN_$c.txt
  for a in FINAL FIN; do awk '$4 == 16' fin_${a}_$c.txt > fins_${a}_$c.txt; done; done
grep -c zold fin_FIN_*.txt     # must be 0
```
Pre-check (as before) on fins_FIN_vetf.txt: face_mmat true, cs2_wbsolve true, am_cs_fold applied, stage_local_be true, face_fct true,
k1_redist 0, NO implicit_physflux_zold / implicit_hr_cmem_rst / fainit lines, dt_retry false, radwave_amp 1e-6, implicit_tol 1e-12.

## Run: smoke (N16, 72 cases, one node, ~1 min), then EVERYTHING AT ONCE
Order inside each list: thick/stiff cells first (tau >= 10 and P 3.33; they found problems before), then the rest.
```
export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena
for c in edd vetf; do for a in FINAL FIN; do RW_RUNS=$PWD/runs_fin${a}_$c python3 ../rw_amp1e6/drive.py run fins_${a}_$c.txt --np 18 > fins_${a}_$c.out 2>&1 & done; done; wait
grep -c "ok w" fins_*.out; grep -h -E "FAILED|FATAL" fins_*.out
```
Full: one job (or node) per list or finer, --np = number of cases in the list (90) so no case waits; e.g. 4 nodes:
`RW_RUNS=$PWD/runs_fin${a}_$c python3 ../rw_amp1e6/drive.py run fin_${a}_$c.txt --np 90 > fin_${a}_$c.out 2>&1` for each of the
4 (a, c). If fewer cores: split the vetf lists by tau (one sub-list per job) and give vetf the cores first (edd lists take ~90 s).
A FATAL is a result: report it, do not add keys or re-enable dt_retry.

## Return (NOTE-2026-11-1x-orion-viper-idortfin-rw.md on this branch; NOTE only; post it AS SOON AS the grid is done)
As 05e5a07c9 with fc -> fin: build md5/rev, host(s)/cores, wall per list; `grep -c "ok w"` / FAILED / FATAL; N256 ok lines;
cmp_wbs.py FIN vs FINAL (edd, vetf); analyze_edd.py / analyze.py through view dirs (PROD = FINAL arm, FINAL = FIN arm), pasted whole.
Also the number of face_mmat-limited faces per FIN case if the log prints it (0 everywhere => FIN's face_mmat is inert on the waves).
Flag every cell where FIN has a lower order (raw or mean-removed, 128 -> 256) or a larger N256 L1 (rho or E) than FINAL by > 1 %.
