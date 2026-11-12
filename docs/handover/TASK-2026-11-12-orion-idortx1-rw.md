# TASK 2026-11-12 viper -> Orion: IDORT linear-wave grid at amp 1e-6 with the X1 radiation keys, arms FINAL and X1 (CPU)

## Why
The IDORT paper suite is being rerun on Raven CPU with the X1 key set (combo-c3-1112 559003d18 = combo-c2-1112 +
vet_gd_rotate_xfade; Raven tree /raven/ptmp/jinma/idort_x1_1112). Its linear-wave part (the amp-1e-6 (P_rad/P_gas, tau)
grid of TASK-2026-11-12-orion-rw-amp1e6 / -rw-edd / -rw-wbsolve) runs here only: 18 cells x N 16..256 x {edd, vetf},
arms FINAL and X1, both on ONE binary (559003d18), so X1 - FINAL is the effect of the X1 keys alone.

## Build (one binary, no PROBLEM: the test uses the built-in pgen rad_m1_beam)
combo-c3-1112 at 559003d18 (pushed to fork). Detached worktree, e.g. /orion/ptmp/jinma/Athenak/wt_x1_c3:
```
git -C /orion/ptmp/jinma/Athenak fetch fork
git -C /orion/ptmp/jinma/Athenak worktree add --detach /orion/ptmp/jinma/Athenak/wt_x1_c3 559003d18
cd /orion/ptmp/jinma/Athenak/wt_x1_c3 && git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_none -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_none -j32
```
(none_cpu = MPI off, Release; -ffp-contract=off as your rw-amp1e6 builds. The cases are one serial process each, no `np`.)
Report `md5sum build_none/src/athena` and `git -C wt_x1_c3 rev-parse HEAD` (must be 559003d18c19...).

## Lists (from the posted rw_wbs/all_FINAL.txt, so the cells, N, tolerances and the edd ny=4 fix are those of rw-wbsolve)
all_FINAL.txt = 180 cases: 18 cells (P 0.0333 / 0.333 / 3.33 x tau 0.01 / 0.1 / 1 / 10 / 100 / 1000) x N 16, 32, 64, 128, 256
(nt 2N, x1, transport implicit, hesdirk2, vimp) x edd (ny=4, ref=edd) and vetf (ref=vetqs_d); amp 1e-6, implicit_tol 1e-12,
lin_tol 1e-13, lin_cnorm 1e-13, tol_rel 0; dt_retry forced false by drive.py.
- FINAL = arm=TIME (PROD19 keys + time2_k1_redist 0, implicit_hr_cmem, cmem_cap 1, time2_vet_x2_fields 15, implicit_physflux,
  physflux_mom fa) + implicit_physflux_zold=fc + implicit_hr_cmem_rst=true (the lines as posted)
- X1 = FINAL + implicit_hr_cmem_fainit=true + time2_stage_local_be=true + time2_face_fct=true + time2_k1_redist=0
  (k1_redist 0 is already in TIME; named again so the X1 key set is explicit). The full X1 rad_m1 additions over TIME are
  `implicit_physflux_zold=fc implicit_hr_cmem_rst=true implicit_hr_cmem_fainit=true time2_stage_local_be=true time2_face_fct=true time2_k1_redist=0`.
  vet_gd_rotate_xfade is NOT added: it applies only to vet_gd with a physical-clock rotation; edd has no rays and vetf is vet_sc.
```
mkdir -p /orion/ptmp/jinma/rwx1_1112 && cd /orion/ptmp/jinma/rwx1_1112
cp -r <bsg-files-1009 checkout>/docs/handover/scripts/{rw_amp1e6,rw_edd,rw_wbs} .
cd rw_wbs
X="rm_implicit_hr_cmem_fainit=true rm_time2_stage_local_be=true rm_time2_face_fct=true rm_time2_k1_redist=0"
for c in edd vetf; do grep " $c " all_FINAL.txt > x1_FINAL_$c.txt; sed "s/\$/ $X/" x1_FINAL_$c.txt > x1_X1_$c.txt
  for a in FINAL X1; do awk '$4 == 16' x1_${a}_$c.txt > x1s_${a}_$c.txt; done; done
wc -l x1_*.txt x1s_*.txt          # 90 per x1_ list, 18 per x1s_ (smoke) list
```

## Pre-check, smoke, full
```
cd /orion/ptmp/jinma/rwx1_1112/rw_wbs
export OMP_NUM_THREADS=1 RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=/orion/ptmp/jinma/Athenak/wt_x1_c3/build_none/src/athena
python3 ../rw_amp1e6/drive.py input x1s_X1_vetf.txt | grep -E "^ *(implicit_physflux_zold|implicit_hr_cmem_rst|implicit_hr_cmem_fainit|time2_stage_local_be|time2_face_fct|time2_k1_redist|dt_retry|radwave_amp|implicit_tol) " | sort | uniq -c
#  must show zold fc, cmem_rst true, cmem_fainit true, stage_local_be true, face_fct true, k1_redist 0, dt_retry false,
#  radwave_amp 9.9999999999999995e-07 (= 1e-6), implicit_tol 1.0e-12, 18 lines each on the N16 smoke list (checked on viper 11-12); x1_FINAL_vetf.txt: no fainit / stage_local_be / face_fct lines
# smoke (N16, 18 cases per list; no FATAL expected):
for c in edd vetf; do for a in FINAL X1; do
  RW_RUNS=$PWD/runs_c3${a}_$c python3 ../rw_amp1e6/drive.py run x1s_${a}_$c.txt --np 18 > x1s_${a}_$c.out 2>&1; done; done
grep -c "ok w" x1s_*.out; grep -h -E "FAILED|FATAL" x1s_*.out
# full (N16 cached from the smoke). The heavy ones are vetf tau 0.01 N 128/256 (~5000 s per core at N256): give vetf more cores.
RW_RUNS=$PWD/runs_c3FINAL_vetf python3 ../rw_amp1e6/drive.py run x1_FINAL_vetf.txt --np 24 > x1_FINAL_vetf.out 2>&1 &
RW_RUNS=$PWD/runs_c3X1_vetf    python3 ../rw_amp1e6/drive.py run x1_X1_vetf.txt    --np 24 > x1_X1_vetf.out 2>&1 &
RW_RUNS=$PWD/runs_c3FINAL_edd  python3 ../rw_amp1e6/drive.py run x1_FINAL_edd.txt  --np 8  > x1_FINAL_edd.out 2>&1 &
RW_RUNS=$PWD/runs_c3X1_edd     python3 ../rw_amp1e6/drive.py run x1_X1_edd.txt     --np 8  > x1_X1_edd.out 2>&1 &
wait
```
A FATAL is a result: report it, do not add keys or re-enable dt_retry.

## Return (NOTE-2026-11-1x-orion-viper-idortx1-rw.md on this branch; NOTE only, no run trees)
1. Build: md5 + rev; hostname, cores; wall per list.
2. `grep -c "ok w" x1_*.out; grep -h -E "FAILED|FATAL" x1_*.out` and the "ok w=... ref=..." line of every N256 case.
3. Per-case X1 vs FINAL (line-aligned lists; ignore cmp_wbs's g^4 prediction column, it belongs to the wbsolve gate):
```
python3 cmp_wbs.py runs_c3FINAL_edd  x1_FINAL_edd.txt  runs_c3X1_edd  x1_X1_edd.txt  > cmp_x1_edd.md
python3 cmp_wbs.py runs_c3FINAL_vetf x1_FINAL_vetf.txt runs_c3X1_vetf x1_X1_vetf.txt > cmp_x1_vetf.md
```
4. Orders and N256 L1 with the existing analyzers, through symlink views (row "PROD" = FINAL arm, row "FINAL" = X1 arm, both
   on 559003d18; analyze.py's "FINAL/PROD L1 at N256" block is then X1/FINAL):
```
mkdir -p view_edd view_vetf
ln -sfn $PWD/runs_c3FINAL_edd view_edd/runs_PROD;   ln -sfn $PWD/runs_c3X1_edd view_edd/runs_FINAL
ln -sfn $PWD/runs_c3FINAL_vetf view_vetf/runs_PROD; ln -sfn $PWD/runs_c3X1_vetf view_vetf/runs_FINAL
RWE_ROOT=$PWD/view_edd  python3 ../rw_edd/analyze_edd.py > ana_x1_edd.md
RWA_ROOT=$PWD/view_vetf RWA_REF=/nonexistent python3 ../rw_amp1e6/analyze.py > ana_x1_vetf.md
```
   Paste cmp_x1_*.md and ana_x1_*.md whole. Flag every cell where X1 has a lower order (raw or mean-removed, 128 -> 256) or a
   larger N256 L1 (rho or E) than FINAL by more than 1 %.
