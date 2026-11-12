# TASK 2026-11-12 viper -> Orion: IDORT linear-wave grid at amp 1e-6 with the FINAL candidate FC, arms FINAL and FC (CPU)

**UPDATE (11-12, supersedes the binary/arms below; same lists, drivers, analyzers, return format):** 7bc19d109 (cmem v1) is
defective, do not run it. There are now TWO FC candidates, each on its own binary (both MPI ON, run as singletons):
- **FC-cmem** on cmemimpl-1112 **aea1cae35** (cmem v2): worktree wt_fcc, binary wt_fcc/build_none/src/athena; arms FINAL and FC
  exactly as below (`implicit_hr_cmem_implicit=true implicit_cs2_wbsolve=true` on top of X1). Run dirs runs_c4cFINAL_* / runs_c4cFC_*.
- **FC-fmm** on facemmat-trans-1112 **c6e42d5f1** (face M-matrix, no cmem_implicit; the two keys FATAL together in the cmem code):
  worktree wt_fcm, arms FINAL and FCM = FINAL + `rm_implicit_hr_cmem_fainit=true rm_time2_stage_local_be=true rm_time2_face_fct=true
  rm_time2_k1_redist=0 rm_implicit_face_mmat=true rm_implicit_cs2_wbsolve=true`. Run dirs runs_c4mFINAL_* / runs_c4mFCM_*.
Build each as below with the sha replaced (aea1cae35 / c6e42d5f1). Run the FC-cmem set first, then FC-fmm. Report per binary:
cmp FC vs FINAL and FCM vs FINAL (cmp_wbs.py + the analyzers through view dirs as below), plus FINAL(aea1cae35) vs FINAL(c6e42d5f1)
per case (max |dev|, should be round-off), and flag cells as below for both FC and FCM.


**Supersedes TASK-2026-11-12-orion-idortx1-rw (05e5a07c9).** USER 11-12: the IDORT suite runs the FINAL method candidate FC,
not X1. Your X1 grid (NOTE 3ac548804: 360/360, X1 keys inert, max dev 6e-6) stands as the X1 record; this TASK repeats the
same grid with FC. Same lists, drivers, analyzers and return format as 05e5a07c9; only the binary, the worktree, the run
dirs and the arm keys change.

## Build: cmemimpl-1112 at 7bc19d109 (= combo-c4-1112 07abb0957 + implicit_hr_cmem_implicit), MPI ON, run as singletons
As in your NOTE 3ac548804 (the none_cpu MPI-off build does not compile at 7bc19d109 either; the fix is 35af2e3bc on
fork/finalc-1112, not in 7bc19d109). Keep the binary at 7bc19d109 so it is the same code as the Raven suite build.
```
git -C /orion/ptmp/jinma/Athenak fetch fork
git -C /orion/ptmp/jinma/Athenak worktree add --detach /orion/ptmp/jinma/Athenak/wt_fc 7bc19d109
cd /orion/ptmp/jinma/Athenak/wt_fc && git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_none -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_none -j32
```
Report `md5sum build_none/src/athena` and `git -C wt_fc rev-parse HEAD` (must be 7bc19d1099e7...). Cases: one serial singleton
process each (no np), OMPI_MCA_hwloc_base_binding_policy=none, as in 3ac548804.

## Arms (both on the ONE 7bc19d109 binary, so FC - FINAL is the effect of the keys alone)
- FINAL = all_FINAL.txt lines as posted (unchanged from 05e5a07c9)
- FC = X1 + `implicit_hr_cmem_implicit=true implicit_cs2_wbsolve=true`, i.e. FINAL +
  `implicit_hr_cmem_fainit=true time2_stage_local_be=true time2_face_fct=true time2_k1_redist=0 implicit_hr_cmem_implicit=true implicit_cs2_wbsolve=true`
  (the cs/sp-only FC keys implicit_cs_e0_metric, vet_col_lat_posfix=lag, implicit_comoving_lat, implicit_cs_closure_metric
  are not added: these waves are Cartesian 1-D; vet_gd_rotate_xfade not added, as in X1)
```
mkdir -p /orion/ptmp/jinma/rwfc_1112 && cd /orion/ptmp/jinma/rwfc_1112
cp -r <bsg-files-1009 checkout>/docs/handover/scripts/{rw_amp1e6,rw_edd,rw_wbs} .
cd rw_wbs
X="rm_implicit_hr_cmem_fainit=true rm_time2_stage_local_be=true rm_time2_face_fct=true rm_time2_k1_redist=0 rm_implicit_hr_cmem_implicit=true rm_implicit_cs2_wbsolve=true"
for c in edd vetf; do grep " $c " all_FINAL.txt > fc_FINAL_$c.txt; sed "s/\$/ $X/" fc_FINAL_$c.txt > fc_FC_$c.txt
  for a in FINAL FC; do awk '$4 == 16' fc_${a}_$c.txt > fcs_${a}_$c.txt; done; done
wc -l fc_*.txt fcs_*.txt          # 90 per fc_ list, 18 per fcs_ (smoke) list
```

## Pre-check, smoke, full
```
cd /orion/ptmp/jinma/rwfc_1112/rw_wbs
export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=/orion/ptmp/jinma/Athenak/wt_fc/build_none/src/athena
python3 ../rw_amp1e6/drive.py input fcs_FC_vetf.txt | grep -E "^ *(implicit_hr_cmem_implicit|implicit_cs2_wbsolve|implicit_hr_cmem|implicit_physflux_zold|implicit_hr_cmem_rst|implicit_hr_cmem_fainit|time2_stage_local_be|time2_face_fct|time2_k1_redist|dt_retry|radwave_amp|implicit_tol) " | sort | uniq -c
#  must show cmem_implicit true, cs2_wbsolve true, implicit_hr_cmem true, zold fc, cmem_rst, fainit, stage_local_be, face_fct true,
#  k1_redist 0, dt_retry false, radwave_amp 1e-6, implicit_tol 1.0e-12 (18 lines each); fcs_FINAL_vetf.txt: no cmem_implicit /
#  wbsolve / fainit / stage_local_be / face_fct lines
for c in edd vetf; do for a in FINAL FC; do
  RW_RUNS=$PWD/runs_c4${a}_$c python3 ../rw_amp1e6/drive.py run fcs_${a}_$c.txt --np 18 > fcs_${a}_$c.out 2>&1; done; done
grep -c "ok w" fcs_*.out; grep -h -E "FAILED|FATAL" fcs_*.out
# FC run logs must print the cmem_implicit banner / not report the two keys as unread; check one:
grep -i -E "cmem_implicit|wbsolve|not read|unread" runs_c4FC_vetf/*/*.log 2>/dev/null | sort | uniq -c | head
RW_RUNS=$PWD/runs_c4FINAL_vetf python3 ../rw_amp1e6/drive.py run fc_FINAL_vetf.txt --np 24 > fc_FINAL_vetf.out 2>&1 &
RW_RUNS=$PWD/runs_c4FC_vetf    python3 ../rw_amp1e6/drive.py run fc_FC_vetf.txt    --np 24 > fc_FC_vetf.out 2>&1 &
RW_RUNS=$PWD/runs_c4FINAL_edd  python3 ../rw_amp1e6/drive.py run fc_FINAL_edd.txt  --np 8  > fc_FINAL_edd.out 2>&1 &
RW_RUNS=$PWD/runs_c4FC_edd     python3 ../rw_amp1e6/drive.py run fc_FC_edd.txt     --np 8  > fc_FC_edd.out 2>&1 &
wait
```
A FATAL is a result: report it, do not add keys or re-enable dt_retry.

## Return (NOTE-2026-11-1x-orion-viper-idortfc-rw.md on this branch; NOTE only, no run trees)
1. Build: md5 + rev; hostname, cores; wall per list.
2. `grep -c "ok w" fc_*.out; grep -h -E "FAILED|FATAL" fc_*.out` and the "ok w=... ref=..." line of every N256 case.
3. Per-case FC vs FINAL:
```
python3 cmp_wbs.py runs_c4FINAL_edd  fc_FINAL_edd.txt  runs_c4FC_edd  fc_FC_edd.txt  > cmp_fc_edd.md
python3 cmp_wbs.py runs_c4FINAL_vetf fc_FINAL_vetf.txt runs_c4FC_vetf fc_FC_vetf.txt > cmp_fc_vetf.md
```
4. Orders and N256 L1 (row "PROD" = FINAL arm, row "FINAL" = FC arm, both on 7bc19d109):
```
mkdir -p view_edd view_vetf
ln -sfn $PWD/runs_c4FINAL_edd view_edd/runs_PROD;   ln -sfn $PWD/runs_c4FC_edd view_edd/runs_FINAL
ln -sfn $PWD/runs_c4FINAL_vetf view_vetf/runs_PROD; ln -sfn $PWD/runs_c4FC_vetf view_vetf/runs_FINAL
RWE_ROOT=$PWD/view_edd  python3 ../rw_edd/analyze_edd.py > ana_fc_edd.md
RWA_ROOT=$PWD/view_vetf RWA_REF=/nonexistent python3 ../rw_amp1e6/analyze.py > ana_fc_vetf.md
```
   Paste cmp_fc_*.md and ana_fc_*.md whole. Also: the FINAL arm here vs your 559003d18 FINAL arm (runs_c3FINAL_*) per case
   (max |dev|; should be round-off: same keys, code differs only by opt-in keys). Flag every cell where FC has a lower order
   (raw or mean-removed, 128 -> 256) or a larger N256 L1 (rho or E) than FINAL by more than 1 %.
