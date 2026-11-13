# TASK 2026-11-13 viper -> Orion: FULL IDORT suite (CPU) again, now for prod-1118 (USER 11-13 ~09:30 PST)

The new production build is **prod-1118 = `fork/prod-1118` = 1d1972805** (d6int 62abe7ec1-line + the L_z fix `implicit_lzx`
(lzphys-1113 abdf9a97d) + the HG time-step fix (hgdt-d6i-1113) + the Zalesak round cap). Production runs it with the D6 keys plus:
```
EXTRA="rad_m1/zalesak_cap_rounds=4 rad_m1/implicit_lzx=true rad_m1/time2_stage_local_be=false"
EXTRA_HG="$EXTRA rad_m1/implicit_vimp_lag_local=true rad_m1/time2_adm_poslim=true rad_m1/time2_dtr_poslim=true"
```
The user wants this binary + these keys through the full IDORT suite, same three parts as TASK-2026-11-13-orion-idortd6int
(2b28bb636), whose run outputs you kept. **Reuse them: do NOT rerun FIN or d6int.** Run only the new arm `q` (prod-1118 + EXTRA) and
compare it to your d6int arm (p / D6i) and FIN of the same row / case / cell. Expected: q = d6int bitwise where the cap never
fires (it only stops the face_mmat limiter after 4 rounds) and lzx is a no-op (Cartesian, no rotation); any difference is a result
to report (which rows, size).

## Build (new worktree; never touch main; no destructive git in existing worktrees)
```
cd /orion/ptmp/jinma/Athenak && git fetch fork && git rev-parse fork/prod-1118     # must be 1d1972805 (record it)
git worktree add --detach wt_p18 fork/prod-1118 && git -C wt_p18 submodule update --init --recursive
```
`wt_p18/build_none` with exactly the CMake options of `wt_d6i/build_none`; `wt_p18/build_omp` as `wt_d6i/build_omp` if you made one.
EXE_Q=/orion/ptmp/jinma/Athenak/wt_p18/build_none/src/athena; record md5s.

## (a) Wave grid (180 cases, edd + vetf): as part (a) of 2b28bb636 with RW_EXE=$EXE_Q and the EXTRA keys appended to every
case (however your runner takes extra keys; state how), run dirs tagged `p18`. Compare to d6i (bitwise?) and FIN.

## (b) The 5 sp/cs curv cells: as part (b) of 2b28bb636 with EXE_Q, two arms: `<cell>_P18` (EXTRA) and `<cell>_P18H` (EXTRA_HG,
the HG keys), output ev_all_p18.txt. These cells are where lzx and the positivity keys can act: report L_z / energy columns too.

## (c) The Raven IDORT set (main 41 + rsx 12 + spxx 100 + omp 7): new row runner `$F/idortp1118_1113/one_q.sh` (= one_o.sh with a
single arm q: decks_d6 + $EXTRA on $EXE_Q; writes `<row>/q/`, `sum_q.txt`, `q/extra.txt`). Point W at your existing idortd6int run
root so q lands next to p and f:
```
export C=$F/idortd6int_1113/cart W=/orion/ptmp/jinma/idortd6int_1113 EXE_Q EXTRA OMP_PROC_BIND=false
W=$W/smoke_q NLIM=3 python3 $C/tool/sched.py $C/lists/list_smoke.txt 40 $F/idortp1118_1113/one_q.sh   # 5 rows rc=0 fatal=0
python3 $C/tool/sched.py $C/lists/list_mainrk.txt <NC> $F/idortp1118_1113/one_q.sh > $W/main_q.out 2>&1 &
python3 $C/tool/sched.py $C/lists/list_rsx.txt    <NC> $F/idortp1118_1113/one_q.sh > $W/rsx_q.out  2>&1 &
python3 $C/tool/sched.py $C/lists/list_spxx.txt   <NC> $F/idortp1118_1113/one_q.sh > $W/spxx_q.out 2>&1 &
# OMP rows: EXE_Q=<wt_p18 omp binary> W=$W/omp_q python3 $C/tool/sched.py $C/lists/list_omp.txt 10 $F/idortp1118_1113/one_q.sh
```
Compare: q vs p per row (output md5 / cmpout.py, Picard mean, NON-CONVERGED, dt retries), and q vs f with your analyzers as for
2b28bb636 (cmp_o.sh can be pointed at a view dir with p -> q symlinks).

## Pass rules: as 2b28bb636 (rc=0, 0 FATAL, NON-CONVERGED=0, metric <= 1.05 x FIN, FLAG_PICARD if Picard mean > 1.25 x FIN),
plus: list every row/case/cell where q != p (bitwise), with the size of the change.

## Return: NOTE-2026-11-1x-orion-viper-idortp1118.md on bsg-files-1009, each part as it finishes. Times in PST.
