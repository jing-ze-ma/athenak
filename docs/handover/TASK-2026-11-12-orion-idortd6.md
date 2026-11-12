# TASK 2026-11-12 viper -> Orion: IDORT D6 (final 2nd-order candidate) — linear-wave grid amp 1e-6 + the 5 sp/cs cells (CPU)

D6 = FINAL + H2 keys MINUS the cmem/physflux family (implicit_hr_cmem, _cap, _rst, _fainit, implicit_physflux, _mom, _zold,
implicit_am_cs_fold, implicit_hr_cmem_implicit) PLUS rad_m1/implicit_blend_dtc=1e-6 (one-flux DTC, FCB_1113). D6 has passed hcyl beam,
HG fresh-start L_z, stiff pulse_cs 2nd order and the accretor 70 cycles on viper; this is its IDORT gate. Same binary, lists, drivers,
analyzers and return format as TASK-2026-11-12-orion-idortfin-rw (bb297426b) and idortfin-more part (2); only the arm changes.
Compare against YOUR FIN results of those TASKs (same binary), so only the D6 arm runs here.

## Binary: the existing 421a3d024 build (wt_fin/build_none, md5 e151dc3b...; it has implicit_blend_dtc). No rebuild.

## (1) Linear-wave grid, arm D6 (180 cases: edd + vetf)
```
cd /orion/ptmp/jinma/rwfin_1112/rw_wbs        # where fin_FIN_*.txt / runs_finFIN_* already are
X="rm_time2_stage_local_be=true rm_time2_face_fct=true rm_time2_k1_redist=0 rm_implicit_face_mmat=true rm_implicit_cs2_wbsolve=true rm_implicit_hr_cmem=false rm_implicit_physflux=false rm_implicit_blend_dtc=1e-6"
for c in edd vetf; do
  sed -e "s/ rm_implicit_physflux_zold=fc//" -e "s/ rm_implicit_hr_cmem_rst=true//" -e "s/\$/ $X/" fin_FINAL_$c.txt > fin_D6_$c.txt
  awk '$4 == 16' fin_D6_$c.txt > fins_D6_$c.txt; done
grep -c "zold\|cmem_rst\|am_cs_fold" fin_D6_*.txt    # must be 0 0
```
(The TIME arm sets implicit_hr_cmem=true, cap, physflux=true, mom=fa; the rm_ lines come after the arm lines in the deck and the
parser keeps the last value, so cmem and physflux are OFF. Pre-check on one fins_D6_vetf deck: implicit_hr_cmem = false (last),
implicit_physflux = false (last), implicit_blend_dtc = 1e-6, face_mmat true, cs2_wbsolve true, stage_local_be true, face_fct true,
k1_redist 0, NO zold / cmem_rst / fainit / am_cs_fold lines, dt_retry false, radwave_amp 1e-6, implicit_tol 1e-12.)
Smoke (N16, 18 cases per list, ~1 min), then everything at once (--np 90, one node per list; vetf first):
```
export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt RW_EXE=/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena
for c in edd vetf; do RW_RUNS=$PWD/runs_finD6_$c python3 ../rw_amp1e6/drive.py run fins_D6_$c.txt --np 18 > fins_D6_$c.out 2>&1 & done; wait
grep -c "ok w" fins_D6_*.out; grep -h -E "FAILED|FATAL" fins_D6_*.out
for c in edd vetf; do RW_RUNS=$PWD/runs_finD6_$c python3 ../rw_amp1e6/drive.py run fin_D6_$c.txt --np 90 > fin_D6_$c.out 2>&1 & done; wait
```
A FATAL is a result: report it, do not add keys or re-enable dt_retry.

## (2) The 5 sp/cs IDORT cells, arm D6 (CPU, as idortfin-more part 2)
Decks: `docs/handover/scripts/idortfin_1112/curv/<cell>_D6/in.athinput` on this branch (= the FIN deck with implicit_hr_cmem = false,
implicit_physflux = false, implicit_hr_cmem_cap / implicit_physflux_mom removed, implicit_am_cs_fold replaced by implicit_blend_dtc = 1e-6).
Cells T2cs_shdC, T2sp_shd, T1sp_star, T3cs_hot5V, T3sp_hot5; ranks: cs cells 6, sp cells 1. Run next to your FIN/X1/P1 runs of part (2) ($F = scripts/idortfin_1112 of a fresh checkout of this branch, $O and $EXE as in idortfin-more part 2):
```
for d in $F/curv/*_D6; do c=$(basename $d); r=$O/curv/$c; mkdir -p $r; cp $d/in.athinput $r/
  np=1; case $c in *cs_*) np=6;; esac
  ( cd $r && mpirun -np $np $EXE -i in.athinput -d $r > run.log 2>&1; echo $? > rc.txt ) &
done; wait
cd $O/curv && ATHENAK_VIS=/orion/ptmp/jinma/Athenak/wt_fin/vis/python python3 $F/curv/ev.py */ > ev_all_d6.txt
```
PASS per cell: D6 rc=0, 0 FATAL, Picard NON-CONVERGED=0, every ev.py metric <= 1.05 x X1 of the same cell (your CPU X1 run).
FLAG any D6 metric > 1.05 x X1 or > P1, and give D6/FIN per metric.

## Return: NOTE-2026-11-1x-orion-viper-idortd6.md on this branch (NOTE only), post each part as it finishes
(1) wall per list; `grep -c "ok w"` / FAILED / FATAL; N256 "ok w" lines (vetf + edd); cmp_wbs.py D6 vs FIN (edd, vetf) and D6 vs FINAL;
analyze_edd.py / analyze.py through view dirs with PROD = FIN arm, FINAL = D6 arm, pasted whole. Flag every cell where D6 has a lower
order (raw or mean-removed, 128 -> 256) or a larger N256 L1 (rho or E) than FIN by > 1 %; state separately whether the P3.33 tau10 fix
(FIN p rho 2.00) and the P3.33 tau0.01 flag (FIN L1 E-mrm 3.57e-4) hold for D6. Also the face_mmat limited-face count (FIN: 0).
(2) ev_all_d6.txt whole + rc / FATAL / NON-CONVERGED of the 5 D6 runs, with the FIN and X1 values beside each D6 metric.
