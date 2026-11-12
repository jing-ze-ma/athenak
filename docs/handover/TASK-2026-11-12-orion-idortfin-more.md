# TASK 2026-11-12 viper -> Orion: more FINAL-method (421a3d024) CPU gates: wave ablation (vetf), 5 sp/cs IDORT cells, pulse_cs ladders

USER 11-12: viper GPUs are full, move work to Orion; run everything at once. Same binary as TASK idortfin-rw:
`/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena` (421a3d024, MPI ON). Files: `docs/handover/scripts/idortfin_1112/` on this
branch (decks inline there). viper runs copies of (2) and (3) too; whichever finishes first counts (post as soon as each part is done).
Common env: `module purge; module load gcc/13 openmpi/4.1; export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
EXE=/orion/ptmp/jinma/Athenak/wt_fin/build_none/src/athena F=<checkout>/docs/handover/scripts/idortfin_1112 O=/orion/ptmp/jinma/idortfin_more`

## (1) Wave key ablation, P 3.33 x tau 0.01 / 0.1, edd (viper: DONE) + vetf (please run)
viper result (edd, N64-256): the ONLY key that changes these cells is implicit_cs2_wbsolve: FIN minus wbsolve = FINAL to 4 digits
(tau 0.01 N256 L1 E 3.089e-4 vs FIN 3.523e-4; tau 0.1 4.339e-4 vs 3.975e-4); stage_local_be, face_fct, face_mmat, am_cs_fold,
+zold, +cmem_rst inert (k1_redist 0 is in TIME; xfade / x2_fields / lat lag do not act on 1-D Cartesian edd). FIN minus wbsolve with
implicit_cs2_vimp_npass=30 = FINAL (0.1 %): the 3-pass fixed point is converged, the wbsolve closed form solves a different relation.
Please repeat for vetf: lists `$F/rw/d_<arm>.txt` are edd lines; make vetf lists from your fin_FINAL_vetf.txt cells P 3.33 tau 0.01/0.1
N 32-256 with the same token edits (arms FINAL, FIN, FINmwbs, FINnp30 = FINmwbs + rm_implicit_cs2_vimp_npass=30, FINmslbe, FINmfct,
FINmfmm, FINmfold; see the edd files for the exact tokens), run with drive.py as in idortfin-rw (one drive per arm, all at once),
and tabulate with `$F/rw/tab.py` (edit B and the edd->vetf reference: refs["vetqs_d"] instead of refs["edd"]). Report the table.

## (2) The 5 sp/cs IDORT FIN cells + references (CPU)
Decks: `$F/curv/<cell>_{FIN,X1,P1}/in.athinput`, cells T2cs_shdC, T2sp_shd, T1sp_star, T3cs_hot5V, T3sp_hot5 (30 cycles each, seconds
on GPU). X1 / P1 = the production / X1b-era references (P1 = PROD keys, X1 = FINAL + X1 keys); their viper GPU evaluations are
`ev_viper_gpu.json` next to the decks. Ranks = MeshBlocks: cs cells (T2cs_shdC, T3cs_hot5V) 6, sp cells 1.
```
for d in $F/curv/*_{FIN,X1,P1}; do c=$(basename $d); r=$O/curv/$c; mkdir -p $r; cp $d/in.athinput $r/
  np=1; case $c in *cs_*) np=6;; esac
  ( cd $r && mpirun -np $np $EXE -i in.athinput -d $r > run.log 2>&1; echo $? > rc.txt ) &
done; wait
cd $O/curv && ATHENAK_VIS=/orion/ptmp/jinma/Athenak/wt_fin/vis/python python3 $F/curv/ev.py */ > ev_all.txt
```
PASS per cell: FIN rc=0, 0 FATAL, Picard NON-CONVERGED=0 (last "Picard iterations" line), and every ev.py error metric of FIN
<= 1.05 x the X1 value of the same cell (your CPU X1 run, not the GPU json; the json is for a sanity check of the CPU references).
FLAG any FIN metric > 1.05 x X1 or > the P1 value. Paste ev_all.txt whole + the rc / FATAL / NON-CONVERGED of every run.

## (3) pulse_cs ladders on FIN vs the TIME base B (same binary)
Decks: `$F/pulse/{B,FIN}_{T_1,T_2,T_4,T_8,C_32,C_64,C_128}/in.athinput` (B = CSP-CM~fA2~norot~wall~x2f15~cap~pf, FIN = B + X1b +
face_mmat + posfix lag + cs_e0_metric + cs2_wbsolve + am_cs_fold applied). Ranks: 6 (C_128: 24). All 14 at once:
```
for d in $F/pulse/*; do c=$(basename $d); r=$O/pulse/$c; mkdir -p $r; cp $d/in.athinput $r/; np=6; [ $c = ${c%C_128}C_128 ] && np=24
  ( cd $r && mpirun -np $np $EXE -i in.athinput -d $r > log 2>&1; echo $? > rc ) &
done; wait
ATHENAK_VIS=/orion/ptmp/jinma/Athenak/wt_fin/vis/python python3 $F/pulse/pulse_ord.py $O/pulse > pulse_ord.txt
```
Report pulse_ord.txt (T-ladder successive differences and orders per field, B vs FIN), rc / FATAL of all 14, and the last
"Picard iterations" line of each. FLAG a FIN T-order below 1.8 where B is >= 1.9 on m1_e or m1_f1. (C ladder: report rc and wall
only; viper evaluates the shell-mean orders from its own copy.)

## Return
NOTE-2026-11-1x-orion-viper-idortfin-more.md on this branch, one section per part, posted as each part finishes.
