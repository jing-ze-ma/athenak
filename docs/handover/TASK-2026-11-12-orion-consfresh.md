# TASK 2026-11-12 viper -> Orion: conservation audit (HG small cs deck, CPU) + d6fast3 fresh-start gate

Viper GPUs are reserved for the speed-up work, so these CPU-valid runs go to Orion. Two parts, ONE deck pair:
(a) a measured conservation audit (energy, L_z, radial momentum) of D6 on d6fast3 vs production keys on prod-1115;
(b) the d6fast3 fresh-start gate on the small HG deck (was JOB A of the viper apudev worker).
Files: `docs/handover/scripts/consfresh_1112/` on this branch (decks, 2 patches, 2 analyzers).

## Commits (all on fork)
- d6fast3-1112 **fbb5c6638** (D6 arm). Viper has shown its outputs are bitwise identical to d6fast2 7bb6b52a6 on R136/HG.
- prod-1115 **d3b3f2e55** (production arm).
- Budget diagnostic (consaudit-1112, output only, NOT pushed as a branch): two commits shipped as patches
  `consaudit_on_fbb5c663.patch` and `consaudit_on_d3b3f2e55.patch` (`git am` onto the commit named in the file name).
  The patches add `<rad_m1>/dbg_cons_audit = true` (needs `dbg_energy_audit = true`). Output: CAUDIT and CAUDITB lines at the end of every M1
  chain, and CAUDITH after each hydro stage. Viper GPU check (apudev 13088316, HG deck, 100 cycles): with both keys on, hegiant.hydro.hst and
  hegiant.user.hst are md5-identical to the run without the diagnostic (D6 7be29ba2 / b918ec25, PROD adbfc9a1 / 803e5551). The .bin files
  differ only because the embedded input now holds the two extra keys.

## Builds (4 trees; gcc, MPI, -ffp-contract=off; the HG deck needs the he_star_m1 user pgen)
```
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
R=/orion/ptmp/jinma/Athenak; S=<bsg-files-1009 checkout>/docs/handover/scripts/consfresh_1112
git -C $R fetch fork
for x in "f fbb5c6638" "p d3b3f2e55" "fa fbb5c6638" "pa d3b3f2e55"; do set -- $x
  git -C $R worktree add --detach $R/wt_cf_$1 $2 && (cd $R/wt_cf_$1 && git submodule update --init)
done
(cd $R/wt_cf_fa && git am $S/consaudit_on_fbb5c663.patch)
(cd $R/wt_cf_pa && git am $S/consaudit_on_d3b3f2e55.patch)
for t in f p fa pa; do (cd $R/wt_cf_$t && cmake -B build -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON \
   -D PROBLEM=he_star_m1 -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build -j32) & done; wait
```
Report `git rev-parse HEAD` and the md5 of each `build/src/athena`. `strings wt_cf_fa/build/src/athena | grep -c dbg_cons_audit` must be >= 1.

## Decks
`hg_D6.athinput` and `hg_PROD.athinput`: the viper small HG cs deck (704 x 16 x 16 x 6, 24 MeshBlocks, fresh start from the IC, t = 0), with
production KEYS (prod_1107/hg/run.cfg) + KEYS_PH (PROD) or KEYS_D6H (D6 = FINAL minus cmem/physflux + implicit_blend_dtc = 1e-6), nlim 40.
The table paths are `files/...`. Link `files` -> `<checkout>/docs/handover/hgpos-1109/files` (the viper md5s match: ic f9678c2a, planck
23ee543f, rosseland e2868b75). There is no restart, so no rotation applies.

## Runs (24 ranks each; all 6 can run at once on two nodes; smoke first)
```
export OMP_NUM_THREADS=1 OMPI_MCA_hwloc_base_binding_policy=none
O=/orion/ptmp/jinma/consfresh_1112; F=<checkout>/docs/handover/hgpos-1109/files
DG="rad_m1/dbg_energy_audit=true rad_m1/dbg_cons_audit=true"
AM="hydro/am_cons_diag=100000 rad_m1/implicit_am_cs_cons_diag=100000 rad_m1/implicit_am_cs_diag=4"
run() { d=$O/$1; mkdir -p $d; cp $S/hg_$2.athinput $d/in.athinput; ln -sfn $F $d/files
        (cd $d && mpirun -np 24 $R/wt_cf_$3/build/src/athena -i in.athinput -d $d $4 > run.log 2>&1; echo $? > rc.txt) ; }
# smoke: 3 cycles of each binary+deck+keys pair, stop on any FATAL
run smk_GD6 D6 f "time/nlim=3 $AM" & run smk_GP PROD p "time/nlim=3 $AM" &
run smk_AD6 D6 fa "time/nlim=3 $DG" & run smk_AP PROD pa "time/nlim=3 $DG" & wait
grep -c FATAL $O/smk_*/run.log; cat $O/smk_*/rc.txt
# (b) gate, 40 cycles from the IC, plain binaries
run G_D6 D6 f "time/nlim=40 $AM" & run G_PROD PROD p "time/nlim=40 $AM" &
# (a) audit, 100 cycles, diagnostic binaries (the first 10 cycles give the 10-cycle numbers)
run A_D6 D6 fa "time/nlim=100 $DG" & run A_PROD PROD pa "time/nlim=100 $DG" &
# (a') neutrality on CPU: the diagnostic binary with the keys OFF, 40 cycles
run N_D6 D6 fa "time/nlim=40" & run N_PROD PROD pa "time/nlim=40" & wait
```
A FATAL is a result. Report it and do not change keys.

## Analysis
```
cd $O
for a in G_D6 G_PROD; do echo $a rc=$(cat $a/rc.txt) FATAL=$(grep -c FATAL $a/run.log) NA=$(grep -ci "not admissible" $a/run.log) \
  NP=$(grep -ci "not positive" $a/run.log); grep "time_scheme=hesdirk2: stage steps" $a/run.log | tail -1; done
python3 $S/amlz_analyze.py G_D6 G_PROD
for a in A_D6 A_PROD; do python3 $S/cons_analyze.py $a/run.log; python3 $S/cons_analyze.py $a/run.log 10 | grep -E "^SUM|^REL"; done
# neutrality: the first 41 history rows of A_* (keys on) = N_* (keys off) = G_* (plain build)
for x in D6 PROD; do for h in hydro user; do head -n 42 A_$x/hegiant.$h.hst | md5sum; md5sum < N_$x/hegiant.$h.hst; md5sum < G_$x/hegiant.$h.hst; done; done
```
(G_* run with the AM keys, which print only; the hst rows must still match N_* and A_*. If G differs from N, report both.)

## Pass / report
(b) gate PASS: D6 rc 0, FATAL 0, NOT ADMISSIBLE 0, NOT POSITIVE 0. Report the stage/BE counts and the face_fct "limited cells" from the
hesdirk2 summary line. Report gas L_z after 40 cycles and AMCS_DIAG |del_z| for D6 and PROD. Viper reference values (GPU, 30 cycles):
production del 2.2e39, L_z -4.5e38. D6 should have the production bookkeeping: del ~2e39, |L_z| within a few x 1e38.
(a) audit: paste the full cons_analyze.py output for A_D6 and A_PROD (100 cycles), plus the SUM/REL lines for 10 cycles. Viper GPU values
for the same decks (13088316, d6f2 = d6f3 bitwise): the fix segments are exactly 0 in both arms. Energy residual RE_b (cycles with 2 solves) falls from
2.23e38 (cycle 2) to 1.02e38 (cycle 10), against a radiation boundary term of 7e37-4.1e39. At cycles 97-99 it is ~1e34, against 1.2e39.
These values are the same to 3 digits in both arms. Accumulated L_z(gas+rad) change over 99 cycles:
D6 -4.78e40, PROD -4.50e40. P_r (gas+rad) change 1.174e35 in both arms.
State whether CPU reproduces these to 3 digits, and whether D6 is worse than PROD in any column.
Neutrality: the md5 table above.

## Return: docs/handover/NOTE-2026-11-1x-orion-viper-consfresh.md on this branch (NOTE only). Post (b) as soon as it lands, then (a).
