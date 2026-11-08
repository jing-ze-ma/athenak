# NOTE 2026-11-08 Orion -> viper: consfix-1109 energy ledger, 4 arms (answer to TASK-2026-11-08-orion-consfix-ledger)

Short answer: AthenaK rad_m1 does not leak like Athena++. The solves' total energy error over 300 s is
2.8e30 erg (A0) and 3.9e29 erg (B0), i.e. 2.3e-10 and 3.2e-11 of L_in x T. In Athena++ the leak was 11-13 % of E_tot.
All three pass criteria hold (details below). The size criterion holds to within a factor of ~3.

## Build and runs
- Code: `consfix-1109` 988c18215, in a detached worktree `/orion/ptmp/jinma/Athenak/wt_consfix1109` (kokkos submodule
  d8e9af03, 4.6.2).
- Build: modules gcc/13 openmpi/4.1 cmake/3.30, Kokkos SERIAL, 108 s with -j32, 0 warnings:
  `cmake -B build_cf -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_cf -j32`
- Inputs: `docs/handover/consfix-1109/*` from bsg-files-1009. All 4 MD5SUMS OK.
- Deck: `inputs/radiation/he_presn_m1_wedge.athinput` from 988c18215, edited with sed only:
  - `HE_IC_FILE`, `HE_TABLES` -> the extracted files; `R_TOP` -> 2.44e11.
  - mesh nx2 = nx3 = 16, meshblock nx2 = nx3 = 8 (4 MeshBlocks).
  - output1 (hst) dt 5.0 -> 1.0.
  - The arm keys go into `<rad_m1>` in the deck, right after `implicit_closure_thin_relax`:
    - all arms: `dbg_cons_ledger = true`
    - A1: `+ implicit_cons_final = true`
    - B0: `implicit_tol_rel = 0.001`, `implicit_tol_rel_max = 1e-6`, and the deck's `implicit_lin_ew_max = 1.0e-2` changed to `0`
    - B1: B0 `+ implicit_cons_final = true`
  - No command-line overrides. tlim = 300 (deck).
- Slurm job 227802: p.shared, `--qos=serial`, `-n 16 --mem=64G -t 06:00:00`, one node (orion001). The four arms ran at
  the same time:
  `srun -n 4 --exact --mem=15G --cpu-bind=cores build_cf/src/athena -i hepresn_<arm>.athinput`
  (A first try, job 227801, ran the steps one after another because the steps had no `--mem`. It was cancelled after 2 min and its outputs deleted.)
- All four arms reached t = 300.000 s (Terminating on time limit), 52 cycles, 103 solves each. No dt_retry, no
  stage fallbacks, 1 backward-Euler step (the first step).
- The hst has 53 rows (one per cycle, since dt ~ 5-7.7 s > 1). All arms compared at the same t_end = 300 s.
- Run dirs: `/orion/ptmp/jinma/consfix_1109/{A0,A1,B0,B1}` (deck, log, hst, wall.txt). Final restart dumps were deleted after the deviation analysis (quota).

## Table

R/L_in = `ledger.py <run>/hepresn.user.hst 6`. The values below are recomputed with the same formula at 9 digits, because
the arm differences are below ledger.py's 5 printed digits.

| | A0 | A1 | B0 | B1 |
|---|---|---|---|---|
| implicit_cons_final | 0 | 1 | 0 | 1 |
| run sum V(c/chat) r [erg] | +2.842332e30 (not applied) | +2.837156e30 (removed) | -3.947164e29 (not applied) | -3.970565e29 (removed) |
| run sum \|r\| [erg] | 2.735545e31 | 2.737195e31 | 7.504679e30 | 7.506187e30 |
| sum r / (L_in T) | +2.32e-10 | (+2.32e-10) | -3.23e-11 | (-3.25e-11) |
| R/L_in total (ALL) | +8.33102704e-05 | +8.33110017e-05 | +8.33110274e-05 | +8.33110078e-05 |
| w0 0-58.1 s | +4.28510115e-04 | +4.28509855e-04 | +4.28509871e-04 | +4.28509881e-04 |
| w1 58.1-115.2 s | -4.71686582e-07 | -4.71663535e-07 | -4.71317799e-07 | -4.71667785e-07 |
| w2 115.2-166.6 s | +1.21575078e-06 | +1.21621966e-06 | +1.21661933e-06 | +1.21635445e-06 |
| w3 166.6-210.1 s | +8.58916628e-07 | +8.63732782e-07 | +8.63451437e-07 | +8.63500126e-07 |
| w4 210.1-257.2 s | -5.41074516e-07 | -5.40888876e-07 | -5.41211064e-07 | -5.40815489e-07 |
| w5 257.2-300 s | +7.85413527e-07 | +7.85201566e-07 | +7.85058792e-07 | +7.85206824e-07 |
| FLR / FLR_DE / FLR_UN / GAS (cell-solves, energy) | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| g0 clipped cell-passes | 626176 | 626176 | 626176 | 626176 |
| NOT POSITIVE / NOT ADMISSIBLE / BE redo | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| Picard mean (max), NON-CONVERGED | 3.83 (9), 0 | 3.83 (9), 0 | 3.35 (9), 0 | 3.35 (9), 0 |
| bicgstab final max\|r\|/max\|b\| mean / max | 2.2e-5 / 8.1e-4 | 2.2e-5 / 8.1e-4 | 2.4e-6 / 3.5e-5 | 2.4e-6 / 3.5e-5 |
| bicgstab inner iterations total | 674 | 674 | 1219 | 1218 |
| sim time reached | 300 s | 300 s | 300 s | 300 s |
| wall (4 ranks, srun step) | 252 s | 253 s | 258 s | 257 s |

Notes on the table:
- `L_in = 4.0775e37` erg/s, T = 300 s.
- The "BE redo" grep hit in the logs is only the time2_auto_be description line. The rollback counter `backward-Euler redo` is 0.
- The total R/L_in (8.3e-5) is the same in all arms and comes almost entirely from w0 (4.3e-4, the start-up transient).
  Windows 1-5 sit at ~1e-6. Inference: this is the hst time-discretisation floor (trapezoid over 5-7.7 s steps). It does not depend on the arm.

### Arm differences (R/L_in; the r column is that window's sum r / (L_in T_w) from the per-solve `CONSR` lines)
| | total | w0 | w1 | w2 | w3 | w4 | w5 |
|---|---|---|---|---|---|---|---|
| A1 - A0 | +7.31e-10 | -2.6e-10 | +2.3e-11 | +4.7e-10 | +4.8e-09 | +1.9e-10 | -2.1e-10 |
| A0 sum r/(L T_w) | +2.32e-10 | -9.0e-10 | -3.0e-11 | +4.5e-10 | +2.2e-09 | +2.0e-10 | -9.9e-11 |
| B1 - B0 | -1.96e-11 | +1.0e-11 | -3.5e-10 | -2.6e-10 | +4.9e-11 | +4.0e-10 | +1.5e-10 |
| B0 sum r/(L T_w) | -3.23e-11 | -2.5e-11 | -2.0e-10 | -1.1e-10 | +4.6e-11 | +1.5e-10 | -9.4e-12 |
| B1 - A1 | +6.1e-12 | +2.6e-11 | -4.3e-12 | +1.3e-10 | -2.3e-10 | +7.3e-11 | +5.3e-12 |
| B0 - A0 | +7.57e-10 | -2.4e-10 | +3.7e-10 | +8.7e-10 | +4.5e-09 | -1.4e-10 | -3.5e-10 |

### Final-state deviation at t = 300 s (optional item 5)
- Max relative difference of E (rad_m1 u0) and of T (general-EOS wtemp cache, T_K = T_code x 1.202724e-8).
- Mask: cells with tau > 1. tau is the radial integral of rho kappa_R(rho, T) dr from the top, using the deck's Rosseland table on A0's state. 8704 of 163840 cells have tau > 1 (tau_max = 84 at the base).
- Check: T(wtemp) = T_rad to 4 digits at tau > 3.

| | tau > 1: max\|dE\|/E | max\|dT\|/T | all cells: max\|dE\|/E | max\|dT\|/T |
|---|---|---|---|---|
| A1 vs A0 | 4.6e-10 | 1.1e-10 | 3.0e-09 | 7.4e-10 |
| B1 vs B0 | 7.3e-10 | 1.8e-10 | 5.0e-09 | 7.2e-10 |
| B0 vs A0 (tolerance noise) | 6.3e-09 | 1.6e-09 | 6.3e-09 | 1.6e-09 |
| B1 vs A1 | 6.2e-09 | 1.6e-09 | 6.2e-09 | 1.6e-09 |

## Pass criteria
1. Run sum r: |A0| = 2.84e30 > |B0| = 3.95e29 (7.2x). In A1 and B1 the same-size r is applied ("removed from the leak").
   **Holds.**
2. R/L_in:
   - A1 vs B1: totals agree to 6e-12; per window to <= 2.3e-10.
   - A0 vs B0: totals differ by 7.6e-10; up to 4.5e-9 in w3.
   - So with cons_final the ledger does not depend on the tolerance; without it, it does. **Holds.**
   - A1 - A0 = +7.3e-10 against sum r/(L_in T) = +2.3e-10: same sign, factor 3.1. Per window, A1 - A0 has the order of
     magnitude of A0's per-window r. The sign matches in w2-w5 (w3: factor 2.2) but not in w0 and w1. **Holds as "of the size of".**
   - Inference, not tested: the factor ~3 is a lag. Energy that r removes or adds is radiated over the next steps, so it shows up in L_top in later windows. The ledger's trapezoid over ~6 s steps also adds noise at ~1e-10 L_in.
3. Deviation at tau > 1: A1 vs A0 = 4.6e-10 (E), below the tolerance noise A0 vs B0 = 6.3e-9. B1 vs B0 = 7.3e-10.
   **Holds.**

Overall:
- The only leak is the linear residual r, as viper predicted.
- At production tolerances it is ~2e-10 of L_in over 300 s on this CPU wedge (16^2 angles). That is far below the ledger's arm-independent floor (~1e-6 per window, 8e-5 total with the start-up window).
- Cost of cons_final: no measurable wall time (253 vs 252 s; 257 vs 258 s).
- The tight tolerance (B) costs 1.8x the inner BiCGStab iterations but fewer Picard passes. Wall is +2 %.

## Commands
```
# worktree + build
git -C /orion/u/jinma/ATHENAK/athenak worktree add --detach /orion/ptmp/jinma/Athenak/wt_consfix1109 988c18215
cd /orion/ptmp/jinma/Athenak/wt_consfix1109 && git submodule update --init
module purge; module load gcc/13 openmpi/4.1 cmake/3.30
cmake -B build_cf -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release \
      -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build_cf -j32
# inputs
for f in MD5SUMS ic_he_presn_m1_gas.txt ledger.py planck_he_x0.0_z0.02_ferg+tops.txt rosseland_he_x0.0_z0.02.txt; do
  git show origin/bsg-files-1009:docs/handover/consfix-1109/$f > inputs/$f; done; (cd inputs && md5sum -c MD5SUMS)
# deck per arm (I=inputs dir, D=deck, extra = the arm keys above)
sed -e "s#HE_IC_FILE#$I/ic_he_presn_m1_gas.txt#" -e "s#HE_TABLES#$I#g" -e "s#R_TOP#2.44e11#" \
  -e '/^<mesh>/,/^<meshblock>/{s/^nx2       = 64/nx2       = 16/;s/^nx3       = 64/nx3       = 16/}' \
  -e '/^<meshblock>/,/^<time>/{s/^nx2       = 32/nx2       = 8/;s/^nx3       = 32/nx3       = 8/}' \
  -e '/^<output1>/,/^<output2>/s/^dt          = 5.0/dt          = 1.0/' \
  -e "s/^implicit_closure_thin_relax = 1.5.*/&\n$extra/" $D > $arm/hepresn_$arm.athinput
# B arms also: sed -i 's/^implicit_lin_ew_max = 1.0e-2/implicit_lin_ew_max = 0/'
# sbatch (p.shared, --qos=serial, -n 16, --mem=64G, -t 06:00:00), one srun per arm in the background, then wait:
srun -n 4 --exact --mem=15G --cpu-bind=cores $EXE -i hepresn_$arm.athinput > log 2>&1
# analysis (python-waterboa/2025.06)
python3 inputs/ledger.py $arm/hepresn.user.hst 6
grep 'consfix-1109:' $arm/log ; grep CONSR $arm/log      # per-solve sumVr, cumulative run
```
- The deviation script reads the restart with a local copy of `docs/handover/scripts/he_remap_rst.py`. Two changes were
  needed: (a) it knows the `M1CSK1S1` / `M1VGDAS1` blocks; (b) the hesdirk2 `t2` block is no longer refused, and the trailing slabs are skipped.
  Only hydro, m1 u0 and wtemp are read; they come before those slabs in the record.
- The scripts are at `/orion/ptmp/jinma/consfix_1109/scripts/{dev.py,rr_local.py}` (not pushed).
