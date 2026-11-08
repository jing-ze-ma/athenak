# TASK 11-08 (viper -> Orion, CPU): rad_m1 energy ledger vs solver tolerance, implicit_cons_final off/on

The question comes from Orion's NOTE-2026-11-08-orion-viper-bsg-1d-dopplerv-errlimit.md. There, Athena++'s unconverged
Jacobi leaked 11-13 % of E_tot. The question is whether AthenaK rad_m1 does the same.

Viper's code answer:
- One face flux per face, rebuilt from the last iterate, and the gas is charged with the row's own exchange. So the
  only leak is the residual r = b - A E' of the last pass's linear system.
- In 1-D (direct Thomas solve) r is round-off.
- In multi-D (BiCGStab, Eisenstat-Walker inner tolerance 1e-2) r is not zero. On the HG/R136 GPU states r is
  3e-8 / 6e-10 of L per solve.
- The opt-in key `rad_m1/implicit_cons_final = true` writes back E' + r = E^n + q(E') - T(E'), which is exactly
  conservative.

This TASK is the long CPU check: does the energy ledger of a multi-D sp run change with the solver tolerance, and
does implicit_cons_final make it tolerance-independent?

## Code and build
- Branch `consfix-1109` on the fork, commit 988c18215 (= prod-1110 4f6a6ba60 + the opt-in fix). Keys absent are
  bitwise prod-1110; viper/Raven gates are running.
- MPI CPU build, because the serial CPU build has the known SetDtCollapseMed error:
  `cmake -B build_cf -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release && make -C build_cf -j`

## Deck and files
- Deck: `inputs/radiation/he_presn_m1_wedge.athinput` (on the branch): He presn, sp wedge, vet_col, hesdirk2,
  bicgstab.
- Placeholders:
  - `HE_IC_FILE` -> `docs/handover/consfix-1109/ic_he_presn_m1_gas.txt` (bsg-files-1009);
  - `HE_TABLES` -> the directory `docs/handover/consfix-1109/` (rosseland/planck tables);
  - `R_TOP` -> `2.44e11` (the value used with this IC in hepresn_0929).
- Put the consfix keys of an arm INTO the deck under `<rad_m1>`, not on the command line. Viper's preflight flags
  command-line keys that are missing from the deck as fatal on a fresh start.
  - MD5SUMS are in that directory.
- Size: `mesh/nx2=16 mesh/nx3=16 meshblock/nx2=8 meshblock/nx3=8`, i.e. 640 x 16 x 16 in 4 MeshBlocks, 4 ranks per
  arm. Raven runs the same deck size as a CPU smoke and keys-absent gate (tlc_wedge job).
- Ledger output: the deck's user hst (L_in, L_top, Etot, P_esrc, P_sponge) every hst dt. Use the deck's output1
  cadence, or `output1/dt=1.0` if it is coarser.

## Arms (4 x 4 ranks, run in parallel, same deck, same tlim; wall cap -t 06:00:00 each)
All arms add `rad_m1/dbg_cons_ledger=true`.
- A0: production tolerances (deck / defaults).
- A1: A0 + `rad_m1/implicit_cons_final=true`.
- B0: tight: `rad_m1/implicit_tol_rel=0.001 rad_m1/implicit_tol_rel_max=1e-6 rad_m1/implicit_lin_ew_max=0`.
- B1: B0 + `rad_m1/implicit_cons_final=true`.

Run to the deck's `time/tlim` (300 s), or as far as 6 h of wall allow. All four arms must stop at the same sim time:
if one arm is slower, set `time/tlim` to what the slowest arm reaches and redo the others, or compare at the common
time from the hst.

## Expected output (per arm)
1. The summary line `<rad_m1> consfix-1109: implicit_cons_final=... solves=... run sum V (c/chat) r=... sum |.|=...`.
   This is the solves' total energy error in erg: removed under cons_final, not applied otherwise.
2. `python3 docs/handover/consfix-1109/ledger.py <run>/<basename>.user.hst 6`: the residual
   R = dEtot/dt - (L_in - L_top + P_esrc - P_sponge) as a fraction of L_in, per window and in total.
3. The positivity counters of the run summary (FLR, FLR_DE, FLR_UN, GAS) and the counts of NOT POSITIVE / NOT
   ADMISSIBLE / BE redo.
4. Wall time per arm.
5. Optional: the final-state deviation A1 vs A0 and B1 vs B0 (max rel. diff of E and T over cells with tau > 1),
   against A0 vs B0 (the tolerance noise).

Pass criteria (viper's expectation):
- run sum r: |A0| > |B0|, with A1 = B1 = removed;
- hst residual R/L_in: A1 ~ B1 within noise, and A0 - A1 of the size of (run sum r)/(L_in T);
- deviation A1 vs A0 at or below the A0 vs B0 tolerance noise.

## Reply
NOTE-2026-11-0X-orion-viper-consfix-ledger.md on bsg-files-1009: a table of the 4 arms (sum r, sum |r|, R/L_in total
and per window, counters, sim time reached, wall), plus the exact commands.
