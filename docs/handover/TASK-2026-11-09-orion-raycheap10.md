# TASK 2026-11-09 viper -> ORION: raycheap10 -- R alone on hesdirk2 (hydro skipped) + R-internal toggles (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap10.md`. raycheap9: dbg_hydro_off is inert on hesdirk2
(driver.cpp needs m1mr). This TASK adds a hesdirk2-path equivalent.

## Binary (REBUILD)
MOD = `fork/raycheap2-1109` head `432b5055c`: new DIAGNOSTIC key `rad_m1/dbg_hydro_skip_t2` (read only when named;
bitwise when absent): on the hesdirk2 path each hydro stage does only its register arithmetic with a zero flux
divergence (stage 1: u1 = u0; stage 2: u0 = gam0 u0 + gam1 u1), so the gas moves only through the radiation exchange
of the stage solves (fluid frozen; gas thermal energy still evolves). viper compile-checked it (incremental MI300A
build, md5 f6e96a8b).

## Runs (your Ar4 deck, t_end 10, 4 levels cfl 0.15 / 0.075 / 0.0375 / 0.01875, 24 ranks)
Part A (first):
| arm | overrides |
|---|---|
| A-PR | rad_m1/time2_vet_col=predict rad_m1/dbg_hydro_skip_t2=true |
| A-X9 | rad_m1/dbg_hydro_skip_t2=true |
Check first that the key acts (the run log of a hydro-skipped run differs from your C0 PR run; hydro task timers ~0).
If early stage fallbacks / BE steps spoil the coarse levels (as C1s did), start ALL levels from one common restart:
run the finest level to t = 1 with an rst output, then restart every level from that file (time/tlim=10, same keys);
the first step after a restart is BE at every level (one step: O(dt^2) globally).

Reading: A-PR p ~ 2 => the 1st-order term is the hydro-radiation coupling (stop after Part A and say so).
A-PR p ~ 1 => the term is inside R: run Part B.

Part B (only if Part A is ~1st order), each on top of A-PR (one toggle each):
| arm | toggle | question |
|---|---|---|
| B1 | rad_m1/implicit_src_stable=false (or =true if the deck runs false) | the stable-source form |
| B2 | rad_m1/implicit_opac_newton=true | opacity re-evaluated at the iterate from pass 0 |
| B3 | rad_m1/implicit_opac_lag_passes=1 | no lagged-opacity passes |
| B4 | rad_m1/time2_k1_stiff=0 rad_m1/time2_k1_redist=0 | the K1 (FSAL) modifications off |
| B5 | rad_m1/implicit_work_fix=false | the work-fix term |
| B6 | rad_m1/implicit_tol=1e-14 rad_m1/implicit_lin_tol=1e-15 | tolerance floor |
(Use the deck's own value to pick the opposite; if a key is refused or stops the run, report and skip it.)

## Analysis
Per arm: `python3 tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> <c0.01875>` (d, p), walls, the hesdirk2 / BE /
fallback counts. Report which toggle (if any) brings p to ~2.

## NOTE back
Part A (and B) conv.py lines with counts; whether the key acted; the common-restart detail if used. No code changes.
