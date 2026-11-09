# NOTE 2026-11-09 ORION -> viper: raycheap10 -- R alone on hesdirk2 (dbg_hydro_skip_t2) + R-internal toggles

Answer to `TASK-2026-11-09-orion-raycheap10.md`. No code changes.

**Summary**
- **The key acts.** With `dbg_hydro_skip_t2 = true` the step stays at its initial value (dt = 0.29563 at c0.15 for the
  whole run, 34 steps to t = 10), while the same deck without it (C0n, rerun on this binary) has dt falling to 0.252
  (39 steps). Final E differs from C0n by 12 % of max E at c0.15 (raycheap9: C1s vs C2 13 %). The log has no hydro
  task timer lines, so this is the evidence.
- **Part A: R alone is still 1st order.** A-PR p = 0.55 / **1.06**, A-X9 p = 0.47 / **0.89** (finest pair in bold).
  All A runs are clean: 1 backward-Euler step (the first), 0 fallbacks, 0 NOT ADMISSIBLE, 0 DT_RETRY, so the
  common-restart scheme was not needed. Per the TASK reading, the 1st-order term is inside R, not the hydro coupling.
- **Part B: no toggle brings p to ~2.** B1 (src_stable false), B3 (opac_lag_passes 1), B6 (tol 1e-14 / lin_tol 1e-15) and B5
  (work_fix false) leave the solution unchanged or nearly so: B3 is bitwise identical to A-PR; B1 and B6 are E-bitwise
  with F differing by <= 1e-9 absolute; B5 differs by 6.4e-8 of max E. Their conv lines equal A-PR's to 4 digits.
  B4 (time2_k1_stiff 0, time2_k1_redist 0) changes the solution (1.7e-3 at c0.15, 2.6e-5 at the finest level) but brings 2-3
  NOT ADMISSIBLE -> BE redos per level (BE steps 6/6/6/11) and stays at p 0.37 / 1.03.
- **B2 refused** at start on all levels: `<rad_m1>/implicit_opac_newton needs implicit_opac_update = true and
  implicit_gas_newton = true` (rc 1, wall 2-3 s). Skipped as instructed (enabling the two prerequisites would be three toggles).
- Part A turned out ~1st order, so Part B was required; both were run together in one job.

## Binary and decks
- MOD `432b5055c` rebuilt (gcc/13, openmpi/4.1, Release, MPI, `-ffp-contract=off`); kept at
  /orion/ptmp/jinma/raycheap10_1109/athena_432b5055c.
- Base: raycheap8 Ar4 deck (`full_decks/evoAr4_c*`). Keys absent from the deck were written into deck copies (sed),
  since fresh starts refuse absent command-line keys; t_end via `time/tlim=10`.
  A-PR = time2_vet_col predict + dbg_hydro_skip_t2 true; A-X9 = deck + dbg_hydro_skip_t2 true; B* = A-PR + one toggle.
- B toggles chosen opposite to the fresh-run default (the deck names none of them except the tolerances):
  implicit_src_stable default true -> false; implicit_opac_newton default false -> true; implicit_opac_lag_passes default 3 -> 1;
  time2_k1_stiff default 0 (already off) and time2_k1_redist default 0.1 -> both 0; implicit_work_fix default true -> false;
  deck implicit_tol 1e-13 / implicit_lin_tol 1e-14 -> 1e-14 / 1e-15.
- One job 228060 (9 exclusive nodes, all 32 runs + C0n at once, 12 min).

## conv.py (t_end 10; cfl 0.15 / 0.075 / 0.0375 / 0.01875)

| arm | d | p | wall (s) | NOT ADMISSIBLE | hesdirk2 stage/BE/fallback | builds |
|---|---|---|---|---|---|---|
| **A-PR** (predict + skip) | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 203/297/443/656 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| **A-X9** (extrap2 + skip) | [1.272e-03, 9.213e-04, 4.959e-04] | [0.47, 0.89] | 203/297/437/663 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| B1 implicit_src_stable=false | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 204/292/431/654 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| B3 implicit_opac_lag_passes=1 | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 203/295/440/656 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| B4 time2_k1_stiff=0, time2_k1_redist=0 | [8.407e-04, 6.489e-04, 3.181e-04] | [0.37, 1.03] | 200/305/451/659 | 2/2/2/3 | 28/6/2 62/6/2 130/6/2 260/11/3 | 36/70/138/274 |
| B5 implicit_work_fix=false | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 210/304/441/664 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| B6 implicit_tol=1e-14, lin_tol=1e-15 | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 226/336/466/737 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| B2 implicit_opac_newton=true | refused (FATAL at start) | - | 2/2/2/3 | - | - | - |

For reference (raycheap9, hydro on): C0 = PR p 0.87 / 0.87; C6 = be 1.03 / 1.07.

Final-E differences vs A-PR (max|dE|/max|E|, c0.15 / c0.01875; cmpbin m1 bitwise count at c0.01875):
- B1 0 / 0 (E bitwise; F 2/4 bitwise, max abs 6e-11); B3 0 / 0 (4/4 bitwise); B6 0 / 0 (2/4, max abs 9e-10);
  B5 6.4e-8 / 6.4e-8; B4 1.7e-3 / 2.6e-5; A-X9 7.6e-4 / 1.9e-4.
- C0n (hydro on) vs A-PR at c0.15: 0.123.

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2 counters)
```
APR_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 203 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
APR_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 297 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
APR_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 443 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
APR_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 656 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
AX9_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 203 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
AX9_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 297 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
AX9_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 437 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
AX9_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 663 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
B1_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 204 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
B1_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 292 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
B1_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 431 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
B1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 654 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
B2_c0.15 rc=1 FATAL=24 NOTADM=0 DTRETRY=0 wall 2 s builds= | 
B2_c0.075 rc=1 FATAL=24 NOTADM=0 DTRETRY=0 wall 2 s builds= | 
B2_c0.0375 rc=1 FATAL=24 NOTADM=0 DTRETRY=0 wall 2 s builds= | 
B2_c0.01875 rc=1 FATAL=24 NOTADM=0 DTRETRY=0 wall 3 s builds= | 
B3_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 203 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
B3_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 295 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
B3_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 440 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
B3_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 656 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
B4_c0.15 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 200 s builds=36 | hesdirk2: stage steps=28 backward-Euler steps=6 stage fallbacks=2
B4_c0.075 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 305 s builds=70 | hesdirk2: stage steps=62 backward-Euler steps=6 stage fallbacks=2
B4_c0.0375 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 451 s builds=138 | hesdirk2: stage steps=130 backward-Euler steps=6 stage fallbacks=2
B4_c0.01875 rc=0 FATAL=0 NOTADM=3 DTRETRY=0 wall 659 s builds=274 | hesdirk2: stage steps=260 backward-Euler steps=11 stage fallbacks=3
B5_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 210 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
B5_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 304 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
B5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 441 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
B5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 664 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
B6_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 226 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
B6_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 336 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
B6_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 466 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
B6_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 737 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
C0n_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 228 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
```
All stage fallbacks other than B4 and all accepted-with-fallback counts are 0.

## Files
/orion/ptmp/jinma/raycheap10_1109: decks/ (APR, AX9, B1-B6, C0n), rc10.sbatch, ana.sh, rc10_ana.out, runs/*/run.log + wall.txt,
athena_432b5055c (kept). bins, rst, sources and the build dir deleted.
