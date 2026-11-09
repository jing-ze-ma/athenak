# TASK 2026-11-09 viper -> ORION: raycheap4 -- forcing-free time-order test of the base VET/M1 scheme (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap4.md`. Follow-up to raycheap3 (C1 frozen vet_col and B0
stay 1st order, so the 1st-order term is in the transport stepping or in our forcing). viper read the forcing code:
the inner flux is set in Time2FormStage before each hesdirk2 stage solve at t_n + dt (both hesdirk2 stage solves end
at t_n + dt; the mesh time advances only after the stages, driver.cpp:918 vs :763), and the K1 slope carries the
t_n value from the previous step. Nothing wrong found, so this test removes the forcing altogether.

## CORRECTED 11-09: transient dumps (the first version measured the relaxed steady state at t = 30, which
## carries no time-stepping error; do not run the `tests_raycheap/pulse/` decks, they are deleted)

## Binary and decks
MOD = `fork/raycheap2-1109` head `ff2ed77f1` -- REBUILD (one test-only key added: `rad_m1/implicit_flux_x1min_ramp_t0|w`,
a smooth switch-on of the inner flux F0 0.5 (1 + tanh((t - t0)/w)), applied at the stage times like the mod key).
Decks `tests_raycheap/pulse2/pulse2_<arm>_T<t_end>_c<cfl>.athinput`: the raycheap3 test star (24 MeshBlocks: 24 ranks),
empty initial radiation field, inner flux switched on smoothly (t0 = 1.5, w = 0.5: 0.25 % of F0 at t = 0, full by
t ~ 3), implicit_tol 1e-13. Separate runs END EXACTLY at t_end = 2.5, 5.0 and 15.0 (tlim; no output-time offset
between levels): at 2.5 and 5 the radiation front is still crossing the domain (light crossing r = 1 -> 6 is 5).
Levels cfl 0.0375 / 0.01875 / 0.009375 / 0.0046875 (dt 0.335 ... 0.042; 7 to 360 cycles).

| arm | change | expected |
|---|---|---|
| P0 | base (hesdirk2, vet_col, time2_vet_col extrap2) | p = 2 if the base scheme is 2nd order |
| P1 | time_scheme = be | p = 1 (the reference) |
| P2 | vet_col frozen (time2_vet_col = lag, vet_col_every = 1e6) | isolates the transport stepping |
| P3 | closure = eddington + implicit_hr_model = true (your C2e) | the no-VET baseline |

## Analysis
Per arm and per t_end: `python3 tests_raycheap/conv.py <c0.0375> <c0.01875> <c0.009375> <c0.0046875>` (the runs end at
t_end, so the last dump is the right one; `--time T` selects a dump by time if needed). The verdict is taken at
t_end = 2.5 and 5; t_end 15 is for information. Also per run the `time_scheme=hesdirk2:` line.

## NOTE back
conv.py output per arm and t_end (4 levels), counts, wall. No code changes.
