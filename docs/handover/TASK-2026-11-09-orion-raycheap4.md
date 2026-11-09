# TASK 2026-11-09 viper -> ORION: raycheap4 -- forcing-free time-order test of the base VET/M1 scheme (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap4.md`. Follow-up to raycheap3 (C1 frozen vet_col and B0
stay 1st order, so the 1st-order term is in the transport stepping or in our forcing). viper read the forcing code:
the inner flux is set in Time2FormStage before each hesdirk2 stage solve at t_n + dt (both hesdirk2 stage solves end
at t_n + dt; the mesh time advances only after the stages, driver.cpp:918 vs :763), and the K1 slope carries the
t_n value from the previous step. Nothing wrong found, so this test removes the forcing altogether.

## Binary and decks
MOD = `fork/raycheap2-1109` head `5e4b46412` (same code as your 9fbd2d3c6 build: inputs only were added; reuse it).
Decks `tests_raycheap/pulse/pulse_<arm>_c<cfl>.athinput`: the raycheap3 test star (24 MeshBlocks: 24 ranks) with a
CONSTANT inner flux (mod_amp 0) and the empty initial radiation field (e_out 1e-8), t_end 30: the radiation fills the
domain from the star (light crossing ~5) and relaxes to the steady state with the hot clump; implicit_tol 1e-13.
Levels cfl 0.15 / 0.075 / 0.0375 / 0.01875 / 0.009375 (dt 1.34 ... 0.084; 22 to 360 cycles, minutes each).

| arm | change | expected |
|---|---|---|
| P0 | base (hesdirk2, vet_col, time2_vet_col extrap2) | p = 2 if the base scheme is 2nd order |
| P1 | time_scheme = be | p = 1 (reference for the order measurement) |
| P2 | vet_col frozen (time2_vet_col = lag, vet_col_every = 1e6) | isolates the transport stepping |
| P3 | closure = eddington + implicit_hr_model = true (your C2e) | the no-VET baseline |

## Analysis
Per arm: `python3 tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> <c0.01875> <c0.009375>` (d and p per pair).
Also per run the `time_scheme=hesdirk2:` line (stage / BE / fallback counts). If the coarsest level fails or
falls back, report it and use the 4 finer levels.

## NOTE back
conv.py output per arm (5 levels), counts, wall. No code changes.
