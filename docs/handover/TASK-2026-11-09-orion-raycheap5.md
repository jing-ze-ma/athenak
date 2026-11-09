# TASK 2026-11-09 viper -> ORION: raycheap5 -- which VET update restores 2nd order in time? (CPU, pulse2 test)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap5.md`. Follow-up to your raycheap4 (pulse2): vet_col frozen
(P2) is 2nd order, the base update (P0: vet_col built once per step at t^n, time2_vet_col = extrap2) loses order.

## Binary and decks
MOD = `fork/raycheap2-1109` head `fa7103f33` (same code as ff2ed77f1, which you built for raycheap4: reuse it).
Decks `tests_raycheap/upd/upd_<arm>_T<t_end>_c<cfl>.athinput` = the pulse2 P0 decks (ramped start, runs ending at
t_end = 2.5 and 5.0, cfl 0.0375 / 0.01875 / 0.009375 / 0.0046875, 24 MeshBlocks = 24 ranks) with ONE change:

| arm | change vs P0 | meaning |
|---|---|---|
| U1 | time2_vet_col = rebuild | a predicted build for stage 1 + a second build at the stage-1 solution |
| U2 | time2_vet_x2_fields = 15 | extrap2 of all hand-over fields (default 9: f_K + q only; + lateral slots, vgd_hr) |
| U3 | time2_vet_col = lag | D(U^n) for both stages (no extrapolation) |
| U4 | time2_vet_col = extrap | D^n + (dt/dt_prev)(D^n - D^{n-1}) (diagnostic linear extrapolation) |
| U5 | time2_vet_col = predict | ONE build per step at the predicted state U^n + dt K1 |
| U6 | vet_col_lat_offdiag = none | DIAGNOSTIC only (lateral off-diagonal coupling off; not a production option) |
| U7 | time2_vet_col = lag, vet_col_every = 2 | build every 2nd step |

(The decks already have vet_gd_tangential = false; the tangential tensor is not part of this test.)

## Analysis
Per arm and t_end: `python3 tests_raycheap/conv.py <c0.0375> <c0.01875> <c0.009375> <c0.0046875>`; compare with your
P0 / P2 numbers from raycheap4. Report d, p (both pairs), wall per run (the cost of rebuild / predict matters), and the
`time_scheme=hesdirk2:` counts (fallbacks would spoil an arm, as in P3).

Also, from your existing logs (no new runs): for raycheap3 C1 and raycheap4 P2 (both "vet_col frozen"), the lines
`<rad_m1> closure = vet_col: N tensor builds` and `<rad_m1> vet_col_lat: N sweeps` / `vet_gd: ... sweeps` -- we need
to know whether "frozen" froze the vet_gd lateral part too, or only the column tensor.

## NOTE back
Per-arm conv.py lines (T 2.5 and 5), walls, counts, the C1/P2 build counts. No code changes.
