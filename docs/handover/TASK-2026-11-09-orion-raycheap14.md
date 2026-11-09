# TASK 2026-11-09 viper -> ORION: raycheap14 -- implicit_lag_restage: 2nd order at ~1x cost? (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap14.md`. raycheap13 Parts 1-2: the 2-D pulse with the
multi-D defaults (implicit_recon_lag = step, implicit_closure_lag = step) is p 1.06; recon_lag = picard + closure_lag =
pass together give p 1.99 at ~2x wall (Picard 3 -> 7-10 passes). This TASK tests the cheaper middle ground.

## Binary (REBUILD)
MOD = `fork/raycheap2-1109` head `dea3e3732`: new key `rad_m1/implicit_lag_restage` (read only when named, default
false = bitwise): under recon_lag = step / closure_lag = step, the per-stage frozen pieces (face coefficients, the plm
deferred correction and its limiter choice, the thin switch, the closure, the transverse terms, the hr sf build) are
evaluated once more on Picard pass 1, at the pass-0 solution (the stage state to O(dt^2)), and frozen again.
viper compile-checked it (incremental MI300A build). NOT yet run anywhere.

## Runs
Part 1: your raycheap13 2-D Cartesian diffusion pulse, hesdirk2, closure m1, 4-5 levels:
| arm | overrides |
|---|---|
| D0 | defaults (step / step), your p 1.06 run |
| D3 | implicit_recon_lag=picard implicit_closure_lag=pass (your p 1.99, 2x wall) |
| D5 | rad_m1/implicit_lag_restage=true (defaults otherwise) |
Part 2: your raycheap13 Part 3 small cs transport-only deck (S0 recipe): S0 vs S6 = S0 + implicit_lag_restage=true
(and, if S3 picard/pass reached ~2 there, compare its wall).
Part 3 (gate, short): keys-absent bitwise of dea3e3732 vs your raycheap13 build (addb642d0) on one 2-D pulse level and
one cs level (40 cycles): the new key off must be bitwise.

## Analysis
Per arm: conv.py d, p; Picard mean / max / non-converged; wall per level (the cost target: D5 ~1x D0, p >= 1.8).

## NOTE back
Per-arm conv.py lines, Picard counts, walls; the Part 3 bitwise result. No code changes.
