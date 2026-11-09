# TASK 2026-11-09 viper -> ORION: raycheap6 -- gas-driven closure change: which VET update is 2nd order? (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap6.md`. Follow-up to raycheap5: there, rebuild/predict were
bitwise = lag because the pulse2 test's only time dependence (the ramped inner-boundary intensity) was not advanced by
the predicted build (gas held, kappa_s = 0: the formal solution saw identical inputs). This TASK (b) moves the time
dependence into the GAS, as in production, and (a) re-tests pulse2 with the predicted build now seeing the boundary at
t_n + dt.

## Binary
MOD = `fork/raycheap2-1109` head `f872899da` (code = a1481b185; later commits are inputs only) -- REBUILD (test-only changes: predicted builds take the test inner flux at
t_n + dt; new key `problem/clump_tmod_amp|period`: the held clump temperature T_c (1 + A sin(2 pi (t+dt)/P))^(1/4)).
viper compile-checked it (incremental MI300A build) and smoked it on one rank.

## (b) gas-driven test: `tests_raycheap/gasmod/gas_<arm>_c<cfl>.athinput`
The raycheap3 forced test star (24 MeshBlocks: 24 ranks), CONSTANT inner flux, clump emission modulated with A = 0.5,
P = 400, t_end 400, implicit_tol 1e-13; cfl 0.15 / 0.075 / 0.0375 / 0.01875 (300 ... 2400 cycles; ~1 h finest).

| arm | time2_vet_col | tangential (full tensor) | note |
|---|---|---|---|
| L | lag | off | D(U^n), the 1st-order reference |
| X9 | extrap2, x2_fields 9 (default) | off | production HG/RSG form |
| X15 | extrap2, x2_fields 15 | off | all hand-over fields extrapolated |
| PR | predict | off | one build per step at U^n + dt K1 |
| RB | rebuild | off | two builds per step (most accurate reference) |
| X9t / PRt / RBt | as X9 / PR / RB | ON (vet_gd_tangential = true) | the user's criterion: the full tensor |

## (a) pulse2 re-test: `tests_raycheap/upd/upd_{U1,U2,U3,U5}_T{2.5,5.0}_c*.athinput` (unchanged decks, new binary)
U1 rebuild, U2 x2_fields 15, U3 lag, U5 predict; t_end 2.5 / 5; cfl 0.0375 ... 0.0046875.

## (c) the VET source relaxation (viper's new prime suspect): `tests_raycheap/srcgas/`
Both the test and production use `vet_source = relaxed_esc` (default): the formal-solution source is relaxed over the
step, so the closure depends on dt itself (an O(dt) term even for a frozen closure, whose single build used dt0 -- this
would explain raycheap3 C1's p = 1.0). The srcgas decks add `vet_source = gas` (unrelaxed B(T), dt-independent):
- `sg_pulseP0_T{2.5,5.0}_c*` (= pulse2 P0 + gas source; 4 levels as (a));
- `sg_gasX9_c*`, `sg_gasX9t_c*`, `sg_gasPRt_c*` (= gasmod X9 / X9t / PRt + gas source; 4 levels as (b));
- `sg_forcedB0_c*` (= raycheap3 B0 + gas source; levels 0.15 ... 0.01875; ~1 h finest).
Expected if the suspect is right: these reach p ~ 2 where their relaxed twins stay at ~1.

## Analysis and PASS
Per arm (and t_end): `python3 tests_raycheap/conv.py <coarse> ... <fine>`. Report d, p, wall per run, the
`closure = vet_col: N tensor builds` counts (cost), and the `time_scheme=hesdirk2:` counts.
PASS (user's VET criterion) for an arm: finest-pair p >= 1.8 in (b) WITH the tangential tensor (X9t / PRt / RBt).
We want the cheapest arm that passes. In (a), predict / rebuild should now differ from lag (report max|dE| vs U3).

## NOTE back
conv.py lines per arm, walls, build counts, counts. No code changes.
