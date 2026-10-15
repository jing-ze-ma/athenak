# NOTE 2026-10-09 (viper -> Delta, DeltaAI, Caltech, Orion): VET closure source bug + pre-production bug hunt

## 1. Bug found (beam worker, branch rad-beam-1008 @ 28d3ea99, not merged)
The VET closure source (vet_col, vet_col_lat, vet_gd) reads the implicit solve's OLD VECTOR
`EN = E^n + dt*(chat/c)*esrc` (rad_m1_implicit.cpp ~8146, m1_impl_i0) instead of the physical E^n. In he_star_m1 the
MLT++ scaffold deposits its flux divergence through esrc; where F_MLT ends just below the photosphere the closure source
S/E reaches 6-11, the vet_gd formal solution over-emits (~3.6 L0 above R_ph) and the Eddington tensor shape is distorted
near the photosphere (a plausible contributor to the thin-atmosphere problems). Fix: opt-in key
`rad_m1/vet_source_noesrc = true` (time_scheme = be; default off = bitwise); at the IC step S/E -> 1.00-1.08 and the FS
luminosity matches M1 below R_ph. Affected: any he_star_m1 run with closure = vet_col/vet_gd while the MLT scaffold
(esrc) is active -- AG Car A/B during the ramp, BSG runs with a scaffold, the He giant during its scaffold phase.
Caltech (He giant N897): please check in your run's input whether the MLT scaffold / esrc is still active at 10.4 d+;
if it is, tell viper (the fix would have to go in at a link boundary; NOTE: switching central->blend on a restart
diverges, see NOTE-2026-10-08-viper-hegiant-blend on hegiant-opn-1007; vet_source_noesrc alone has not been tested
on a restart yet).

## 2. Pre-production bug hunt (viper, running now)
Before the AG Car productions start (~18-20 h on Raven), two viper workers audit the production path: a code review for
the same class of bug (physical state read from solver work arrays with folded-in terms; double-applied sources; silent
fallbacks -- e.g. the vet_col_lat positivity fallback drops the lateral VET term for the whole step and is not printed
before 9d13ed76; restart key/cache inconsistencies; units/frames; sp metric factors) and dynamic consistency tests
(restart = continuous, energy budget by term, formal-solution consistency, symmetry, determinism, key echo).
Results will be pushed as NOTE-2026-10-09-viper-audit-*.md on rt-integration.

## 3. Delta
Keep running the radiation test batches (TASK-2026-10-08-delta-test-runner.md). Audit test batches may be added on
branches named in a later update of that TASK; same budget rules. No production on Delta (user).
