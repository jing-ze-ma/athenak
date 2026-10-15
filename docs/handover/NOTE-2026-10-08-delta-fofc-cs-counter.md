# NOTE 2026-10-08 (Delta -> viper): fofc_cs test -- what it counts is not the floor the cubed sphere applies

Follow-up to NOTE-2026-10-08-delta-bringup §3. This supersedes the "FMA mismatch between the floor-TEST and
RaiseVel" reading there. It also bears on viper's "shared non-inlined floor-test/RaiseVel arithmetic" option.
Read-only analysis at rt-integration 3272e5c4; **no code changed**. The fix is viper's decision.

## Findings
1. **The counter the test asserts on records nothing that is applied on the cubed sphere.**
   - `test_hydro_fofc_cs_cpu.py:81` asserts that `eos_efloor` == 0.
   - On the cubed sphere `eos_efloor` is incremented only by `IdealHydro::ConsToPrim` (`eos/ideal_hyd.cpp:91`).
   - Under `defer_cons_floors` (always on for the cubed sphere, `eos/eos.cpp:165`), that check uses the
     ORTHONORMAL kinetic energy and writes nothing (`eos/ideal_c2p_hyd.hpp:57`:
     `if (!eos.defer_cons_floors) u.e = ...`). `cons(IEN)` is rewritten with its own value, and `w0(IEN)` is
     overwritten later by RaiseVel.
   - So the 309 events (Cray FMA build) and the 232917 events (control run) are cells where the wrong-KE check
     dips below the floor. They are not floors applied to the solution.
2. **The floor that IS applied is in a kernel FOFC does not mirror, and that kernel counts nothing.**
   - With no floor switch set (`floors_legacy`, the default; this test), the real pass runs THE LEGACY KERNEL,
     `coordinates.cpp:840-911`. The floor is at `:905`, `ApplyEnergyFloor`.
   - The FOFC floor-TEST calls `GnomonicRaiseVelFloors` instead (`hydro_fofc.cpp:187`), a separate copy. The
     comment at `coordinates.cpp:840` itself says the two are "algebraically the same ... but not bitwise" under
     contraction.
   - The legacy kernel accumulates neither `efloor_de` nor a count. `efloor_de` is 0 in every run, including
     a control run forced to floor with `hydro/pfloor=0.3`, where the true floor certainly fires.
3. **Consequence:**
   - The FMA failure is the ConsToPrim wrong-KE check on the real post-update state disagreeing with the same
     check on the FOFC trial state. The trial state is pre-source-term and FMA-contracted differently. Without
     FMA the two happen to agree.
   - Nothing in the test, and no event-log column, measures whether FOFC removes the cubed-sphere floors that
     actually change the solution on the default (legacy) path.
   - Also: FOFC's own test-mode `ConsToPrim(utest_)` (`hydro_fofc.cpp:141`) flags cells on that same wrong-KE
     check. Some of the ~1e5 FOFC firings in this test may be spurious first-order downgrades. Not measured.

## Evidence (Delta CPU, rt-integration 081dd60b; logs in /work/nvme/bivj/jma20/delta_1008/trace_fofc_cs/)
| build | fofc=true eos_efloor | control eos_efloor | efloor_de (all runs) |
|---|---|---|---|
| Cray CC (-march=znver3, FMA) | 309 (264 by cycle 5 at fofc count 0; 45 at cycle 46) | 232917 | 0 |
| Cray CC + -ffp-contract=off | 0 | 239660 | 0 |
| plain g++ (no -march) | 0 | – | 0 |
| no-FMA, control, pfloor=0.3 | – | 181501 by cycle 10 | 0 |

## Suggested fix, in order (viper decides)
1. **Count the real cubed-sphere floor.**
   - Add a count to the legacy RaiseVel kernel: a new counter, or `efloor_de`.
   - To keep the legacy kernel bitwise (its stated purpose), do the count in a separate read-only pass run only
     when the event log is written, not by turning that par_for into a reduction.
2. **Fix what is reported and tested.**
   - On the cubed sphere, stop ConsToPrim adding its deferred-mode checks to `eos_efloor`. Either drop them or
     move them to a separate diagnostic column.
   - Point `test_hydro_fofc_cs_cpu` at the real count from (1).
   - Add an FMA build (e.g. `-march=native` or `-ffp-contract=fast`) to the CPU test matrix.
3. **Decide the floor-TEST mirror.**
   - Either route the legacy path through the shared `GnomonicRaiseVelFloors`, which costs bitwise reproduction
     of legacy runs, or accept a rounding-level residual.
   - Separately: should FOFC's test-mode ConsToPrim still flag on the wrong-KE check on the cubed sphere?
   - A floor-test margin, which I proposed earlier on Delta, is NOT recommended: it would only add first-order
     cells against the wrong-KE check.

Delta keeps the CPU test build at the Cray defaults until viper decides. Not a blocker for the AG Car runs:
they are spherical-polar, not cubed-sphere.
