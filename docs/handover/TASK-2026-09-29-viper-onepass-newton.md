# TASK for viper: validate `implicit_one_pass` ON TOP OF `implicit_opac_newton` (He box) before any default change

**User 09-29:** DeltaAI measured a further ~10 % from `implicit_one_pass = 8` now that Newton is on. It stays
**opt-in** until viper has re-run the ke-dt gates with it. The user remembers the reasons it was switched off, and
they were measured WITHOUT Newton:
- the ke-dt recipe (09-27) had `one_pass = 0` next to the converged opacity update;
- it changes results by 5-22x the round-off spread;
- at cfl 0.9 and on stiff wedges nothing was accepted, which cost +4.6 %;
- audit-m1-0928 found -2.8..0 % on MI300A.

Without Newton, q = res_1/res_0 ~ 0.17, so the test `res_0 qe/(1-qe) < tol` with `qe = 30 q` could almost never
pass. With Newton q ~ 4e-5.

## DeltaAI result (GH200; branch m1-perf-0928 6432a25b, section "implicit_one_pass on top of Newton" of NOTE-2026-09-29-m1-perf.md)

He box 84x104x104, 1 GPU, cfl 0.3, 300 cycles, fresh start:

| setting | ms/cycle | Picard/solve |
|---|---|---|
| Newton (default) | 47.71 | 2.02 |
| + one_pass 8, time2 safety 30 (the input's value) | 44.55 (-7 %) | 1.33 |
| + time2 safety 10 | 42.93 (-10 %) | 1.18 |

Order test from the evolved t = 40 restart (t → 80, cfl 0.6/0.3/0.15/0.075, safety 30):
- the orders equal Newton-only to 3 digits: dens 1.98/1.89, ener 1.97/1.86, E 1.47/2.02;
- ||newton - onepass|| is 1e-4..5e-4 of the time error;
- no ringing (R <= 0.053).

## Asked of viper (the ke-dt run tree and gates, `docs/dev/ke_dt_0926.md` sections 7-8)

Arms: Newton default vs Newton + `implicit_one_pass = 8` (auto on). Use `time2_one_pass_safety` 30, and 10 if 30
passes. Build: rt-integration at or after d26b7364.
1. Gate (a): gamma(KE_h) over 4000-5000 s from the t = 3800 restart at cfl 0.075 / 0.15 / 0.3 / 0.6 / 0.9. It must
   stay dt-independent (2.15e-4 to 3 digits, as in section 7).
2. Gate (b): the order table 0.8→0.4 / 0.4→0.2 (runs_ord), all variables, both arms.
3. Cost on MI300A at cfl 0.3 and 0.9, 1 and 2 GPUs: wall per simulated second, Picard/solve, and the one_pass line
   (accepted / measurements / auto switch-offs). At cfl 0.9, check whether anything is accepted and whether
   auto switches it off.
4. Optional: the moving wedge (gateC) cost, both arms.
5. Push a NOTE with the table and a verdict. If gates (a) and (b) pass and the gain is >= ~3 % on MI300A, propose
   the default: `implicit_one_pass = 8` wherever `implicit_opac_newton` is on. The user decides.
