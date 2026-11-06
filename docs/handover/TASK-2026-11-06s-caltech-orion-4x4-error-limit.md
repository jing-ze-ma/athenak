# TASK 2026-11-06 Caltech -> Orion (USER standing OK): 4 x 4 column with a tighter error_limit

## Why
With t = 0 identical and the force law the same, the 1-D gap between your column and our `Lr` looks like an energy leak in Athena++: at r 24-35 Rsun, your Er falls behind ours by about 7e-7 per solve, steady at 0.05 and 0.1 d (ln(Er ours/yours) per step, dt 10.55 s), and it's about 0 at r 35-45 early. Shell budgets at r 30-35 show -0.25 to -0.41 L* in yours not accounted for by fluxes (noisy, snapshot-based). Your Jacobi takes a mean of 1.0-1.26 iterations with error_limit 1e-6. **Hypothesis, inferred, not tested:** stopping at error_limit leaves a biased residual of about 7e-7 per solve. Your total dE at 1 d is about 5x ours.

## Please run
**Same 4 x 4 deck as job 227389, with only `radiation/error_limit` changed**, in a new directory. Keep the old one:
1. `error_limit = 1e-9`, to at least 2 d, with `nlimit` high enough that it isn't reached (say so if it is; report iterations).
2. If cheap, also `1e-8` or `1e-7`, for the scaling: does the deficit go as error_limit?

## Please report (`repro3d/col4x4/orion/el<value>/`)
- `*_ke.txt` at hst cadence.
- `col4x4_prof_<t>.txt` as before, at t = 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2 d, nearest dumps. Full 256 cells, double.
- dE, dU, dEgrav against time, Jacobi iterations mean and max.
- Cost per cycle against 1e-6.
- Anything the code reports as a per-solve energy residual, if it's available without code changes.

The question: does v at 22-27 Rsun, dE and KE collapse toward our `Lr` (v at 22, 25, 27 Rsun at 1 d: -29, -109, -199 m/s; KE at 3 d 1.6e39; yours -92, -222, -287 m/s and 4.3e40)? Does the 7e-7 go away? Don't disturb `prod_seeded`. Short NOTE when first results exist.
