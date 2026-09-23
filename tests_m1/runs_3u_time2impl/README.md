# runs_3u_time2impl: H-ESDIRK2 for implicit rad_m1 -- STOPPED at the pre-implementation model check

Branch `m1-time2` (worktree /viper/ptmp2/jinma/wt_time2, base 51a9adb2). No AthenaK source changed.

## What was checked

The brief: keep the `(4/3) v E` enthalpy term implicit inside the solve, and check with the
model that H-ESDIRK2 is still order 2 that way.

What the code does now (`rad_m1_implicit.cpp`):
- `M1_IW_V1..V3` are set ONCE from the gas momentum at the start of `ImplicitSolve` (4533-4536).
- They are never updated in the Picard loop.
- The enthalpy flux `a E'` (`a = v + v.D`, 4841) is therefore implicit in `E` only. The
  velocity is lagged at the solve's "old" gas state.

In the linear model the enthalpy row is `-(4/3) E0 i k v`, which depends on `v` alone. So the
code's in-solve form is EXPLICIT in the model's sense. Inside an H-ESDIRK2 stage the old state
is `rhs_i`, and that form is `v* = v(rhs_i)`.

`enth_lag_model.py` is H-ESDIRK2 with the FSAL slope carried exactly (10x10 step map). It
covers five treatments:
- `imp`: `v* = v(Y)`, fully implicit;
- `rhs`: the code's lagged `V1`;
- `pred`: `v(rhs) + g dt K_prev`;
- `exp`: Heun-weighted explicit, design option 1;
- `code`: the present Lie-BE.

`picard_v_contraction.py` measures the gain of a Picard pass that re-evaluates `v` from the
previous pass's force.

## Results

Source files: `RESULTS_enth_lag_model.txt` and `RESULTS_picard_v_contraction.txt`.

**Order**, p(1024->2048):

| variant | order |
|---|---|
| imp | 1.97-2.00 in all 12 cases |
| rhs (the code's form) | 1.00 at (1,1e3), (10,1e3), (100,10), (100,1e3); 1.22 at (0.1,1e3); 1.57 at (10,10) |
| pred | 1.97-2.00 |
| exp | 1.98-2.00 |

**Stability, gas pressure off (as C3)**: the maximum over k up to Nyquist (N = 64) and
tau = 1e2-1e6 of rho(G) - 1.

| cfl | P | imp | rhs | pred | exp | code |
|---|---|---|---|---|---|---|
| 0.15 | 1 | 0 | 1.2e-2 | 0 | 0 | 0 |
| 0.15 | 10 | 0 | 0.22 | 0 | 0 | 0 |
| 0.15 | 100 | 0 | 1.5 | 0.75 | 0.71 | 14 |
| 0.30 | 10 | 0 | 0.72 | 1.4 | 1.5 | 3.5 |
| 0.30 | 100 | 0 | 38 | 200 | 140 | 65 |

**Picard pass with `v` re-evaluated per pass**: the largest gain over k, and over tau = 0.1-1e6,
at H-ESDIRK2's `g dt`:

| cfl | P = 1 | P = 10 | P = 100 |
|---|---|---|---|
| 0.15 | 0.014 | 0.15 | 1.5 |
| 0.30 | 0.054 | 0.59 | 6.0 |

At P = 100 the pass diverges even at cfl 0.15.

## Consequence

- **Only `imp` is both order 2 and bounded.** That needs `v'` inside the linear system: the
  force-dependent `v'` enters the upwinded `a(v') E'`. The E row then couples to cells two away,
  which is beyond the present tridiagonal / 7-point operator. That is a change to the flux
  assembly and the solver, which is exactly the area the m1-space2 branch is editing.
- **Keeping the code's lagged `V1` gives G1 < 1.9** in 4-6 of the 12 Eddington cases (p = 1.00).
  In the continuum symbol it is also less stable than the present code at P <= 10.
- **Plain Picard on `v` is not an option.** It diverges at P >= 10 (cfl 0.3) and at P = 100
  (cfl 0.15).

A decision on the enthalpy treatment is needed before implementation.
