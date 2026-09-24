# runs_5k_h2acc: the remaining accuracy defects of time_scheme = hesdirk2

- **Branch.** `m1-h2acc` from rt-integration d2572b4f; worktree `/viper/ptmp2/jinma/wt_h2acc`.
- **Work dir.** `/viper/ptmp2/jinma/h2acc_0924`: `bin/` (CPU, gcc 14 + openmpi 5, Release,
  MPI), `cpu/`, `rw/`, `runs/` (shadow), `lists/`, `scripts/` (copied to `scripts/` here).
- **Binaries.** `base` = d2572b4f; `new` = d2572b4f + `vet_gas_start.patch` (the tensor-timing
  change of sect. 2); `exp` = `new` + `exp_gamma3.patch` (the ESDIRK family of sect. 1,
  `<rad_m1>/time2_dbg_g3`). Both patches are kept here as files; **neither is committed**.
- **Verdict.** No code change. Both candidate fixes were measured and are worse than
  today's hesdirk2. The branch carries this README, the analysis script and the results only.

## 1. Stage-1 stability: a structural limit, not a tableau choice

Today's tableau (rad_m1_time2.cpp): explicit Heun hydro paired with the stiffly accurate
ESDIRK `A = [[0,0,0],[1-g,g,0],[1/2,b2,g]]`, g = 1 - 1/sqrt 2, b2 = 1/2 - g, c = (0, 1, 1),
FSAL K1. That makes it an H-ESDIRK: the first stage is explicit, and K1 = f_I(U^n).

For y' = lambda y, z = dt lambda:
- stage: R1(z) = (1 + (1-g) z)/(1 - g z). R1(-inf) = -(1-g)/g = **-2.414**, with a pole at z = 3.41;
- step: R(-inf) = 0 (L-stable), and |R(iy)| <= 1 (A-stable).

**Why no 2-solve tableau removes it.** The explicit part is fixed to Heun (c~ = (0,1,1),
b~ = (1/2,1/2,0)), and the scheme must be stiffly accurate with 2 implicit stages. The IMEX
order-2 coupling conditions (b~.c = 1/2, b.c~ = 1/2, b.c = 1/2, sum b = 1) then force:
- c2 = 1 and b1 = a31 = 1/2 (the explicit K1 always carries weight 1/2);
- a21 = 1 - g2 and a32 = 1/2 - g3, with two free diagonals g2 and g3.

The limits are then:
- stage: R1(-inf) = -(1-g2)/g2;
- step: R(-inf) = -(1/2 + (1/2 - g3) R1(-inf))/g3.

L-stability needs **g2 = (1 - 2 g3)/(2 (1 - g3)) < 1/2**, so |R1(-inf)| > 1 for every
L-stable member. The stage value of a stiff mode always overshoots, and its sign flips.
- A positive stage (g2 = 1, backward Euler) gives R(-inf) = -1/(2 g3): not L-stable.
- A stage-order-2 stage (g2 = 1/2) gives R(-inf) = -1 for any g3.

Choosing g3 only trades the overshoot for accuracy (`stab.py`, `stab.txt`):

| g3 | g2 | R1(-inf) | R(-inf) | max\|R(iy)\| | PR order at lambda = -1e4 | phase error at y = 1 |
|---|---|---|---|---|---|---|
| 0.2929 (today) | 0.2929 | -2.41 | 0 | 1 | 1.00-1.06 | -3.75e-2 |
| 0.20 | 0.375 | -1.67 | 0 | 1 | 1.00-1.07 | -4.20e-2 |
| 0.15 | 0.412 | -1.43 | 0 | 1 | 1.01-1.08 | -4.75e-2 |
| 0.10 | 0.444 | -1.25 | 0 | 1 | | |
| SDIRK2 (not pairable with Heun: b.c~ = g) | | 0 | 0 | 1 | 0.97-1.05 | |

**Order reduction is not a stage-1 artefact.** On the Prothero-Robinson problem every member
drops to order 1 at intermediate stiffness (lambda = -1e4 and -1e6), and so does SDIRK2, which
has a positive stage. The cause is stage order 1, which all of them share. Stage order 2 at
c2 = 1 needs g2 = 1/2, which gives R(-inf) = -1. The 1.2-1.3 orders of the runs_5h thin-top
relaxation are this effect, with an error ~ dt/|lambda|. TV2 below gives the same errors for
vet_col (1.24 / 1.22 / 1.30) as in runs_5h.

**Measured, g3 = 0.15 (`exp`).**
- Radwave G1, Eddington (`RESULTS_g1.txt`): medians 1.87-2.07.
  - (1, 1e3) falls from 1.94 (min 1.82) to 1.87 (min 1.23), and its e1024 doubles (6.6e-6 ->
    1.2e-5).
  - (10, 1e3) e1024 goes from 1.1e-5 to 1.7e-5.
  - (100, 10) improves from 1.99 to 2.07.
- Stiff corner P = 100, tau 1e5, 40 periods (`RESULTS_stiff.txt`): the same as today. Bounded;
  end 0.715 / 0.799 / 0.802 vs today's 0.719 / 0.800 / 0.801 (exact 0.799).
- Shadow (runs_5j reproducer, all binaries): 0 stage fallbacks.

It is not at least as accurate everywhere, so it is not adopted and no `hesdirk2_old` is
needed. (The LE / Marshak exp arms were not run: runlist overrides of a key absent from the
input are fatal. They were not needed after G1.)

## 2. The vet_col held-gas steady-state mismatch: hesdirk2 is right, be carries the O(dt) error

`new` built the stage-1 formal solution the way be does: E^n with the stage-start gas (T and
opacities after the Heun predictor) instead of T^n. This is `time2_vet_gas = start`,
patch here.

**Steady states** (`RESULTS_steady.txt`, held-gas atmospheres, 3000 / 12000 steps).
max|dE|/max E:

| case | h2 base vs be, cfl 0.3 / 0.075 | h2 new vs be, cfl 0.3 / 0.075 | own cfl 0.3 vs 0.075: be / h2 base / h2 new |
|---|---|---|---|
| vatm (sp wedge, vet_col) | 1.93e-4 / 5.25e-5 | 1.17e-5 / 8.1e-7 | 1.36e-4 / **4.4e-6** / 1.47e-4 |
| vstr (stretched r) | 5.88e-5 / 1.54e-5 | 1.33e-6 / 1.6e-7 | 3.95e-5 / **4.2e-6** / 4.06e-5 |
| vpp (Cartesian Milne, vet_col) | 3.45e-5 / 8.75e-6 | 2.16e-6 / 1.3e-7 | 2.57e-5 / **1.7e-7** / 2.78e-5 |
| vsc (the same, vet_sc) | 3.45e-5 / 8.75e-6 | 2.16e-6 / 1.3e-7 | |

- `start` makes hesdirk2 match be, but only because it inherits be's O(dt) error.
- Today's hesdirk2 steady state is almost dt-independent (4e-6 or less between cfl 0.3 and
  0.075). be and `start` move by 1.4e-4.
- The held-gas test is not a hydro steady state: the hydro step heats the gas and the hold
  undoes it. So the post-hydro T that be reads is one step of that heating away from the
  held T. T^n (hesdirk2) is the right source temperature, and the 1.9e-4 was be's error.

**Transient and order** (`RESULTS_relax.txt`, `RESULTS_g1.txt`):
- TV (vet_col with absorption relaxing, t = 6.4): errors 8.9e-6 .. 8.9e-7 (base) vs 4.1e-5 ..
  4.5e-6 (`start`, the same as be's 3.7e-5 .. 4.2e-6).
- Radwave G1 vet_sc: `start` falls to first order: medians 1.04-1.08 at tau = 10 (base
  1.87-1.97) and 1.61-1.70 at tau = 0.1.
- T1 vet_col transient: 1.98-1.99 with `start` (no gas there).

So `start` was rejected; the timing stays T^n. The only way to make be agree would be to give
be the T^n tensor, and be is first order anyway.

## 3. Other stiff-limit points

- The enthalpy velocity at P = 100, tau >= 1e5 is bounded with `time2_enth_vel = central`
  (the default since runs_3x). `exp` changes nothing there (sect. 1).
- The remaining stiff-corner inaccuracies are:
  - NC or stage fallbacks at tol 1e-11 when the Picard loop stalls (runs_3x);
  - the P = 10, tau = 1e5 end amplitude at nt 48.

  Both are solver tolerance and stage-order-1 effects, not stage-1 positivity. Stage
  fallbacks stay a clean backward-Euler redo (runs_5j).

## 4. Gates on the `new` binary (they apply to today's code for everything `start` did not touch)

`RESULTS_gates_cpu.txt`:
- be-named inputs (box, vet_sc slab, nd slab, cart_sym, marshak, rw_cart, vet_col atm) vs base:
  7/7 BITWISE.
- hesdirk2 Eddington (nd slab, 3-D box, cart_sym, marshak, rw_cart): 5/5 BITWISE.
  - He slab and box: NON-CONVERGED 0 and fallbacks 0.
  - cart_sym: 1 fallback, as before.
- Restart, slab Eddington and vet_sc: BITWISE.
- `tests_m1/gates/gates.py`: PASS.
- `tst/test_suite/rad_m1`: 3 passed.
- (The hesdirk2 vet arms of gate B aborted on the command-line key; sect. 2 compares them
  directly.)

No GPU timing: nothing is changed.
