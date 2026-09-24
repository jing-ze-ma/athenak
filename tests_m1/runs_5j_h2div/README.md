# runs_5j_h2div: why hesdirk2 diverged in the thin diffuse-wall shadow, and the fix

- **Branch.** `m1-h2div` from rt-integration 7c16ff86; worktree `/viper/ptmp2/jinma/wt_h2div`.
- **Run tree.** `/viper/ptmp2/jinma/h2div_0924`: `bin/`, `runs/` (CPU, 1 rank), `cpu/` (gates),
  `gpu/` (GPU job), `inp/`, `run.sh` (writes the overrides into a copy of the input, so keys
  absent from the input work), `ts.py` (max/min E per dump), `loc.py`, `l1.py` (L1 against
  the exact shadow solution, vetdo_0924/ana/shadow.py).
- **Binaries.** `base` = 7c16ff86, `old` = 9c1a12d7 (hesdirk2 before m1-h2fast),
  `fix` = b67e736a (change 1 only), `scb` = b67e736a with the SC bottom intensity hard-wired
  to the bath (the diagnostic that found change 2; not committed), `f2` = the branch (changes
  1 + 2). gcc 14 + openmpi 5, Release, MPI; GPU rocm 6.3, gfx942.
- **The case.** `/viper/ptmp2/jinma/vetdo_0924/inp/m1_shadow.athinput` (vet_sc 4x8, c = 1,
  thin medium kappa rho = 0.1, clump chi = 100, Marshak bath E_b = 1 at x1min), n = 64 x 40,
  t = 20, dt = 7.5 dx/c (hydro-CFL limited, implicit_cfl 10), 175 steps.

## 1. Reproduction (CPU, base)

| run | E_max(t) | Picard NON-CONVERGED | stage fallbacks |
|---|---|---|---|
| hesdirk2 (default) | 0.50 (t=4), 1.04 (5.5), 1.59 (9), 2.5 (12), 1.4e6 (16), **6.6e12 (20)** | 120 | 120 (every step from cycle 55) |
| be | 0.52 (4), 1.3 (5.5), 1.58-1.60 from t = 7 on, bounded | 174 of 175 | - |
| exact | 0.5 at the wall | | |

- **Onset.** Both schemes leave the exact solution at t = 4-5, in the wall cell in front of
  the clump (x = 0.02, y = 0): E rises above the bath value 1 with f = 1. From t = 6.5 the
  minimum hits the floor 1e-30. hesdirk2 and be are identical up to here (E_max 0.498 vs 0.524
  at t = 4).
- **What fails first.** The linear solve: the very first (backward-Euler) step already ends
  with a line-Jacobi residual of 5e9 after 100 Picard passes (E jumps from 1e-12 to 1), and
  under be every step is NON-CONVERGED. Under hesdirk2 the stage solves converge up to cycle
  54; from cycle 55 on every stage 1 is converged but has E < 0 somewhere, so the step is
  redone with backward Euler, and that redo does not converge.
- **Where it diverges.** Only hesdirk2 grows without bound, from t = 11: x1.58 per step (fit
  of log E_max over t = 13-20), every step a backward-Euler redo.

## 2. Bisection (CPU, base unless noted; E_max at t = 20)

| arm | E_max(20) | verdict |
|---|---|---|
| hesdirk2 default | 6.6e12 | diverges |
| implicit_predictor_order = 1 | 1.57 | bounded (like be) |
| implicit_predictor = none | 1.57 | bounded |
| all h2fast levers off (order 1, one_pass 0, stage safety 0, lin_tol_fac 1, fast_kernels off) | 1.57 | bounded |
| 9c1a12d7 (before h2fast; predictor order 1) | 1.57 | bounded |
| implicit_one_pass = 0 / time2_one_pass_safety = 0 / time2_lin_tol_fac = 1 | 6.6e12 / 6.6e12 / 6.8e12 | diverge |
| time2_enth_vel = old / implicit_enthalpy = upwind | 6.6e12 / 6.6e12 | diverge (v = 0: no-ops) |
| time2_vet_extrap = true | 2.5e14 | diverges |
| implicit_solver = bicgstab | 0.91 (be 0.97) | bounded, 75 fallbacks, but the final linear residual is 37 (mean): not solved either |
| implicit_cfl 2 (dt 2 dx/c) | 1.58 (be 1.57) | bounded, 452 fallbacks |
| implicit_cfl 0.5 (every solve converged) | 1.42 (be 1.44) | bounded but WRONG: the same f = 1 wall state |
| closure = eddington | 0.358, both schemes | stable |
| scb (SC bottom rays = the bath) | 0.504, both schemes, every dt, both solvers | stable and right |

## 3. Mechanism

Two separate faults, one in the closure and one in the time scheme.

**(a) The vet_sc boundary closes a lagged loop (every scheme, every dt).** The short-
characteristics formal solution gave the rays entering at x1min the "diffusion intensity" of
the bottom cell's M1 state, eps = E + 3 F.Omega/c, instead of the bath of the Marshak BC
(rad_m1_vet.cpp, four copies of the boundary code). In a thin wall cell that intensity is
strongly forward-peaked (eps(mu = 1) = E (1 + 3 f), up to 4 E), so the tensor from the formal
solution points along +x1 with chi -> 1, the moment solve with that tensor pushes f toward 1
and E up, and the next formal solution sees an even more peaked boundary intensity. The loop
saturates at the chi = 1 clamp: E = 1.4-1.6 E_bath at the wall (exact 0.5), f = 1, and a
floor-valued E elsewhere. Converged solves at dt = 0.5 dx/c (be or hesdirk2) reach the same
state, so this is not a time-integration or solver effect. With the bath intensity the steady
E_max is 0.504, every scheme, dt and solver agree to 5e-6, and L1(E) at t = 20 falls from
0.234 (be) / diverged (hesdirk2) to 0.055.

**(b) hesdirk2's backward-Euler redo started from the failed stage.** Stage 1 of the ESDIRK has
the internal amplification Y1/U^n = (1 + (1-g) z)/(1 - g z), which tends to -(1-g)/g = -2.41
for a stiff mode (z = dt lambda -> -inf) and has a pole at z = 1/g = 3.41 (table below). In the
wrong wall state of (a) this stage overshoots to E < 0 in some cell every step, so the step is
redone with backward Euler (correct), but:
- the stage-1 solve stores its increment in the implicit predictor BEFORE the admissibility
  test (rad_m1_implicit.cpp, the pstore block);
- under implicit_predictor_order = 2 (the hesdirk2 default since m1-h2fast) a backward-Euler
  step does not store (pskip), so the redo starts its Picard loop from
  U^n + (dt / (g dt)) x (the failed stage increment) + dt dt1 h, i.e. 3.41 x the overshoot plus
  a rate term built from two failed stages;
- the redo does not converge here (line-Jacobi, 100 passes), so its result stays near that
  start, and the next failed stage extrapolates it again: x1.58 per step to 1e12.
Predictor order 1, no predictor, or 9c1a12d7 all stay bounded, which isolates the predictor;
none of one_pass, the stage safety, lin_tol_fac, enth_vel, enthalpy or vet_extrap matters.

Amplification of one step for dy/dt = lambda y, z = dt lambda (FSAL K1 exact):

| z | stage 1 | hesdirk2 | be |
|---|---|---|---|
| -100 | -2.30 | -0.044 | 0.010 |
| -10 | -1.55 | -0.204 | 0.091 |
| -3 | -0.60 | -0.069 | 0.250 |
| -1 | 0.23 | 0.350 | 0.500 |

(The step itself is L-stable; the internal stage is not, which is why stage values of a stiff
mode can be negative and the admissibility test fires. For growing modes, Re z > 0,
|R_hesdirk2| > 1 up to z = 11.7 with a pole at 3.41, while |R_be| < 1 for z > 2: backward Euler
damps a spurious growing mode at large dt, hesdirk2 follows it.)

## 4. The fix (branch m1-h2div)

1. **A failed stage invalidates the predictor** (`Time2Restore`, rad_m1_time2.cpp:
   `pred_ok = pred2_ok = false`). The backward-Euler redo starts cold from U^n, and the next
   stage solves rebuild the history. Bitwise wherever no stage fails.
2. **`<rad_m1>/vet_bc_bath`** (rad_m1_vet.cpp, rad_m1.hpp; default true, a restart whose file
   lacks the key keeps false, echoed): at an x1 end whose `implicit_bc` is marshak the rays of
   the vet_sc formal solution entering the box carry the incident bath, eps =
   `implicit_ebath_x1min/max` (0 = the vacuum of a free surface), instead of the diffusion
   intensity of the end cell (bottom) or vacuum (top). The Milne diagnostic keeps its own
   intensity. Bitwise for flux / reflect / periodic ends and for a marshak top with ebath 0.

## 5. Results of the fix

**Shadow, n = 64, t = 20 (CPU; `l1.py` vs the exact steady-state solution).**

| run | E_max(20) | NON-CONVERGED | stage fallbacks | L1(E) | L1b(E) behind the clump |
|---|---|---|---|---|---|
| base hesdirk2 | 6.6e12 | 120 | 120 | diverged | |
| base be | 1.55 | 174 | - | 0.234 | 0.268 |
| change 1 only, hesdirk2 | 1.57 | 120 | 119 | 0.237 | 0.273 |
| **f2 hesdirk2** | **0.504** | **1** (the first step) | **0** | **0.0549** | **0.145** |
| f2 be | 0.504 | 101 | - | 0.0549 | 0.145 |
| f2 hesdirk2, implicit_cfl 2 | 0.504 | 1 | 0 | 0.0549 | 0.145 |
| f2 hesdirk2, vet_nmu 8 x nphi 16 | 0.502 | 1 | 0 | 0.0411 | 0.127 |

(The 101 NON-CONVERGED be steps of f2 are line-Jacobi at c dt/dx = 7.5, with the Picard
residual at 1e-9 to 1e-1. The state agrees with hesdirk2 and with a dt = 0.5 dx/c run to
5e-6. bicgstab converges them. That is a solver-cost issue, not this bug.)
At t = 3 (the vetdo_0924 table), L1 goes from 0.117 to 0.094 (n = 64) and from 0.110 to 0.088
(n = 128); L1b goes from 0.31 to 0.21 and from 0.29 to 0.19. Marshak (vetdo_0924 m1_marshak,
vet_sc, n = 128, t = 3): L1(E) goes from 5.04 to 4.07 % (cfl 0.4) and from 4.65 to 3.83 % (cfl 10).

**GPU (job 11966336, apudev, 1 GCD, same node, interleaved, 2 repeats; `RESULTS_gpu.txt`).**
Main-loop seconds, shadow n = 64 to t = 20:
- f2 hesdirk2: 8.16 / 7.93 (E_max 0.504, 0 fallbacks);
- base hesdirk2: 14.7 / 14.2 (E_max 5.5e12, 120 fallbacks);
- f2 be: 8.27 / 7.75;
- base be: 8.78 / 8.81.

**Gates (CPU; `RESULTS_cpu_gates.txt`, `scripts/gate_cpu.sh`).**

| gate | result |
|---|---|
| A: be named, bitwise vs base (box 1 rank; slab vet_sc 2 ranks; slab2d_nd 2 ranks; cart_sym, marshak_cart, rw_cart) | 6/6 BITWISE |
| B: hesdirk2 by default (He slab2d_nd 2 ranks, box3d_nd 2 ranks Eddington and vet_sc, marshak_cart, rw_cart) | 5/5 BITWISE, He slab and box NON-CONVERGED 0, stage fallbacks 0 |
| B: cart_sym (one natural stage fallback at cycle 1) | E and F1 equal to 2e-13 (the redo now starts cold); F2/F3 differ only at round-off (1e-13) in a 1-D-symmetric problem. The same at implicit_tol 1e-12 |
| F: a forced stage-1 failure (time2_dbg_fail), slab2d_nd / box3d_nd vet_sc / box3d_nd Eddington | fix vs base within solver tolerance (hst_cons 3.8e-9 / 4.7e-10 / 5.7e-10); at tight tolerance baseT vs fixT 8.9e-12 (box); fix is as close to the tight run as base is |
| D: restart across a forced fallback (slab2d_nd 2 ranks, rst at t = 100, failure at cycle 900) | RESTART BITWISE (rst + hst rows) |
| `tests_m1/gates/gates.py` (f2 box binary) | GATES: PASS, 14 pairs (`RESULTS_gates_py.txt`) |
| radwave G1 (runs_5f_h2fast lists/g1.txt, hesdirk2, 12 cases x 2 closures) | medians identical to base: Eddington >= 1.94, vet_sc >= 1.87 (at (100,10)); lowest single 1.58 at (100,1e3), as in base (`RESULTS_g1.txt`) |
| `tst/test_suite/rad_m1` (ATHENAK_M1_DATA = faces_0924/m1data), f2 snapshot | 3 passed (354 s) |
| cpplint on the 3 changed files | 0 errors |

## 6. Remaining limits

- Line-Jacobi does not converge at c dt/dx >= 2 in this thin 2-D problem (every be step at
  7.5). Stage failures are harmless now (the redo is cold), but the be result then sits at the
  Picard tolerance reached in 100 passes. bicgstab solves it.
- vet_sc behind the clump is still 15 % off (the uniaxial projection and the 4 x 8 rays). That
  is the closure, not the time scheme.
- Stage 1 of the ESDIRK is not positivity preserving for stiff modes (limit -2.41), so an
  inadmissible stage and a backward-Euler redo remain possible. The redo is now a clean
  fallback.

## Files

- `scripts/gate_cpu.sh` (A/B/F/D), `gcmp.py`, `rstcmp.py`, `shadow_gpu.sh`, `run.sh`, `ts.py`,
  `l1.py`, `build2.sh`.
- `RESULTS_cpu_gates.txt`, `RESULTS_gates_py.txt`, `RESULTS_g1.txt`, `RESULTS_gpu.txt`,
  `RESULTS_pytest.txt`.
