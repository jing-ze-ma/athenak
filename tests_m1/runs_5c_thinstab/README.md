# runs_5c_thinstab: the multi-D instability of implicit M1 (chi(f) closures) in thin cells

- **Date:** 2026-09-24, viper. **Branch:** `m1-thinstab` (from rt-integration 304fb38f),
  worktree `/viper/ptmp2/jinma/wt_thinstab`.
- **Run tree:** `/viper/ptmp2/jinma/thinstab_0924` (`runs/` CPU, `gpu/` apudev,
  `gate_thick/`, `he/`; binaries in `bin/`: `s0` = 304fb38f + the pgen seed, `s1` = this
  branch, `thin` = m1-thin cb7ea901 + the pgen seed, `sp2ctr` = m1-sp2 b5f3afea + the
  switch). Scripts: `scripts/`, `spread.py`, `top.py`, `vs1d.py`, model in `model/`.
- **Source changes:**
  - `src/pgen/tests/rad_m1_tests2.cpp`: `m1_test = atmosphere` gets `atm_seed`,
    `atm_seed_k` (a transverse E seed; default 0 = IC untouched).
  - `src/rad_m1/rad_m1.hpp`, `rad_m1_implicit.cpp`: the switch
    `<rad_m1>/implicit_closure_thin_relax = C` (default 0 = off; read only when named;
    nothing allocated). About 80 lines. No Cartesian expression is touched.

## Verdict

1. **Mechanism: the scheme, not the PDE.** With a lagged chi(f) closure the F-dependent
   part of div P is explicit:
   - dP11/dF1 = chi' and dP12/dF2 = b/f;
   - this is an explicit advection of divergence-free flux perturbations at c dt/dx,
     damped only by the implicit c dt rho kappa;
   - so for c dt >> dx, |g| ~ max(chi', b/f)/tau_cell per step, independent of dt;
   - 1-D is stable because div F = 0 pins dF1 there.

   Measured on a uniform scattering slab: g = 1.18-1.84 per step at tau_cell 0.125-0.03.
   It decays at tau_cell >= 0.25. The model's local bound is 2.6-18 (section 2), and the
   converged-closure (implicit) scheme gives |g| <= 0.2.
2. **Existing switches do not cure it.**
   - `trans_limit lp` (+ offdiag none, the 09-24 choice): 1.10-1.18.
   - m1-thin `implicit_thin = hll`: off at c dt >> dx, and worse when forced on.
   - `closure_lag = pass` with relax or Anderson: only down to tau_cell 0.125, at 5-7x
     the outer passes, and on the GPU 4.1x the time of (d).
3. **Cure: `implicit_closure_thin_relax = 1.5`** (new, default off). It relaxes the
   step-lagged closure from step to step with a tau-dependent weight.
   - The seed decays or stays <= 1e-5 over 500 steps at tau_cell 0.006-0.5 (CPU and GPU).
   - 0 positivity fallbacks.
   - The fixed point is unchanged: the run converges to the 1-D M1 solution.
   - Cost: no extra passes; GPU 12.4 s vs 51.4 s for pass + Anderson (128x64 slab).
   - Limits: cells with tau << 1 adapt their closure only slowly (a quasi-frozen tensor,
     3 % from M1 at tau_cell 0.006). Not cured for f > 0.69 in thin cells.
4. **The 09-22 He-slab thin-top damage: yes, the same mechanism.** Seeded He slab
   (runs_3b7 A0 input, 50 s, current code):
   - pass: 1848/1848 steps NON-CONVERGED at 200 passes;
   - step: 332 positivity fallbacks;
   - step + switch: 0 fallbacks, 4.1 passes/step, 0 NON-CONVERGED.

## 1. Reproducer and growth table

`inp/uni.athinput` (`m1_test = atmosphere` with `atm_scale_h = 1e8`, i.e. uniform):
- 32 x 16 cells, x1 in [0, 4], x2 in [0, 2], so dx1 = dx2 = 0.125. One MeshBlock, serial.
- Pure scattering, rho kappa_s = 1 (tau_cell = 0.125). Imposed flux F = 1 at x1min
  (`implicit_bc_x1min = flux`), Marshak q = 0.5 at x1max, periodic x2.
- c = 100, `implicit_cfl = 300` (c dt/dx1 = 300). Closure m1, bicgstab,
  offdiag auto (= operator), `implicit_closure_lag = step`.
- Seed: E *= 1 + 1e-5 cos(2 pi y/L2). The IC's E is the Eddington guess for F = 100, so
  E relaxes down by 100x in the first ~10 steps. f runs from 0.07 at the bottom to 0.5
  at the top.
- Measure (`spread.py`): s = max over x1 of (max_x2 E - min_x2 E)/mean_x2 E. g = the
  per-step factor fitted between s = 1e-5 and 1e-2. The bin dumps are float32, so
  "0" means below 1e-7.

**Growth table** (CPU, `s0`, 60-100 steps; "decays" = s at the float32 floor by step 15):

| arm | g/step | note |
|---|---|---|
| base, tau_cell 0.125 | 1.18 (1.46 with a 1e-6 seed) | O(1) by step 40; 21-61 positivity fallbacks |
| tau_cell 0.031 / 0.0625 | 1.80 / 1.84 | O(1) by step 10 / 15 |
| tau_cell 0.25 / 0.5 / 1.0 | decays | |
| seed k = 2 / 4 / 8 (Nyquist) | 1.22 / 1.42 / 1.41 | 4 transverse cells: 1.42 |
| c dt/dx1 = 1 / 3 | decays | |
| c dt/dx1 = 30 / 300 / 3000 / 30000 | 1.11 / 1.18 / 1.22 / 1.22 | 3000 and 30000 identical: dt-independent |
| closure kershaw / eddington | 1.29 / decays | |
| offdiag none / lagged | 1.11 / 1.18 | lagged: 4 NON-CONVERGED steps |
| dbg_trans_memory = 0 | 1.20 | |
| 1-D replica (nx2 = 1), q = 0.5 and 0.9, tau_cell 0.125 and 0.031 | steady | max change between dumps 0 after step 80 |
| 2-D, Marshak q = 0.9, tau_cell 0.25 | unstable | the f > 0.62 channel below (tau-independent) |

The S2 reproducer (`/viper/ptmp2/jinma/s2_0924/cpu/zc_base`: sph_atm on Cartesian,
rho kappa = 25/r^2, surface tau_cell 0.125) is the same case: 2e-6 -> 6e-3 in 10 steps.

## 2. Linear analysis (`model/fourier.py`, `scan1.py`, `scan2.py`)

Von Neumann analysis of the discrete step about a uniform state (E0 = 1, F = f0 c x^,
frozen coefficients, rho kappa_t = K, v = 0).
- Unknowns: cell E, x1-face F1, x2-face F2, exactly as the code has them:
  - face elimination with theta = 1/(1 + c dt K) and the F^n memory;
  - P_dd = D_dd E' in the matrix;
  - off-diagonal d_e P_de as the face mean of centred cell differences (M1OffDiv);
  - cell F = mean of the two faces (m1_impl_f1).
- The closure perturbations, linearised:
  - dP11 = chi dE + chi' df;
  - dP22 = (1-chi)/2 dE - chi'/2 df;
  - dP12 = (b/f0) dF2 (from n2 = F2/|F|), with b = (3 chi - 1)/2;
  - df = dF1/(cE0) - f0 dE/E0.
- These are taken from the OLD state (`implicit_closure_lag = step`, equivalently one
  Picard pass under `pass`) or from the NEW state (the converged Picard fixed point =
  the nonlinear backward-Euler step).

**Mechanism.** The lagged tensor D(F/E) makes the F-dependent part of the pressure
divergence explicit:
- the part that depends on F1 is dP11/dF1 = chi'(f);
- the part that depends on F2 is dP12/dF2 = b/f.

In the momentum equation these terms ADVECT a flux perturbation along n at a speed of
order c. The scheme treats that advection explicitly, at Courant number c dt/dx, and the
only thing that damps it implicitly is c dt rho kappa. For c dt >> dx the ratio is
dt-independent, and per step

    |g| ~ max(chi', b/f) sin(kx dx) / tau_cell      (tau_cell = rho kappa dx along n).

With the E part added there is a tau-independent term f chi'/min(chi, 1-chi). It exceeds
1 only for f > 0.62 (the transverse P22) and f > 0.69 (P11).

**Why 1-D is stable.** In 1-D, div F = 0 at c dt >> dx pins dF1 to the boundary flux, so
the F loop has nothing to act on. What is left is the E part of the lag,
g = f chi'/chi = 0.25 / 0.60 at f = 0.3 / 0.5. In multi-D the divergence-free flux
perturbations are unconstrained and grow:
- dF1(y), column-to-column, through chi';
- dF2(x), through the flux direction n.

**Why the other settings behave as they do:**
- Eddington and vet do not depend on F, so there is no loop.
- offdiag none removes only the b/f channel (the model gives the same max |g|; measured
  1.11).
- trans_limit lp adds an opacity proportional to |G_f|, which is second order about a
  state with F2 = 0. It bounds the saturation (2 instead of 4) but not the linear growth.

Model numbers, c dt/dx = 300, dy = dx (max over all modes):

| tau_cell | f0 = 0.3 | f0 = 0.5 | measured (f from 0.07 at the bottom to 0.5 at the top) |
|---|---|---|---|
| 0.031 | 9.7 | 18.1 | 1.80 |
| 0.0625 | 5.1 | 9.5 | 1.84 |
| 0.125 | 2.6 | 4.9 | 1.18 |
| 0.25 | 1.33 | 2.5 | decays |
| 0.5 | 0.67 | 1.25 | decays |
| 1.0 | 0.61 | 0.77 | decays |

- c dt/dx = 1 / 3 / 30 / 300 / 3000: 1.02 / 1.46 / 3.95 / 4.88 / 5.00 (f0 0.5, tau 0.125).
- Kershaw 6.0, Eddington 0.08, offdiag none 4.88, memory 0 4.88.
- Converged Picard (implicit closure): **0.195**.

The model is a local, frozen-coefficient bound at the top-cell f. The measured global
rates are lower because only the top few cells have f near 0.5. The model gets the
signs, the 1/tau scaling, the dt independence and the closure and offdiag dependence
right. It puts the threshold at tau 0.35-0.6, where the code has it between 0.125
and 0.25.

**PDE or scheme: the scheme.**
- The M1 flux Jacobian has real eigenvalues in every direction at f = 0.5, 0.8 and
  0.95. It is hyperbolic.
- The frozen-coefficient relaxation system (the PDE with rho kappa F) has max Re(lambda)
  < 0 at f <= 0.65 for all k.
- At f >= 0.7 it has a bounded, k-independent growth of 0.09 c K (f 0.7) to 4.3 c K
  (f 0.95) for oblique modes. That is a frozen-coefficient artefact about a
  non-equilibrium background, and backward Euler at c dt K >> 1 maps it to |g| -> 0.
- The same discretisation with the closure taken implicitly (converged Picard) has
  |g| <= 0.2 in the model and decays in the code (section 3, pass + anderson).

## 3. Cures tried (all on the reproducer, CPU)

| cure | tau_cell 0.125 | tau_cell 0.0625 | cost |
|---|---|---|---|
| (a) m1-thin `implicit_thin = hll`, default tau0 (binary `thin`) | 1.18 (unchanged) | - | the temporal tau (rho kappa c dt = 37) sets w = 0: HLL is never on at c dt >> dx |
| (a) `implicit_thin = hll`, tau0 = 100 | grows, 1.2 by step 90 | - | |
| (a) `implicit_thin = hll`, tau0 = 1e4 (full HLL) | O(1) by step 15 (worse) | tau 0.031: O(1) by step 10 | |
| (a) `trans_limit = lp`, fmax 1 / 0.5 | 1.18 / 1.18 | - | saturates lower |
| (a) `lp 0.5 + offdiag none` (the 09-24 decision) | 1.10 | - | **does not cure** |
| (b) `closure_lag = pass` | O(1) by step 5 | - | Picard NON-CONVERGED 100/100 |
| (b) pass + `closure_relax` 0.3 / 0.1 | decays | 0.3: O(1), NON-CONVERGED 100/100 | 25.8 / 31.7 passes (base 3.6) |
| (b) pass + `implicit_accel = anderson` (m 5 / 10) | decays | m 5: O(1), NON-CONVERGED 100/100 | 20.2 / 15.9 passes, 1 / 0 NC; inner its 11209 / 10938 vs 10380 |
| (c) HLL/upwind diffusion on the faces | = m1-thin above: diffusion of E does not reach the divergence-free flux modes | | |
| **(d) `implicit_closure_thin_relax = 1.5`** | **decays** | **decays** | **+1 kernel branch, 0 extra passes** |

**(d), the cure.**
- Under `closure_lag = step`, the (chi, n) a cell starts its step with is relaxed from
  the (chi, n) it used last step: c_used = (1-w) c_prev + w c(F^n, E^n), with
  w = min(1, 2/(1+G^2)).
- G = C max(chi', b/f)/tau_c + f chi'/min(chi, 1-chi) is the bound of the lagged gain.
  tau_c = rho kappa_t / sum_d |n_d|/dx_d is the cell optical depth along n (sp wedge:
  rho kappa_t dx1).
- The relaxed step map has eigenvalues 1 - w + w g. |1 - w + w g|^2 < 1 iff
  w |1-g|^2 < 2 (1 - Re g).
  - With w = 2/(1+G^2) this holds for every eigenvalue the model finds at f <= 0.65:
    max |eig| <= 1.0000 at tau 0.003-3.
  - It does not hold in general.
  - The fixed point is unchanged (c_used = c(F) there).
- Cells with G <= 1 (every optically thick cell) are left bitwise untouched.
- It does NOT cure f > 0.69 in thin cells: there Re g = f chi'/chi > 1, and relaxation
  cannot help.

Results over 500 steps (`s1`, `spread.py`, s at steps 25 -> 500; the reference
without the switch reaches O(1) with 461/500 positivity fallbacks):

| arm | s(25) -> s(500) | positivity fallbacks |
|---|---|---|
| tau_cell 0.00625 | 9.0e-6 -> 7.1e-6 (bounded, slow) | 0 |
| tau_cell 0.031 | 4.6e-6 -> 2.6e-6 | 0 |
| tau_cell 0.0625 | 1.2e-6 -> 9.3e-8 | 0 |
| tau_cell 0.125 (C = 1.5 and C = 1.0) | 3.4e-6 -> 0 by step 225 | 0 |
| tau_cell 0.25 / 0.5 | -> 0 | 0 |
| Nyquist seed / kershaw / offdiag none / c dt/dx 30 | -> 0 | 0 |
| seed 1e-2 | 3.3e-3 -> 0 by step 375 (0.97/step) | 0 |
| Marshak q 0.9, tau_cell 0.25 / 0.031 | 3.9e-6 -> 8e-8 / 6.4e-6 -> 3.4e-6 | 0 |

**Consistency** (`vs1d.py`): the relaxed runs converge to the 1-D M1 solution at the same
parameters, since the fixed point is unchanged. max |<E>_x2 - E_1D|/E_1D at steps
100 -> 500:
- tau_cell 0.125: 3.4e-3 -> 0;
- tau_cell 0.031: 2.8e-2 -> 1.2e-2;
- tau_cell 0.00625: 3.2e-2 -> 3.0e-2.

The price:
- in cells with tau_c << C max(chi', b/f), the closure adapts only over ~G^2/2 steps
  (~8000 at tau_cell 0.006);
- there it behaves as a slowly moving FROZEN tensor, and the state sits 3 % from the M1
  state (here the closure started from the Eddington-like IC);
- the step-restart loses the memory (not in the restart file), so the first step after
  a restart is unrelaxed.

## 4. Gates

| gate | result |
|---|---|
| CPU bitwise off, reproducer (`s0` vs `s1`, key absent) | 22/22 bin + log identical |
| CPU bitwise off, `tests_m1/gates/bitwise_off.sh` (box1, box2o, box4o, slab1, slab2o, slabv2, slabnd2; box_convection builds `s0bc`/`s1bc`) | all 7 rst bitwise, hst 0 |
| GPU bitwise off, apudev job 11958205 (HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1): reproducer + radwave xy | 17/17 and 33/33 files identical |
| GPU cure: reproducer tau 0.125 / 0.031, 500 steps | 6.5e-7 -> 0 / bounded 2.5e-6..4.7e-6, 0 fallbacks |
| tst rad_m1: slab_cpu, opcheck_mpicpu, restart_mpicpu | 3 passed (ATHENAK_M1_DATA set to symlinks) |
| thick limit, switch ON vs OFF (`scripts/gate_thick.sh`, CPU): radwave x2, radwave xy | 33/33 and 33/33 files identical in data (`cmpdata.py`). The byte cmp reports all files different only because the bin header carries the input file (1.5 vs 0.0). |
| Marshak (marshak_cart, closure m1, lag step), ON vs OFF | bin identical; the final tab differs by max dE = 5.3e-5 E_max at the front (x = 0.605, E there 3.7e-3 E_max). This is expected: the cold medium ahead of the front is thin, and the switch acts there. |
| 2-D diffusion pulse | NOT RUN: the arm failed on an out-of-range slice_x2 in the input. |
| implicit_op_check | not needed: the operator is untouched (only the lagged closure values change) |

GPU timing (job 11958205, 128x64, tau_cell 0.125, 200 steps, same binary, 2 interleaved
repeats):
- off: 31.4 / 30.0 s, unstable (162 positivity fallbacks), 4.5 passes;
- switch: 12.4 / 12.4 s, 3.0 passes;
- pass + Anderson: 51.4 / 51.4 s, 37.9 passes, 7 NON-CONVERGED.

S2 check (m1-sp2 b5f3afea + this switch, binary `sp2ctr`, `sp2/run_sp2.sh`, 300 steps):
- zc_base (Cartesian): spread O(1) by step 50 with 284 fallbacks; with the switch, 0 by
  step 50 and 0 fallbacks;
- the sp wedge at the same parameters: O(1) without, 0 by step 150 with.

The patch applied to m1-sp2 cleanly.

## 5. The 09-22 He slab

`he/run_he.sh`, box_convection `s1bc`, serial, tlim 50 s, vpert 1e-2:

| arm | positivity fallbacks | Picard mean | NON-CONVERGED |
|---|---|---|---|
| H0 closure_lag pass (the 3b7 setting) | 1848 | 200 | 1848 |
| H1 closure_lag step | 332 | 11.9 | 0 |
| H2 step + thin_relax 1.5 | 0 | 4.08 | 0 |
| H3 pass + anderson | crashed (rc=1) | - | - |

KE and the tau distribution of v^2 were NOT measured (no time). The fallback and Picard
counts point at the same thin-cell closure loop.

## HANDOVER

- Code: the switch lives in `rad_m1_implicit.cpp`, closure kernel (b), marked runs_5c.
  It needs multi-D implicit + `implicit_closure_lag = step`.
- To do:
  1. The He slab KE / F1top comparison of H1 vs H2 over 200 s.
  2. Store `ctr_mem` in the restart. The first step after a restart is unrelaxed now.
  3. The 2-D pulse gate: fix slice_x2 in the gate input.
  4. f > 0.69 in thin cells needs F implicit (Newton on P(F)). The switch does not
     cover it.
- The vet_col agent: the switch applies on the wedge unchanged (tau_c = rho kappa dx1
  there).
- Cleanup: build dirs `b_*` in `/viper/ptmp2/jinma/thinstab_0924` were deleted; the
  runs are kept.

## 6. After merging rt-integration 08a82d86 (merge fc2b0b9c)

- `bitwise_off.sh` with the switch off, new vs 08a82d86 (box_convection CPU): all 7 cases
  have rst bitwise, hst 0.
- tst rad_m1: 3/3 passed (ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data).
- Reproducer under `time_scheme = hesdirk2` (`inp/uni_ctr_h2`), 500 steps, spread s at
  steps 75 -> 475:

| tau_cell | switch off | switch 1.5 |
|---|---|---|
| 0.125 | O(1), 933 positivity fallbacks | 1.1e-7 -> 0, 1 fallback |
| 0.031 | O(1), 987 | 3.6e-6 -> 1.1e-6, 3 |
| 0.006 | O(1), 990 | 1.3e-5 -> 1.4e-5 (bounded), 2 |

The cure holds under hesdirk2. The relaxation acts on the first pass of each stage solve.
