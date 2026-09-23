# runs_3x_hesdirk2: `<rad_m1>/time_scheme = hesdirk2` (H-ESDIRK2 for implicit M1)

- **Branch and base.** `m1-time2b` is branched from `m1-vimp` 2c3c4178; worktree `/viper/ptmp2/jinma/wt_time2b`.
- **Build and run tree.** `/viper/ptmp2/jinma/h2_3x`:
  - `new/` is the snapshot of the worktree; `base/` is the `git archive` of 2c3c4178.
  - `runs/` holds the radwave runs, `cpu/` the He slab, `gpu/` the GPU jobs.
- **References.** Design: `docs/dev/rad_m1_time2_design.md`. Model: `tests_m1/runs_3t_time2`, `tests_m1/runs_3u_time2impl`.
- **Spatial and coupling setting of every gate.** `implicit_enthalpy = plm` + `implicit_vimp = true`.

## STATUS: stopped at G1, stiff-limit boundedness

**Passed:**
- **G0.** Bitwise when off, on CPU and on GPU.
- **G1, Eddington time order.** p >= 1.9 in all 12 cases.
- **G1, listed stiff arms.** Every stiff arm of the vimp gate is bounded.

**Failed: an extended stiff probe.** At P_rad/P_gas = 100 and tau_lambda = 1e5 (tau_cell = 1.6e3), with `implicit_enthalpy = plm`, the scheme is linearly unstable:
- It blows up at N64 nt32, N64 nt48 and N32 nt24.
- `time_scheme = be` with the same spatial scheme is bounded.
- `hesdirk2` with `implicit_enthalpy = central` is bounded and matches the exact decay.

**Not run** (scripts are ready, see the end): vet_sc G1, the He slab gates, restart G8, the fallback runs, GPU cost G7 and the He-box CFL pair G5.

## What the code does

The new file is `src/rad_m1/rad_m1_time2.cpp`. The hooks are in these files:
- `rad_m1_implicit.cpp`: ImplicitSolve;
- `driver.cpp`: the stage loop;
- `restart.cpp` and `pgen.cpp`: the slope block;
- `rad_m1.hpp` and `rad_m1_implicit.hpp`: state and constants.

**The switch.** `<rad_m1>/time_scheme = be | hesdirk2`:
- The default is be. The parameter is read only when named, so the parameter dump and the restart files of a be run are unchanged.
- It requires `<time>/integrator = rk2` and implicit transport; `closure = tau` is refused.
- `time2_dbg_fail = <cycle>` is a DEBUG option: it declares stage 1 of that cycle not admissible, to exercise the fallback.

**One step.** Stage i is the Heun stage i of the hydro, unchanged, followed by the stage solve i. The stage solve is the backward-Euler ImplicitSolve, changed as follows:
- **Time step.** `dt_sub = g dt`, with `g = 1 - 1/sqrt 2`.
- **First iterate.** The first iterate and the only state that the EOS and the opacities see is the **stage start state**:
  - the radiation: U^n at stage 1, and the Heun average `(U^n + Y2)/2` at stage 2, for E and for all face fluxes F0;
  - the gas: whatever the hydro stage produced.
- **Old vector.** The old vector is `start + t2inc`, with:
  - stage 1: `t2inc = (1-g) dt K1`;
  - stage 2: `t2inc = dt (g/2 K1 + (b2 - g/2) K2)`.
- **Where t2inc is added.** It is added to `EN` and to `f0x*n` inside the solve. It is added to the hydro `u0` only after the start T has been computed.
- **Slope.** The write-back stores `K = (Y - old)/(g dt)`. This covers E, the gas momentum, the gas total energy and the x1/x2/x3 face F0. The stage-2 slope is K1, the FSAL slope.
- **vet_sc.** Stage 1 runs the formal solution at U^n: E^n, T^n and the opacities of the hydro `u1`. Both stages then use `D* = D^n + (dt/dt_prev)(D^n - D^{n-1})`.
  - A cell whose D* is not realizable keeps D^n. Such cells are counted as "vet clips".
  - After a backward-Euler step, the next stage step uses D^n.
- **Predictor.** The stage-2 solve has its own increment, `ipred2`. Increments are measured from the stage start.

**Admissibility.** After each stage solve, the stage is not admissible if any of these holds:
- the Picard loop did not converge;
- the od or vimp positivity fallback fired;
- min(E', T', e_gas') <= 0.

The Driver then restores U^n (the radiation `u1` and the face copies, the hydro `u1`, then an exchange and ConToPrim), reruns both Heun stages, and takes the backward-Euler step.

**Backward-Euler steps with slope storage** (`K1 = (Y - rhs)/dt`) are taken:
- at the first step;
- after a restart from a file without the slope;
- after a fallback;
- after a pack-size change.

**Restart.** Block `M1TIME2A` sits behind `M1PRED01`: an int32 flag, an int32 nch, `pred2_dt`, `dt_prev` and `vprev`. It is followed by nch cell slabs: K1 (8 channels), `ipred2` (3) and `vet_prev` (10). It is written only when a slope is stored.

**Model check of the in-solve enthalpy.** The data are in `tests_m1/runs_3u_time2impl`. H-ESDIRK2 with v' implicit (`imp`) has order 1.97-2.00 in all 12 cases and is bounded. This is what `implicit_vimp` implements.

## Gates run

### G0: switch off, bitwise vs 2c3c4178 (CPU)

- **Radwave.** Four cases (`lists/g0.txt`: Eddington default, Eddington plm+vimp at P=1/tau=10 and at P=100/tau=1e3, and vet_sc full plm+vimp). All 34 tab dumps of each case are identical with `cmp`. Location: `h2_3x/g0rw/{new,base}`.
- **He slab 2-D** (`slab2d_plm_vimp`, 200 s). All hst, bin, rst and `rt_profile.bin` files are identical. Location: `h2_3x/cpu/g0_{new,base}`.
- **GPU.** Job 11949915 on apudev, script `g0_gpu.sh`. It ran the 3-D He box (plm + vimp) for 40 cycles, base (md5 998ab6bb) vs new (md5 6efbbf8f). All 5+5 bin files, both hst files and `rt_profile.bin` are identical. Location: `h2_3x/gpu/g0_*`, log `h2_3x/gpu/log.out.11949915`.

### G1: radwave time order, Eddington (`RESULTS_g1_edd.txt`)

The setup is an x1 wave on a 64x4 mesh with bicgstab, plm + vimp, tol 1e-11, nt = 128 ... 2048. Each entry is the error at nt = 128 / the error at nt = 1024 / the median successive-difference order:

| P \ tau | 0.1 | 10 | 1e3 |
|---|---|---|---|
| 0.1 | 5.0e-4 / 7.8e-6 / 2.00 | 4.1e-4 / 6.7e-6 / 1.99 | 4.2e-4 / 6.5e-6 / 2.00 |
| 1 | 4.4e-4 / 7.1e-6 / 1.99 | 3.9e-4 / 6.2e-6 / 1.99 | 4.0e-4 / 6.9e-6 / 1.95 |
| 10 | 3.9e-4 / 6.5e-6 / 1.97 | 4.0e-4 / 6.4e-6 / 1.99 | 1.3e-3 / 1.8e-5 / 1.98 |
| 100 | 3.8e-4 / 6.1e-6 / 1.99 | 4.9e-4 / 7.7e-6 / 1.91 | 7.0e-4 / 9.4e-6 / 2.00 |

- **Individual orders.** The lowest single successive-difference order is 1.78, at (100, 10), 128/256/512. The median over each case is >= 1.91.
- **Error at nt = 1024.** It is at most 1.8e-5, against the design's criterion of <= 2e-5.
- **Model agreement.** The model gives e128 = 3.8e-4 to 5.0e-4. The measured values match it except at (10, 1e3), which is 1.3e-3.
- **Comparison with be.** The same runs with `time_scheme = be` give p = 0.91-1.11 and e128 = 3.5e-3 to 5.0e-2.
- **Solver behaviour.** No run was NON-CONVERGED and there were no stage fallbacks: each run has exactly one BE step, the first.

### G1: stiff limits (`RESULTS_stab.txt`, `RESULTS_stiff_tau1e5.txt`)

**The vimp gate's stiff arms, 10 periods: all bounded.** They cover P=100 with tau 1e3 and 10, nt 24/48, with and without a drift v0 = c_s. Findings:
- max |Z|/A0 = 1.000 in every arm.
- NON-CONVERGED 0 and 0 fallbacks.
- **Damping of the fundamental at P=100, tau=1e3.** After 10 periods the fundamental keeps 5.8e-2 (nt 24) and 1.8e-2 (nt 48) of its amplitude, against 3.6e-3 exact. vimp under be keeps 0.22 and 0.115 (runs_3v).

**Extended probe: P=100, tau_lambda=1e5, 40 periods.**

| arm | result |
|---|---|
| hesdirk2 + plm, N64, nt 24 / 64 / 96 | bounded. Grid noise is 7 / 4 / 3 % of A0, against 0.1-0.5 % under be |
| hesdirk2 + plm, N64, nt 32 / 48 | **unbounded**: \|Z\| up to 280x / 396x, noise up to 1e3 x A0 |
| hesdirk2 + plm, N32, nt 24 | **unbounded** |
| hesdirk2 + plm, N64, nt 48, amp 1e-8 | **unbounded**: exponential, about 3 %/step, so it is not amplitude-triggered |
| hesdirk2 + plm, tol 1e-8 | **unbounded** |
| be + plm (N64 nt48, N32 nt24) | bounded |
| hesdirk2 + central, N64 nt48 | bounded; end amplitude 0.801 vs exact 0.799 |
| be + central | bounded |
| hesdirk2 without vimp | NaN by cycle 100, as the model's `rhs` form predicts |
| upwind enthalpy, be or hesdirk2 | NaN: the known anti-diffusion of runs_3v |
| hesdirk2 + plm at P=100 tau=1e3 / 1e4, and at P=10 tau=1e5 | bounded |

- **P=1, tau=1e5, nt 24, 10 periods.** The Picard loop stalls at resid 1e-10 against tol 1e-11 in both schemes: be has 440 NON-CONVERGED solves. Under hesdirk2 this causes 265 stage fallbacks (266 of 612 steps backward Euler). The run stays bounded.
- **Reading.**
  - The instability needs H-ESDIRK2 and the plm (van Leer limited) enthalpy together.
  - It sits only at the stiffest corner, tau_cell ~ 1e3 with P_rad = 100 P_gas, in a window of dt.
  - It is not in the linear continuum model, which has neither the limiter nor the discrete plm face values.
- **Suspected mechanism (not tested).** Inside a stage, the plm `a_f` is reconstructed from the OLD-vector velocity `v_rhs`. That vector carries 2.41x the previous stage's force increment through the FSAL term. vimp then adds the implicit increment as a face MEAN, so the limited part sees a non-physical, grid-scale velocity.
  - Candidate fix: reconstruct `a_f` from the stage-start velocity and move `(v_rhs - v_start)` into the linear face-mean increment.
  - Alternative: use `implicit_enthalpy = central` under hesdirk2, which is bounded in this probe.

### Supplementary (not a gate): He slab 2-D, 200 s, hesdirk2 (`h2_3x/cpu/h_200`)

- **Solver.** NON-CONVERGED 0. There are 2481 solves: 1240 stage steps and 1 backward-Euler step.
- **Admissibility.** Stage fallbacks 0; od and vimp positivity fallbacks 0. The Picard count is 2.01 per solve.
- **Timing.** CPU time is not reported: cost is measured on the GPU only (G7, not run).

## Not run (stopped): scripts ready

| script | purpose |
|---|---|
| `h2_3x/cpu_gate.sh` | slab G0/NC, hesdirk2 slab 200 s / 1000 s, 2 ranks, vet_sc, `time2_dbg_fail`, restart with slope (G8), restart from a be file without slope |
| `lists/g1_vet.txt` | vet_sc G1 |

Input variants are in `h2_3x/inp`: `*_h2` (hesdirk2) and `*_be` (be named explicitly). A command-line override of a parameter that is not in the input file is fatal, so hesdirk2 has to be named in the input.

## Files

| file | content |
|---|---|
| `run_radwave.py` | runner from runs_3v, plus the keys `ts`, `dbgfail` and `amp` |
| `time_table.py`, `stab_table.py` | analysis scripts |
| `lists/` | case lists |
| `RESULTS_g1_edd.txt`, `RESULTS_stab.txt`, `RESULTS_stiff_tau1e5.txt` | results |

## Follow-up (09-23): fixing the P=100, tau_lambda=1e5 instability (`time2_enth_vel`)

`<rad_m1>/time2_enth_vel` sets how a stage solve under `implicit_vimp` builds the enthalpy
face coefficient a_f:

| value | a_f is built as |
|---|---|
| `old` | the plm a_f of the old-vector velocity (the behaviour of commit 58598ee6) |
| `start` | the plm a_f of the **stage-start** a (a - DA), with DA = the old-vector velocity increment (`t2inc/rho`) times (1 + D), added back as its face **mean** (`M1EnthCorrT`, `M1EnthEfT`) |
| `central` | a_f = the mean of the two cells; E_f stays plm (van Leer, upwind of a_f) |

- **Default is `central`.**
- **Implementation.** The donor-cell matrix part stays the implicit `a(v_old) E'`; only the lagged correction changes. `DA` is 3 iw components behind the vimp block, allocated only when `time_scheme = hesdirk2` is named.
- **be is untouched.** time_scheme = be is BITWISE unchanged with this binary: the radwave G0 list (4 cases x 34 dumps) and the 2-D He slab (hst, bin, rst: 12 files) are identical to 2c3c4178.
- **Binaries.** `h2_3x/athena_mpi_v4` (md5 02b2096f) runs `central` and `old` via `ev=`. `start` was measured with `athena_mpi_v3` (a274425), when `start` was the default.

### Stability (40 periods; `RESULTS_fix_central_af.txt`, `RESULTS_fix_start_af.txt`)

| arm | `start` | `central` | be |
|---|---|---|---|
| P=100 tau=1e5 N64 nt48, amp 1e-5 / 1e-8 | NaN / NaN | bounded: end 0.801 / 0.797 (exact 0.799), noise 7.8e-4 / 1.9e-2 | 0.112 / 0.068, noise 1.4e-3 / 0.10 |
| P=100 tau=1e5 N64 nt32, amp 1e-5 / 1e-8 | NaN / NaN | 0.800 / 0.792, noise 1.2e-3 / 3.5e-2 | - |
| P=100 tau=1e5 N32 nt24, amp 1e-5 / 1e-8 | NaN / NaN | 0.719 / 0.731, noise 2.2e-3 / 4.2e-2 | 0.192 / 0.077, noise 4.9e-3 / 0.18 |
| P=100 tau=1e5 N64, nt 16 / 24 / 64 / 96 / 128 | nt 24: 114x growth; nt 96: 444x; nt 128: 408x | all bounded (max \|Z\|/A0 <= 1.02), noise <= 1.3e-3 | - |
| P=10 tau=1e5 N64, nt 16 ... 128 | bounded, noise up to 4.8e-2 | bounded, noise <= 1.3e-3; end 0.39-0.54 (exact 0.533) | - |
| P=10 and P=100, tau=1e6, nt 16 ... 128, tol 1e-9 | NaN in every arm | bounded, NC 0, 0 fallbacks; end 0.93-0.94 (exact 0.939) at P=10 and 0.80-0.98 (exact 0.978) at P=100; noise 7e-4 to 7e-3 | nt32/96: 0.031/0.031 (P=10) and 0.131/0.053 (P=100) |

Notes on the table:
- **Run status.** Some `start` arms were killed early once the variant was known to fail. Each `start` arm's status is in its `log.txt` under `h2_3x/runs_v3`.
- **tau=1e6 needs `tol=1e-9`.** At the default tol 1e-11 the Picard loop stalls at resid 2e-10 to 1e-9 under both schemes: be has 1090 (P=100) and 1280 (P=10) non-converged solves at nt32. Under hesdirk2 those solves become stage fallbacks, 191-4610 per run. The runs stay bounded but are backward-Euler damped.
- **P=10 tau=1e5 nt48 at tol 1e-11:** 29 NC and 29 fallbacks, and the end amplitude is 0.39.
- **Amplitude-1e-8 noise.** At amplitude 1e-8 the relative grid noise of `central` (2-4 %) is below that of be (10-18 %).

### G1 Eddington order (`RESULTS_g1_edd_central_af.txt`, `RESULTS_g1_edd_start_af.txt`)

The metric is the same as before: nt 128 ... 2048, amplitude 1e-5, tol 1e-11, median successive-difference order.

| variant | lowest medians | other 10 cases |
|---|---|---|
| central | (100,10): 1.89; (1,1e3): 1.92 | 1.96-2.00 |
| start | (100,10): 1.90; (1,1e3): 1.91 | same as central |
| old | (100,10): 1.91 | - |

- **Errors.** e128 = 3.8e-4 to 1.3e-3 and e1024 <= 1.7e-5 in every variant.
- **Finer dt** (`RESULTS_g1_fine_dt.txt`). At (100,10) and (1,1e3) the successive-difference orders fall below 2 beyond nt 1024:
  - `central`: p(2048-4096-8192) = 1.49 / 1.75 and 1.80 / 1.65.
  - `old`: p = 1.21 / 1.29 at (100,10), and 1.86 / 1.77 at (1,1e3).
- **What does not change it.** The tail does not move with the Picard tolerance (1e-11 and 1e-12 give identical differences). It is not proportional to the amplitude: amp 1e-4 still gives 1.35 at 4096/8192, and amp 1e-7 is noisier (0.17 ... 1.53).
- **What does.** At amp 1e-4 the 128 ... 2048 medians rise to 1.93 at (100,10) and 1.90 at (1,1e3).
- **Cause.** Not identified. It predates the fix: `old` shows it too.

## STATUS after the follow-up

- **`start` fails.** It is unstable in every tau >= 1e5 arm at P=100, and in all tau=1e6 arms. Its G1 order is 1.90-2.00.
- **`central`** (a_f central inside stage solves only, E_f plm):
  - bounded in every probe;
  - grid noise at or below be;
  - G1 order >= 1.9 in 11 of 12 cases, with (100,10) at 1.89 by the gate metric (1.93 at amp 1e-4).
- **Stopped as instructed.** Per the follow-up brief (`start` failed, so try `central`, report both and stop), the remaining gates were not run: vet_sc G1, He slab, G8, fallback, G7, G5.
