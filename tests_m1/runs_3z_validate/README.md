# runs_3z_validate: nonlinear validation of implicit M1 (plm + vimp + hesdirk2) vs be

## Verdict (2026-09-23)

**1. Optically thick radiative shocks (T7): PASS in every arm.**
- The arms cover M0 = 2, 5 and 7 (radiation-dominated), under both be and h2.
- L1 is 1e-4 to 2e-3 in rho, T_gas and T_rad, against a gate of 2e-2.
- NON-CONVERGED 0, stage fallbacks 0, one BE step per h2 run.
- **h2 and be are indistinguishable (within 1-20 %).** The error does not depend on dt
  (cfl 0.4 vs 0.2), so no time order can be measured: the time error is below the
  spatial/model error.
- **No clean space order either.**
  - M2: L1 is flat at about 4e-4 from N256 to N1024. It is the same for the explicit
    control, so it is a model/reference floor.
  - M5 and M7: rho converges at about 0.5-order.
  - Implicit transport is worse than explicit Eddington in two cases. M5 rho is 1.5-1.8e-3
    against 4.5e-4, and the shock drifts 0.9-1.6 cells. M7 T is 7.7e-4 against 7.7e-5.
  - Both of these are independent of N, dt and scheme. The cause is not identified.
    Candidates are the spatial operator (central flux, dc recon) and the Marshak ends.

**2. Free streaming (packet, shadow): FAIL for implicit transport, under both be and hesdirk2.**
- **Packet** (kappa = 0, f = 1 - 1e-6): E becomes non-finite within the first 0.125 in
  every implicit arm. This covers 64/128/256, implicit_cfl 0.5 and 2, and background E
  from 1e-15 up to 1e-2. Explicit M1 is fine in all of them.
- **Shadow:** a positivity fallback on essentially every solve.
  - h2 falls back to BE on every step (stage fallbacks = number of steps).
  - The depth ratio is erratic, 0.3 to 20, where explicit M1 gives 6.5e-2 at 280x80.
- **Arrival at 280x80 is fine.** Error -1.3 % (h2, k = 1), -2.0 % (be, k = 1), -2.0 %
  (explicit). It is -8 / -12 % at k = 8.
- **The limitation is in the implicit transport (M1 closure, f -> 1), not in the time
  scheme.**

**3. He box, vet_sc, from a common restart: both schemes are dt-independent at the level
  the chaotic flow allows.**
- KE_h(0.3)/KE_h(0.15): h2 0.999-1.061, be 0.953-1.071.
- T(z) max |dT/T|: h2 1.8-8.7e-5, be 4.9-8.9e-5.
- hesdirk2 is marginally closer up to t = 400 (KE_h within 2 % to t = 200, against 5 % for
  be). At t = 500 the two are equal within the noise.
- h2 vs be at equal cfl differ by as much as cfl 0.3 vs 0.15 does (KE_h up to 16 %).
- The earlier x4 cfl dependence is absent in all four arms, confirming it was the start-up
  kick.

**4. `implicit_bc = efix` breaks T7 (separate finding, see Constraints 4).**

## Setup

- **Code.** Branch `m1-validate` = d0c59f7c + 61c9b67a (`problem/beam_packet`, default off,
  read only when named). Bitwise check: `rad_m1_beam.athinput`, 3/3 bin dumps `cmp`-equal
  against d0c59f7c (`lists/bitwise_run.sh`).
- **Worktree.** `/viper/ptmp2/jinma/validate_0923/wt`.
- **Binaries.**
  - CPU `build_mpi` (d0c59f7c).
  - GPU `build_gpu_gen` (built-in pgens) and `build_gpu` (box_convection).
- **Runs.** `/viper/ptmp2/jinma/validate_0923/{t7,gpu2d,he3d}`.
- **Tables.** `RESULTS_all.txt` (`analyse_all.sh`) and `he_pairs.txt` (`he_pairs.py`).
- **Implicit block (`inp/impl_block.txt`).** transport implicit, bicgstab, offdiag auto,
  closure_lag step, tol 1e-8, lin_tol 1e-10, implicit_flux central, recon dc, enthalpy plm,
  vimp true, time_scheme hesdirk2 | be.

## Constraints found

1. **`implicit_vimp` fatals on a 1-D mesh.** The T7 shocks therefore run on nx2 = 4, x2
   periodic, dx2 = dx1.
2. **Implicit physical x2/x3 faces are always reflecting**, and the x1 boundary condition
   is one value per face. The T1 patch beam cannot be posed, so a periodic free-streaming
   packet (`beam_packet`) is used instead. The shadow test uses `ox2_bc = reflect` in all
   arms.
3. **vimp needs `gas_feedback = true`.** The shadow gas is re-imposed by the pgen each step.
4. **`implicit_bc = efix` (Dirichlet E in the end cell) breaks the T7 outflow end.**
   - Runs: `/viper/ptmp2/jinma/validate_0923/t7_efix_failed`.
   - What happens: after about 1-1.5e-10 s the x1max cell collapses to e_floor with
     F/cE -> 0.999, and the downstream state then drains out.
   - Where: M2 N1024 at cfl 0.4 and 0.2, and M5 N1024 at cfl 0.4 and 0.2. It happens under
     both be and hesdirk2.
   - Examples at t = 2e-10: M2 N1024 c0.4 L1 T_rad = 1.9e-2 (be) and 8.5e-2 (h2); M5 be
     N1024 c0.4 L1 T_gas = 2.8.
   - N256/N512 had not broken out by 2e-10.
   - **Replacement BC:** marshak at both ends, with bath = a T^4 of each equilibrium.
5. **Radiation-dominated reference.**
   - Lowrie-Edwards M0 = 27 / 50 (P0 = 1e-4) are not solvable with our tools.
     - `t7_radshock.py` hangs (LSODA underflow).
     - `t7x_radshock.py` (x-space, relative eigensteps): M27 has no (Theta, F_r)
       crossing, and M50 hits a step underflow.
   - Used instead: **M0 = 7, P0 = 1** (downstream P_rad/P_gas = 2.95).
     - A continuous shock, with the two branches joined at M = 1 (mismatch 6.1e-4).
     - Conservation residual 5.6e-16.
     - It is our own ODE reference, not a published profile.
     - arad x 1e4; domain 1.825 cm; tlim 4e-9 s.
     - `--max-shift 400`; shifts are reported relative to the t = 0 dump.

## 1+2. T7 (CPU; L1 rho / T_gas / T_rad at t_end = 2e-10 for M2/M5 and 4e-9 for M7)

| case | N | h2 cfl0.4 | h2 cfl0.2 | be cfl0.4 | be cfl0.2 | explicit Edd |
|---|---|---|---|---|---|---|
| M2 | 256 | 3.1e-4/3.7e-4/7.7e-5 | 3.1e-4/3.8e-4/7.7e-5 | 3.8e-4/4.7e-4/1.6e-4 | 3.3e-4/4.2e-4/1.2e-4 | |
| M2 | 512 | 4.5e-4/4.9e-4/7.7e-5 | 4.5e-4/4.9e-4/7.6e-5 | 4.8e-4/5.2e-4/1.2e-4 | 4.5e-4/5.0e-4/9.7e-5 | 4.4e-4/5.0e-4/7.1e-5 |
| M2 | 1024 | 4.5e-4/3.8e-4/9.2e-5 | 4.4e-4/3.9e-4/9.2e-5 | 4.5e-4/3.9e-4/1.1e-4 | 4.4e-4/3.9e-4/9.6e-5 | |
| M5 | 512 | 1.7e-3/3.9e-4/9.1e-5 | 1.7e-3/4.1e-4/8.8e-5 | 2.0e-3/5.1e-4/1.7e-4 | 1.9e-3/4.4e-4/1.2e-4 | |
| M5 | 1024 | 1.8e-3/4.4e-4/6.6e-5 | 1.8e-3/4.5e-4/6.6e-5 | 1.8e-3/4.5e-4/9.7e-5 | 1.8e-3/4.3e-4/7.0e-5 | 4.5e-4/7.6e-4/3.0e-4 |
| M5 | 2048 | 1.5e-3/7.3e-4/1.8e-4 | 1.5e-3/7.3e-4/1.7e-4 | 1.5e-3/7.3e-4/2.5e-4 | 1.5e-3/7.3e-4/2.1e-4 | |
| M7 | 512 | 1.3e-3/8.2e-4/8.1e-4 | 1.3e-3/8.2e-4/8.1e-4 | 1.4e-3/7.9e-4/8.2e-4 | 1.2e-3/8.0e-4/8.0e-4 | |
| M7 | 1024 | 9.4e-4/7.9e-4/7.8e-4 | 9.5e-4/7.8e-4/7.8e-4 | 8.9e-4/7.7e-4/7.8e-4 | 9.1e-4/7.7e-4/7.7e-4 | 1.1e-3/9.8e-5/9.7e-5 |
| M7 | 2048 | 4.6e-4/7.8e-4/7.7e-4 | 4.6e-4/7.8e-4/7.7e-4 | 5.0e-4/7.6e-4/7.5e-4 | 4.7e-4/7.7e-4/7.6e-4 | 7.5e-4/7.7e-5/7.6e-5 |

- **Zel'dovich spike.**
  - M2: +1.6 % (N256) -> +0.3 % (N1024).
  - M5: -1 % to +2.6 %. The explicit N1024 control is -8.5 %.
  - M7: within 0.4 %.
- **Precursor length.** Within 0.6 %.
- **Shock drift over the run.**
  - M5: -0.5 (N512), -0.9 (N1024) and -1.6 (N2048) cells, the same in all four implicit
    arms. Explicit: -0.07.
  - M7: up to 0.5 cells.
- **dt.** Hydro-set: at M2 N512, dt = 2.28e-13 (cfl 0.4) and 1.14e-13 (cfl 0.2), i.e.
  chat*dt/dx of about 240 and 120.
- **Solver.** Picard mean 1.7-2.7; NC 0; h2 stage fallbacks 0.

## 3. Free streaming (GPU apudev)

**Packet** (periodic box, 45 degrees, Gaussian half-width 0.05, t = 0.5). Explicit M1 is
exact to within:

| N | sig_perp / sig_0 | sig_par / sig_0 | peak | travel/(ct) |
|---|---|---|---|---|
| 64 | 1.208 | 1.210 | 0.567 | 0.999 |
| 128 | 1.032 | 1.032 | 0.801 | 0.999 |
| 256 | 1.004 | 1.004 | 0.924 | 1.000 |

- **Every implicit arm diverges (E = inf by t = 0.125).** This covers be and h2, 64/128/256,
  and implicit_cfl 0.5 and 2.
- **What the logs show.** In `pk_be_64_k2`, for example:
  - BiCGStab: 4377 breakdowns and 1459 line_jacobi fallbacks.
  - Positivity fallbacks on 14 of 16 solves; min E from the solve -3.8e242.
  - h2 fails every stage.
- **Background does not help.** Raising the background E to 1e-8 / 1e-4 / 1e-2 (packet
  contrast 1e8 down to 1e2) still diverges.

**Shadow** (depth E(0.8, 0)/E(0.8, 0.12) at 10 crossings; arrival at x = 0.98 from the
280x80 fine series):

| arm | depth 70x20 | depth 140x40 | depth 280x80 | arrival 280 | positivity fallbacks (280) |
|---|---|---|---|---|---|
| explicit M1 | 8.7e-3 | 5.6e-2 | 6.5e-2 | -2.0 % | - |
| be k=1 | 1.5 | 0.73 | 1.01 | -2.0 % | 6481 |
| be k=8 | 0.30 | 14.4 | 1.28 | -8.3 % | 820 |
| h2 k=1 | 0.72 | 0.45 | 0.48 | -1.3 % | 12993 (stage fb 6538) |
| h2 k=8 | 5.6 | 0.57 | 2.8 | -11.9 % | 1638 (stage fb 820) |

- The implicit depth time series swing between 0.05 and 20, so the lit/shadow contrast is
  not preserved.
- Arrival is not measured at 70 and 140 (no fine output series).
- The 560x160 arm was replaced by 70x20. A single 560-wide MeshBlock aborted with "HIP
  could not find a valid team size" (explicit arm; the implicit solve needs a whole x1 line
  per MeshBlock).

## 4. 3-D He box, vet_sc (GPU apudev, jobs 11953494-8)

- **Restart.** `kehdt_0924/R0/rst/m1slab.00001.rst`, t = 20. It was written with
  eddington + be + dc enthalpy at cfl 0.15.
- **All arms.** Switched to vet_sc + plm + vimp at t = 20 (`-r rst -i inp/he_vet_*.athinput`)
  and run to t = 500.
- **Solver.** NC 0 in all four arms. h2: 0 stage fallbacks, 0 vet clips, 1 BE step.
- **Cost.** The table gives `cpu time used` from the run logs. It covers only one run per
  arm; the arms were not interleaved or repeated.

| pair | KE_h ratio t = 25/100/200/300/400/500 | KE_v ratio (range) | max \|dT/T\| t = 45 ... 495 |
|---|---|---|---|
| h2 0.3 / 0.15 | 0.999 / 1.018 / 1.021 / 1.049 / 1.061 / 1.042 | 1.001-1.046 | 1.8e-5 ... 8.7e-5 |
| be 0.3 / 0.15 | 0.990 / 0.957 / 0.953 / 1.028 / 0.994 / 1.071 | 1.000-1.061 | 4.9e-5 ... 8.9e-5 |
| h2 / be at 0.15 | 1.010 / 1.067 / 1.086 / 1.038 / 1.094 / 1.013 | 0.990-1.014 | 2.9e-5 ... 5.1e-5 |
| h2 / be at 0.3 | 1.019 / 1.136 / 1.163 / 1.059 / 1.168 / 0.986 | 0.955-1.029 | 6.1e-5 ... 1.5e-4 |

| arm | `cpu time used` (s) |
|---|---|
| h2, cfl 0.15 | 408 |
| h2, cfl 0.3 | 232 |
| be, cfl 0.15 | 292 |
| be, cfl 0.3 | 161 |
