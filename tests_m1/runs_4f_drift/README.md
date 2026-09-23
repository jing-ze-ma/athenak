# runs_4f_drift: the implicit-M1 T7 shock drift (M5) and T error (M7) are a Marshak-end defect

## Verdict (2026-09-23)

**Cause.** The implicit x1 Marshak end face (`M1_IBC_MARSHAK`) carries only the comoving
flux `F0_f = c q (E - E_bath)`; the advected enthalpy flux `A E = v (1 + chi) E`, which
every interior face carries (upwinded), is left out at the end faces
(`src/rad_m1/rad_m1_implicit.cpp`, rows at HEAD 145e7ae8 lines 6191-6193 (x1max) and
6247-6249 (x1min); the stored face flux at 6483-6484). The lab-frame radiation energy the
gas carries in at the inflow end and out at the outflow end cannot cross the boundary, so
the solve moves the end cell's E off the bath until `c q (E - E_bath)` carries it:
`E_end/E_bath - 1 = -(A/(c q))`, i.e. about `-(8/3) v/c` at the inflow end and `+(8/3) v/c`
at the outflow end. Measured (`cons.py`), M7 N512: -1.27e-2 / +2.51e-3 (formula -1.08e-2 /
+1.92e-3); M5 N512: -3.4e-3 / +1.8e-3. The spurious non-zero boundary F0 also pushes on the
end cell through the momentum write-back. The result is a wrong incoming energy flux, a
precursor that is too short (M7: -0.62 %) and an N-independent T error in the
radiation-dominated M7, and a shock drift that grows with N in M5.

**Not the cause.** The interior coupling, the enthalpy/vimp terms and the time scheme are
not responsible. The same binary with `implicit_bc = efix` at both ends (Dirichlet E; needs
the efix write-back fix 1b1dd460, cherry-picked here) gives the same numbers as the fix, to
within 1-5 % in L1. Both reach the explicit-Eddington level.

**Fix.** Set `<rad_m1>/implicit_bc_advect = true` (default false, bitwise off). It adds the
enthalpy flux to a Marshak end face, upwinded with the end cell's velocity:
- outflow: `A E'` of the end cell, in the diagonal;
- inflow: `A E_bath`, on the right-hand side.

The stored face flux `f0` stays comoving (`A` is added separately for the cell flux, as in
the interior). Under vimp the end face's `A` uses the lagged `v^k`, the same as without
vimp. Recommended for every implicit run with a Marshak x1 end where the gas flows through
the end. For the He slabs the change is O(v/c) at the top and is not measured here.

## Results (CPU serial `build_cpu`, all cfl 0.4; L1 at t_end = 2e-10 (M5) / 4e-9 (M7); drift = shock shift(last) - shift(t=0) in cells; `table.py`)

| run | L1 rho | L1 T_gas | L1 T_rad | drift |
|---|---|---|---|---|
| M5 N512 be, marshak (before) | 1.98e-3 | 5.11e-4 | 1.69e-4 | -0.53 |
| M5 N512 be, marshak + bc_advect | 1.17e-3 | 7.50e-4 | 3.72e-4 | -0.20 |
| M5 N512 be, efix ends | 1.19e-3 | 7.24e-4 | 3.52e-4 | -0.20 |
| M5 N1024 be, marshak (before) | 1.81e-3 | 4.48e-4 | 9.65e-5 | -0.89 |
| M5 N1024 be, marshak + bc_advect | 5.31e-4 | 8.00e-4 | 4.07e-4 | -0.09 |
| M5 N1024 be, efix ends | 5.32e-4 | 8.14e-4 | 4.18e-4 | -0.09 |
| M5 N1024 h2, marshak + bc_advect | 4.49e-4 | 7.59e-4 | 3.27e-4 | -0.07 |
| M5 N2048 be, marshak (before, runs_3z) | 1.48e-3 | 7.26e-4 | 2.48e-4 | -1.64 |
| M5 N2048 be, marshak + bc_advect | 4.36e-4 | 4.22e-4 | 2.81e-4 | -0.46 |
| M5 N2048 be, efix ends | 4.39e-4 | 4.32e-4 | 2.86e-4 | -0.46 |
| M5 N1024 explicit Eddington (runs_3z) | 4.46e-4 | 7.55e-4 | 3.03e-4 | -0.07 |
| M7 N512 be, marshak (before) | 1.43e-3 | 7.85e-4 | 8.17e-4 | -0.10 |
| M7 N512 be, marshak + bc_advect | 1.68e-3 | 1.08e-4 | 2.11e-4 | -0.04 |
| M7 N512 h2, marshak + bc_advect | 1.56e-3 | 1.49e-4 | 1.32e-4 | -0.06 |
| M7 N512 be, efix ends | 1.69e-3 | 4.27e-5 | 1.89e-4 | -0.03 |
| M7 N1024 be, marshak (before) | 8.89e-4 | 7.70e-4 | 7.80e-4 | -0.17 |
| M7 N1024 be, marshak + bc_advect | 1.04e-3 | 8.28e-5 | 1.15e-4 | -0.07 |
| M7 N1024 h2, marshak + bc_advect | 1.07e-3 | 1.03e-4 | 9.73e-5 | -0.09 |
| M7 N1024 be, efix ends | 1.05e-3 | 4.36e-5 | 8.13e-5 | -0.06 |
| M7 N1024 explicit Eddington (runs_3z) | 1.09e-3 | 9.83e-5 | 9.69e-5 | -0.08 |

- **Residual M5 drift.** At N2048 there is a drift of -0.46 cells, which is 0.23 N1024 cells,
  against -1.64 before the fix. The efix ends give the same value, so it is not a boundary
  effect. Its time series, shift(t) at t = 0 / 0.5 / 1 / 1.5 / 2 (units of 1e-10 s), is
  0.14 / 0.11 / -0.01 / -0.19 / -0.32. There is no explicit M5 N2048 control to compare
  with. Not pursued.
- **M7 precursor length:** -0.62 % before, +0.00 % after (N512, `t7_table.py`).
- **M5 T_gas / T_rad:** L1 goes up a little with the fix (4.5e-4 -> 8.0e-4 at N1024).
  This is the explicit level (7.6e-4 / 3.0e-4), with the shock now in the right place. The
  lower "before" T error came from the misplaced shock offsetting the profile error.
- **End cells, `cons.py`.** E_rad(end cell)/E_bath - 1:

  | run | before (lo / hi) | after (lo / hi) |
  |---|---|---|
  | M5 N1024 | -3.4e-3 / +2.0e-3 | +1.2e-4 / -4.3e-5 |
  | M7 N1024 | -1.17e-2 / +2.2e-3 | -4.0e-4 / +2.5e-4 |

- **Domain momentum budget.** Net end-flux imbalance per `t * rho v^2` upstream:
  - M5 N1024: 8.7e-5 before, 1.8e-4 after (explicit 1.6e-4).
  - M7 N1024: -1.9e-3 before, -1.3e-3 after (explicit N1024 -1.1e-3).

  The "s/dx" columns of `cons.py` (shock shift implied by the mass and energy totals) do not
  isolate a leak: they disagree in the explicit control too (-0.05 / -0.23), because the
  precursor also relaxes. The write-back is conservative by construction (gas energy from
  the assembled row, face-shared momentum), and nothing here points at it.

## Gates

- **Bitwise, default off.** `m5_be_n512_c0.4` and `m7_be_n512_c0.4` (all tab dumps) are
  `cmp`-identical to `/viper/ptmp2/jinma/validate_0923/t7` (runs_3z binary d0c59f7c). The
  new code is only reached under `implicit_bc_advect = true`.
- **Solver.** NC 0 and no stage fallbacks in any fix arm (see the run.log tails).
- **Style.** The cpplint count is unchanged against HEAD.

## Files

- `run_t7.sh <name> <M5|M7> <N> <be|h2> [overrides]`: one run in
  `/viper/ptmp2/jinma/drift_0923/<name>`, with the runs_3z inputs plus the new parameter
  (`inp/`).
- `jobs1.txt` / `jobs2.txt` / `jobs3.txt`: the run lists.
- `table.py`: L1 and drift.
- `cons.py`: end-cell E and domain budgets.
- `go1.sh`: build and run.
