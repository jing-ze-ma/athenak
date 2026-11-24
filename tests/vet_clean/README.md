# Clean VET radiation test suite

A self-contained CPU test suite for the clean VET radiation scheme: one command, a PASS / FAIL
line per test, about 5 min (`--quick`) or 10-15 min (full) on a 16-core workstation.

```bash
# build a CPU binary (into build_m1test/) and run the quick subset
tests/vet_clean/run_tests.sh --quick
# or test a binary you built yourself
tests/vet_clean/run_tests.sh --exe build_cpu/src/athena            # full suite
tests/vet_clean/run_tests.sh --exe build_cpu/src/athena --only restart_sp,dt_retry
```

Requirements: python3 with numpy (nothing else), and for `--build` CMake, MPI (`mpicxx`) and
a C++17 compiler. The binary needs only the built-in problem generators, so the default
build (no `-D PROBLEM`) is enough. Every test runs on one MPI rank with OpenMP threads.

Options (`run_tests.py --help`): `--jobs N` tests at the same time (default: cores / 4),
`--threads T` OpenMP threads per run (default 4), `--workdir DIR` (default
`./vet_clean_runs`), `--keep` (keep the run directories of passed tests; failed tests are
always kept), `--record FILE` (write every measured number to a JSON file), `--keys
'block/key=value ...'` (extra input keys for every run, e.g. to test a new key on the whole
suite; a test's own settings win).

Exit status: 0 if all tests pass, 1 otherwise. A test whose run fails (non-zero return code)
is a FAIL with the path of its `run.log`.

## What the tests check

| test | deck | what | criterion |
|---|---|---|---|
| `conservation_sp` | `shell_sp` | closed rotating shell, 3-D spherical polar, 40 cycles (16 quick): total energy (gas + radiation), mass, total L_z | dE/E, dM/M, dLz <= 1e-13 |
| `conservation_cs` | `shell_cs` | the same on a cubed sphere, 12 cycles | as above |
| `conservation_lz` | `shell_sp`, `shell_cs` | the rotating shells with every L_z-exact cap / cut on (`clean_stag_fcone_lz`, `clean_lz_spin`, `clean_stag_fcone_emax`, `clean_wb_mcap`, `floor_from_ekin` with `vceil = 5e-3`; sp also `floorcap_lz` on every cell), 40 / 12 cycles; restart 8 + 8 and a forced dt_retry (factor 1) with the same keys | dE, dM, dLz <= 1e-13 on both, dLz_gr <= 1e-13 on sp (floorcap_lz included; cs: reported, as conservation_cs); restart and retry bitwise; no `CLEAN WARNING` |
| `restart_sp` | `shell_sp` | 16 cycles straight vs 8 + restart + 8 | every dump and history row identical (difference 0) |
| `restart_cs` | `shell_cs` | 8 straight vs 4 + restart + 4 | identical |
| `dt_retry` | `shell_sp` | a forced step retry at cycle 4 with factor 1 (the step is repeated at the same dt) and with factor 0.5 | factor 1: identical to the run without retry; factor 0.5: runs on; one retry each |
| `layout_tile` | `shell_sp` (16 x 16 lateral) | `<mesh>/layout_invariant_mode = tile`, `layout_tile = 8`, 8 cycles: 1 MeshBlock vs 1 x 2 MeshBlocks (along phi) and vs 2 x 2 (theta x phi, vacuum ray walls) | 1 x 2: identical (difference 0); 2 x 2: largest difference / largest value <= 1e-8 (see below) |
| `order_time` | `deckD` | time step halved at a fixed grid (3 levels quick, 5 full); successive-difference L1 of E, F_r, T | order of the finest pair >= 1.7; errors within 5 % of the reference |
| `order_space` | `wave_cart` | linear radiation-hydro wave, one period, N = 16 ... 128 cells with dt following the CFL | order of the finest pair >= 1.7 (full only) |
| `nlte_eps1e-1`, `nlte_eps1e-4` | `nlte_eps*` | two-level-atom non-LTE atmosphere (Jiang 2021, ApJS 253, 49, sec. 4.5): J vs the exact transport solution (`ref/nlte_exact_eps*.txt`) | L1 within 5 % of the reference value |
| `thermal_eq` | `thermal_eq` | uniform gas + radiation relaxing to equilibrium (Jiang 2021, sec. 4.1) vs the exact ODE | max relative error of E(t), T(t) <= bound |
| `wb_sp`, `wb_cs` | `rhse_sp`, `rhse_cs` | grey scattering envelope in hydrostatic and radiative equilibrium (spherical polar 2-D, cubed sphere) | max over the run of max |v| / c_s <= 1e-8 |
| `radwave` | `wave_cart` (64 x 4) | linear radiation-acoustic wave (amplitude 1e-6, the exact eigenmode, 2 periods, cfl 0.3) at P_rad/P_gas, x = (10, 0.01), (100, 0.1), (10, 1), (1, 100), x = tau_cell c_eq / c the grid-scale trapping parameter; phase speed and damping of a k = 1 fit vs the analytic root (`vett.radwave_root`); the x = 0.01 case also without `hydro/rad_signal_speed_diff` (full only, ~1 min) | phase speed within 3e-3, damping within 5e-3 per radian, no growth; dt gain of the key within 2 % of 2.77 |
| `wedge_sym` | `shell_sp` | the shell at rest, spherically symmetric, GL 7 x 18 rays, 8 cycles, on a 20 deg wedge and (with `vet_gd_wall_sym = true`) 3.7 and 7.3 deg wedges: fphi = \|Σ F_phi\| / Σ \|F_r\|; the 3.7 deg wedge without the key; the 20 deg shell made theta-periodic with `vet_gd_theta_seam_sym = true`: fth = \|Σ F_theta\| / Σ \|F_r\|, and without the key | fphi <= 1e-8 (fth reported); both runs without their key stop with the FATAL |
| `floor_pay` | `shell_sp` | `hydro/eos = general`, `tfloor` just above the gas temperature (the floor acts in most cells every cycle), `floor_pay = rad_column`, `floor_pay_skip`, `dfloor_flux_limit`, `floor_audit = 1`, 40 cycles: total-energy drift from t = 0; then 16 cycles with `floor_from_ekin = true` and `vceil = 1.2e-2`: energy, mass and L_z drift | dE <= 1e-13; floor energy paid > 0, unpaid 0; ekin dE, dM, dLz <= 1e-13; no `CLEAN WARNING` |
| `ceiling_metric` | `shell_sp`, `shell_cs` | the rotating shells (`gas_rot = 0.3`, `gas_vr = 0.2`, r 1-3), `hydro/eos = general`, `floor_from_ekin = true`, `vceil = 0.5` (below \|v\| in 1408 / 2928 cells) against the same run without the ceiling, cell by cell at t = 0 (cs: the velocity rebuilt from the contravariant v^a on the unit tangents, \|V\| the metric norm): v_phi kept (L_z), \|V\| = max(vceil, \|v_phi\|) where capped, else unchanged, v_theta 0, the radial part scaled in [0, 1]. A ceiling on Σ v_a² of mixed components or a split along the covariant phi component fails it on cs (dvphi 0.31). The cell-level C++ test `unit/ceiling_cell.cpp` (build line in its header) checks the same, the tfloor payment and the heat = metric KE change, at corner and ghost-extended points of all six panels, on sp, Cartesian and cs | dvphi, dcap, vtheta <= 1e-13; scale_ok |
| `momentum_box` | `wave_cart` (8^3, periodic 3-D) | total linear momentum, gas (`hydro_u`) + radiation (face momenta of `vet_face` / c^2), every cycle for 8 cycles; radiation wave along the xyz diagonal, amplitude 1e-3, `radwave_eig = false`, `nghost = 3`, tile 8: `dP_cyc` = largest component of \|P - P_0\| / Σ\|p_gas\| per cycle, `dP_max` its maximum (quick and full, ~10 s) | **a report, not a criterion**: passes when the run completes, the numbers are finite and there is no `CLEAN WARNING` |

Notes on the results of this release (465aff485):
- **L_z.** `dLz` is the drift of `Lz_tot` (gas + radiation + the L_z ledger of the clean
  solve). It is round-off (2e-16) by construction: the ledger books what the step does not
  conserve, so `dLz` checks the bookkeeping, not the conservation. The conservation error
  is the ledger's growth, `dLz_gr` (gas + radiation alone): 7e-7 on `shell_sp`, from the
  face-flux cone cut `clean_stag_fcone` (with `clean_stag_fcone_gas = true` and
  `clean_stag_fwb_ring = 1`, 1.6e-15 with no ledger; without the cut 2e-15). On the cubed
  sphere the ledger stays at round-off; there `dLz_gr` (6e-5) is not a conservation measure,
  because the history's `Lrz` is the cell vector u0 F while the scheme's radiation momentum
  is its face momenta (`CleanLzTotal`; ~0.6 % of the radiation's L_z apart).
- **Layout.** Tile mode gives bit-identical results for every layout of whole-tile
  MeshBlocks with one MeshBlock along x1 and x2 (theta): 1 x 2 (along phi) is identical.
  2 x 2 differs after 8 cycles by 1.6e-10 of the largest value: splitting theta puts ray
  walls (mirror or vacuum) inside the mesh, which is outside the tile guarantee, so the test
  bounds that difference instead of requiring 0.
- **Floors.** The floor energy of the initial conversion is paid before the t = 0 output, so
  `floor_pay` measures dE from t = 0 (`dE0`, the first step alone, is printed).
- **Narrow wedges.** On a phi-periodic wedge every ray set needs the width to be a multiple
  of its azimuthal period (GL n_theta x n_phi: 360/n_phi deg; LQ10: 90 deg; the default
  icosahedral set: 180 deg), else the run stops with a FATAL; `vet_gd_wall_sym = true` (with
  `vet_gd_wall_interp = true`) lifts it on narrow wedges. All suite decks meet the rule. A
  theta-periodic mesh needs `vet_gd_theta_seam_sym = true` (FATAL otherwise; `wedge_sym`).
- **Linear momentum.** The scheme does not conserve linear momentum (gas + radiation) in
  general, and keeps no ledger for it: the face-flux cone cut (`clean_stag_fcone` with
  `clean_stag_fcone_gas = false`), the well-balancing residual (`clean_stag_fwb`) and the
  floor velocity rules (`floorcap`, `dfloor_keep_velocity`, the velocity ceiling) change it
  and are not booked. Where none of them acts the momentum exchange is exact; at amplitude
  1e-3 in `momentum_box` none acts and dP is round-off (5.8e-17; per cycle 5.4-5.8e-17, measured at 69bfe5d04).
  In a non-linear box (amplitude 0.1) the cut changes P by up to ~4e-3 of Σ|p_gas| per
  cycle. `momentum_box` therefore reports dP (printed, and kept with `--record`) and does
  not bound it; on spherical polar and cubed-sphere meshes Cartesian P is not a discrete
  invariant at all (local bases), so the test is Cartesian only.
- **Order.** Deck D reaches second order only as dt -> 0: successive orders 1.57, 1.61,
  1.84 (E) over the five time steps. The criterion is on the finest pair.

The reference values and tolerances are in `ref/references.json` (measured with commit
512cf768a, CPU, GCC 13, `-ffp-contract=off`; reproduced to all printed digits on 465aff485 with
the IA12.2 scheme keys of the decks; `radwave` measured there). Bitwise tests (restart, dt_retry, layout) need
no reference: they compare two runs of the same binary. `momentum_box` has only
`"equal": {"warnings": 0}`: its dP is recorded, not bounded.

`references.json` format, per test: `"max": {name: bound}`, `"min": {name: bound}`,
`"equal": {name: value}`, `"near": {name: [reference, relative tolerance]}`; `"last": true`
applies list-valued bounds (orders) to the last element only.

## The decks

`decks/*.athinput` are complete input files. Each starts with a comment on the problem; its
`<rad_vet>` block begins with the scheme keys that production pins (the values of the
production star decks in `inputs/rad_vet_clean/`), then the test's physics and solver
set-up. A test changes only a few keys on top (for example `time/nlim`, `time/dt_fixed`,
`mesh/nx1`); `run_tests.py` writes the deck it actually ran as `in.athinput` in each run
directory. To run a deck by hand:

```bash
mkdir run && cd run && OMP_NUM_THREADS=4 ../build_cpu/src/athena -i ../tests/vet_clean/decks/shell_sp.athinput
```

## GPU builds

The same suite runs a GPU binary built WITHOUT `-D PROBLEM` (the built-in problem
generators; a `he_star_vet` binary cannot run these decks): `--exe build_tests/src/athena
--threads 1 --jobs 1` (on MI300A also `export HSA_XNACK=1`). The bitwise tests compare two
runs of one binary and must pass on a GPU too. The reference tests (`order_*`, `nlte_*`,
`thermal_eq`, `wb_*`) must agree with the CPU references within their tolerances; a much
larger difference points at the build (see the CUDA 12.2 warning in
`docs/vet_clean/BUILD.md`). The quick subset passed on MI300A (8/8) before `radwave`,
`wedge_sym` and `floor_pay` were added; those three have not been run on a GPU yet.

## Notebooks

`notebooks/` has one Jupyter notebook per test (conservation, thermal equilibrium, order,
well-balancing, non-LTE, momentum box, restart, radiation-acoustic waves, narrow wedges,
floors): each runs its test through this suite (or loads
an existing run) and plots the result against the exact solution or the criterion. Only numpy,
matplotlib and jupyter are needed; see `notebooks/README.md`.

```bash
cd tests/vet_clean/notebooks && ATHENA_EXE=$PWD/../../../build_m1test/src/athena \
  jupyter nbconvert --to notebook --execute --inplace thermal_eq.ipynb
```

## Files

- `run_tests.sh`: builds a CPU binary if needed, then runs `run_tests.py`.
- `run_tests.py`: the tests and the PASS / FAIL logic.
- `vett.py`: readers for AthenaK `.bin` dumps and `.hst` files, the deck writer and the run
  helper (Python standard library + numpy).
- `decks/`: the input decks. `ref/`: reference values and exact solutions.
