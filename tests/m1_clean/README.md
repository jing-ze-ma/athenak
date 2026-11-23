# Clean M1 radiation test suite

A self-contained CPU test suite for the clean M1 radiation scheme: one command, a PASS / FAIL
line per test, about 5 min (`--quick`) or 10-15 min (full) on a 16-core workstation.

```bash
# build a CPU binary (into build_m1test/) and run the quick subset
tests/m1_clean/run_tests.sh --quick
# or test a binary you built yourself
tests/m1_clean/run_tests.sh --exe build_cpu/src/athena            # full suite
tests/m1_clean/run_tests.sh --exe build_cpu/src/athena --only restart_sp,dt_retry
```

Requirements: python3 with numpy (nothing else), and for `--build` CMake, MPI (`mpicxx`) and
a C++17 compiler. The binary needs only the built-in problem generators, so the default
build (no `-D PROBLEM`) is enough. Every test runs on one MPI rank with OpenMP threads.

Options (`run_tests.py --help`): `--jobs N` tests at the same time (default: cores / 4),
`--threads T` OpenMP threads per run (default 4), `--workdir DIR` (default
`./m1_clean_runs`), `--keep` (keep the run directories of passed tests; failed tests are
always kept), `--record FILE` (write every measured number to a JSON file).

Exit status: 0 if all tests pass, 1 otherwise. A test whose run fails (non-zero return code)
is a FAIL with the path of its `run.log`.

## What the tests check

| test | deck | what | criterion |
|---|---|---|---|
| `conservation_sp` | `shell_sp` | closed rotating shell, 3-D spherical polar, 40 cycles (16 quick): total energy (gas + radiation), mass, total L_z | dE/E, dM/M, dLz <= 1e-13 |
| `conservation_cs` | `shell_cs` | the same on a cubed sphere, 12 cycles | as above |
| `restart_sp` | `shell_sp` | 16 cycles straight vs 8 + restart + 8 | every dump and history row identical (difference 0) |
| `restart_cs` | `shell_cs` | 8 straight vs 4 + restart + 4 | identical |
| `dt_retry` | `shell_sp` | a forced step retry at cycle 4 with factor 1 (the step is repeated at the same dt) and with factor 0.5 | factor 1: identical to the run without retry; factor 0.5: runs on; one retry each |
| `layout_tile` | `shell_sp` (16 x 16 lateral) | `<mesh>/layout_invariant_mode = tile`, `layout_tile = 8`, 8 cycles: 1 MeshBlock vs 1 x 2 MeshBlocks (along phi) and vs 2 x 2 (theta x phi, vacuum ray walls) | largest difference / largest value <= 1e-8 (see below) |
| `order_time` | `deckD` | time step halved at a fixed grid (3 levels quick, 5 full); successive-difference L1 of E, F_r, T | order of the finest pair >= 1.7; errors within 5 % of the reference |
| `order_space` | `wave_cart` | linear radiation-hydro wave, one period, N = 16 ... 128 cells with dt following the CFL | order of the finest pair >= 1.7 (full only) |
| `nlte_eps1e-1`, `nlte_eps1e-4` | `nlte_eps*` | two-level-atom non-LTE atmosphere (Jiang 2021, ApJS 253, 49, sec. 4.5): J vs the exact transport solution (`ref/nlte_exact_eps*.txt`) | L1 within 5 % of the reference value |
| `thermal_eq` | `thermal_eq` | uniform gas + radiation relaxing to equilibrium (Jiang 2021, sec. 4.1) vs the exact ODE | max relative error of E(t), T(t) <= bound |
| `wb_sp`, `wb_cs` | `rhse_sp`, `rhse_cs` | grey scattering envelope in hydrostatic and radiative equilibrium (spherical polar 2-D, cubed sphere) | max over the run of max |v| / c_s <= 1e-8 |

Notes on the results of this release (512cf768a):
- **L_z.** `dLz` is the drift of the scheme's total angular momentum, gas + radiation + the
  L_z ledger of the clean solve (history column `Lz_tot`); it is exact to round-off
  (2e-16). `dLz_gr`, gas + radiation alone, is printed for information: 7e-7 on
  `shell_sp` (the face-flux cone limiter `clean_stag_fcone` books the momentum it cuts) and
  6e-5 on `shell_cs` (the cubed-sphere torque that the gas does not take is booked).
  Without `clean_stag_fcone`, `shell_sp`'s `dLz_gr` is 2e-15.
- **Layout.** Tile mode is designed to give bit-identical results for every layout of
  whole-tile MeshBlocks. On this release the layouts differ after 8 cycles by 1.5e-10
  (1 x 2) and 4.5e-10 (2 x 2) of the largest value, so the test bounds the difference
  instead of requiring 0. The 1 x 2 difference is the round-off of the join of the
  per-MeshBlock F line solves; a development fix of that join makes 1 x 2 identical, while
  2 x 2 still differs by 7e-11 with it (open).
- **Order.** Deck D reaches second order only as dt -> 0: successive orders 1.57, 1.61,
  1.84 (E) over the five time steps. The criterion is on the finest pair.

The reference values and tolerances are in `ref/references.json` (measured with commit
512cf768a, CPU, GCC 13, `-ffp-contract=off`). Bitwise tests (restart, dt_retry, layout) need
no reference: they compare two runs of the same binary.

`references.json` format, per test: `"max": {name: bound}`, `"min": {name: bound}`,
`"equal": {name: value}`, `"near": {name: [reference, relative tolerance]}`; `"last": true`
applies list-valued bounds (orders) to the last element only.

## The decks

`decks/*.athinput` are complete input files. Each starts with a comment on the problem; its
`<rad_m1>` block begins with the scheme keys that production pins (the values of the
production star decks in `inputs/rad_m1_clean/`), then the test's physics and solver
set-up. A test changes only a few keys on top (for example `time/nlim`, `time/dt_fixed`,
`mesh/nx1`); `run_tests.py` writes the deck it actually ran as `in.athinput` in each run
directory. To run a deck by hand:

```bash
mkdir run && cd run && OMP_NUM_THREADS=4 ../build_cpu/src/athena -i ../tests/m1_clean/decks/shell_sp.athinput
```

## GPU builds

The same suite runs a GPU binary (`--exe build_cuda/src/athena --threads 1 --jobs 1`). The
bitwise tests compare two runs of one binary and must pass on a GPU too. The reference
tests (`order_*`, `nlte_*`, `thermal_eq`, `wb_*`) must agree with the CPU references within
their tolerances; a much larger difference points at the build (see the CUDA 12.2 warning in
`docs/m1_clean/BUILD.md`). The suite has been run on CPU only so far.

## Files

- `run_tests.sh`: builds a CPU binary if needed, then runs `run_tests.py`.
- `run_tests.py`: the tests and the PASS / FAIL logic.
- `m1t.py`: readers for AthenaK `.bin` dumps and `.hst` files, the deck writer and the run
  helper (Python standard library + numpy).
- `decks/`: the input decks. `ref/`: reference values and exact solutions.
