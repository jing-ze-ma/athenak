# `tests_m1` — verification gates for the grey M1 module (`<rad_m1>`)

Analysis / reference-solution scripts for the test plan of
`docs/dev/rad_m1_design.md` section 8 (T1-T7).  Each script takes AthenaK
`file_type = bin` dumps (read through `vis/python/bin_convert.py`) or a history
file, prints one verdict line starting with `PASS` or `FAIL` plus the measured
numbers, and exits non-zero on `FAIL`.  Every script has a `--selftest` that
runs the same analysis path on synthetic analytic data, and
`--selftest --selftest-fail` which must FAIL; that is how the scripts are
verified before the code exists.

Everything is in **code units**.  `c`, `a_r`, `rho`, `kappa`, `v` are never
inferred from a dump: pass them on the command line and keep them equal to the
athinput of the run.  `--quiet` prints only the verdict line; `--all-ranks`
takes the rank0 file of a one-file-per-rank output; `--bin-convert-dir`
overrides the location of `bin_convert.py`.

The scripts are uniform-mesh only (they refuse a multi-level dump); the T1-T6
inputs are uniform by construction.

---

## `common.py`

Shared helpers, no command line of its own:

* `load_dump` / `load_series` — `bin_convert.read_binary` (or
  `read_all_ranks_binary`), MeshBlocks stitched into global `(nx3, nx2, nx1)`
  arrays plus cell-centre coordinates; a `.npz` with keys `time, x1, x2, x3`
  and one array per variable is accepted as well, which is what the selftests
  feed in, and so is a `.tab` file (see below).

  MeshBlocks are placed from `mb_geometry` (each block's own `x1min/x1max/…`),
  NOT from `mb_index`: that field is the block's index range **inside its own
  output slab**, and for a collapsed or sliced direction it comes back as `-2`.
  The milestone-1b version of this loader used it as a global index range, so
  for every 1-D and 2-D dump it wrote nothing at all and each gate silently saw
  an array of zeros.  1-, 2- and 3-D, many blocks and sliced outputs all work
  now; overlapping or non-tiling blocks raise instead of returning zeros.
* `load_tab` — an AthenaK `file_type = tab` output.  **The `bin` writer is
  single precision**; wherever a gate needs more than ~7 digits (the T3b thin
  side changes `E` by 3e-6 per cell, the T4b drift target is 1e-12) the run
  writes `file_type = tab` with `data_format = %26.17e` and the gate reads
  that.  Two traps in `formatted_table.cpp` are handled here: the header names
  a coordinate pair for every direction that is not *sliced* while the rows
  carry one only for directions with more than one cell, and with several
  MeshBlocks the rows come out in block order, not coordinate order.
  `file_type = tab` writes one file per `<output>` block, so the radiation and
  hydro variables of one snapshot are in two files; `t4_advect.py --hydro`,
  `t5_equil.py --hydro` and `t6_marshak.py --hydro` take the second one.
* `python3 common.py --selftest` — reads the REAL dumps in `data/`
  (`loader_1d.bin`, `loader_1d.tab`, `loader_2d.bin`, produced by
  `data/loader_1d.athinput`, 32 cells in 2 MeshBlocks and 16 x 8 in 4) and
  checks them against the analytic initial condition of that input.
  `--selftest-fail` must FAIL.
* `extract_1d` (1-D cut along x1/x2/x3, `--reduce mid|mean`), `sample_2d`
  (bilinear), `moments`, `linfit`, `l1_rel`, `nyquist_amplitude`,
  `nyquist_power_fraction`, `read_history` / `hist_column`,
  `base_parser` / `report` / `verdict`.

Variable names expected in the dumps: `m1_e`, `m1_f1`, `m1_f2`, `m1_f3`, and
`dens`, `velx`, `eint` or `press` where hydro is present.

---

## `t1_beam.py` — T1, beam

```
python3 t1_beam.py beam.out1.00010.bin --c <c> [--max-fwhm 25]
python3 t1_beam.py --selftest
```

Input: one 2-D dump of the 128^2, 45-degree beam run (`kappa = 0`,
`thick_flux = none`).  `--c` must equal `<rad_m1>` code-unit light speed;
`--angle`, `--x0`, `--y0` set the beam axis (default: 45 degrees from the lower
left corner of the domain).

The profile is cut **perpendicular to the beam axis** at the point where the
axis crosses mid-domain, bilinearly sampled at `--oversample` points per cell,
and the FWHM is taken from the two half-maximum crossings by linear
interpolation; it is reported in cells.  Also reported: `min(E)`, `max(E)` and
the largest reduced flux `f = |F|/(cE)` over cells above `--f-thresh`.

Pass (design note section 8, T1): FWHM <= 25 cells (24 expected with the
computed eigenvalues, 30 with `lam = +-c`; Gonzalez et al. 2007),
`min(E) >= -1e-12 max(E)` and `f <= 1 + 1e-10` (the M1 admissibility
conditions, section 4).

## `t1_pulse.py` — T1, free-streaming translation

```
python3 t1_pulse.py --c <c> \
    --run 128 pulse128.out1.00000.bin pulse128.out1.00010.bin \
    --run 256 pulse256.out1.00000.bin pulse256.out1.00010.bin \
    --run 512 pulse512.out1.00000.bin pulse512.out1.00010.bin
python3 t1_pulse.py --selftest
```

Input: for each resolution the **initial** and the **final** dump of the same
1-D `kappa = 0` periodic run; all runs must end at the same time.  The exact
solution is the initial profile translated by `c (t_f - t_0)` (periodic linear
interpolation), and the relative L1 error is computed per resolution, then the
order between consecutive pairs.

Pass: every measured order >= `--min-order` (default 1.8, i.e. the PLM second
order of section 3).

## `t3_pulse.py` — T3, static thick pulse (the AP gate)

```
python3 t3_pulse.py pulse.out1.000*.bin --c <c> --rho <rho> --kappa <kappa_s>
python3 t3_pulse.py pulse.out1.000*.bin ... --expect-excess 0.04   # scaled
python3 t3_pulse.py pulse.out1.000*.bin ... --nyquist
python3 t3_pulse.py --selftest [--nyquist]
```

Input: a time series (>= 3) of 1-D dumps of the scattering-dominated Gaussian,
no hydro; `--rho`, `--kappa` and `--c` must reproduce the run's
`D = c / (3 rho kappa)`.  The variance of `E` is fitted linearly in time
(`--skip` drops early dumps, `--subtract-min`/`--background` removes a floor)
and compared with `2 D`.

Pass: the ratio is within `--tol` (default 2 %) of `1 + --expect-excess`.  Use
the default target 1 for `thick_flux = ap_hll` at every `tau_cell`
(0.1, 10, 1e3, 1e6); for `thick_flux = scaled` pass
`--expect-excess 0.04` (the predicted `3<|mu|>/(2 * scaled_prefactor)` at the
default prefactor 20); `thick_flux = none` is the failing control at
`tau_cell >= 10` — if it passes, the test is not testing.

`--nyquist` seeds-mode variant: the amplitude of the odd-even mode
(projection on `(-1)^i`) is fitted exponentially and compared with the decay
rate of the compact 3-point diffusion operator at the Nyquist wavenumber,
`4 D / dx^2`.  Same `--tol` / `--expect-excess` logic (`scaled` is predicted to
reach only ~4 % of that rate).

## `t3b_jump.py` — T3b, opacity jump (REDEFINED in milestone 1c)

```
python3 t3b_jump.py jump.out2.00002.tab --c 1 --kappa 64 --flux 1e-6
python3 t3b_jump.py --selftest
```

The original gate ("the cell-centred `F` must be uniform to 1e-3") was wrong;
design note section 10 says why.  What is gated now, on one steady-state dump:

* `max |E/E_exact - 1|` against the analytic two-slope solution
  `dE/dx = -3 rho kappa F/c`;
* the two slopes, fitted away from the jump (`--fit-cells`), and their RATIO,
  which is independent of the overall flux the two Dirichlet ends settle on and
  so isolates the treatment of the jump itself;
* the FACE energy flux reconstructed from `E` with the scheme's own compact
  two-point diffusion flux and the arithmetic-mean face opacity: its
  uniformity away from the jump, its offset from the imposed flux, and, as a
  separate number, its peak deviation inside `--jump-cells` of the jump.

The cell-centred `m1_f1` error is printed as a **diagnostic only**.  Read the
run from the `tab` output: on the thin side of a 1e3 jump `E` changes by
~3e-6 per cell and single precision is not enough.

## `t4_advect.py` — T4, advected thick pulse

```
python3 t4_advect.py adv.out1.000*.bin --c <c> --v <v> --rho <rho> \
    --kappa <kappa> [--sigma0 <w>] [--drift-tol 1e-10] [--nyq-tol 1e-8]
python3 t4_advect.py --selftest --v 0.1
```

Input: a time series of 1-D periodic dumps with prescribed velocity `v`
(QUOKKA AP-paper parameters, design note T4), including `eint` or `press` so
that the gas energy drift can be formed; `--gamma` is only used in the
`press` case.  Moments are circular, so a pulse that has wrapped is handled.
`--x0` and `--sigma0` default to the measured centre/width of the first dump.

Reported and checked: the centre against `x0 + v t` (in cells), the width
against the static solution `sqrt(sigma0^2 + 2 D t)`, the relative drift of the
total gas energy, and the fraction of fluctuation power sitting in the Nyquist
mode of `E`.

Pass: centre within `--centre-tol` (1 cell), width within `--width-tol` (2 %),
gas drift below `--drift-tol`.  The odd-even diagnostic is report-only unless
`--nyq-tol` is given; the design note requires no odd-even mode at 512 cells
with `ap_hll`.

Added in 1c: `--hydro` (the parallel hydro tab series, matched by time),
`--temp-tol` (relative drift of the mean gas temperature — this is the T4b
gate) and `--static DUMP` (the final dump of the STATIC run at the same time;
the advected profile is compared with it shifted by `v t`, which is QUOKKA's
"advected vs static" number).  Note that the gas-energy drift of the PULSE is
physical — a radiating pulse exchanges energy with the gas — so the 1e-10
target of the design note belongs to T4b (uniform medium), not to T4.

## `t5_equil.py` — T5, single-zone equilibration

```
python3 t5_equil.py equil.out1.000*.bin --rho <rho> --c <c> [--chat <chat>] \
    --kappa-p <kP> [--kappa-e <kE>] --a-rad <a> --cv <c_v> [--eq-tol 1e-10]
python3 t5_equil.py --hist equil.hst --hist-erad m1_e --hist-egas eint ...
python3 t5_equil.py --selftest
```

Input: a single-zone time series, either as dumps (cell means are taken) or as
a history file (columns selected by name substring or 1-based index).  The
parameters must reproduce the run's source term.

The reference is built here, not read from a file: with `v = 0` the design
note's section 1 / 4(a) system reduces, using the invariant
`E_gas + (c/chat) E`, to one stiff ODE for `E`, integrated with SciPy `BDF`
(a bisection-based backward-Euler fallback if SciPy is missing), plus the exact
fixed point `kappa_E E = kappa_P a T^4` from a bracketed root find.  The EOS is
pluggable: `IdealEOS` supplies `e(rho, T)` and `T(rho, e)` from `--cv` or
`--gamma`/`--mu`, and a tabulated EOS only has to provide the same two methods.

Pass: final state within `--eq-tol` (1e-10 relative) of the exact equilibrium
(skip with `--no-equilibrium` for runs stopped short), invariant conserved to
`--cons-tol` (1e-13), and no sign change of `T_gas - T_rad` (differences below
`--sign-eps` of `T` count as converged).  `--ode-tol` additionally requires the
whole trajectory to follow the stiff reference.

## `t6_marshak.py` — T6, Marshak wave / Su-Olson

```
python3 t6_marshak.py marshak.out1.00010.bin --c <c> --sigma <rho*kappa> \
    --eps <4a/alpha> --x0 0.5 --t0 <t_off> --q0 <q> --tend <t_dump> \
    --lz 10 --nx 4000 --nmu 32 [--xmax 5] [--tol 0.02]
python3 t6_marshak.py --selftest [--ref-convergence]
```

**Provenance (what could not be sourced).**  The published Su & Olson (1996,
JQSRT 56, 337) benchmark table could not be accessed from this machine, and the
paper's semi-analytic quadrature could not be transcribed either.  No table
values are invented.  Instead the script computes its **own** reference for the
same linearised problem with a discrete-ordinates transport solver:

```
(1/c) dI/dt + mu dI/dx + sigma I = (sigma c V + Q)/2
dV/dt = eps sigma c (E - V),   E = (1/c) int I dmu,  V = a T^4 = eps U_material
```

(grey, constant opacity, `c_v = alpha T^3`, `eps = 4a/alpha`, unit source in
`0 < x < x0` switched off at `t0`, cold start, reflecting at `x = 0`).

**`--marshak` (milestone 1c).**  `c_v = alpha T^3` means `e ~ T^4`, which no
EOS in this tree can represent, so what the code runs is the CONSTANT-`c_v`
non-equilibrium Marshak wave and the same switch puts the script's reference on
the same problem: `U = rho c_v T` with `dU/dt = sigma c (E - a T^4)`
(nonlinear emission) and an incident isotropic bath `I(mu>0) = c E_b/2` at
`x = 0` instead of a reflection.  `--cv --rho --a-rad --t-init --e-bath` set
it, and the material profile of the dump is converted to `a T^4` with the same
`c_v`.  Use `--cfl 0.15`: at 0.4 the explicit reference sits on the stability
edge of its own emission term and overflows.  This is NOT Su-Olson, and the
published Su-Olson table does not apply to it.  `--hydro` takes the matching
hydro tab file.


Inputs from the run: `m1_e` and the material energy density
(`--material-var`, default `eint`, converted with `V = eps * U_material`;
`--material-is-v` if the dump already holds `a T^4`).  `--centre` is the
symmetry point of the run in x1, `--xmax` restricts the comparison window.

Pass (design note section 8, T6): relative L1 of `E` and of the material
energy <= `--tol` (2 %; Bloch et al.'s AP scheme reaches 1.1 %, uncorrected
HLL 84 %).  QUOKKA's nonlinear Marshak variant (4.5 %) is not covered by this
script.

## `t7_radshock.py` — T7, Lowrie & Edwards radiative shocks

```
python3 t7_radshock.py --m0 2                       # reference + IC numbers
python3 t7_radshock.py --m0 5 --athinput --write-ref ref_m5.txt
python3 t7_radshock.py shock.out1.00100.bin --m0 2 [--tol 0.02]
python3 t7_radshock.py --m0 2 --nondim shock.out1.00100.bin
python3 t7_radshock.py --m0 2 --selftest [--selftest-fail]
```

Semi-analytic reference for the **steady planar grey two-temperature radiative
shock** of Lowrie & Edwards (2008, Shock Waves 18, 129): Mach 2 subcritical and
Mach 5 supercritical, `L1 < 2 %` in `rho`, `T_gas`, `T_rad` (design note
section 8, T7).

**Provenance (what could not be sourced).**  The Lowrie & Edwards paper itself
could not be read from this machine.  The ODE system solved here is **derived
in the script's docstring** from the steady grey non-equilibrium-diffusion
equations as written out in Skinner & Ostriker 2013 (arXiv:1306.0010, Eqs.
94-97) and Ferguson, Morel & Lowrie 2017 (arXiv:1612.06346, Eqs. 9-11), which
restate the Lowrie-Edwards non-dimensionalisation.  The connection algorithm
(matching the two branches at equal `Theta` and `F_r`) is our own.  No
reference number is taken on faith; see the validation table below.

**Non-dimensionalisation.**  `rho`, `T`, velocity in the upstream `rho0`, `T0`,
`a0`; `p = rho T / gamma`; sound speed `sqrt(T)`; `P0 = a_r T0^4/(rho0 a0^2)`,
`C = c/a0`, `sigma` = total cross section per unit reference length,
`kappa = C/(3 sigma)` the radiative diffusivity, and `sigma_a` the
Lowrie-Edwards emission coupling (`sigma_a = sigma~_a L~ C`).  The standard
case `--sigma-a 1e6 --kappa 1` implies `C = sqrt(3 sigma_a) = 1732.05` and
`sigma~_t L~ = 577.35`, which is the value Ferguson et al. quote for the same
problem.  Only **constant** `sigma_a` / `kappa` are implemented; the power-law
(Bremsstrahlung) opacities of the paper are **not**.

**Method.**  (1) The two far-field equilibria follow from the three steady
integrals (mass; momentum including `P0 Theta/3`; energy including the
radiation enthalpy flux `4/3 P0 v Theta` and `F_r = -P0 kappa dTheta/dx`) by a
2-D root find.  (2) Both equilibria are saddles; the precursor is launched from
the upstream point along its unstable eigendirection and the relaxation region
from the downstream point along its stable one (`--eps`, default `1e-6`), and
both are integrated in **Mach space** with SciPy `LSODA` (in `x` the right-hand
sides are singular at `M = 1`; `dx/dM` and `dT/dM` are not).  (3) The embedded
hydrodynamic shock conserves `Theta` and `F_r` as well as the three integrals,
so it maps the supersonic branch onto the subsonic one at equal `(Theta, F_r)`:
it is located as the intersection of the two trajectories in that plane.
Subcritical, supercritical and the Zel'dovich spike all come out of the same
construction.

**Caveat.**  Lowrie & Edwards use Eddington closure, `P_r = E_r/3`.  These
shocks are optically thick, so an M1 code is in the diffusion regime and closes
at 1/3 as well — which is why everyone gates on this solution — but the
comparison is a diffusion-limit comparison and not a test of the M1 closure.
Skinner & Ostriker attribute their own M1 residual to exactly this, and report
1.7 / 6.1 / 7.8 % (`rho` / `T_gas` / `T_rad`) with computed eigenvalues against
0.42 / 0.49 / 0.42 % with `lambda = +-c/sqrt3`; record both, as the design note
asks.

**Validation.**  Cross-checks actually run, not quoted from memory:

| check | this script | reference |
| --- | --- | --- |
| RHS vs the published Lowrie-Edwards closed form (`cross_check`, run in every invocation) | max rel. diff `6.7e-16` | QUOKKA `extern/LowrieEdwards/radshock.py` |
| `M0 = 3` downstream `rho1/rho0` | `3.0021677` | `17.08233/5.69 = 3.0021670` (QUOKKA `testRadhydroShockCGS.cpp`) |
| `M0 = 3` downstream `T1/T0` | `3.6619127` | `7.98297e6/2.18e6 = 3.6619128` (same) |
| `M0 = 3` downstream `v1` (cgs, `a0 = 1.73e7`) | `1.728751e7` | `1.72875e7` (same) |
| `M0 = 3` domain length at `dT/T0 = 1e-4` | `1.57361e-2 cm` | `Lx = 0.01575 cm` (same) |
| `M0 = 3` shock position | `1.31828e-2 cm` | `0.0132 cm` (same) |
| `M0 = 3` character | subcritical (`T_p = 3.5536 < T1 = 3.6619`) | "sub-critical" (Skinner & Ostriker Sect. 4.3.5) |
| `M0 = 2` character | subcritical, Zel'dovich spike, `T_rad,s = 1.872 < T1 = 2.078` | Ferguson et al. Fig. 8 caption |
| `M0 = 5` character | supercritical (`T_p/T1 = 0.99938`), `rho1 = 3.598` | Ferguson et al. Fig. 5 and Fig. 11 (`rho` axis to 3.6) |

(The Ferguson figures use `P0 = 8.53e-5`, `sigma_t = 577.35`, so those rows are
qualitative.)

**Parameter sets.**  Non-dimensional: `gamma = 5/3`, `P0 = 1e-4`,
`sigma_a = 1e6`, `kappa = 1`, `M0 = 2` or `5`.  The cgs scaling is the Skinner
& Ostriker / QUOKKA one (`--rho0 5.69 --t0 2.18e6 --mu 1.67353e-24
--sigma-cgs 577`, giving `a0 = 1.7313e7 cm/s`, `L~ = 1.0003 cm`,
`P0 = 1.0019e-4`):

| | `M0 = 2` (subcritical) | `M0 = 5` (supercritical) |
| --- | --- | --- |
| `rho_l`, `v_l`, `T_l` | 5.69, 3.46264e7, 2.18e6 | 5.69, 8.65660e7, 2.18e6 |
| `rho_r`, `v_r`, `T_r` | 13.0078, 1.51467e7, 4.52910e6 | 20.4721, 2.40601e7, 1.86547e7 |
| nondim `rho1`, `T1` | 2.286075, 2.077570 | 3.597911, 8.557199 |
| domain, shock at | 1.52978e-2 cm, 1.00251e-2 | 3.35316e-2 cm, 3.31574e-2 |

`--athinput` prints exactly this as a pasteable `<problem>` block
(`m1_shock_rho_l/_v_l/_t_l/_rho_r/_v_r/_t_r`, plus `m1_shock_xs` and the mesh
line).  Start from the two equilibrium states separated at `m1_shock_xs`,
Dirichlet on both x1 faces, and run several shock-crossing times to steady
state.  At 512 cells the `M0 = 2` domain gives 335 cells of precursor and 177
of relaxation; the `M0 = 5` domain is precursor-dominated (506 / 6), so raise
`--ncells` (or `--dt-domain`) for that case — the script prints both cell
counts.

**Comparison mode** takes one dump plus the parameters, forms
`T_gas = eint/(dens c_v)` (or `press/(gamma-1)` if `eint` is absent) and
`T_rad = (m1_e/a_r)^{1/4}`, then slides the reference over the data
(bounded search, `--max-shift` cells, initialised at the steepest `dens`
gradient) minimising the summed L1 — a numerical shock drifts by a cell or two.
The applied shift is reported in cells.  Errors are `common.l1_rel`, i.e.
normalised by the L1 norm of the reference.  Use `--cv` / `--a-rad` if the run
is not in the cgs units above, or `--nondim` if it is in the Lowrie-Edwards
units (then `a_r -> P0`, `c_v -> 1/(gamma(gamma-1))`, `sigma -> 577.35`).

Pass: all three relative L1 below `--tol` (2 %), **and** the internal checks —
the three conservation integrals satisfied along the profile to `--cons-tol`
(1e-8; they come out at 1e-16), the two profile ends within `--end-tol` of the
equilibrium states, and the right-hand side agreeing with the published closed
form to `--cross-tol`.

`--selftest` resamples the reference onto `--selftest-nx` (512) cells and runs
the whole comparison path on it (L1 ~ 1e-8).  `--selftest --selftest-fail`
stretches the precursor by 20 % and puts a 5 % error on the Zel'dovich
spike/relaxation region, giving L1(`T_gas`) = 3.1 % at `M0 = 2` and 13 % at
`M0 = 5` — a FAIL.

---

## Selftest and lint

```
cd tests_m1
for s in t1_beam t1_pulse t3_pulse t3b_jump t4_advect t5_equil t6_marshak \
         t7_radshock; do
  python3 $s.py --selftest --quiet; done
python3 -m flake8 --max-line-length 90 .
```

(`t3_pulse.py` additionally has `--selftest --nyquist`, `t4_advect.py` needs
`--v`, `t7_radshock.py` takes `--m0 2` / `--m0 5`, and each script has the
`--selftest-fail` counterpart that must exit 1.)

## Not covered here

T2 (shadow), T3c (clipped extremum), T8 (conservation/restart/rank-invariance,
which is a bitwise diff, not an analysis) and T9 (radiation-supported
atmosphere) have no script in this directory yet.
