# `tests_m1` — verification gates for the grey M1 module (`<rad_m1>`)

Analysis / reference-solution scripts for the test plan of
`docs/dev/rad_m1_design.md` section 8 (T1-T6).  Each script takes AthenaK
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
  feed in.
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

## `t3b_jump.py` — T3b, opacity jump

```
python3 t3b_jump.py jump.out1.00050.bin [--x-jump 0.5] [--tol 1e-3]
python3 t3b_jump.py --selftest
```

Input: one dump of the steady constant-flux state across a 1e3 jump in
`rho kappa_F` (Bloch et al. 2021 Sect. 5.2).  The jump position is taken from
`--x-jump`, or found automatically as the largest step in `log dens` when the
dump carries `dens`.  `--trim` excludes cells at the two ends (boundary
layers), `--jump-cells` sets the half-width of the window around the jump.

Pass: `max |F/<F> - 1| < --tol` (default 1e-3) both over the whole profile and
inside the jump window (no peak at the jump — the arithmetic-mean face opacity
of section 3).

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
gas drift below `--drift-tol` (1e-10 — this is the test that the `beta^2` terms
of section 1 cancel; for T4b use 1e-12).  The odd-even diagnostic is
report-only unless `--nyq-tol` is given; the design note requires no odd-even
mode at 512 cells with `ap_hll`.

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
`0 < x < x0` switched off at `t0`, cold start, reflecting at `x = 0`).  The
solver is minmod-limited upwind per ordinate (2nd order in space), SSP-RK2 in
time, Gauss-Legendre in `mu`.  `--ref-convergence` doubles `nx` and `nmu` and
prints the self-error: **0.42 % in E at the default `nx=2000, nmu=16` (3 s) and
0.13 % at `nx=4000, nmu=32` (40 s)** for the selftest configuration — so run
the gate at 4000/32, where the reference is ~15x below the 2 % threshold.
If the published table is obtained later, feed it with
`--table FILE` (whitespace columns `x E V`, converted to code units) and the
solver is bypassed; nothing else changes.

Inputs from the run: `m1_e` and the material energy density
(`--material-var`, default `eint`, converted with `V = eps * U_material`;
`--material-is-v` if the dump already holds `a T^4`).  `--centre` is the
symmetry point of the run in x1, `--xmax` restricts the comparison window.

Pass (design note section 8, T6): relative L1 of `E` and of the material
energy <= `--tol` (2 %; Bloch et al.'s AP scheme reaches 1.1 %, uncorrected
HLL 84 %).  QUOKKA's nonlinear Marshak variant (4.5 %) is not covered by this
script.

---

## Selftest and lint

```
cd tests_m1
for s in t1_beam t1_pulse t3_pulse t3b_jump t4_advect t5_equil t6_marshak; do
  python3 $s.py --selftest --quiet; done
python3 -m flake8 --max-line-length 90 .
```

(`t3_pulse.py` additionally has `--selftest --nyquist`, `t4_advect.py` needs
`--v`, and each script has the `--selftest-fail` counterpart that must exit 1.)

## Not covered here

T2 (shadow), T3c (clipped extremum), T7 (Lowrie & Edwards radiative shocks),
T8 (conservation/restart/rank-invariance, which is a bitwise diff, not an
analysis) and T9 (radiation-supported atmosphere) have no script in this
directory yet.
