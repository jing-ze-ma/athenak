# Magnetized Kelvin-Helmholtz growth rate: `hlld` vs `lhlld`

Branch `lhlld` worktree `/viper/u2/jinma/ATHENAK/bench/wt_lhlld`, HEAD `ae9f5336`.
Date 2026-09-13. Serial CPU (Kokkos SERIAL), double precision.

## Literature expectation (Minoshima & Miyoshi 2021, JCP 446, 110639, sec. 5)

The low-dissipation HLLD solver (LHLLD) should recover a KH growth rate that is
**insensitive to grid resolution and to the Mach number**, whereas plain HLLD's
numerical dissipation scales with the fast-magnetosonic speed and therefore
**suppresses the growth rate at low resolution and at low Mach number**.

## Setup

Problem generator: `src/pgen/kh.cpp`, `iprob=4` (Lecoanet et al. KH), MHD, ideal EOS,
gamma=5/3, `reconstruct=wenoz`, `integrator=rk3`, `cfl=0.4`, uniform `bx = b0`,
1 passive scalar. Domain x1 in [-0.5,0.5], x2 in [-1,1], periodic; single MeshBlock
covering the whole mesh. "64^2" = `nx1=64, nx2=128`; "128^2" = `nx1=128, nx2=256`.

**Mach number.** `iprob=4` **hard-codes `pres = 10.0`** (kh.cpp:168) and rho=1, so the
sound speed cannot be set from the input file. Per the instructions I therefore took the
second route: the low-Mach cases **divide `vshear` and `amp` by 10** (1.0 -> 0.1,
0.01 -> 0.001) and stretch `tlim` by ~10. Note this also lowers the shear Alfven Mach
number by 10 (`b0` was left at 0.01, i.e. v_A = 0.01); the field is weak in both cases
(v_A/v_shear = 0.01 and 0.1) so the instability is hydrodynamically dominated throughout.
Nominal shear Mach numbers: M = v_shear/c_s = 1.0/4.082 = 0.245 and 0.0245.

## Commands

Binary: `.../scratchpad/kh/build/src/athena` — see "Why a rebuild" below.
Base input `.../scratchpad/kh/kh.athinput` is `inputs/mhd/kh2d-lecoanet-mhd.athinput`
with `output1/dt = 0.05`, the VTK output removed, an `output2 file_type=rst` block added,
and `nlim=1e7`. Each run:

```
athena -i kh.athinput -d . \
   mhd/rsolver={hlld|lhlld} \
   mesh/nx1={64|128} mesh/nx2={128|256} \
   meshblock/nx1={64|128} meshblock/nx2={128|256} \
   problem/vshear={1.0|0.1} problem/amp={0.01|0.001} \
   time/tlim={4.0|32.0}
```

The two 128^2 low-Mach runs were long, so they were advanced in wall-clock-limited
chunks (`-t hh:mm:ss`) and continued with `athena -r <last>.rst -d . -t hh:mm:ss
time/tlim=22.0`; the history file is the concatenation of all segments.

## Metric

`2-KE` (column 9) from `KH.mhd.hst`. Reported quantity is the **energy** growth rate
`sigma_E = d ln(2-KE)/dt` (so the amplitude rate is sigma_E/2), from a least-squares
straight-line fit of ln(2-KE) vs t. Two windows are given:

* **auto** — between 2-KE = 10 x its initial value and 2-KE = 0.30 x its maximum,
  restricted to times before the maximum;
* **matched** — one window per (Mach, resolution) row, used for the ratio column so the
  two solvers are compared over identical times.

## Results

### Matched-window fits (the comparison)

| Mach | resolution | window in t | sigma_E (hlld) | R^2 | sigma_E (lhlld) | R^2 | lhlld/hlld |
|---|---|---|---|---|---|---|---|
| M = 0.245 (default) | 64^2  | 0.70 - 1.50 | 6.4135 | 0.99947 | 6.4158 | 0.99947 | **1.0004** |
| M = 0.245 (default) | 128^2 | 0.70 - 1.50 | 6.4725 | 0.99947 | 6.4726 | 0.99947 | **1.0000** |
| M = 0.0245 (10x lower) | 64^2  | 6.25 - 14.50 | 0.6701 | 0.99912 | 0.6723 | 0.99908 | **1.0034** |
| M = 0.0245 (10x lower) | 128^2 | 5.00 - 12.00 | 0.7043 | 0.99847 | 0.7041 | 0.99849 | **0.9996** |

Rescaling the low-Mach rows by the 10x slower shear (sigma_E x 10) gives 6.70 / 6.72
(64^2) and 7.04 / 7.04 (128^2), against 6.41 / 6.42 and 6.47 / 6.47 at the default Mach.

### Per-run auto-window fits

| run | t_end reached | 2-KE(0) | 2-KE max | auto window | sigma_E | R^2 | n pts |
|---|---|---|---|---|---|---|---|
| hlld  64^2  M=0.245  | 4.00  | 1.2533e-05 | 1.0623e-01 | 0.700 - 1.501 | 6.3728 | 0.99935 | 17 |
| lhlld 64^2  M=0.245  | 4.00  | 1.2533e-05 | 1.0595e-01 | 0.700 - 1.501 | 6.3749 | 0.99935 | 17 |
| hlld  128^2 M=0.245  | 4.00  | 1.2533e-05 | 1.1328e-01 | 0.700 - 1.550 | 6.3788 | 0.99909 | 18 |
| lhlld 128^2 M=0.245  | 4.00  | 1.2533e-05 | 1.1154e-01 | 0.700 - 1.500 | 6.4295 | 0.99933 | 17 |
| hlld  64^2  M=0.0245 | 32.00 | 1.2533e-07 | 1.1035e-03 | 6.251 - 14.650 | 0.6673 | 0.99909 | 168 |
| lhlld 64^2  M=0.0245 | 32.00 | 1.2533e-07 | 1.0423e-03 | 6.251 - 14.501 | 0.6711 | 0.99910 | 165 |
| hlld  128^2 M=0.0245 | 17.06 | 1.2533e-07 | 8.3746e-04 | 6.201 - 14.001 | 0.6816 | 0.99926 | 156 |
| lhlld 128^2 M=0.0245 | 12.34 (partial) | 1.2533e-07 | 8.9549e-05 | 6.201 - 10.501 | 0.7016 | 0.99804 | 86 |

All 8 matrix cells ran; none went non-finite and every run exited 0. The two 128^2
low-Mach runs were stopped inside/just after the linear phase on a time budget:
`hlld` reached t = 17.06 and `lhlld` t = 12.34; neither had saturated (the 64^2 low-Mach
runs saturate at t ~ 19.5), so the quoted "max" for those two is their last point.
Their auto windows are therefore not directly comparable, which is exactly why the
matched window 5.0 - 12.0 is used for that row.

## Fit quality

Every fit has R^2 >= 0.998 over 16-168 points, and the residuals are the usual gentle
curvature at the two ends of the linear phase. Sensitivity check on the one row whose
runs have unequal length (128^2, low Mach): window 5-10 gives 0.7167 (hlld) / 0.7153
(lhlld), ratio 0.998; window 5-12 gives 0.7043 / 0.7041, ratio 0.9996. So the window
choice moves sigma_E by ~2%, far more than the hlld-lhlld difference.

## Verdict

**Inconclusive with respect to the literature claim — because the two solvers agree.**

In this configuration `lhlld` and `hlld` give the same KH energy growth rate to within
0.04% at the default Mach number and within 0.35% at the 10x lower Mach number, at both
resolutions. There is no sign of the HLLD suppression that Minoshima & Miyoshi report:
between 64^2 and 128^2 the rate changes by +0.9% for hlld and +0.9% for lhlld, i.e.
**both** solvers are already resolution-insensitive here, and after removing the trivial
10x shear scaling the low-Mach rate is ~5-9% *higher* than the default-Mach rate for both
solvers rather than lower for hlld only.

This does not refute the paper. The likely reason the test cannot see the effect is that
it is not dissipation-limited:

* `reconstruct = wenoz` at 5th order makes the reconstruction, not the Riemann solver,
  the leading truncation error; MM21's sec. 5 comparisons are far more sensitive at low
  reconstruction order. A repeat at `reconstruct = plm` (or `dc`) is the obvious next step.
* The lowest Mach number reached is 0.0245 in *shear* Mach but the mode is resolved by
  64-128 cells per wavelength in x1; MM21's suppression appears when the mode is
  marginally resolved (their low-resolution cases) and at Mach numbers where the HLLD
  dissipation term ~ c_f |du| dominates the physical flux.
* The initial field is weak (plasma beta ~ 2e5), so the specifically *magnetic* part of
  the LHLLD modification barely engages.

So: the measurement **neither supports nor contradicts** MM21 here; it shows that under
`rk3 + wenoz` on this KH problem the choice of hlld vs lhlld is immaterial for the linear
growth rate. To test the claim, rerun this matrix with `plm` reconstruction and with a
coarser grid (e.g. 32^2 and 16^2).

## Why a rebuild was needed

`/viper/u2/jinma/ATHENAK/bench/wt_lhlld/build/src/athena` is configured with
`PROBLEM=built_in_pgens`, and `kh` is **not** one of the built-in problem generators
dispatched in `src/pgen/pgen.cpp` (the chain has advection, cpaw, gr_bondi, cshock,
linear_wave, implode, gr_monopole, mri3d, gresho_mhd, orszag_tang, rad_linear_wave,
rad_beam, shock_tube, shwave, wb_atm, z4c_*, spherical_collapse, diffusion). `kh.cpp`
lives in `src/pgen/` as a *user* problem, so it is only compiled with
`-D PROBLEM=kh`. That binary was built into the scratch directory
(`.../scratchpad/kh/build`) from the same unmodified worktree source; the existing
`build/` tree and everything under `src/` were left untouched.
