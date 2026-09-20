# tests_r22 — low-resolution relaxation route for the 3-D He4 envelope wedge

The production wedge (`tests_r11/r11_wedge23.sh`, argument list of job 11865437 =
`wedge11`) spends its first several turnovers getting through the initial transient: the
1-D `ic_profile` stratification is not an equilibrium of the 3-D code, the envelope
inflates, `L_out/L` drifts, and nothing of that is physics we want at 192². Standard
practice is to pay for the transient on a coarse grid and start the expensive run from
the relaxed state. This directory is that route.

Everything here is **preparation only** — no long run has been launched.

    r22_lowres.sh      the low-resolution arm (1 node, 2 GPUs, self-chaining)
    mk_ic_from_run.py  turns a run's shell-mean state into a new problem/ic_profile
    README.md          this file

---

## 1. The grid: 144 x 96 x 96

| | production | this arm |
|---|---|---|
| mesh | 144 x 192 x 192 | 144 x **96 x 96** |
| meshblock | 144 x 24 x 24 | 144 x 24 x 24 |
| blocks | 8 x 8 = 64 | 4 x 4 = 16 (ADI limit NADIB = 8: fine) |
| cells | 5.31e6 | 1.33e6 |
| nodes / GPUs | 4 / 8 | **1 / 2** |
| cells per rank | 6.64e5 | 6.64e5 (identical) |

Radial: unchanged, 144 cells in a **single** MeshBlock, as the implicit column solve and
the shell-mean dump both require.

**What 96² can and cannot resolve.** At R = 2.3717e11 cm, one angular cell is

    192 cells / 90 deg -> 0.469 deg = 1.94e9 cm
     96 cells / 90 deg -> 0.938 deg = 3.88e9 cm
     64 cells / 90 deg -> 1.41  deg = 5.82e9 cm

* The **entropy seed** (`vpert_sp_rand`, `vpert_sp_kmax = 8`) has a shortest wavelength of
  90/8 = 11.25 deg — 12 cells at 96, 8 cells at 64. Both grids carry the whole seed
  spectrum; the seed is not what sets the resolution.
* The **porous channels** measured at 3.5 turnovers have autocorrelation lengths
  5–14e9 cm = 0.5–1.2 deg, i.e. **3–9 cells at 192**, **1.3–3.6 cells at 96**,
  **0.9–2.4 cells at 64**.

So, honestly: at 96² only the *largest* channels exist at all, and they are 3–4 cells
across — marginally advected, heavily numerically diffused, with a convective flux that
will be biased low and a porosity (the low-density channels that let radiation leak) that
will be too smooth. Everything below ~8e9 cm is absent. At 64² *no* channel is resolved
(1–2 cells), so the mean state would relax under a flux the grid cannot generate — that
is a different star, not a cheaper one. 96² is the coarsest grid on which the relaxation
is still driven by (under-resolved) convection instead of by numerical diffusion, which
is why it is the choice here. Treat the relaxed profile as a **mean stratification and
energy balance**, never as a converged convective structure.

96² also keeps the *per-rank* load identical to production (8 blocks of 144x24x24 per
GPU), so nothing about kernel occupancy or the ADI/column workspace changes — only the
number of ranks and the number of blocks per angular ADI line (4 instead of 8).

## 2. Smoke measurement (job 11865564 / 11865620, apudev, 1 node, 2 GPUs)

100 cycles from cold start with all outputs off (11865564), then 200 cycles with the full
production output cadence on (11865620: `hst`/`log`/`rt_profile` every 47 s, `rt_surface`
every 94 s, `rst` every 940 s) — the same rate, so I/O is not a cost at this size.

| quantity | value |
|---|---|
| dt | 6.73 s at cycle 0, 8.35 s at cycle 100, **8.42 s** at cycle 200 (t = 1657 s), still rising |
| s / cycle | **0.238 s** (outputs off), **0.229 s** over 200 cycles (outputs on) |
| zone-cycles / s | 5.6–5.8e6 (2 GPUs) |
| host memory | MaxRSS 1.52 GB / rank, 3.19 GB / node (of 220 GB) |
| GPU memory | RT column workspace 466.6 + 177.0 MB, mode-3 exact block-tridiagonal |

Wall-clock estimates (1 turnover = 4705 s):

| mean dt | 6 turnovers (28230 s) | 10 turnovers (47050 s) |
|---|---|---|
| 8.2 s (as measured, no erosion) | 3.4e3 cyc, **0.23 h** | 5.7e3 cyc, **0.38 h** |
| 5.0 s (dt erodes as wedge10's did: 6.35 -> 4.63 by 3.1 turnovers) | 5.6e3 cyc, **0.37 h** | 9.4e3 cyc, **0.62 h** |
| 3.2 s (pessimistic, the wedge2 value) | 8.8e3 cyc, **0.58 h** | 1.5e4 cyc, **0.97 h** |

So the whole relaxation is **well under one node-hour** on a single node, against ~11
node-hours for the same 10 turnovers at 192² (wedge10 ran at 1.036 s/cycle on 4 nodes).
Note the arm is 4.4x faster per cycle than wedge10 at the *same per-rank cell count*:
part of that is wedge10's `mlt_split_sync = 4` passes (not used here, and not in wedge11
either), part is the angular ADI coupling over 4 blocks per line instead of 8. Caveat:
if dt collapses the way the earlier `d1`/`d2` arms did (hot voids above the swept shell
erode dt to ~0.07 s), these numbers are meaningless — watch dt in the first turnover.

## 3. The two-step recipe

### Step 1 — run the low-resolution arm

```bash
cd /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r22
sbatch --parsable --nodes=1 --time-min=01:00:00 -o lr1.log r22_lowres.sh lr1
```

Resubmitting exactly the same line continues from `tests_r22/lr1/rst/`'s newest restart
(same chaining as `r11_wedge23.sh`; `--export=ALL,OLDJOB=<id>` takes over from a running
job of the same arm). `tlim` in the script is 47050 s = 10 turnovers.

**When to stop it.** Both of:

* the face budget in the log (`grep "face budget" lr1.log`) has `L_rad,out/L` within a few
  per cent of 1 and no longer trending — this is the thermal relaxation of the envelope;
* the total energy in `lr1/he4.hydro.hst` is flat (no systematic drift) over at least one
  full turnover, and the shell-mean profiles stop moving: run
  `python3 mk_ic_from_run.py tests_r22/lr1 --tmin T-1 --tmax T` at two times one turnover
  apart and check the band table changes by less than the transient did.

Expect this at 5–10 turnovers; do not stop before the swept shell has settled.

### Step 2 — build the new IC and launch the high-resolution run

```bash
python3 mk_ic_from_run.py tests_r22/lr1 --tmin 5.0 --tmax 6.0 --smooth 3 \
        --out /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r22/ic_he4_relaxed.txt
```

`--tmin/--tmax` average over a window of `rt_profile.bin` records (in turnovers; default
is the last record alone), `--smooth N` is a radial boxcar of N run cells, `--blend N`
(default 4) is the smoothstep join width where the run's radial range meets the original
file. Rows outside the run's range are copied unchanged from
`tests_r14/ic_he4_tall_neww.txt`, and the output keeps that file's row grid, row density
and 3-column `r[cm] rho[g/cm^3] eint[erg/cm^3]` format, so it satisfies the coverage check
in `red_giant.cpp` (the file must span the grid *plus its ghosts*).

Then the production line of job 11865437 with `problem/ic_profile` repointed:

```bash
cd /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11
sbatch --parsable --nodes=1-4 --time-min=01:00:00 -o wedge12.log r11_wedge23.sh wedge12 \
  time/cfl_number=0.3 problem/rt_top_vacuum=true \
  mesh/nx2=192 mesh/nx3=192 meshblock/nx2=24 meshblock/nx3=24 \
  problem/vpert_var=eint problem/vpert=3.0e-2 problem/vpert_sp_rand=true \
  problem/vdamp_all_mean_only=true problem/vdamp_all_until=18820.0 \
  problem/vdamp_all_time=500.0 time/tlim=28230 \
  hydro/eos_rad_rho_hi=1.0e-10 hydro/eos_rad_rho_lo=1.0e-11 \
  hydro/eos_rad_t_hi=0.0 hydro/eos_rad_t_lo=0.0 \
  problem/rt_rad_force=false problem/rt_force_tau_gate=false hydro/rad_gate_rho=0.0 \
  problem/ic_profile=/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r22/ic_he4_relaxed.txt \
  mesh/x1max=2.964625e11 mesh/nx1=144 meshblock/nx1=144 \
  mesh/f_stretch_r_c1=0.992525 mesh/f_stretch_r_c2=-0.404821 \
  mesh/f_stretch_r_c3=-0.372255 mesh/f_stretch_r_c4=-1.499759 problem/mlt_alpha=0.0
```

(The radial grid, `x1max` and the stretch coefficients must stay *exactly* the production
values — the new file is only the stratification, not a new grid.)

## 4. Verification of `mk_ic_from_run.py`

**Is slot 6 the IC file's `eint`?** Yes. `python3 mk_ic_from_run.py --check
tests_r11/wedge10` round-trips the t = 0 record against
`tests_r14/ic_he4_tall_neww.txt` interpolated (linear in r, log in value — the code's own
interpolant) to the cell centres:

| band r/R | rho, median / max | eint, median / max |
|---|---|---|
| 0.50–0.90 | 5.3e-7 / 6.8e-6 | 7.4e-4 / 3.8e-3 |
| 0.90–0.97 | 2.9e-6 / 1.8e-5 | 2.7e-3 / 5.2e-3 |
| 0.97–1.01 | 2.2e-5 / 5.9e-5 | 4.5e-2 / 2.6e-1 |
| 1.01–1.30 | 3.5e-8 / 2.1e-3 | 2.1e-2 / 2.9e-1 |

So slot 6 *is* the internal energy density in erg/cm^3, the same quantity and units as the
file's third column: below 0.97 R the shell mean reproduces the file to 1e-6 in rho and
5e-3 in eint. Above 0.97 R eint disagrees by up to 29 % — that is not a unit error, it is
(a) the ±3 % random entropy seed, which does not average to zero once the floors act
(57600 `eos_efloor` events already at cycle 0), and (b) the shell mean of a quantity that
falls by decades per cell in the near-surface/vacuum transition, which is not the
log-interpolated value. **Everything the pipeline writes above ~0.97 R carries a
few-to-30 % uncertainty; below it, sub-per-cent.**

**Does it reproduce the IC at t = 0?** `--tmin 0 --tmax 0` on `wedge10` writes a file
that, compared row by row with the original, agrees to median 2.8e-4 / max 3.4e-3 in rho
and 7.2e-4 / 4.1e-3 in eint below 0.90 R, and ~1e-3 to 1.4e-1 above 1.01 R. Rows outside
the mesh are bit-identical. The residual inside the mesh is the unavoidable cost of
representing a 19489-row file on the run's 144 radial cells — the new IC is only as
detailed as the run grid, which matters only in the vacuum transition.

**How far has the state moved?** At 2.4–2.6 turnovers of `wedge10` (`--smooth 3`),
relative to the original IC:

| band r/R | mean d(rho)/rho | mean d(eint)/eint |
|---|---|---|
| 0.50–0.70 | -1.5e-3 | +8.6e-3 |
| 0.70–0.90 | -1.1e-2 | +1.3e-2 |
| 0.90–0.97 | -1.1e-2 | +6.0e-2 |
| 0.97–1.01 | +1.1e-1 | +9.5e-2 |
| 1.01–1.10 | **+34** | **+3.1e2** |
| 1.10–1.30 | +2.1e-1 | +1.8e-1 |

i.e. the deep envelope has barely moved (~1 %), the photospheric layer by ~10 %, and the
region just above the photosphere by **two orders of magnitude** — the swept/inflating
shell filling what the IC left as near-vacuum. That band is exactly what the relaxation is
for, and exactly where the numbers above are least trustworthy.

## 5. Caveats — read before trusting the relaxed IC

**(a) A shell-mean restart throws the 3-D state away.** Only `<rho>` and `<eint>` survive;
the velocity field, the channels and the density contrast are gone and the high-resolution
run re-seeds from the same `vpert` entropy noise. You pay the *convective* onset time
again (in the box runs, ~10 turnovers to saturation) — what you save is the *thermal /
hydrostatic* transient, which is the expensive part here.

**(b) Turbulent pressure is not in the 1-D balance.** Measured on `wedge10` at 2.5
turnovers, `rho <v_r^2> / p` per shell peaks at **1.3 %** (at 0.917 R), is ≤ 0.5 % below
0.86 R and ≤ 0.6 % above 0.95 R; the transverse term `rho <v_perp^2> / p` peaks at 1.7 %
near 0.96 R. So writing a purely thermal (rho, eint) profile and letting the code rebuild
its own hydrostatic column mis-states the support by ~1 % where convection is strongest.
That is small — but convection was still growing at 2.5 turnovers, and at saturation
(v ~ v_MLT) this could reach 5–10 %, so re-measure it on the relaxed arm before quoting it.

**(c) The mean state of a porous envelope is not a 1-D equilibrium for a homogeneous
shell.** This is the real limitation. In the 3-D run, radiation escapes preferentially
through the low-density channels: the *mean* (rho, eint) profile is in balance only
*given* that porosity. Re-imposed as a homogeneous shell, the same profile has a different
(smaller) effective escape rate, so the high-resolution run starts out of thermal balance
— too optically thick for its own luminosity — and will have its own secondary transient
while its porosity develops, in the opposite direction to the low-res run's. Worse, the
96² porosity is itself under-resolved (§1), so the relaxed mean is the equilibrium of a
*less* porous envelope than the 192² run will build. Expect the high-resolution run to
re-expand somewhat, and monitor `L_out/L` over its first few turnovers rather than
assuming it starts relaxed.

**(d) Carrying the 3-D state instead would be better — and the code cannot do it.**
Interpolating the low-res 3-D state onto the fine grid would keep the velocity field, the
channels and the porosity, and would avoid (a) and most of (c). There is **no supported
path for it in this code**:

* Restart does not re-grid. `BuildTreeFromRestart` (`src/mesh/build_tree.cpp`:470-490)
  reads `mesh_size`, `mesh_indcs` and `mb_indcs` out of the restart-file header and
  derives the root-grid block counts from them; the `<mesh>`/`<meshblock>` values in the
  input file are ignored on a restart, and the block list must match the file exactly
  (`"Incorrect number of MeshBlocks in restart file"`, line 518). There is no
  prolongating-restart option anywhere in `src/mesh`.
* `MeshRefinement`'s prolongation exists only for SMR/AMR *within* a running mesh. One
  could in principle start the fine run as an AMR/SMR job that refines the relaxed coarse
  mesh, but that produces a refined *AMR* mesh (and the ADI column/NADIB constraints and
  the shell-mean dump both assume every block spans the full radius), not the uniform
  144x192x192 wedge the production setup needs.
* `red_giant.cpp` has no 3-D initial-condition reader: `problem/ic_profile` is the only
  external-state input and it is strictly a 1-D three-column (r, rho, eint) table
  (`red_giant.cpp`:1593-1620, 2261-2340). No pgen in `src/pgen` reads a `bin`/`cbin` dump
  back in.

Adding a 3-D IC reader (read a `hydro_w` dump, trilinearly interpolate onto the new grid
in the problem generator) is a contained piece of work and would be the right thing to do
if this route is used more than once — but it is a source change, out of scope here.
