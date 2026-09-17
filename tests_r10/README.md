# tests_r10: the spherical two-stream's TOP FACE / THIN-LAYER behaviour under mode 3

2026-09-18, viper, branch `he4-presn-global` at `e73c6add`, binary `build_gpu_rg` (HIP,
`PROBLEM=red_giant`, unchanged; no source edits in this round). Measurement only. Read
`tests_r2/thick/README.md` (the diffusion-limit unit test and its 8-case matrix) and
`tests_r9/README.md` sections 0, 2.1, 3, 8 (the He4 star's own 0.956-0.981 emergent
shortfall and the coordinator's caveat pointing at the top face / thin layers) first.

## 0. What was measured, in short

* The 8 existing `tests_r2/thick` cases at `t = 0` (their own `mltfaces.txt`, **mode 0**,
  since `rt_col3_skip_sweep = false` in that input) plus 2 new optically-thin-top cases
  (`n8`, `n9`), all pushed through **mode 3** (`rt_col3_skip_sweep = true`) to 200 cycles
  and, separately, to each case's own **0.1 x top-cell thermal time**.
* The He4 star's own 1-D baseline (`tests_r9/base9` configuration, `he4_presn_cs.athinput`,
  mode 3 by default) run fresh for 200 cycles with a face dump at both ends.
* The plane-parallel box (`tests_r8/g1_new`, `box_w8`/`he_box_w8.athinput`, mode 3, already
  on disk from 09-17) re-read for `F_top/F_imposed`.
* Legality of `rt_top_re` / `rt_top_vacuum` as boundary-condition knobs on the sphere under
  mode 3, checked live.

**Headline finding, ahead of the tables**: the "0.97 R pile-up" defect (`tests_r9`) is
correctly attributed to the **last few cells / top face**, not the interior, but the
picture is time-scale dependent in a way `tests_r9` did not have the tool to see. At
**~0.1 of the top cell's own radiative thermal time** every thick-shell case, thin or
thick, sits at a modest **+20% to +35%** excess at the top face -- same sign and same
order as the existing mode-0 measurements and as the He4 star's own 1.28 -> 0.98 startup
transient. Only when a case is pushed **many** top-cell thermal times (which the fixed
200-cycle recipe silently does whenever the top cell is optically thin, because the
domain-wide hydrodynamic `dt` does not scale down with it) does a real divergence appear:
thin-top cases run away to 5x-77x the requested flux, while thick-top cases instead decay
smoothly to a 0.4x-0.7x deficit. The He4 star itself, run to the same 200 cycles, does
**not** run away -- it decays smoothly to 0.98, matching `tests_r9`. So the operator's
near-term defect (the thing that actually matters for a star, whose thermal time is
always >> 0.1 turnover) is the modest ~20-30% one; the catastrophic blow-up is a genuine,
separate long-time instability of mode 3 in a thin top layer that the He4 production grid
happens not to sit close enough to trigger over a few hundred cycles.

## 1. The harness

* `two_stream_sph_thick_relax.athinput` = `inputs/tests/two_stream_sph_thick.athinput` +
  `<output2> file_type = rst, dt = 1.0e30` (forced at `Driver::Finalize` regardless of
  `dt`, so a run with `time/nlim = N` and no other output cadence still leaves a restart
  at cycle `N`). `problem/mlt_dump` writes **once**, at the run's first two-stream call
  (`red_giant.cpp:3636`), so a single run can only ever report `t = 0`-ish; getting a
  face dump at a later cycle requires a two-stage run: stage 1 to `nlim = N` (no
  `mlt_dump`, keeps the restart), stage 2 restarted from that `.rst` with a **fresh**
  `mlt_dump` filename and `nlim = N+1`, so the dump fires on the first RT call after the
  restart -- i.e. on the relaxed state. Every run here used that pattern;
  `run_relax.sh` does it for all 10 thick cases and `he4_relax/run_he4_relax.sh` for the
  star. `.rst`/`.bin`/`.cbin*` were deleted immediately after each dump; only the small
  text dumps, `.hst` and the two run logs are kept (5.4 MB total in this directory).
* Cases 1-8 reuse `tests_r2/thick`'s own `ic_thick.txt` and `cases.txt` row (n, tau_tot,
  r_out/r_in, kappa, lstar, rad_flux_inner unchanged). Cases 9-10 (`n8_t100_r2.0`,
  `n9_t100_r2.0`) are new: same generator (`make_ic_thick.py --n 8|9 --tau 100 --ratio
  2.0`), chosen so the **last active cell's own optical depth** (`dtau_cell`, computed
  from the generator's own `kappa`/`rho(r)` -- not a sim output) is thin:

  | tag | n | tau_tot | dtau_top1 (last cell) |
  | --- | - | ------- | ---------------------- |
  | n0_t100_r2.0 | 0 | 100 | 0.781 |
  | n1_t100_r2.0 | 1 | 100 | 0.565 |
  | n3_t100_r2.0 | 3 | 100 | 0.262 |
  | n3_t300_r2.0 | 3 | 300 | 0.786 |
  | n3_t30_r2.0  | 3 | 30  | 0.079 |
  | n3_t100_r1.1 | 3 | 100 (r=1.1) | 0.677 |
  | n3_t100_r1.5 | 3 | 100 (r=1.5) | 0.418 |
  | n7_t100_r2.0 | 7 | 100 | 0.038 |
  | **n8_t100_r2.0** | 8 | 100 | **0.022** |
  | **n9_t100_r2.0** | 9 | 100 | **0.013** |

  `n8`/`n9` sit exactly in the requested "optically thin outer layer, `dtau_cell` at the
  top ~0.02-0.05" band (`n8` a little above, `n9` a little below); `n7` (already in the
  original 8) and `n3_t30` also land in the thin regime and are folded into the same
  analysis rather than treated as a separate "thin" set, since the generator's density
  power law puts every case's outer few cells at *some* `dtau_cell`, and it is that
  number, not a label, that turns out to control the outcome.
* `analyze_top.py <tag> <r_in> <lstar> <t0 file> <relaxed file>` prints the last 8 faces
  of `F_2s/F_req`, `w_blend` and `L_out/L = 4 pi r^2 F_2s / L` for a given pair of dumps;
  `topface_tables.txt` is its output for all 10 cases (t=0 vs. 200-cycle-relaxed).

## 2. t = 0 (existing mode-0 dump) vs. 200-cycle mode-3 relaxation, last 8 faces

`L_out/L` at the true top face (`r/r_in` in the last column of each block):

| case | dtau_top1 | t_final [s] | t_therm(top cell) [s] | t_final / 0.1 t_therm | `L_out/L` @ t=0 | `L_out/L` @ 200 cyc (mode 3) |
| --- | --- | --- | --- | --- | --- | --- |
| n0_t100_r2.0 | 0.781 | 1.02e5 | 1.22e6 | 0.84 | 0.672 | 0.199 |
| n1_t100_r2.0 | 0.565 | 1.04e5 | 7.11e5 | 1.46 | 0.671 | 0.344 |
| n3_t100_r2.0 | 0.262 | 1.07e5 | 1.07e6 | 1.00 | 0.850 | 0.627 |
| n3_t300_r2.0 | 0.786 | 1.06e5 | 1.53e5 | 6.95 | 0.675 | 0.650 |
| n3_t100_r1.1 | 0.677 | 2.29e4 | 8.23e4 | 2.79 | 0.680 | 0.018 |
| n3_t100_r1.5 | 0.418 | 6.82e4 | 2.55e5 | 2.68 | 0.702 | 0.352 |
| n3_t30_r2.0  | 0.079 | 1.08e5 | 4.63e5 | 2.32 | 1.131 | **4.689** |
| n7_t100_r2.0 | 0.038 | 1.13e5 | 3.30e4 | 34.3 | 1.205 | **76.58** |
| n8_t100_r2.0 | 0.022 | 1.15e5 | 3.62e4 | 31.8 | 1.355 (mode 3, own t=0) | **64.54** |
| n9_t100_r2.0 | 0.013 | 1.17e5 | 2.77e4 | 42.1 | 1.370 (mode 3, own t=0) | **51.15** |

(`n8`/`n9` have no pre-existing mode-0 dump; their "t=0" column is mode 3's own first-call
value instead, from the same run that produced the relaxed one -- still essentially
un-relaxed, at cycle 0-1.)

**The fixed-200-cycle recipe does not hold `t_final/t_therm` fixed across cases.** The
hydrodynamic timestep is set by the bulk of the domain (`dt` fell from ~1.4e3 s to
~3.5-4.0e2 s over 200 cycles in every case, regardless of `n`), so `t_final` lands within
a factor 2 of `1.1e5` s for every `r_out/r_in = 2` case -- but the **top cell's own
thermal time** varies by two orders of magnitude with `n` (`1.2e6` s at `n=0` down to
`2.8e4` s at `n=9`, because a steeper density power law puts less mass, and hence less
heat capacity, in the last cell). The intended "0.1 relaxation time" therefore lands
close to the mark only for `n<=3` at `tau=100` (`ratio` 0.84-1.46); every thin-top case is
run 7-42 top-cell thermal times instead of 0.1.

## 3. What the top cell actually looks like at ~0.1 of ITS OWN thermal time

Because case 2's ratio column shows the 200-cycle recipe overshoots badly for the thin
cases, `n7`, `n8` and `n9` were re-run with a properly scaled duration: stage 1 to
`nlim = 2` (~2.8e3 s, matching their own `0.1 t_therm` of 2.8-3.6e3 s), stage 2 restarted
with a fresh dump at `nlim = 3`. (Commands: same two-stage pattern as `run_relax.sh`,
`time/nlim=2` then `time/nlim=3` from the cycle-2 restart; not kept as separate scripts,
folded into this file since each is a 3-line variant.)

| case | t [s] | 0.1 t_therm [s] | `F_2s/F_req` top face | `L_out/L` top face |
| --- | --- | --- | --- | --- |
| n7_t100_r2.0 | 2821 | 3297 | 1.278 | 1.278 |
| n8_t100_r2.0 | 2814 | 3617 | 1.279 | 1.279 |
| n9_t100_r2.0 | 2809 | 2771 | 1.265 | 1.265 |

**At the correct time scale every thin-top case sits at +26% to +28% at the top face --
the same order as the existing t=0 mode-0 numbers (1.21-1.37) and as the He4 star's own
1.28 startup value (`tests_r9` section 3).** None of the three shows any sign of the
50-77x runaway yet at this point; the runaway is something that develops **later**, well
past the top cell's own thermal time, not an instantaneous mode-3 defect. It is real
(section 2 shows it lands the same whether reached by `n3_t30`'s milder 2.3x overshoot or
`n7-n9`'s 30-42x overshoot) but it is a **long-time secular growth**, distinct from the
near-term `1 - H_T/(2r)`-type error the earlier rounds characterised.

## 4. Direction of the long-time drift depends on `dtau_top1`, not on the overshoot ratio

Reading section 2's table by `dtau_top1` rather than by overshoot ratio: every case with
`dtau_top1 >~ 0.1` (`n0, n1, n3_t100_r2/r1.1/r1.5, n3_t300`) decays smoothly to a
**deficit**, `L_out/L = 0.018-0.65`, monotonically falling as the run proceeds (spot
checked on `n3_t100_r2.0`'s full radial profile at 200 cycles: faces 20-64, `r/r_in =
1.16-1.51`, are themselves grossly wrong there too -- 170, 114, 9.4 -- settling back to
0.88-0.99 by `r/r_in = 1.79-1.87` and only then declining to 0.63-0.72 at the top; i.e.
by 200 cycles the whole outer half of even a "thick" case has been disturbed, not just
the last 5 cells, because the star has genuinely evolved away from the exact-diffusion
IC under the wrong flux divergence -- `dt` dropping 4x over the run is the hydrodynamic
signature of that). Every case with `dtau_top1 < ~0.08` (`n3_t30, n7, n8, n9`) instead
**grows** to `L_out/L = 4.7-76.6`, and by 200 cycles the growth has already saturated to
a near-flat value across the last 8 faces (76.62 down to 76.58, `n7`; compare the smooth
monotonic profile at `t=0`), i.e. it is not actively diverging further at that instant,
it has run away to a different, still-wrong, plateau. `~0.1` in `dtau_top1` is a
reasonably sharp line between the two long-time fates in this matrix; it is not a
overshoot-ratio effect, since `n3_t30` (only 2.3x overshoot, `dtau_top1 = 0.079`) already
blows up while `n3_t100_r1.1` (2.8x overshoot, `dtau_top1 = 0.677`) instead decays.

## 5. The He4 star's own baseline, run fresh (not the runaway regime)

`tests_r9/base9`'s own `he4_presn_cs.athinput` configuration (mode 3 by default -- no
`rt_col3_skip_sweep` override needed, the input already carries `rt_implicit_column = 3`),
`1-D` (`nx2=nx3=4`), `inner_bc=wall`, `rt_bottom_flux=true`,
`rad_flux_inner=1.305278e15`, `mlt_alpha=1.5`, run fresh for 200 cycles (stage 1) then
restarted for a face dump (stage 2). `r_in = x1min = 1.18585e11` cm,
`L = 2.3066e38` erg/s. Last 8 faces (`r/r_in = 1.985-2.024`, i.e. the outermost ~2% of
the domain, `tau` there ~0.02-0.03 by the run's own column dump -- comparable
`dtau_cell` to `n8`/`n9`):

| `r/r_in` | `F_2s/F_req` @ t=0 | `F_2s/F_req` @ 200 cyc |
| --- | --- | --- |
| 1.985 | 1.072 | 0.997 |
| 1.991 | 1.111 | 0.988 |
| 1.996 | 1.159 | 0.986 |
| 2.002 | 1.213 | 0.979 |
| 2.007 | 1.257 | 0.977 |
| 2.013 | 1.276 | 0.978 |
| 2.018 | 1.278 | 0.978 |
| 2.024 (top face) | 1.278 | 0.978 |

**No runaway.** The He4 star decays smoothly from the same 1.28 startup value the box and
the idealized cases share, down to 0.978 -- the `0.956-0.981` `tests_r9` already
measured, confirmed here with an independent fresh run. Its outer `dtau_cell` (~0.02-0.03)
sits inside the idealized matrix's "runaway" band (`< 0.08`), yet it behaves like the
"decays to a deficit" branch instead. The opacity there is a real tabulated Rosseland
opacity with a hard floor (`opac_floor = 1e-5 cm^2/g`, `red_giant.cpp` startup print),
not the idealized tests' exactly-held `kappa_const` all the way to `r_out`; a floored,
weakly `rho`-dependent opacity does not reproduce the idealized power-law profile's very
steep *local* density/opacity falloff right at the boundary even at matched `dtau_cell`,
so `dtau_cell` alone is not a sufficient predictor once a real EOS/opacity table is
involved -- it is a necessary condition in the idealized matrix, not a universal one.
**Practically: the He4 production grid is not in the regime that runs away over the
~0.1-turnover timescales these campaigns actually run (200 cycles here is 200/6.2e4 =~
0.003 turnover already reaches the plateau seen at 0.978); the runaway found in section 2
is a real defect of the operator but not the one presently biting the star.**

## 6. The plane-parallel box: no such defect at the same thin layer

`tests_r8/g1_new` (already on disk from 09-17, re-read here, no new run): `box_w8`
production configuration (`he_box_w8.athinput`), mode 3
(`rt_col3_skip_sweep = true` from the box's own printed startup line), `nx2=nx3=16`,
`meshblock 134x8x8`, 50 cycles, top boundary VACUUM (`rt_top_vacuum` default true in
`box_convection.cpp`, see section 7). `rt_surface.bin` carries one record (`rt_surface_dt
= 75.22` s is longer than the 50-cycle run, whose own turnover is only 7.98 s, so only the
`t = 0` record exists):

```
mean F_top = 2.485458e15 erg/cm^2/s,  F_bot (imposed) = 2.47520e15,  F_top/F_imposed = 1.0041
```

**1.004, not the 1.2-1.35 the sphere shows at the same point in its own relaxation.** The
plane-parallel top boundary carries the imposed flux to within 0.4% from the very first
cycle; there is no equivalent of the sphere's startup transient or its later runaway. This
isolates the defect to the spherical geometry specifically -- consistent with `tests_r2/thick`'s
own diffusion-limit analysis (`F_2s/F_req = 1 - H_T/(2r)`, which is identically 1 in the
plane-parallel limit `r -> infinity`) but now demonstrated for the *thin-layer / top-face*
regime as well as the diffusion regime.

## 7. `rt_top_re` / `rt_top_vacuum` as boundary-condition variants under mode 3

* **`rt_top_re = true` is illegal under mode 3 and FATALS**, live-verified:
  `problem/rt_top_re = true` added to `two_stream_sph_thick_relax.athinput`'s `<problem>`
  block (a command-line override on a parameter absent from the file is itself rejected by
  `parameter_input.cpp`, so the parameter has to be *in the file* to override it) and run
  one cycle: `two_stream_rt.hpp`'s own mode-3 gate fatals immediately with
  `rt_top_re = false` listed as a hard requirement alongside `rt_grey`, `rt_use_cons`,
  `rt_src_direct`, `!rt_layer_legacy`, `rt_outer_iter = 1`. There is therefore no way,
  input-file only, to test the "radiate at the ghost's own temperature" vs. "radiative
  equilibrium slab" boundary variants against each other on the sphere while mode 3 is
  active -- only the code's compiled default (`rt_top_re = false`, "the ghost's own Planck
  function") is reachable.
* **`rt_top_vacuum` is not wired into `red_giant.cpp` at all.** `grep` of the pgen finds no
  `GetOrAddBoolean("problem", "rt_top_vacuum", ...)` call; the flag exists only in
  `box_convection.cpp` (default `true` there). So on the sphere it is permanently at its
  namespace-scope compiled default, `false` (`two_stream_rt.hpp:421`) -- i.e. the
  background text's premise that "`rt_top_vacuum` is off under mode 3" is right about the
  *value* but not about *why*: it isn't a switch the red-giant/He4 problem generator ever
  reads, so there is no `problem/rt_top_vacuum = true` override available to test the
  vacuum boundary on the sphere without editing `red_giant.cpp` to add the read (a source
  change, out of scope for this measurement round).
* **Net: neither boundary-condition lever the background text proposed is legally
  reachable from the input file on the sphere under mode 3.** The only boundary model in
  play for every measurement above is the compiled default, "the ghost mirrors the top
  cell's own Planck function" (ghost = active cell's `T`).

## 8. Verdict

1. **Confirmed**: the emergent-luminosity defect lives in the last few cells / the top
   face, not the interior, at the time scale that matters for a star (a few tenths of a
   top-cell thermal time) -- `+20%` to `+35%` there, matching the pre-existing 1.007-1.28
   mode-0 numbers and the He4 star's own 1.278 startup value, and it decays toward
   `L_out/L ~ 0.65-0.98` (deficit) within ~1 top-cell thermal time for every case whose
   `dtau_top1 >~ 0.1`, the He4 star included.
2. **New**: mode 3 (`rt_col3_skip_sweep = true`, the production solver) has a **separate,
   long-time runaway** in the optically thin top layer, distinct from the near-term
   diffusion-limit or startup-transient defects already documented. It appears only for
   `dtau_top1 < ~0.1` (bottom-up cell optical depth of the outermost active cell) and only
   after many (`> ~7`, saturating by `~30-40`) of that cell's own thermal times; it settles
   to a wrong plateau (5x-77x the requested flux) rather than diverging further, and it
   does not appear at all in the plane-parallel box or (over 200 cycles / 0.003 turnover)
   in the He4 star's own tabulated-opacity grid, whose outer `dtau_cell` (~0.02-0.03) is
   in the idealized matrix's runaway band by that number alone -- so `dtau_cell` is a
   necessary but not sufficient trigger once a real floored opacity table replaces the
   idealized cases' held-constant `kappa`.
3. **Not tested, and not testable without a source change**: `rt_top_re` (illegal under
   mode 3, fatals) and `rt_top_vacuum` (not read by `red_giant.cpp`) as alternate top
   boundary models on the sphere. Whatever is driving either the near-term excess or the
   long-time runaway, it cannot be bisected against the boundary condition from the input
   file as things stand.
4. **Recommendation for the next round** (no code touched here): if the He4 production
   grid is ever pushed to a coarser outer cell, a lower opacity floor, or run for enough
   turnovers to approach `~10-40` of its own outer cell's thermal time (its 200-cycle /
   0.003-turnover check here is nowhere close), re-run this same 3-stage recipe
   (`run_relax.sh` + the `nlim=2/3` short-restart pattern of section 3) on the star itself
   before trusting a long production run's top boundary; the idealized matrix here shows
   the transition from "modest deficit" to "50-77x runaway" is sharp in `dtau_top1`, not
   gradual.

## 9. Files

* `two_stream_sph_thick_relax.athinput` -- the thick-test input plus a forced restart
  output.
* `cases.txt` -- the 10-case matrix (8 original + `n8_t100_r2.0`, `n9_t100_r2.0`).
* `n{0,1,3,7,8,9}_t*_r*/ic_thick.txt` -- ICs (8 copied from `tests_r2/thick`, 2 new from
  `make_ic_thick.py --n 8|9 --tau 100 --ratio 2.0`).
* `run_relax.sh` -- the two-stage (stage1 `nlim=200` + stage2 restart) driver for all 10
  cases; each case directory keeps `mltfaces_t0.txt`/`mltfaces_relax.txt`, `.hst` and the
  two stage logs; `.rst`/`.bin` deleted after each dump.
* `he4_relax/` -- the same two-stage recipe on `he4_presn_cs.athinput` (`base9`
  configuration), 200 cycles; `mltfaces_t0.txt`, `mltfaces_relax.txt`, `.hst`,
  `column_he4_ic.txt` kept, `.rst`/`bin`/`cbin` deleted.
* `analyze_top.py`, `topface_tables.txt` -- the last-8-face `F_2s/F_req` and `L_out/L`
  reduction for all 10 thick cases, t=0 vs. 200-cycle-relaxed.
* Section 3's `n7`/`n8`/`n9` short (`nlim=2`+restart `nlim=3`) runs and section 6's box
  read used the same binary and the two-stage pattern but were not kept as separate
  scripts (each a 3-line variant of `run_relax.sh`/already-existing `tests_r8/g1.sh`); the
  commands are given inline in sections 3 and 6.
