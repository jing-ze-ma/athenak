# `rad_angular`: an on/off switch for the HORIZONTAL radiative conduction

`<hydro>/rad_angular` and `<mhd>/rad_angular` (boolean, **default `true`**) switch the
horizontal (x2/x3, "angular"/"transverse") half of the radiative conduction on and off.
`true` is bitwise today's arithmetic. `false` removes the horizontal operator:

* `Conduction::AddIsotropicHeatFluxRadiative` returns after the x1 faces, so **no
  radiative heat flux is added to any x2/x3 face** — on every mesh type, because the
  spherical-polar and cubed-sphere seam terms live inside those same two kernels
  (`radcond2`, `radcond3`);
* `Conduction::BuildAngularCoeffs` returns immediately, so no frozen transverse
  conductances, no `cap_x`/`cap_c2`/`cap_c3`, no `rad_cap_ang` report — this also covers
  the problem generators that call it themselves under `problem/rt_split_transverse`;
* `Conduction::NewTimeStep` drops the transverse conduction limits `dt2` and `dt3`.

Untouched: the x1 face flux, `rad_implicit_x1`, the tau blend (`rad_tau_lo/hi`,
`rad_blend_radial`), `rad_flux_inner`, the two-stream hand-over (`rad_blend_use_2s`,
`rad_f2s`) — and ordinary (`constant` / `spitzer`) thermal conduction, which never reads
the switch (the dt hunk is guarded on `radiative`).

Combinations:

| with `rad_angular = false` | behaviour |
| --- | --- |
| `rad_implicit_ang`, `rad_sts_all`, `rad_sts_split`, `rad_ang_solver = adi` | **FATAL** at start-up (they *are* the transverse treatment) |
| `rad_cap_ang > 0` | accepted **silently**, becomes a no-op — the production inputs carry `rad_cap_ang = 0.5` and need no edit |

One line is printed at start-up:

```
Conduction: rad_angular = false -- no radiative heat flux through any x2/x3 face, no
angular coefficients, and the transverse conduction dt (dt2, dt3) is dropped
```

Caveat (not fixed here, it is in `mesh.cpp`): the dt-collapse report keys its `"capped"`
label on `rad_cap_ang` alone, so with `rad_angular = false` **and** `rad_cap_ang = 0` it
prints the sentinel `dt2 = dt3 = -1` instead of a word. `-1` means "no constraint".

Touched files: `src/diffusion/conduction.hpp`, `src/diffusion/conduction.cpp`. Nothing
else. cpplint (`--filter=-build/include_subdir`) and the custom checks (90 columns, tabs,
`}}`, trailing whitespace, mode 644) are clean on both.

## The test

`base.athinput` = `inputs/tests/dhj_ck_spherical.athinput` (3-D **spherical polar**,
`deep_hot_jupiter_rt`, correlated-k RT, tabulated general EOS, `nx1 = 128`,
`nx2 = nx3 = 8`, 2 MeshBlocks, `isotropic_conduction = radiative`, `rad_implicit_x1 =
true`, `rad_cap_ang = 0.5`, tau blend 30/300) with three keys appended to `<hydro>` at
their default, inert values — `rad_angular = true`, `rad_tr_tau_lo = 0.0`,
`rad_tr_tau_hi = 0.0` — because AthenaK refuses a command-line override of a key that is
not in the file. A 1-D column has no x2/x3 faces, which is why this 3-D few-column
variant is used. Driver: `./go.sh <dir> <binary> [overrides]`, 20 cycles serial CPU,
hst every cycle, one `hydro_w` binary dump every cycle. Payloads are compared with
`cmpbin.py` (the `.bin` *header* embeds the parameter dump and therefore differs whenever
an input key differs; only the decoded variable arrays are physics).

Binaries: `athena_ref` = the unmodified tree at HEAD `131eabe8`; `athena_new` = the same
tree plus these two files. Both `-DPROBLEM=deep_hot_jupiter_rt`, Release, CPU serial.

## Gates

**(a) default is bitwise the reference.** `ref` (athena_ref) vs `new` (athena_new, no
overrides), a run on which the angular radiative conduction is active
(`sum|delta flx.x2f(IEN)| = 1.4e3 … 2.6e3` per call, see (c)):

```
dhj.hydro.hst          bitwise identical (cmp)
dhj.hydro_w.00021.bin  BITWISE IDENTICAL payload (40960 values)
```

**(b) `rad_angular = false` == the taper route.** `off` (`hydro/rad_angular=false`) vs
`tap` (`hydro/rad_tr_tau_lo=1e30 hydro/rad_tr_tau_hi=1e31`), same run:

```
dhj.hydro.hst          bitwise identical (cmp)
dhj.hydro_w.00021.bin  BITWISE IDENTICAL payload (40960 values)
```

The taper route **is honoured on this spherical-polar mesh** — no fatal, no silent
no-op; it needs the tau blend (`rad_tau_hi > 0`), which this input has, and it prints
`Conduction: transverse tau taper on, w = 0 below tau = 1e+30`. A GPU production arm can
rely on it. (It is the more expensive way: it still builds and runs the angular kernels
with a weight of exactly `0.0`; see (e).)

For reference, the switch is not vacuous on this problem — `ref` vs `off` after 20
cycles: `max|d| = 1.22e-04`, `max rel = 1.11e-07`.

**(c) the x2/x3 radiative face fluxes are exactly zero.** Measured, not inferred, with a
throw-away instrumented build (`-DRADANG_PROBE`, since reverted): `flx.x2f` / `flx.x3f`
were deep-copied before `AddIsotropicHeatFluxRadiative` and
`sum|flx - flx_before|` reduced after it, i.e. exactly what conduction added to the
transverse faces.

| run | `sum|delta x2f|` | `sum|delta x3f|` |
| --- | --- | --- |
| `rad_angular = true`, `rad_cap_ang = 0` | 1.40e3, 1.58e3, 2.63e3 (cycles 2–3) | 3.11e-1, 3.19e-1, 4.59e-1 |
| `rad_angular = false`, `rad_cap_ang = 0` | **0.000000e+00** every call | **0.000000e+00** |
| `rad_angular = false`, `rad_cap_ang = 0.5` | **0.000000e+00** every call | **0.000000e+00** |

The x1 side is unaffected: the pre-conduction `flx.x2f/x3f` sums (which the hydro
fluxes set) agree digit for digit between the on and off runs at the same cycle
(`5.202026e+13`, `1.648084e+11` at cycle 3), i.e. the radial/hydro state is identical
while the transverse contrast is or is not being relaxed. Over 20 cycles the two states
then separate by the `1.1e-7` above, which is the transverse transport that was removed.
(The stronger form asked for — the radial result identical to a run with no horizontal
contrast at all — is not constructible here without editing the problem generator, which
was out of scope.)

**(d) the conduction dt no longer sees dt2/dt3.** Same instrumented build, printing
`Conduction::dtnew` and the winning cell's `dt1/dt2/dt3` from `dt_diag`. With
`rad_cap_ang = 0` (so the transverse limits are live; with the production
`rad_cap_ang = 0.5` they are already dropped):

```
rad_angular = true   ### RADANG_DT cycle 0 conduction dtnew = 6.20192e+07
                         | dt1 = -1  dt2 = 6.20192e+07  dt3 = 2.29159e+09
rad_angular = false  ### RADANG_DT cycle 0 conduction dtnew = 2.99616e+307
                         | no winning cell: the operator constrains nothing
rad_tr_tau_lo = 1e30 ### RADANG_DT cycle 0 conduction dtnew = 2.99616e+307
                         | no winning cell: the operator constrains nothing
```

`dt1 = -1` is `rad_implicit_x1`, so `dt2` *was* the conduction timestep; with the switch
off radiative conduction constrains nothing at all. (The run's own dt is hydro-limited at
2.9e1 s either way, so the wall-clock win here is (e), not a larger step.)

**(e) cost.** 60 cycles, serial CPU, same node, best of two:

| | cpu time | zone-cycles/cpu_second |
| --- | --- | --- |
| `rad_angular = true` | 10.99 s | 4.472e4 |
| `rad_angular = false` | 10.16 s | 4.836e4 |

**+8.2 %** zone-cycles/s (−7.5 % cpu). The taper route `rad_tr_tau_lo = 1e30` measures
4.487e4, i.e. the same as `true` — it buys the physics but none of the cost, which is the
reason this switch exists alongside it.

**Validation.** `rad_angular = false` plus `rad_implicit_ang = true`, `rad_sts_all =
true` or `rad_ang_solver = adi` each exit 1 with the conduction.cpp fatal (the check sits
right after `rad_ang_solver` is parsed, so it fires before the unrelated mesh guards of
`rad_implicit_ang`); `rad_angular = false` with the input's `rad_cap_ang = 0.5` starts
normally.

## Files kept here

`base.athinput`, `fatal.athinput`, `go.sh`, `gof.sh`, `cmpbin.py`, and each run's
`run.log` + `dhj.hydro.hst`. The binaries, the `bin/` and `rst/` dumps, `build_radang`
and `build_radang_probe` were deleted after the gates were read.
