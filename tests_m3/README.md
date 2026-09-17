# Mode-3 column solve on a SPHERICAL mesh (2026-09-17, viper, `he4-presn-global`)

Follow-up to `tests_r2/` (the r^2 dilution fix, commit `1159a8f3`) and to `tests_1d/`
gate D(iii), which reported that the mode-3 exact column solve
(`problem/rt_implicit_column = 3`) was six orders out on the 1-D He column while the
explicit sweep (`rt_implicit_column = 0`) was correct.

Everything here is CPU, `cmake -B build_cpu -D PROBLEM=red_giant` (Serial Kokkos, no
MPI) and `cmake -B build_cpu_box -D PROBLEM=box_convection` for the box regression.
All `.bin`/`.rst`/`.cbin` and the 0.5 MB per-run `column*.txt` IC dumps have been
deleted, and the few logs above 200-400 kB were truncated (every number in the tables
below was read off the full log before truncation; the commands regenerate them).

    cmake -B build_cpu     -D PROBLEM=red_giant     && (cd build_cpu     && make -j32)
    cmake -B build_cpu_box -D PROBLEM=box_convection && (cd build_cpu_box && make -j32)

## Verdict

**The radial arithmetic of the mode-3 residual is CORRECT.**  Cell by cell, at the entry
state, the mode-3 formal solution reproduces the sweep's source to round-off on a
spherical mesh, in the transparent limit AND at tau/cell = 18 (tables below); the
converged solution telescopes to the faces to 5e-12; the serial (`thomas`) and
partitioned (`pcr`) paths agree to every printed digit; and the Newton residual reaches
2e-12 in four passes.  No missing or spurious `Ac`/`Av`/`W` factor was found in the
residual, the cut/interface injection, the `Dtop`/`Ucut`/`fbot`/`FL` scaling, the
`wrflux` writer or `stat(16)`.

**One real regression of `1159a8f3` was found and fixed: the areas were ABSOLUTE.**  The
solve's unknowns `D_q`, `U_q` became `J = A I ~ 1e22 x b`, while the fifth unknown of the
same 5x5 block is `b` itself.  The transport rows' `b`-column entries (`t^2 A`, and every
source coefficient, all O(A)) therefore sat 22 decades above their unit diagonals, the
row equilibration in `RTCol3Inv5X` drove the transport pivots to ~1/A, and the
single-precision pivot test fired on **every** block: `problem/rt_impl_mixed = 2` --
the production setting -- silently degenerated into the double path plus a wasted float
factorisation.  Measured on the 1-D He column, per RT call:

| `rt_impl_mixed = 2`, pcr | `nmixfb` (float blocks redone in double) |
| --- | --- |
| before | **10464** (= every cell of every column) |
| after  | **0** |

The fix carries all three geometry helpers RELATIVE to the column's top face,
`Ar = A(ie+1)`.  The scaling cancels identically -- every source endpoint, boundary datum
and face flux carries one power of the area and every divergence one power of `1/V`, so
`J`, `Src`, `Fb`, the Jacobian and every ratio-valued reduction slot are unchanged -- but
`A(i)/A(ie+1)` is O(1), so the block is commensurate again.  Under
`problem/rt_plane_parallel` `Ar` is not read at all and the box is bitwise.

**The 1-D He 5-turnover gate still fails, and it is not the geometry** -- see
*He 1-D: the remaining failure* below.

## A. Transparent spherical gate -- PASS

`thin_ex.athinput` = `inputs/hydro/red_giant_1d_dilution.athinput` (cubed sphere,
`nx2 = nx3 = 4`, 128 uniform radial cells, `r` 1.6e12 .. 4.0e12 so the area ratio is
6.25, `kappa_const = 1e-5`, `vpert = 0`) plus `rt_col3_ex_iter = true`, i.e. the
production handover.  `L = 4 pi r^2 F` from the solver's own face flux
(`problem/rt_apply_debug`), read with `../tests_r2/rg1d/lprof.py` / `lprof3.py`.

    cd tests_m3/g_m3p && $A -i ../thin_ex.athinput problem/rt_implicit_column=3 \
        problem/rt_impl_solver=pcr problem/rt_impl_mixed=0 problem/rt_impl_warm=0 \
        time/nlim=4 problem/rt_apply_debug=1000

| arm | `L` max/min, cycle 0 | cycle 4 |
| --- | --- | --- |
| mode 0 (`g_m0`) | 1.000305 | 1.000705 |
| mode 3, `thomas`, mixed 0, warm 0 (`g_m3s`) | 1.000298 | 0.998674 |
| mode 3, `pcr`, mixed 0, warm 0 (`g_m3p`) | 1.000298 | 0.998674 |
| mode 3, `pcr`, mixed 2, warm 1 (`g_m3pm`) | 1.000298 | 0.998676 |

Gate was `< 1.005`.  `thomas` and `pcr` agree to all printed digits.

## B. Optically THICK spherical gate -- PASS

`thick.athinput` is the same column with `ptop = 3.0e6` (density x 1e6, so the column
optical depth is 200 instead of 4e-4) and `rad_tau_lo/hi = 1e5/1e6`, i.e. the whole
column is carried by the two-stream as in the He production.  And, for the real thing,
the 1-D He column itself (`he1d.athinput`, `tau/cell = 18` at the base, column `tau` 362).

The gate is the operator identity, which is what an area error in the diffusive limit
would break:

* **entry state** (`problem/rt_impl_maxit = 1`): mode 3's own per-cell source `srcdx`
  must equal the sweep's `sw_srcdx`, because the mode-3 formal solution at `b = b_entry`
  IS the sweep;
* **converged** (`maxit = 8`): `srcdx` must equal `divF` built from the same converged
  intensities -- the energy row and the transport rows telescoping to the faces.

Read off `problem/rt_outer_verbose = true`'s `### rt_col3_flux` dump with `dumpstat.py`.

| column | tau/cell (base .. top) | max rel `srcdx - sw_srcdx` (entry) | max rel `srcdx - divF` (converged) |
| --- | --- | --- | --- |
| thick synthetic (`k_m3s`) | 0.77 .. 1.9e-13 | **4.6e-14** | **4.6e-14** |
| He 1-D (`he_dump`) | 18 .. 4.7e-05 | **2.7e-12** | **4.5e-12** |

Global slots on the thick synthetic, mode 3 converged: `telescope_rel = 2.6e-15`,
`budget_rel = -3.5e-09`, `Ftop_impl/Ftop_sweep = 0.99564`.

Cycle-0 solve of the He column is identical across the three solver configurations:

| arm | `resid_max` | `Ftop_impl/Ftop_sweep` | `nmixfb` |
| --- | --- | --- | --- |
| `thomas`, mixed 0, warm 0 (`he_f_m0w0`) | 2.384996e-12 | 0.7998022 | 0 |
| `pcr`, mixed 0, warm 0 (`he_f_m0p`) | 2.385005e-12 | 0.7998022 | 0 |
| `pcr`, mixed 2, warm 1 (`he_f_m2w1`) | 2.394806e-12 | 0.7998022 | 0 |

(the same three arms BEFORE the fix, `he_m0w0` / `he_m0w0p` / `he_m2w1`, carry
`nmixfb = 0 / 0 / 10464`.)

## C. Box regression (plane-parallel) -- BYTE-IDENTICAL

`bench/hestar_fecz/box_w8/he_box_w8.athinput` in the production mode-3 configuration
(`rt_implicit_column = 3`, `rt_impl_solver = pcr`, `rt_impl_mixed = 2`,
`rt_impl_warm = 1`, `rt_col3_skip_sweep = true`), shrunk to `nx2 = nx3 = 16`,
`meshblock 134x8x8`, 20 cycles, serial CPU.  Reference binary built from HEAD
(`1512b882`), new binary from this working tree.

    cmp box_ref/feczrt.hydro.hst box_new/feczrt.hydro.hst   -> IDENTICAL
    cmp box_ref/feczrt.user.hst  box_new/feczrt.user.hst    -> IDENTICAL
    cmp box_ref/feczrt.log       box_new/feczrt.log         -> IDENTICAL
    cmp box_ref/column_used.txt  box_new/column_used.txt    -> IDENTICAL
    cmp box_ref/rt_surface.bin   box_new/rt_surface.bin     -> IDENTICAL
    cmp box_ref/rt_profile.bin   box_new/rt_profile.bin     -> IDENTICAL
    diff -r box_ref/bin box_new/bin                         -> IDENTICAL
    diff -r box_ref/cbin_hydro_w_2 box_new/cbin_hydro_w_2   -> IDENTICAL
    diff -r box_ref/rst box_new/rst                         -> IDENTICAL

(the `.bin`/`.rst` pairs were compared and then deleted for quota; `column_used.txt` is
kept only as its md5 pair in `box_column_used_md5.txt`, the `.hst` and the event log are
kept in `box_ref/`, `box_new/`.)

## D. He 1-D: the remaining failure -- NOT the mode-3 arithmetic

`he1d.athinput` = `inputs/hydro/he4_presn_cs.athinput` with the debug knobs added, run as
`mesh/nx2=4 mesh/nx3=4 meshblock/nx2=4 meshblock/nx3=4 problem/vpert=0.0`.

| arm | t at cycle 500 | `L_rad,cut/L` | `L_rad,out/L` | `it_mean` |
| --- | --- | --- | --- | --- |
| mode 0 (`he_g500m0`) | 6392.2 s | 0.6986 | -- | -- |
| mode 3, pcr, mixed 2, warm 1, `rt_rad_force=false` (`he_g500`) | 552.4 s | 8.168 | 4.21e+03 | 4.48 |

Still a failure (the `tests_1d` D(iii) number was `L_out/L = 8.6e5`; the conditioning fix
buys two orders but not the gate).  With the production `rt_rad_force = true`
(`he_gA`) the run collapses exactly as `tests_1d/gA.log` reported: `dt` falls to 2e-13,
`t` sticks at 58.5 s, and `rt_use_cons gave a non-positive internal energy in 36 cell
reads`.  5 turnovers is unreachable, so the `eos_fail`/mass-drift gates cannot be read.

**Why it is not the solver.**  Backward Euler must converge in `dt`, and it does.  Mode 3,
`thomas`, mixed 0, warm 0, `rt_rad_force = false`, all run to `t = 250 s`:

| `time/cfl_number` | `L_rad,cut/L` at t = 250 s |
| --- | --- |
| 0.3 (`he_c0.3`) | 2.419 |
| 0.1 (`he_c0.1`) | 2.070 |
| 0.03 (`he_c0.03`) | 1.945 |

First order in `dt`, converging on ~1.9.  Mode 3 is solving the right problem; the
problem itself moves.

**What moves it.**  The imported 1-D structure is nowhere near two-stream radiative
equilibrium at the base.  At cycle 0, from the sweep's own face fluxes (`he_p0/run.log`,
`rt_apply`): `4 pi r^2 F` = 1.96e38 at the cut face and 1.12e38 eight cells above it --
the deep radiative luminosity loses 43 % over 8 cells, and the bottom cell's
`e/(src) = 2.4e7/9.97e4 = 240 s` thermal time against that imbalance.  In the real star
the missing flux is convective; in a 1-D column with `vpert = 0` there is nothing to carry
it, so the base must heat.

**Why mode 0 looks healthy anyway.**  It is frozen, not balanced.  Both modes compute the
SAME source -- at `i = 3`, cycle 0, `src = 9.9704e+04` in both, bit for bit -- but mode 0
divides the applied `de` by its ALI/semi-implicit factor `lamdt = 4 Em/e bdt = 8977`:

    mode 0:  src = 9.9704e+04  lamdt = 8.9767e+03  de = 7.1802e+01  de/e = 3.0e-06
    mode 3:  src = 9.9704e+04  (skip_de)           de = 6.4454e+05  de/e = 2.7e-02  (= src*bdt)

`lamdt` is the LOCAL emission rate; at `tau = 340` emission and absorption cancel and the
net source is the flux divergence, which that factor has no business damping.  So mode 0's
stable 1-D He column is the ALI diagonal holding the base at its initial entropy, and the
`tests_1d` D(iii) comparison was between a solved column and a frozen one.

**What would actually gate.**  Either seed/allow convection in the 1-D column (impossible
by construction), or give the base the flux the star carries convectively -- i.e. an inner
boundary that injects `L` rather than `rad_flux_inner = 0` with `inner_bc = open` --
before asking for a 5-turnover 1-D run.  That is a setup decision, not a code fix, so it
is left to the user.

## E. Two diagnostics that were wrong, now fixed

* `problem/rt_apply_debug`'s `divF` still printed the PLANE-PARALLEL `-(Ft-Fb)/dx1` on a
  radial mesh; it now prints `-(A_top Ft - A_bot Fb)/V` and agrees with `srcraw` to 13
  digits (`t_area/run.log`).  The line also gained `Ab/At/Ac/V` so the mesh geometry can
  be read off it; `tests_r2/rg1d/lprof.py` still parses it.
* `problem/rt_force_verbose` normalised by `problem/rt_force_grav = GM/r_in^2`, a constant
  that is 4x the local `g` at the He photosphere, and subtracted that same constant in
  `resid`.  It now uses `EffGravAt(r)`, the gravity the force actually works against
  (`he_fv/run.log`: `i = 91`, `a_p/g = 1.549`, `a_f/g = -0.565`, `resid/g = -0.0154`,
  against the hand conversion 1.58 / -0.575 / -0.013 in `tests_1d/README.md`).
* `### rt_col3_flux`'s `F3lo` omitted the lower face's area weight at `i = icut` and
  printed a J-flux next to a per-area `Fb`, so its `divF`/`dif` columns were meaningless
  in the cut cell (`-3.1e38` against a source of `2.6e35`); it now carries the weight and
  prints `F3lo/A`.  `### rt_col3_deep`'s `srcdx` used `Dx` where the solve uses `Vc`.
* `red_giant`'s opacity-table window refusal now skips when `problem/kappa_const > 0`
  replaces the table -- without that the thick gate above cannot be built.

## Files

* `thin_ex.athinput`, `thick.athinput`, `he1d.athinput` -- the three test inputs.
* `lprof3.py` (L(r) from `rt_apply`, first or last RT call), `dumpstat.py`
  (`rt_col3_flux` operator/telescoping metrics), `lprof2.py`, `dumpcmp.py`.
* `g_*` transparent gate, `k_m3s` thick synthetic, `he_dump` / `he_f_*` / `he_m*` the He
  operator and conditioning arms, `he_c*` the dt convergence, `he_g500*` / `he_gA` the
  gate attempts, `box_ref` / `box_new` the plane-parallel regression, `t_area` the
  `divF` check.
