# The NIGHT-SIDE TOP SLAB under the centre-to-centre correlated-k layers

What `d8f00d49` (centre-to-centre layers, `problem/rt_layer_legacy = false`) does to the
tenuous night-side top of a deep-hot-Jupiter column, reproduced on a CPU, and what the
fix is.  Everything below was run with `build_cktop` (`-D PROBLEM=deep_hot_jupiter_rt`,
serial CPU) on `inputs/tests/dhj_ck_spherical.athinput` — the production radial grid
(`nx1 = 128`, `x1 = 9.44e9 .. 2.0556e10`, the polynomial stretch, the production floors
`pfloor = 1e-3`, `dfloor = 5e-14`, `tfloor_kelvin = 200`) shrunk to `8 x 8` columns.
`base.athinput` is that file plus four lines that only make existing switches reachable
from the command line (`rt_layer_legacy`, `rt_floor_consistent`, `rt_outer_verbose`,
`rt_apply_debug`).

## 1. The reproducer

From the analytic initial state the night side cools and drains; by cycle ~1400 its top
5-10 cells sit on `dfloor = 5e-14` at `T ~ 1.4-1.6e3 K` and `p` within a factor of a few
of `pfloor`.  From there on, **4000 cycles, identical input, identical final time
`t = 1.087e5 s`**:

| run | layers | `pfloor` | `eos_efloor` | `eos_dfloor` | cells clipped by `rt_de_max` | total E |
| --- | --- | --- | --- | --- | --- | --- |
| `leg` | legacy | 1e-3 | **0** | 2.089e6 | 0 | 6.38769e38 |
| `c2c` | centre-to-centre | 1e-3 | **27048** | 2.088e6 | 2 | 6.38758e38 |
| `fix` | c2c, `rt_floor_consistent = true` | 1e-3 | 27048 | 2.088e6 | 2 | 6.38758e38 |
| `pf5` | c2c | **1e-5** | **0** | 2.087e6 | 2 | 6.38758e38 |

`c2c` against `leg`: the energy-floor rate goes from nothing to a steady **16 cells per
cycle** and the total energy is **1.7e-5 lower**.  That is the same signature the
production cubed-sphere restart shows (efloor x10, `dE/E = -4.5e-5 .. -8.9e-5`).

The A/B also runs **from one shared state** (`f15n/rst`, cycle 2000, written by the run
itself): 560 `eos_efloor` per 35 cycles with the c2c layers against 0 with the legacy
ones, the difference appearing within the first interval.  That is the cleanest form of
the reproducer and costs two 300-cycle restarts.

## 2. Which term

`problem/rt_apply_debug` (wired to the input by this work; it prints `src`, `Em`, `de`,
`de/e` for the top cells of one column) on the same state, night column `mu0 < 0`
(`m = 0, k = 2, j = 5`), top cell `i = ie`:

| | `src` | `Em` | `A = src + Em` | `de/e` |
| --- | --- | --- | --- | --- |
| centre-to-centre | 4.889e-6 | **1.271e-5** | 1.760e-5 | +1.71e-2 |
| legacy | 9.161e-6 | **8.724e-6** | 1.788e-5 | +3.30e-2 |

**The absorption is the same to 1.6 %; the whole difference is the top cell's OWN
EMISSION, which the c2c layers raise by 46 %.**  It is the c2c value that is right: the
closing half layers give the cell `kappa rho B_ie` over its own two halves, once per ray,
while the staggered layer wrote `kappa rho (B_ie + B_ghost)/2` — an average with the
GHOST cell above the domain, which here has `B_ghost ~ 0.37 B_ie`.  The legacy top cell
was therefore under-emitting and sitting artificially hot, with its pressure
artificially above `pfloor`.

Checked and **not** guilty: the top closure does not double count (the down sweep takes
the upper half of `ie` and, on the next iteration, its lower half; the up sweep takes the
lower half in the `i = ie-1` iteration and the upper half in the closing block, so each
half is crossed exactly once per ray), and `Em_g` is accumulated once per cell.

Also **refuted, on this reproducer**: the mechanism of `bench/bisect_me/README.md` round 2
("the net source grows until `e/|src| < dt`, the `rt_de_max` clip fires, the residual
over-cooling is caught by the energy floor").  A scan of 24 columns at the reproducer
state gives `max |de/e| = 3.3e-2`, fifteen times below the `rt_de_max = 0.5` clip; the
clip fired on 2 cells in 4000 cycles, and never in the floored columns.  The `eos_efloor`
cells are top cells whose pressure sits AT `pfloor`, trimmed once per cycle.

## 3. The fix

### 3a. `problem/rt_floor_consistent` (new, default `false` = bitwise today)

Bounds the RT decrement by the state the EOS floors would restore the cell to:
`e^{n+1} >= max(e(rho,tfloor), e(rho,pfloor))`, cooling steps only, one EOS evaluation per
cooling cell.  It is the right safety net for a stiff floored cell — and on this
reproducer it is a **measured no-op**: 4000 cycles with it on reproduce `c2c` to every
digit of the history file and to the last `eos_efloor` count (27048).  The sweep never
takes a cell across the floor here.  It is kept, off, for the case the production state
does reach (the clip warning did fire there).

### 3b. What actually fixes it: `<hydro>/pfloor`

The corrected emission makes the night-side top slab genuinely colder, so its pressure
falls to `1e-3 barye` — which is the floor.  Lowering the floor by 100x removes **all**
27048 events and changes nothing else:

```
                         i=122     i=125     i=127     i=129 (top)
T [K]   pfloor 1e-3      1575.5    1537.7    1495.5    1405.9
        pfloor 1e-5      1575.5    1537.8    1495.5    1405.9
e       pfloor 1e-3    3.4820e-2 2.6577e-2 1.9396e-2 1.0416e-2
        pfloor 1e-5    3.4826e-2 2.6579e-2 1.9397e-2 1.0417e-2
```

identical to five digits, same `dt` (both reach `t = 1.087e5` in 4000 cycles), same
`eos_dfloor`, same total energy to six digits.  `pfloor = 1e-3` was chosen (see the note
in the input) because at that value "the floor never fires at all" — that calibration was
made with the staggered layers and does not survive the correction.

### 3c. The energy

`pf5` (no floor events at all) has the same total energy as `c2c` to six digits, so the
`-1.7e-5` against the legacy layers is **not** the floor laundering energy: it is the
extra radiative loss of a top slab that is now emitting at its own Planck function.  On
this evidence the production `dE/E = -4.5e-5 .. -8.9e-5` is physical, not an artefact.

## 4. Gates

* **Bitwise default.**  `inputs/tests/dhj_ck_spherical.athinput`, 50 cycles, both
  `ck_spherical = false` and `true`, 52 dumps each: every data array bit-identical to the
  reference binary built from the same tree before this work.  (The files themselves
  differ by the length of the embedded `PAR_DUMP`, which now carries the new switches.)
  That is the dhj gate (c) of `tests_cleanup_0922/gates.sh`.
* Reference binary reruns are bit-identical to themselves (determinism control).
* `red_giant`, which shares `two_stream_rt.hpp`, compiles clean.
* `cpplint` on the two touched files reports nothing on any added line; no line over 90
  columns, no tabs.
* The `tests_ck_sph` thermal gates are untouched by construction: the default arithmetic
  is bitwise unchanged, and 3b is an input value, not code.

## 5. GPU confirmation (2026-09-22, apudev job 11941092)

Done on the real 3-D cubed-sphere state.  Full numbers in `bench/pfloor_ab/README.md`.
Two arms restarting `bench/cs_mhd_prod3/rst/dhj.00567.rst` (rot ~283) with the production
`ck_spherical + ck_beam_sph` input, the **same binary**, 2 GPUs, 13 min each, differing in
one line: `bench/ck_sph_ab/sphbeam` (`pfloor = 1e-3`, the control) and
`bench/pfloor_ab/p1e5` (`pfloor = 1e-5`).  8200 cycles, `t = 8.64096e7 -> 8.6515e7`.

**It holds.**  `eos_efloor` falls **7.4x** (3391 -> 457 per cycle) and in the written
state the pressure floor is gone completely: 1 active cell of 786432 sat exactly on
`pfloor` in the control, **0** in the `1e-5` arm; the lowest active pressure goes from
`1.0000e-3` barye (on the floor) to `8.135e-4` barye, 81x above its floor.
`dt` is **identical**
(same plateau `1.285230e+01`, same set of values), throughput costs 0.4 %, `eos_dfloor` is
unchanged (+0.4 %), and `eos_fail`, `fofc`, `c2p_it`, `vceil` are 0 in both.  Night-side
T(p) at 1e-6 .. 1e-3 bar agrees to <= 15 K (0.9 %), and the history agrees on mass and
total energy to all printed digits, with the same `dE/E = -2.0273e-4` — so section 3c's
conclusion (the c2c energy loss is physical, not the floor laundering energy) survives on
the production state.

Two corrections to the text above:

* the key is **`<mhd>/pfloor`** for this MHD production run, not `<hydro>/pfloor`: it is
  read from the physics block (`src/eos/eos.cpp:109`).  The units are code pressure =
  cgs barye, so `1e-3` means `1e-9` bar, not `1e-3` bar.
* `eos_efloor` does **not** go to zero here as it did on the CPU reproducer.  The residual
  457/cycle, and a compensating **4.25x rise in `eos_tfloor`** (751 -> 3192 per cycle),
  are
  ghost-zone/intermediate work at the top radial boundary and the cubed-sphere seam halos:
  the counters cover every C2P call over the full index range, ghosts included
  (`mhd_tasks.cpp:694`).  **No active cell is near the 200 K floor in either arm** — the
  coldest is 456 K against 561 K in the control, and zero cells sit on `tfloor` — so this
  is a bookkeeping swap, not new floor activity in the solution, and it carries none of
  the
  50 K-floor failure signature (`dt` 10 s -> 1e9 s in ~1500 cycles); `dt` here is
  unchanged
  to every digit.

## 6. Left open

* Unrelated, found on the way: `src/eos/eos.cpp:87-91` fatals on a **hydro restart**
  (`<hydro>/hlld_bx_zero_tol is implemented only for MHD`) because the parameter is
  recorded into the restart by its own `GetOrAddReal`.  Every restart here had to be run
  through a byte-preserving rename of that key.  Owner: whoever is holding `src/eos/`.
