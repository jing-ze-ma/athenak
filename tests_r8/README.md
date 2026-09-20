# he4_presn: THE TAPER NOW HAS A TEMPERATURE GATE, AND IT REMOVES THE 1.18-TURNOVER
# KILLER -- BUT A SECOND, FORCE-INDEPENDENT DRAIN IS LEFT

2026-09-17, viper, branch `he4-presn-global`, against `65cd5bd1` (= the `tests_r7`
commit).  Read `tests_r7/README.md` sections 5 and 5.1 first: this round executes its
proposed remedy 1.

Binaries: `tests_r8/athena_box_ref` / `athena_box_new` (serial CPU,
`PROBLEM=box_convection`, gate 1) and `tests_r8/athena_v4` (HIP,
`Kokkos_ARCH_AMD_GFX942_APU`, MPI, `PROBLEM=red_giant`, gates 2-3).  Every `.bin`,
`.rst` and `.cbin` of the star runs has been deleted after measuring.  Turnover = 4705 s.

---

## 1. WHAT WAS BUILT, AND WHY IT IS A TEMPERATURE GATE AND NOT A RADIUS GATE

`tests_r7` section 5.1 asked for `w = max(w_rho, w_tau)` -- "gate `w` on optical depth or
radius as well as density, so that an interior cell cannot be handed the surface treatment
because its own density dropped".  What is implemented is

    w_eff = max( w_rho(log10 rho) , w_T(log10 T) )

with `w_T = 1` at and above `<hydro>/eos_rad_t_hi`, `0` at and below `eos_rad_t_lo`, the
same C1 cubic `3s^2 - 2s^3` in between, and the two temperatures read off the SAME two
optical depths the density window came from (tau = 3.0 and 0.30).  Both new parameters
default to zero; with `eos_rad_t_hi = 0` every expression is on its original branch.

| file:line | what |
| --- | --- |
| `src/utils/rad_taper.hpp:76-124` | `rad_taper::WeightGated(x, xlo, xhi, y, ylo, yhi, gate, w, dwdx, dwdy)` -- `Weight()` then `max` with the T smoothstep.  Exactly ONE of `dwdx`, `dwdy` is non-zero, because `max()` picks a branch; `gate = false` returns `Weight()`'s own answer with `dwdy = 0`, bit for bit |
| `src/eos/eos_table.hpp:107-113` | the two new fields `rad_tgate`, `rad_lt_lo/hi` |
| `src/eos/eos_table.hpp:289-321` | `EvalFromLogs` -- the gated weight, and the `dw/dlnT` term the gate adds to `chi_T`, `c_v` and `dln e/dln T`, so each is still the derivative of the `p` and `e` returned beside it |
| `src/eos/eos_table.hpp:420-476` | the same in `EvalEOnly` / `EvalPOnly` (both currently unreferenced, kept in step) |
| `src/eos/eos_table.hpp:517-542` | `EvalResidual` -- the gate's `dw/dlnT` in the root find's derivative for MODE 0 and 1.  It is NON-NEGATIVE, so `e` and `p` stay monotonically increasing in `T` and `SolveLog()`'s bracket update stays valid |
| `src/eos/eos_table.cpp:131-166` | the parameter read and its two guards (`eos_rad_t_hi` without a density window, and `t_lo >= t_hi`, are both fatal) |
| `src/eos/eos_table.cpp:295-345` | the pressure-floor bound `efbnd` is taken over BOTH ends of the gate (the weight at a given density is no longer known at build time), and `eminspec` picks up the gate's own weight at the table's lowest row, where it is exact |
| `src/utils/two_stream_rt.hpp:4082-4088, 4776-4830` | the FORCE: the same gated weight, and the `P_rad dw/dlnT grad T/T` partner of the existing `P_rad dw/dlnrho grad rho/rho`, radial and both transverse, on the same arc-length derivatives |
| `src/pgen/box_convection.cpp:652-676` | the box's `E_rad` diagnostic uses the gated weight too |
| `bench/hestar_presn/make_ic_sph.py:57-58, 146-160` | the ic builder's `w` is `max(w_rho, w_T)` on the SAME table `T` the run will solve for |
| `inputs/hydro/he4_presn_cs.athinput:125-144` | `eos_rad_t_hi = 6.3281e4`, `eos_rad_t_lo = 4.5213e4` |

### 1.1 WHY NOT RADIUS OR tau, HONESTLY

The brief asked for the radius gate first and called it "the minimal, always-consistent
option".  It is not minimal in this code, and that is a property of the EOS interface
rather than a choice: `EOSTable` is evaluated per cell from `(rho, T)` and NOTHING else --
no `(m,k,j,i)`, no coordinates, no column.  Threading a radius to the point where `w` is
formed means an extra argument on `EvalFromLogs`, `EvalImpl`, `Eval`, `EvalNoMu`,
`EvalEOnly`, `EvalPOnly`, `EvalResidual`, `SolveLog`, `ClampedSolveLogT`, the four
`SolveLogTemperature`/`SolveTemperature` overloads, `SolveTemperatureFromP`,
`SolveDensity` and `EvalTMin`, then on the fourteen `EOS_Data` wrappers over them, and
then at every call site -- and the radiation-bearing entry points are called from
`wb_background.hpp` (19 sites), `red_giant.cpp` (17), `coordinates.cpp` (11),
`pgen_eos_utils.hpp` (10), `box_convection.cpp` (9), `solar_convection.cpp` (8),
`conduction.cpp` (7), the two `gnomonic_raisevel` headers (11), `ppm.hpp`/`wenoz.hpp`
(12), four Riemann solvers, `two_stream_rt.hpp`, `runaway_scan.hpp`, `hydro_newdt.cpp`,
`srcterms.cpp` -- about 150 sites over 28 files.  A gate that only some of them apply is
NOT a gate: the EOS would return one `p` and the reconstruction another, and the taper's
complementarity with the force (the thing that makes the total force `-grad P_rad` deep
and `kappa F/c` thin with no gap) would be broken in the middle of the window.  Optical
depth has the same problem and one more: the EOS is called before the sweep in every
stage, so a per-cell `tau` array would be a stage stale.

Temperature is the one stratification variable that is already in hand at EVERY one of
those sites, because it is an ARGUMENT of the EOS.  On this star it is also the better
variable: `T(r)` is strictly monotone over the whole column (checked, 15999 nodes), while
`rho(r)` is NOT -- the FeCZ density inversion is exactly what broke `w(rho)`.  And the
margin is larger where it matters: the FeCZ ambient sits **1.21x** above `eos_rad_rho_hi`
in density but **2.59x** above `eos_rad_t_hi` in temperature, so it takes a factor 2.6 in
`T` rather than 1.2 in `rho` to walk a convection-zone cell into the window.  The cell
that actually killed `Gnr` was at `T = 1.3-2.3e6 K`, twenty times `T_hi`: the gate pins
it at `w = 1`.

**The radius gate is still the better statement of the physics** and is worth doing when
someone is willing to thread a coordinate through the EOS interface; `tau` is better
still.  Neither is a one-round change, and neither is needed to test the diagnosis.

### 1.2 THE ic REBUILD

`make_ic_sph.py` now applies the same `max(w_rho, w_T)`, on the table `T` it solves for
(not the column's ideal-gas `T`), and prints the two thresholds it wants:

    tau = 3.00 at r = 2.338661e+11 (r/R 0.9861): rho = 1.9475e-09  T_col = 6.32457e+04  T_table = 6.32813e+04
    tau = 0.30 at r = 2.382576e+11 (r/R 1.0046): rho = 4.8608e-10  T_col = 4.51922e+04  T_table = 4.52131e+04

which is where `eos_rad_t_hi = 6.3281e4` and `eos_rad_t_lo = 4.5213e4` come from.  Because
the two windows are read off the SAME two optical depths, the gate is very nearly a NO-OP
on the initial column and a large change only on cells that leave it:

| check | result |
| --- | --- |
| `r`, `rho` columns of the ic | byte-identical to the pre-gate file |
| nodes where `w != w_rho` | **67 of 15999**, all between r/R 0.986 and 1.005 |
| `max |w - w_rho|` | **6e-4** |
| `max |eint_new/eint_old - 1|` | **1.8e-3**, at r/R 1.0042 (where `eint` is 5e3 erg/cm^3, 1e-8 of the base) |
| cycle-0 read-back of `p` | see gate 2's `column_x_gate.txt`: the run's own `p` at cycle 0 against the column's, max `|dp/p|` **4.0e-6**, well inside the 1e-4 the ungated file was held to and inside the 9.3e-7-to-1.8e-5 accuracy of the offline bilinear solve itself |

## 2. ALSO IN THIS COMMIT: `RegionIndcs` IS ZERO INITIALISED

`tests_r7` B.2a found that `restart.cpp` STEP 1 writes `mesh_indcs` and `mb_indcs` as raw
structs while their coarse-cell fields `cnx1..cke` are only ever assigned on a multilevel
mesh, so every restart file of a uniform grid carried 36 bytes of indeterminate stack.
`src/mesh/mesh.hpp:37-51` now gives every member a `= 0` initialiser.  Read straight out
of the files (the nine `int`s after `mb_indcs`):

    OLD: [538976288, 589307936, 541206560, 1918989427, 774971450, 540304688, ...]
    NEW: [0, 0, 0, 0, 0, 0, 0, 0, 0]

-- the old values are ASCII fragments of the input dump that happened to be on the stack.

---

## 3. GATE 1 -- THE PLANE-PARALLEL BOX REGRESSION: BITWISE

`g1.sh` (a copy of `tests_r4/g1.sh` pointed at this directory), both modes, serial CPU,
`athena_box_ref` = `65cd5bd1` unchanged against `athena_box_new` = this commit.  The gate
is OFF there (`eos_rad_t_hi` unset), which is what this checks.

| check | mode 3 (production) | mode 0 (explicit sweep) |
| --- | --- | --- |
| `feczrt.hydro.hst` | **IDENTICAL** | **IDENTICAL** |
| `feczrt.user.hst` | **IDENTICAL** | **IDENTICAL** |
| `column_used.txt` | **IDENTICAL** | **IDENTICAL** |
| `rt_surface.bin` | **IDENTICAL** | **IDENTICAL** |
| `rt_profile.bin` | **IDENTICAL** | **IDENTICAL** |
| `bin/ cbin/ rst/` | differ by **exactly 152 bytes** | same |

The 152 bytes are the two lines

    eos_rad_t_hi               = 0            # Default value added at run time
    eos_rad_t_lo               = 0            # Default value added at run time

that `GetOrAddReal` records in the parameter dump every output embeds, which shifts every
byte after them.  `diff` of the two dumps as text shows those two lines and nothing else.
Adding any parameter to this code does that; the five physics files are byte-identical.

## 4. GATE 2 -- THE 1-D COLUMN: THE KILLER IS GONE

`r8_1d.sh`.  The configuration is `tests_r7`'s `x_base` exactly -- arm `Gnr` with
`nx2 = nx3 = 4`, `vpert = 0`, one rank -- with the gated taper and
`problem/rt_rad_force = true`.  `x_gate` is that; `x_gate_noforce` is the same with the
force off, as the control.  Both were cut off by a 5 min 30 s wall, not by physics.

| arm | `dt COLLAPSE` | `eos_fail` | `vceil` | reached | turnovers |
| --- | --- | --- | --- | --- | --- |
| `tests_r7` `x_base` (ungated, force on) | **YES**, cell (5,5,6,36) r/R 0.791, `v1` = `vceil`, T at the table edge | 0 | fires inside the collapse | t = 5541 s | **1.18, DEAD** |
| `x_gate` (gated, force on) | **NONE** | **0** | 8432 (t/turn 1.5) | t = 10645 s | **2.26, alive at the wall** |
| `x_gate_noforce` (gated, force off) | **NONE** | **0** | **0** | t = 7606 s | **1.62, alive at the wall** |

**The single-cell runaway is removed.**  `tests_r7` established that `rt_rad_force = false`
was the ONLY switch of fifteen that removed it and that every other one moved the death by
0.03-0.22 turnover; with the gate, `rt_rad_force = true` runs past that death by a factor
1.9 in time with no collapse print, no `eos_fail` and a `dt` that falls smoothly rather
than by three decades in 90 s.  That is the diagnosis of `tests_r7` section 5 confirmed by
its own remedy: the killer was a cell being handed the optically-thin force because its
DENSITY had dropped, and closing that door closes the loop.

### 4.1 AND THE SECOND PROBLEM, WHICH IS NOT THE FORCE

`gate.py` / `an6.py` on the two arms.  The star does not settle, and **the control does the
same thing faster**, so what is left is not the radiative force:

| at 1.62 turnover | `x_gate` (force ON) | `x_gate_noforce` (force OFF) |
| --- | --- | --- |
| `d rho` at 0.97 R | +0.69 | **-0.87** |
| `d T` at 0.97 R | +0.50 | **+1.39** |
| `vr_rms/v_MLT` | 0.03 | **0.80** |
| `dt` [s] | 6.0 | 0.54 |
| `vceil` / `fofc` | 8432 / 892451 | **0 / 0** |

and `x_gate`'s own history, which is quiet for 1.75 turnovers and then turns over:

| t/turn | 0.25 | 1.00 | 1.50 | 1.75 | 2.00 | 2.25 |
| --- | --- | --- | --- | --- | --- | --- |
| `dt` [s] | 12.5 | 10.8 | 2.55 | 6.01 | 2.48 | 0.109 |
| `Gamma` at the kappa peak (0.76 R) | 1.156 | 1.166 | 1.146 | -- | -- | -- |
| `Gamma` at the base (0.63 R) | 0.982 | 0.979 | 0.986 | 0.998 | **0.791** | **0.477** |
| T at the base [K] | 1.875e5 | 1.873e5 | 1.867e5 | 1.856e5 | **2.038e5** | **2.481e5** |
| `Mdot(0.9 R)` [g/s] | +5.5e20 | +2.8e21 | +1.0e21 | +3.4e20 | **-5.0e21** | **-6.8e22** |
| domain mass drift | -2.5e-5 | +8.6e-5 | +4.9e-5 | -5.6e-4 | **-1.6e-2** | -1.6e-2 |
| `vr_rms/v_MLT` | 0.009 | 0.023 | 0.018 | 0.036 | 0.338 | **1.219** |

So from t = 0 to 1.75 turnovers the 1-D column is doing what a settling star should do --
`Gamma` at the kappa peak flat to 2 %, base T flat to 1 %, shell drifts of a few per cent,
`Mdot` outward and falling from its 1.0-turnover peak -- and then between 1.75 and 2.0
turnovers the whole column turns round: `Mdot` becomes strongly INWARD at every radius,
1.6 % of the domain mass leaves, the base heats 33 % and `Gamma` there falls to 0.48.
`vceil` and `fofc` fire hard at 1.5 turnovers, just before it.  **This is a new event, and
the force-off control shows the same 0.97 R shell being drained and heated (by 87 % and
+139 % at 1.6 turnovers) without any of those counters**, i.e. there are two distinct
mechanisms left and neither is `rt_rad_force`.  The obvious suspects, in order: the MLT
closure's deposition at the top of its own window, the bottom flux / wall, and the
two-stream blend at the 0.92-0.96 R shell that `tests_r6` already called "the swept
shell".  `mlt_alpha = 0` and `rt_bottom_flux = false` arms are the cheap way in and were
not run here.

**Gate 2 therefore passes on its stated test -- survive past 2.4 turnovers with the force
on -- and FAILS the 5-turnover settling test**, which it never reached: `dt` falls to
0.1 s, so 5 turnovers is ~1.3e5 cycles and about three GPU-hours for the 1-D column.
A 3 h `apu` slot was run for it (`r8_1d_long.sh`, arm `x_gate5`); see section 8.

### 4.2 THE ic READ-BACK, THROUGH THE CODE

`icchk.py <rt_profile.bin> <ic_file>` on `x_gate`'s cycle-0 record, 96 radial nodes:

| quantity | agreement |
| --- | --- |
| shell-mean `rho` vs the ic file | max `|d rho/rho|` **6.8e-5**, rms 1.3e-5 (linear interpolation of the ic onto the stretched grid) |
| shell-mean `T` vs the column's own `T` | max ratio **1.0003** over 0.50-1.01 R |

`red_giant.cpp:2219-2230` converts the file's `(rho, eint)` to `(p, T)` with the run's own
EOS and everything downstream reads that, so a `T` that reproduces the column to 3e-4 IS
the statement that the rebuilt `eint` and the gated EOS agree: feed the wrong `eint` to the
tapered EOS and `T` moves, which is the failure mode `make_ic_sph.py`'s header warns about.
**One thing did NOT check out and is left open**: `rt_profile.bin`'s slot 6 (`eint`) is 1.02x
the ic file's `eint` at 0.50 R and 80x at 1.01 R, with `rho` matching to 1e-5 and `T` to
3e-4 at the same nodes -- so the slot is not the quantity its comment says, or is scaled.
It is a DIAGNOSTIC DUMP and pre-existing (nothing in this commit touches it), but `an6.py`
reads it, so it needs settling before that reader's energy columns are trusted.

## 5. GATE 2b -- THE RESTART FILE IS NOW REPRODUCIBLE

Two independent 8-cycle runs of the same 1-D configuration, new binary:

    rst 00000: BYTE-IDENTICAL
    rst 00001: BYTE-IDENTICAL

and the 36 bytes that used to be stack garbage now read as nine zeros (section 2).

## 6. GATE 3 -- THE 3-D SMOKE ARM: IT PASSES THE WINDOW EVERY tests_r6 ARM DIED IN

`r8_3d_dev.sh`, arm `Gnr8dev`: the `tests_r6` arm `Gnr` configuration unchanged
(`inner_bc = wall`, `rt_bottom_flux = true`, `rad_flux_inner = 1.305278e15`,
`mlt_alpha = 1.5`, `vpert = 1e-3`, the input's own `dfloor_keep_velocity`, `vceil = 1e8`,
`efloor_as_tfloor`), with the gated taper, 2 ranks, ONE piece, no chaining.  It was cut by
a 12-minute `apudev` wall at **t = 6239 s = 1.326 turnovers**.

| gate | `tests_r6` `Gnr` (ungated) | `Gnr8dev` (gated) |
| --- | --- | --- |
| alive through 1.04-1.08 turnover | **NO -- died at 1.077** | **YES, no `dt COLLAPSE` print at all** |
| `eos_fail` | 0 | **0** |
| `vceil` count | 0 until the collapse | 146 (t/turn 1.25-1.33) |
| `fofc` count | -- | **4** |
| mass drift over the run | -- | **+9.5e-5** |
| `Gamma` at the kappa peak (0.761 R) | 1.154 -> 1.247 over one turnover | **1.154 -> 1.168**, flat to 1.2 % |
| T at the base (0.633 R) | -- | 1.871e5 -> 1.871e5 K, flat to 3e-4 |
| `Gamma` at the base | -- | 0.989 -> 0.982 |
| `Mdot(0.9 R)` | +5.4e20 -> +2.8e21 monotone up | +5.5e20 -> +2.8e21 (1.0 turn) -> **+1.5e21 (1.32), turning over** |
| FeCZ (0.68-0.92 R) drift | -9.7 %/turn | -9.7 (1.0) -> **-10.7 %/turn** |
| swept shell (0.92-0.96 R) | +13.8 %/turn | +13.5 (1.0) -> **-6.7 %/turn** |
| top atmosphere (0.96-1.01 R) | +5.4 %/turn | +5.8 (1.0) -> **+47.0 %/turn** |
| `vr_rms/v_MLT` | 0.034 at 1 turnover | 0.021 (1.0), 0.024 (1.32) |
| `F_res/F_req` at 0.70/0.80/0.90 R | 0.000 | **0.000** (only the t = 0 `mltfaces` dump exists in a one-piece run) |
| `dt` | 10.7 -> collapse | 10.7 (1.0) -> 5.2 (1.25) -> **0.27 (1.33)** |

**The 3-D result is the 1-D result: the killer is gone, and the second problem is at
0.97 R.**  Everything about the interior improves or holds -- the kappa-peak `Gamma` that
`tests_r6` watched climb 1.154 -> 1.247 in one turnover is now flat at 1.168, the base is
flat to 3e-4 in T, the mass is conserved to 1e-4, and `Mdot(0.9 R)` has stopped growing
for the first time in this campaign.  What takes the timestep at 1.32 turnovers is the
0.969 R shell: `rho` there goes +7 % (1.0 turnover) -> +33 % (1.25) -> **+96 % (1.32)**,
`v1` reaches 1.7e6 cm/s and its local `Gamma` crosses **1.001**.  Mass is piling into the
0.96-1.01 R shell at 47 %/turnover and leaving the swept shell below it.  That is a
LAYER, not a cell, it is in both 1-D arms with and without the force, and it sits just
BELOW the taper window (0.986-1.005 R) -- i.e. exactly where the two-stream, the MLT
closure and the taper all hand over to one another.

### 6.1 THE CHAINED ARMS, AND THE RESTART

`chain8.sh`, six legs of 0.5 turnover on `apu`, `athena_v4`: `Gnr8` = the same
configuration chained, `Gv2_8` = the same with `vpert = 1e-2`.  Both ran out of the 4-min
per-leg wall rather than reaching their `tlim`, because `dt` erodes (below), so they got
to 1.53 and 1.29 turnovers in six legs:

| | `Gnr8` (`vpert` 1e-3, chained) | `Gv2_8` (`vpert` 1e-2) |
| --- | --- | --- |
| reached | **1.530 turnovers** | **1.285 turnovers** |
| `dt COLLAPSE` | **none** | **none** |
| `eos_fail` | **0** | **0** |
| `vceil` / `fofc` | 5 / 551 | 2834 / 7646 |
| mass drift | **-1.9e-4** | **+8.9e-5** |
| `dt` at the end [s] | 0.67 | 0.31 |
| `d rho` at 0.97 R | +1.34 | +0.60 |
| `vr_rms/v_MLT` | 0.016 | 0.011 |

`Gamma_rad = Gamma x F_2s/F_req` from the per-leg `mltfaces` dumps, which is the gate
`tests_r6` section 4 set and could not meet:

| t/turn | 0.00 | 0.50 | 1.00 | 1.50 |
| --- | --- | --- | --- | --- |
| 0.70 R | 1.004 | 1.018 | 1.021 | 1.028 |
| **0.80 R** | 1.007 | 0.996 | 0.995 | **1.000** |
| 0.90 R | 1.036 | 1.083 | 1.086 | **1.020** |

**`Gamma_rad -> 1 within 0.02 at 0.8 R` now PASSES** (1.007 -> 1.000, inside 0.02 at every
dump) and 0.90 R, which `tests_r6` watched go the WRONG way (1.036 -> 1.084), now comes
back to 1.020.  `F_res/F_req` is still **0.000** at all three radii: the resolved flow is
not taking any of the flux, so the hand-over question still cannot be asked.

**The restart kick is gone in 3-D too**, which is `tests_r7` part B measured on the real
arm rather than on a 1-D chain.  `d ln KE1/dt` across the two leg boundaries that landed
exactly on a restart:

| boundary | `d lnKE1/dt` at it | its two neighbours |
| --- | --- | --- |
| t = 2352.5 s | 7.59e-4 | 8.32e-4, 7.73e-4 |
| t = 4705.0 s | 6.60e-4 | 6.57e-4, 6.76e-4 |

against `tests_r6`'s 2.5x jump in `1-KE` in the interval containing a restart.  The
horizontal kinetic energy is equally smooth.  A `mlt_alpha > 0` run can be chained.

## 8. THE 5-TURNOVER ARM IS STILL RUNNING

`r8_1d_long.sh` (arm `x_gate5`, 1-D, gated, force ON, a 3 h `apu` slot, `tlim` = 5
turnovers) was **at t = 11703 s = 2.487 turnovers with dt = 0.075 s and no collapse** when
this was written, with about 2.5 h of its slot left; at that dt it will reach roughly 4
turnovers and not 5.  So the 5-turnover settling gate is NOT answered here, and what
answers it is not more wall time -- it is section 4.1's 0.97 R layer, which is what took
dt from 12.9 to 0.08 in the first place.  Check `tests_r8/x_gate5/` for where it stopped.

## 7. THE VERDICT

1. **`tests_r7`'s diagnosis is confirmed by its own remedy.**  The 0.79 R killer was the
   optically-thin radiative force being handed to an interior cell because its DENSITY had
   dropped into the taper window.  Gating the taper weight on temperature removes the
   event in 1-D (no collapse to 2.26 turnovers against a death at 1.18) and in 3-D (past
   1.33 turnovers against 1.04-1.08 for four ungated arms), with `rt_rad_force = true`.
2. **The box is untouched**, bitwise, in both sweep modes; the gate is off unless
   `eos_rad_t_hi` is set.
3. **The smoke model is NOT settling yet, and it is no longer the force's fault.**  The
   interior is now quiet and conservative over 1.3 turnovers on the sphere -- kappa-peak
   `Gamma` flat to 1.2 %, base T flat to 3e-4, mass to 1e-4, `Mdot(0.9 R)` turning over --
   but a 0.97 R layer is piling up at 47 %/turnover, reaching `Gamma = 1.0` locally and
   taking `dt` to 0.27 s.  The 1-D control with `rt_rad_force = false` drains and heats the
   same shell (-87 % in `rho`, +139 % in T by 1.6 turnovers), so this is a THIRD mechanism
   and the next round's target.
4. **No production run is justified yet** and none has been started.  The next round is
   cheap and 1-D: `mlt_alpha = 0` and `rt_bottom_flux = false` against the gated baseline,
   to say whether the 0.97 R layer is the closure, the bottom flux or the two-stream
   blend; then the `rt_profile.bin` `eint` slot (section 4.2); then 5 turnovers.

## Files

* `README.md` -- this file.
* `g1.sh` -- gate 1; `g1_ref/ g1_new/ g1m0_ref/ g1m0_new/` (+ `g1_ref2/ g1_new2/`, the
  determinism pair) keep the `.hst`, the event log and `column_used.txt`.
* `r8_1d.sh` -- gate 2 and 2b; `x_gate/ x_gate_noforce/` keep the `.hst`, the event log,
  the column and `mltfaces`.  The determinism pair `rdet_d1/ rdet_d2/` was 58 MB of
  restart file and has been deleted; `r8_1d.log` carries its verdict.
  `r8_1d_long.sh` -- the 3 h `apu` version, arm `x_gate5/`.
* `r8_3d_dev.sh` -- gate 3, the `apudev` slot that decided it; `Gnr8dev/`.
  `chain8.sh` -- the chained `apu` version; `Gnr8/ Gv2_8/` and their `.log`s.
* `icchk.py` -- the cycle-0 ic read-back (section 4.2).  `tests_r6/gate.py` and
  `tests_r6/an6.py` are the gate tables, unchanged.
* `make_ic_sph.py` -- a COPY of the ic generator for the record; see `NOTE-ic.md`.
* `build_gpu.log`, `build_gpu2.log`, `g1_redo.log` -- the build logs.
* `athena_box_ref`, `athena_box_new`, `athena_v4` -- the binaries.  NOT committed;
  rebuild `athena_v4` from this commit in `build_gpu_rg`.
