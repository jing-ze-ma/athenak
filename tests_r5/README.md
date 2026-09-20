# he4_presn: WHERE THE MASS GOES, THE FLOOR MADE CONSISTENT, AND ARMS B''/D'/E''/F''

2026-09-17, viper, branch `he4-presn-global`, against `9e0b6f41` (= the `tests_r4`
commit).  Binaries: `build_gpu_rg` (MI300A, `PROBLEM=red_giant`, copied to
`tests_r4/athena_v2`) for the arms, `build_cpu_box` (serial CPU,
`PROBLEM=box_convection`) for the G1 bitwise gate.  Every `.bin`, `.rst` and `.cbin` has
been deleted after measuring and the four job logs truncated to the startup block, the
`dt COLLAPSE` blocks and the cycle/face-budget prints; every number below was read off
the full output first.

Read `tests_r4/README.md` (the fixed spherical two-stream, G1-G5) and
`tests_3d/arms/README.md` (arms A/B) first.  `an5.py <armdir>` builds the mass-budget and
Eddington tables, `gate.py <armdir>` the per-0.02/0.25-turnover gate tables (it is
`tests_3d/arms/gate.py` plus the two new event-log columns).

---

## 1. WHERE DOES THE MASS GO?  **VERDICT (a): RADIATIVELY DRIVEN INFLATION OF A
## SUPER-EDDINGTON LAYER.  Nothing escapes and nothing is numerical.**

### 1.1 The layer is super-Eddington on radiation alone, at t = 0, from the ic

`g(r) = GM/r^2` (the pgen prints this) and the star is a point mass, so the local
Eddington ratio collapses to a function of the opacity alone:

    Gamma(r) = kappa F_rad/(c g) = kappa L/(4 pi c G M) = 1.4642 kappa(r)

with `L = 2.3066e38`, `M = 6.2651e33`.  Read straight off the run's own initial column
(`Bpp/column_he4_Bpp.txt`, which carries `kappa_R`):

| r/R | kappa_R | **Gamma** | v_MLT |
| --- | --- | --- | --- |
| 0.60 | 0.665 | 0.973 | 0 |
| 0.635 (FeCZ base) | 0.677 | 0.991 | 1.8e6 |
| 0.70 | 0.738 | **1.080** | 1.28e7 |
| 0.744 | 0.772 | **1.131** | 1.41e7 |
| 0.76 (kappa peak) | 0.787 | **1.153** | 1.45e7 |
| 0.80 | 0.805 | **1.178** | 1.43e7 |
| 0.82 | 0.817 | **1.196** (max) | 1.40e7 |
| 0.90 | 0.788 | **1.154** | 1.06e7 |
| 0.95 | 0.638 | 0.934 | 2.6e6 |
| 0.97 (FeCZ top) | 0.565 | 0.828 | 0 |

**0.66 R to 0.94 R is unbound if the radiation carries all of L**, by up to 20 %.  With
`mlt_alpha = 0` and no resolved convection there is nothing else to carry it, so the layer
must inflate.  That is not a defect of the run; it is the statement that 15-25 % of `L`
*has* to be convective, which is what `tests_3d/arms/README.md` opened with.

### 1.2 The mass budget of B'' (= arm B' with the floor fix; identical to B')

`an5.py Bpp` -> `massBpp.txt`.  Shell mass flux `Mdot(r) = 4 pi r^2 <rho v1>` from
`rt_profile.bin` slot 2 (`<rho v1>`, area weighted), + = outward, [g/s]:

| t/turn | 0.60 R | 0.70 R | 0.80 R | 0.90 R | 0.97 R | top face (1.013 R) |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0 | 0 | +4e14 | -7e14 | +1e15 | 0 | +4.9e15 |
| 0.2 | -1.1e19 | +2.9e20 | +1.4e21 | +9.4e20 | -6.8e20 | -9.4e17 |
| 0.5 | -3.8e18 | +4.5e20 | **+3.0e21** | +5.7e21 | -1.3e21 | -1.3e19 |
| 0.8 | +6.7e17 | +2.2e20 | +2.0e21 | **+9.2e21** | -2.3e21 | -2.9e19 |
| 1.0 | -3.2e18 | -2.0e21 | +1.7e21 | +7.5e21 | **-2.6e21** | **-3.3e19** |

and the mass that actually crosses the domain faces, from the pgen's own `face budget`
print (`g in/out`, cumulative grams) and from the `hst` total:

* **top face: -3.3e19 g/s, i.e. 1.5e23 g over the turnover = 8e-4 of the domain.**
  Two to three ORDERS OF MAGNITUDE below the internal fluxes above.
* `hst` total mass **1.912482e26 -> 1.913006e26 g, drift +2.7e-4**.
  (The `face budget` `g in/out` numbers are 1e11 too large to be masses -- their unit
  conversion does not reproduce the `hst` drift and they should not be quoted as grams;
  the `hst` column and the profile flux above agree with each other and are what is used.)

**So the mass does not leave.  It moves.** Shell-integrated:

| zone | t=0 | 0.25 | 0.5 | 0.75 | 1.0 | drift |
| --- | --- | --- | --- | --- | --- | --- |
| 0.50-0.68 R, below the FeCZ | 5.055e25 | 5.049e25 | 5.038e25 | 5.028e25 | 5.016e25 | **-0.8 %** |
| 0.68-0.92 R, the Gamma>1 layer | 7.459e25 | 7.544e25 | 7.500e25 | 6.527e25 | 4.751e25 | **-36.3 %** |
| 0.92-0.96 R, the swept shell | 4.489e25 | 4.504e25 | 4.778e25 | 6.096e25 | 8.313e25 | **+85.2 %** |
| 0.96-1.01 R, the top atmosphere | 2.142e25 | 2.047e25 | 1.828e25 | 1.496e25 | 1.071e25 | **-50.0 %** |

The super-Eddington layer expands outward at `v1 = +2.6e6 .. +5.1e6` cm/s, piles the
overlying envelope into a dense shell at 0.92-0.96 R (the shell-mean `rho` at 0.951 R goes
8.5e-9 -> **2.7e-8**, x3.1), and the atmosphere above that shell falls inward onto it at
`v1 = -1.7e6`.  **The "-54 % at 0.97 R" of `tests_r4`'s G5 table is the layer boundary
moving outward past a FIXED diagnostic radius, plus the collapse of the atmosphere onto
the shell -- it is not a drain out of the star.**

### 1.3 The velocities are slow, subsonic and far below escape

At 1 turnover (B''), from `massBpp.txt` TABLE 4:

| r/R | v1 | c_s | v_esc | v1/c_s | v1/v_esc | v1/v_MLT |
| --- | --- | --- | --- | --- | --- | --- |
| 0.70 | -5.5e6 | 7.6e7 | 7.1e7 | -0.073 | -0.078 | -0.38 |
| 0.80 | +2.6e6 | 7.7e7 | 6.6e7 | +0.034 | +0.039 | +0.18 |
| 0.90 | +5.1e6 | 6.4e7 | 6.3e7 | +0.081 | +0.082 | +0.35 |
| 0.97 | -1.7e6 | 4.5e7 | 6.0e7 | -0.038 | -0.028 | -0.12 |

7 % of the sound speed and 8 % of escape.  **Nothing is escaping and nothing is
supersonic** -- this is a slow, subsonic readjustment, exactly what inflation at
`Gamma - 1 ~ 0.2` should look like.  Hypothesis (b) is a *consequence*, not the driver:
the top atmosphere does fall inward, but only above the expanding layer that took its
support away.  Hypothesis (c) is out: the sponge and the top face pass 1e-3 of the
internal flux and mass is conserved to 3e-4.

### 1.4 Does it settle?  **Yes, once the convective flux exists -- and that is the test.**

The Eddington ratio is self-regulating: as the layer drains, `kappa` falls with it.
`massBpp.txt` TABLE 3, the FeCZ body at 0.69-0.86 R, `Gamma` from the opacity table on
the shell-mean `(rho, T)`:

| r/R | Gamma(t=0) | Gamma(0.5) | Gamma(1.0) |
| --- | --- | --- | --- |
| 0.744 | 1.131 | 1.070 | 1.009 |
| 0.791 | 1.172 | 1.126 | 1.014 |
| 0.831 | 1.195 | 1.172 | 1.041 |
| 0.877 | 1.182 | 1.220 | 1.122 |
| 0.928 | 1.088 | 1.135 | **1.422** |
| 0.944 | 0.983 | 0.944 | **1.610** |

The body relaxes toward `Gamma = 1` -- but the swept shell becomes strongly
super-Eddington in its turn (1.6 at 0.944 R), because piling mass up at falling `T` puts
it back on the rising side of the Fe opacity bump.  **With `mlt_alpha = 0` the front keeps
propagating; it is a runaway, not a settled inflated state.**  The closure test is
section 3.

---

## 2. THE FLOOR: A RETRACTION AND A FIX

### 2.1 RETRACTION -- the "EOS-inversion garbage T = 5.2e11 K" is NOT garbage

`tests_r4/README.md` G5 and `tests_3d/arms/README.md` both read the `dt COLLAPSE` print's
`T` as kelvin and concluded "EOS-inversion garbage, the
`red-giant-seam-floor-eos-garbage` signature".  **That reading is wrong.**
`hydro_newdt.cpp:222` prints `wtemp`, which under a general EOS is `p/rho` at
`mu_ref = 1`, i.e. `8.314e7 x T[K]`.  So

    T_print = 5.24607e11  ==  5.24607e11/8.314e7  =  6310 K  =  10^3.8 K

which is **exactly `eos_logt_min`, the lowest row of the tabulated EOS**, and it repeats
to six digits at every collapse because it is a CLAMP (`EOSTable::ClampLogT`), not a root
find.  And the state is thermodynamically consistent at that temperature:
`p = rho R T/mu` with `rho = 3.52449e-9`, `T = 6310 K` gives `mu = 4.03` -- **neutral
helium**, `p = 459` against the printed 459.013.  The same check on arm B's cell
(`rho = 5.37e-10`, `p = 69.43`) gives `mu = 4.04`.  Nothing failed to invert;
`eos_fail = 0` and `max c2p it = 0` in every arm, as the tables always said.

**What IS wrong** is that the clamp leaves the cell inconsistent in the other direction:
`T`, `p`, `Gamma_1` and the sound speed are the clamp row's, while the conserved energy
is still that of a gas colder than the table can represent.  The Riemann solver then sees
a pressure the cell's own energy does not support, the timestep sees the clamp row's
sound speed, and the two-stream reads a Planck function the energy cannot pay for.

**And the temperature floor could not fire.**  `tfloor = 5.0e3` K is BELOW
`10^3.8 = 6310` K, so every state that would have tripped it is clamped to 6310 K first:
`tfloor = 0` in every arm's event log, for the whole run, by construction.  The comment in
`gnomonic_raisevel.hpp` ("it is the TEMPERATURE floor that actually keeps the lookup in
range") was describing a floor that was never reachable.

### 2.2 THE FIX -- `<block>/efloor_as_tfloor`, default OFF

New switch, added to the `floors_legacy` list, meaningful only for the tabulated EOS:

* **`src/eos/eos.cpp`, `BuildGeneralEOS`** -- with the switch on, `tfloor` is RAISED to
  the table's lowest tabulated temperature when it is set below it, and the run says so:
  `General EOS: <hydro>/efloor_as_tfloor raised tfloor 6.01e-05 K -> 6309.57 K, the
  table's lowest tabulated temperature; a floor below it is unreachable`.
  (6.01e-05 is `tfloor = 5e3` K in code units.)  This also repairs the two-stream: its
  own `rt_use_cons` guard clamps a non-positive internal energy to
  `e(rho, eos.tfloor)` (`two_stream_rt.hpp:1614`), which until now asked the table for a
  state three decades below its own edge; every other temperature the two-stream forms
  goes through `eos.Temperature`, which `ClampLogT` already holds in range.
* **`src/eos/general_c2p_hyd.hpp`, `SingleC2P_GeneralHyd`** and
  **`src/coordinates/gnomonic_raisevel.hpp`, `GnomonicRaiseVelFloors`** (the latter is
  where the floors are actually applied on the cubed sphere -- `defer_cons_floors`) --
  when a floor fired OR the inversion was clamped, `e` is set to `e(rho, T)` for the
  temperature that was settled on (`max(temp, tfloor)`, and the clamp row where the clamp
  bound) and the conserved energy is corrected to match, so **re-inverting the cell
  returns that same T by construction**.  At the low edge this raises `e` -- the
  temperature floor doing what the pressure floor could not reach; at the high edge it
  lowers it, removing energy the table cannot represent instead of leaving a saturated
  state.  `RaiseVelFloors` now also calls the clamp-REPORTING overload of
  `TemperaturePressureGamma1` (the six-argument form wraps it, so the arithmetic is
  unchanged).
* **`floored` is NOT set by the repair.**  FOFC flags cells on `floored`, and the table
  clamp is not a floor -- it fires wherever the trial state leaves the table, which in
  the draining top atmosphere is a fifth of the grid.  The first build that folded it in
  put **5.3e7 FOFC firings** into arm B''; with `floored = floor_fired` FOFC is back to
  **0**, as in B'.  That one line is the whole difference between the two B'' runs.
* **A new event counter, `EventCounters::neos_tset`** -- "floor set by T" -- APPENDED as
  the last column of the event log (`eos_tset`), so every column index an existing reader
  uses is unchanged.  Plumbed through `general_hyd_floors.cpp` (7th reducer) and
  `coordinates.cpp` (3rd reducer of the RaiseVel kernel); `hydro_fofc.cpp` and
  `prolong_prims.cpp` take the extra argument and ignore it.
* `inputs/hydro/he4_presn_cs.athinput` sets `efloor_as_tfloor = true` and carries
  `efloor_from_ekin = false` so the latter can be overridden from the command line.

### 2.3 THE BITWISE CHECK -- G1, both modes, IDENTICAL

`g1.sh` (= `tests_r4/g1.sh` with the output directory moved), the box_w8 production
mode-3 configuration at `nx2 = nx3 = 16`, 50 cycles, serial CPU, and the mode-0 variant.
The box input sets no floor switch, so it stays on the `floors_legacy` path and every
change above is in code it does not execute.

    cmp tests_r4/g1_new/feczrt.hydro.hst  g1_r5/...     IDENTICAL
    cmp tests_r4/g1_new/feczrt.user.hst   g1_r5/...     IDENTICAL
    cmp tests_r4/g1_new/column_used.txt   g1_r5/...     IDENTICAL
    cmp tests_r4/g1_new/rt_surface.bin    g1_r5/...     IDENTICAL
    cmp tests_r4/g1_new/rt_profile.bin    g1_r5/...     IDENTICAL

and the same five lines for `g1m0_new` vs `g1m0_r5` (`rt_implicit_column = 0`).  The only
difference anywhere is the event log's HEADER, which gained ` eos_tset`; there are no
data rows, because no floor fires in the box.

---

## 3. THE ARMS

`chain.sh` (`--export=ALL,ARM=,WALL=,TLIM=,VPERT=,MLT=,EXTRA=`), the `tests_r4/g5` grid
and cadences: `nx1 = 96` stretched from 0.50 R, `nx2 = nx3 = 32` per panel,
`meshblock 96x16x16` = 24 blocks, 2 MI300A, turnover = 4705 s, `apudev`.  All four arms
are the arm-B boundary (`inner_bc = wall`, `rt_bottom_flux = true`,
`rad_flux_inner = 1.305278e15`) with the section-2 fix ON.  A turnover costs **~80 s** of
wall time, which is why four arms fit in one sitting.

| arm | `mlt_alpha` | `vpert` | extra | lifetime | `dt` at the end | **verdict** |
| --- | --- | --- | --- | --- | --- | --- |
| B' (`tests_r4/g5/Bp`, no fix) | 0 | 1e-3 | -- | 1.000 (tlim) | 4.09 s | drains |
| **B''** | 0 | 1e-3 | -- | 1.000 (tlim) | 4.09 s | **identical to B'** |
| **D'** | 1.5 | 1e-3 | -- | 1.000 (tlim) | **10.7 s** | **drain CLOSED** |
| E' (`tests_r4/g5/Ep`, no fix) | 0 | 1e-2 | -- | **0.782** | died | dies |
| **E''** | 1.5 | 1e-2 | -- | **1.034** | died | dies at the clamp |
| **F''** | 1.5 | 1e-2 | `efloor_from_ekin` | **1.039** | died | **null** |

Full tables: `tableBpp.txt`, `tableDp.txt`, `tableEpp.txt`, `tableFpp.txt`;
mass budgets `massBpp.txt`, `massDp.txt`, `massEpp.txt`.

### 3.1 B'' -- the floor fix ALONE changes nothing about the star

Gate columns at the quarter turnovers (`tableBpp.txt`):

| t/turn | dt[s] | eos_fail | floors | fofc | tclamp | tset | `L_out/L` | d_rho base | peak | top | `vr/vMLT` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 12.9 | 0 | 0 | 0 | 0 | 0 | 0.996 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.26 | 12.6 | 0 | 1.0e6 | 0 | 0 | 3.7e5 | 0.983 | 0.000 | -0.039 | -0.050 | 0.102 |
| 0.50 | 11.2 | 0 | 3.3e6 | 0 | 616 | 2.2e6 | 0.988 | -0.001 | -0.082 | -0.158 | 0.099 |
| 0.76 | 9.15 | 0 | 3.6e6 | 0 | 1716 | 2.5e6 | 0.970 | -0.002 | -0.164 | -0.323 | 0.188 |
| 1.00 | 4.09 | 0 | 3.6e6 | 0 | 2467 | 2.5e6 | **0.917** | -0.003 | **-0.292** | **-0.544** | 0.123 |

`eos_fail = 0`, `fofc = 0`, `tfloor = 0`, `vceil = 0`; `dfloor = 2.29e6`,
`efloor = 1.29e6`, `tclamp = 2467`, `tset = 2.50e6` cell-cycles.  **Every structural
number is B' to three digits** (`L_out/L` 0.917, `d_rho` -0.003/-0.292/-0.544, and the
mass tables of `massBpp.txt` and `tests_r4`'s B' agree to the last printed digit).  So
the fix is a genuine null on the mean structure, which is what a floor repair should be:
it removes an inconsistency, it does not remove the physics.  **The drain is physical.**

### 3.2 D' -- `mlt_alpha = 1.5` CLOSES THE DRAIN

| t/turn | dt[s] | eos_fail | floors | fofc | tclamp | tset | `L_out/L` | d_rho base | peak | top | `vr/vMLT` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 12.9 | 0 | 0 | 0 | 0 | 0 | 0.999 | 0.000 | 0.000 | 0.000 | 0.000 |
| 0.26 | 12.5 | 0 | 9.8e5 | 0 | 0 | 3.0e5 | 0.961 | -0.008 | +0.029 | +0.007 | 0.010 |
| 0.50 | 12.4 | 0 | 3.3e6 | 0 | 643 | 2.1e6 | 0.959 | -0.017 | +0.010 | -0.035 | 0.006 |
| 0.76 | 11.7 | 0 | 3.4e6 | 0 | 2409 | 2.1e6 | 0.951 | -0.025 | +0.169 | +0.131 | 0.019 |
| 1.00 | **10.7** | 0 | 3.4e6 | 0 | 4268 | 2.1e6 | **0.944** | -0.029 | **+0.390** | **+0.091** | 0.033 |

`eos_fail = 0`, `fofc = 0`, `tfloor = 0`, `vceil = 0`, and no `dt COLLAPSE` anywhere.
Against B'' at the same time:

| | B'' (mlt 0) | **D' (mlt 1.5)** |
| --- | --- | --- |
| `dt` at 1 turnover | 4.09 s | **10.7 s** (2.6x) |
| `L_out/L` range | 0.917 .. 0.996 | **0.944 .. 0.999** |
| `d_rho` FeCZ top | **-0.544** | **+0.091** |
| `d_rho` kappa peak | -0.292 | +0.390 |
| FeCZ mass (0.68-0.92 R) | **-36.3 %** | **-9.3 %** |
| swept shell (0.92-0.96 R) | **+85.2 %** | **+13.6 %** |
| top atmosphere (0.96-1.01 R) | **-50.0 %** | **+5.3 %** |
| `Mdot` at 0.90 R, 1 turnover | +7.5e21 g/s | **+2.8e21 g/s** |
| `Mdot` at 0.97 R, 1 turnover | -2.6e21 g/s (inward) | **+1.5e21 g/s (outward)** |
| `lnKE1` at 1 turnover | 89.78 | **86.91** (KE_radial /18) |

**and the reason is exactly the Eddington ratio.**  `Dp/mltfaces.txt` at `t = 0`:

| r/R | `F_mlt/F_req` | `F_2s/F_req` | Gamma_total | **Gamma_rad = Gamma x F_2s/F_req** |
| --- | --- | --- | --- | --- |
| 0.70 | 0.102 | 0.921 | 1.080 | **0.995** |
| 0.744 | 0.151 | 0.881 | 1.131 | **0.996** |
| 0.80 | 0.194 | 0.852 | 1.178 | **1.004** |
| 0.85 | 0.222 | 0.848 | 1.193 | **1.012** |
| 0.90 | 0.232 | 0.899 | 1.167 | **1.049** |

The MLT closure takes 10-23 % of `L` off the radiation and brings the radiative Eddington
ratio from **1.08-1.20 down to 1.00-1.05**.  The layer stops being unbound, the collision
at 0.92-0.96 R stops (the flux at 0.97 R reverses from inward to outward), and the drift
becomes a slow, coherent, monotone expansion at 1/4 the rate -- `d_rho` at the peak
+39 % and at the top +9 % after a turnover, with `Gamma_rad` still 1-5 % above 1, which
is exactly the residual that is still driving it.  **This is a settling inflated state,
not a runaway** -- but one turnover does not prove it settles, because the residual
`Gamma_rad - 1 ~ 0.03` still has a 30-turnover timescale in it.

### 3.3 E'' and F'' -- the seed still kills it, and the death is now a VACUUM HOLE

E'' = D' with `vpert = 1e-2`, 3 turnovers requested.  The mean structure tracks D' face
for face (`L_out/L` 0.945-0.999, `d_rho` -0.029/+0.399/+0.081 at 1 turnover, the same
mass table to 0.4 %) -- only `lnKEh` differs, by the 100x the seed puts in.  It dies at
**t = 4867 s = 1.034 turnover** (E' died at 0.782, so 1.32x), on

    cell (4,16,18,33) r = 1.8226e11 (r/R = 0.768)
    rho = 4.389e-9  T_print = 5.246071e11 (= 6310 K, the clamp)  p = 570.6
    v = (-1.34e9, 1.7e4, -9.8e4)      then  v1 = +2.08e10 one cycle later

`p/(rho T_print) = 0.2475`, i.e. `mu = 4.04`: the SAME consistent neutral-helium clamp
state as before, only now it is consistent BY CONSTRUCTION rather than by accident.  So
the fix removed the inconsistency and the star still dies -- because the cell is not
being mis-inverted, it is being **emptied**: `|v1|` is already 1e9-2e10 cm/s (60-1400
v_MLT) before the floor is anywhere near, so `e = E - e_kin` goes to nothing and the cell
lands on the floor as a *consequence* of a velocity runaway, not as its cause.

F'' = E'' + `efloor_from_ekin = true`, which holds the conserved total fixed and rescales
the momentum instead of donating the cell's whole kinetic energy back into its own
acceleration.  **Null: 1.039 turnover against 1.034.**  Its collapsing cell says why:

    cell (5,17,3,37) r = 1.8936e11 (r/R = 0.798)
    rho = 1.000000e-13  ( = dfloor EXACTLY )   p = 1.000000e-01 ( = pfloor EXACTLY )
    v = (6.37e12, -6.3e9, +6.5e9)     T_print = 1.857e12 (= 22335 K, on the table)

**The cell has been evacuated by 4.4 orders of magnitude** (local `rho` is 2.3e-9) and the
velocity is `m/dfloor` = 6.4e12 cm/s = 212 c.  `dfloor = 1e-13` is 5.7e-3 of the initial
TOP cell but 4.3e-5 of the FeCZ, so inside the layer a cell can be emptied by four
decades before any floor notices, and when the density floor finally fires the default
(`dfloor_keep_velocity = false`) raises `rho` and leaves the momentum alone -- which is
what turns an empty cell into a 212 c one.

---

## 4. RECOMMENDED NEXT STEP

**Run D' (not E'') for 3-5 turnovers, with a density floor and a velocity ceiling that
are sized to the FeCZ rather than to the top cell.**  Three things, in this order:

1. **`dfloor_keep_velocity = true` and a radius-appropriate `dfloor`/`vceil`.**  This is
   the only thing F'' identified and it is a two-line input change:
   `dfloor_keep_velocity = true` (scale the momentum with the density, so a floored cell
   keeps its velocity instead of acquiring `m/dfloor`), and `vceil ~ 1e8` cm/s = 1.4 c_s
   = 7 v_MLT, which is above anything physical in this star and 60x below the 6.4e12 the
   run reached.  Both are already implemented for cubed-sphere hydro and both are
   deferred to `GnomonicEquiangleRaiseVel`, so nothing new is needed.  Raising `dfloor`
   itself is the cruder alternative and would change the top atmosphere.
2. **Then D' for 3-5 turnovers with `vpert = 1e-3`, to decide whether the inflated state
   SETTLES.**  The gate is `Gamma_rad -> 1` and `Mdot(0.90 R) -> 0`: `Mdot` there already
   peaked at 2.9e21 g/s near 0.9 turnover and the `Gamma_rad` residual is only 1-5 %, so
   a settled state is the expected outcome and 5 turnovers is 7 minutes of `apudev`.
   `vpert = 1e-2` should NOT be used again until (1) is in: the seed does not change the
   mean structure at all (E''/F'' = D' face for face) and it buys nothing but the
   velocity runaway.
3. **Only then the seed and the resolution.**  `mlt_alpha = 1.5` is a shell-mean closure
   standing in for convection that is not resolved; the real test of this star is whether
   resolved convection can take over the 10-23 % of `L` the closure is carrying, and
   `vr_rms/v_MLT` is 0.03 in D' after one turnover.  That needs turnovers, not cells.

## Files

* `README.md` -- this file.
* `an5.py <armdir> [label]` -- the mass-budget, Eddington-ratio and velocity tables.
* `gate.py <armdir> [label]` -- `tests_3d/arms/gate.py` plus the `tclamp`/`tset` columns.
* `chain.sh`, `g1.sh` -- the submit script and the G1 driver.
* `g1_r5/`, `g1m0_r5/` -- G1, the `.hst`, `column_used.txt`, `rt_{surface,profile}.bin`
  and event log kept; `bin/ rst/ cbin/` deleted.
* `Bpp/ Dp/ Epp/ Fpp/` -- `he4.hydro.hst`, `he4.log`, `rt_profile.bin`, `rt_surface.bin`,
  `column_he4_*.txt`, and `mltfaces.txt` for the `mlt_alpha > 0` arms.
* `table{Bpp,Dp,Epp,Fpp}.txt`, `mass{Bpp,Dp,Epp}.txt` -- the tables as generated.
* `{Bpp,Dp,Epp,Fpp}.*.log` -- the job logs, truncated.
