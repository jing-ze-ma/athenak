# he4_presn: THE FLOORS WORK, THE STAR STILL DIES AT 1.07 TURNOVER -- AND CHAINED
# RESTARTS KICK IT

2026-09-17, viper, branch `he4-presn-global`, against `e247210d` (= the `tests_r5`
commit).  Binary `tests_r4/athena_v2` (= `build_gpu_rg`, MI300A, `PROBLEM=red_giant`),
unchanged: **this round changes only the input file.**  Grid, cadences and boundary
configuration are `tests_r5`'s arm D': `nx1 = 96` stretched from 0.50 R, `nx2 = nx3 = 32`
per panel, `meshblock 96x16x16` = 24 blocks, 2 MI300A on `apudev`, `inner_bc = wall`,
`rt_bottom_flux = true`, `rad_flux_inner = 1.305278e15`, `mlt_alpha = 1.5`,
`efloor_as_tfloor = true`.  Turnover = 4705 s.  Every `.bin`, `.rst` and `.cbin` has been
deleted after measuring; `hst`, `rt_profile.bin`, `rt_surface.bin`, the event logs, the
`column_*` dumps and the `mltfaces_*` dumps are kept.

Read `tests_r5/README.md` first (the inflation diagnosis and arms B''/D'/E''/F''); this
round executes its section-4 recommendation.  `gate.py <armdir>` and `an6.py <armdir>` are
the analysers (`gate.py` = `tests_r5/gate.py` plus an `eos_vceil` column; `an6.py` =
`tests_r5/an5.py` rescaled to a multi-turnover run, with a new TABLE 4 that reads the
hand-over budget out of the `mltfaces_*` dumps).  Tables as generated:
`table_<arm>.txt`, `mass_<arm>.txt`.

---

## 0. THE INPUT CHANGE

`inputs/hydro/he4_presn_cs.athinput`, `<hydro>`:

    dfloor_keep_velocity = true      # scale the momentum with the density
    vceil                = 1.0e8     # cm/s = 1.4 c_s at the FeCZ = 7 v_MLT

Both are read in the `<hydro>` block (`src/eos/eos.cpp:39,66`), both are refused outside
hydrodynamics, and on the cubed sphere both are **deferred** to
`Coordinates::GnomonicEquiangleRaiseVel`, which owns the metric `|v|` and kinetic energy
(`src/coordinates/gnomonic_raisevel.hpp:75-99`; the `ConsToPrim` block only DECIDES,
`general_c2p_hyd.hpp:186`), and which increments `EventCounters::neos_vceil`
(`coordinates.cpp:977`).  Verified honoured on this path: the effective parameter list
dumped into the restart file carries `dfloor_keep_velocity = true` and `vceil = 1.0e8`,
and the counter fires (section 2).  cgs throughout (`<units>` = 1,1,1), so `vceil` is
cm/s.  The reasons are documented in the input file itself.

## 1. THE ARMS

| arm | floors | `vpert` | restarts | lifetime | died on |
| --- | --- | --- | --- | --- | --- |
| **G** | on | 1e-3 | every 0.5 turnover | **0.63** | the restart kick (section 3) |
| **H** | on | 1e-2 | every 0.5 turnover | **0.59** | the restart kick |
| **Gnr** | on (`vceil` 1e8) | 1e-3 | none | **1.077** | cell 0.776 R, `v1` at `vceil` |
| **Gv5** | on (`vceil` 5e7) | 1e-3 | none | **1.078** | cell 0.791 R, `v1 = vceil` exactly |
| **Hnr** | on | 1e-2 | none | **1.037** | cell 0.791 R, `v1 = vceil` exactly |
| **Dnr** | **off** (= D' continued) | 1e-3 | none | **1.075** | cell 0.791 R, `v1 = 5.2 c` |

`Gnr` is arm G as it was meant to be run, `Dnr` the control that tells the floors from the
star, `Gv5` the tighter ceiling, `Hnr` the seed test.  **All four no-restart arms are the
same run to three digits through 1.00 turnover** -- and all four die within
1.037-1.078 turnover.

## 2. THE FLOORS: A NULL ON THE STRUCTURE, A WIN ON THE PATHOLOGY, NO WIN ON LIFETIME

### 2.1 Null on the mean structure

`table_{Dnr,Gnr,Gv5}.txt`, identical to each other AND to `tests_r5`'s D' to three digits:

| t/turn | dt[s] | `L_out/L` | d_rho base | peak | top | `vr_rms/vMLT` | `eos_vceil` | fofc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.00 | 12.9 | 0.999 | 0.000 | 0.000 | 0.000 | 0.000 | 0 | 0 |
| 0.25 | 12.5 | 0.968 | -0.008 | +0.025 | +0.008 | 0.009 | 0 | 0 |
| 0.50 | 12.4 | 0.959 | -0.017 | +0.010 | -0.035 | 0.006 | 0 | 0 |
| 0.75 | 11.8 | 0.951 | -0.024 | +0.159 | +0.110 | 0.018 | 0 | 0 |
| 1.00 | 10.7 | 0.944 | -0.029 | +0.415 | +0.081 | 0.034 | **0** | **0** |

(D' at 1.00 turnover: `dt` 10.7, `L_out/L` 0.944, `d_rho` -0.029/+0.390/+0.091.)
`eos_fail = 0`, `tfloor = 0`, mass drift +2.9e-4, **and the velocity ceiling does not fire
once in the healthy phase** -- `eos_vceil = 0` through 1.00 turnover in every arm, in
`Gv5` too.  So the ceiling is NOT doing physics: it is catching a numerical runaway, which
is what it was sized for.

### 2.2 Win on the pathology: the superluminal evacuated cell is gone

The collapsing cell, from each arm's own `dt COLLAPSE` print (code `T` / 8.314e7 = T[K]):

| arm | t[s] | t/turn | cell (m,k,j,i) | r/R | rho | `v1` [cm/s] | T[K] | `c_s` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **Dnr** (no floors) | 5056 | 1.075 | (11,18,17,36) | 0.791 | 4.9e-9 | **+1.57e11 = 5.2 c** | 6310 (clamp) | 4.6e5 |
| **Gnr** (`vceil` 1e8) | 5067 | 1.077 | (11,18,18,34) | 0.776 | 5.1e-10 | -5.78e7 | 2.3e6 | 3.3e9 |
| **Gv5** (`vceil` 5e7) | 5074 | 1.078 | (6,18,6,36) | 0.791 | 5.1e-9 | **+5.000000e7 = `vceil`** | 1.3e6 | 1.4e9 |
| **Hnr** (`vpert` 1e-2) | 4879 | 1.037 | (7,3,4,36) | 0.791 | 3.6e-9 | **+1.000000e8 = `vceil`** | 1.9e6 | 3.6e9 |

`tests_r5`'s F'' reached 6.4e12 cm/s = 212 c on an evacuated `rho = dfloor` cell; `Dnr`
reaches 1.57e11 = 5.2 c on the same clamp state.  **With the floors on, `v1` sits exactly
on `vceil` and `rho` is the local density, not `dfloor`** -- no evacuation, no
superluminal cell, and `Gnr`'s later top-cell collapse (i = 92) likewise sits at
-9.9999e7.  That half of the `tests_r5` recommendation is delivered.

### 2.3 No win on lifetime, because the runaway is THERMAL, not kinetic

`Dnr` dies with `c_s = 4.6e5` (a cold clamped cell, `dt` set by `|v|`); every floored arm
dies with `c_s = 1.4e9 .. 3.6e9`, i.e. **`dt` is now set by the sound speed of a cell
heated to 1.3-2.3e6 K** at r/R = 0.78-0.79, where the ambient is 1.7e5 K.  Capping the
velocity removes the velocity from the collapse and leaves the heating, so the lifetime
moves by 0.003 turnover (1.075 -> 1.077/1.078).  A tighter ceiling (`Gv5`, 5e7 = 3.4
v_MLT) buys nothing either.

**And the seed is now provably irrelevant**: `Hnr` (`vpert` 1e-2) dies at 1.037, which is
`tests_r5`'s E''/F'' (1.034/1.039) to the third digit, on the SAME cell index i = 36 and
the same radius 0.791 R as the 1e-3 arms.  The "velocity runaway of the 1e-2 seed" of
`tests_r5` section 3.3 was never about the seed: **there is a deterministic event at
r/R = 0.79 at t = 4.9-5.1e3 s in this configuration**, and every arm hits it.

## 3. THE NEW FINDING: A CHAINED RESTART KICKS THE STAR (mlt_alpha > 0)

Arms G and H were run as 0.5-turnover legs, because `problem/mlt_dump` is one-shot per
process and a leg boundary is the only way to get a `Gamma_rad` time series.  Both died at
0.6 turnover -- 0.45 turnover EARLIER than the same configuration run in one piece.  The
`hst` across the single restart at t = 2352.5 s:

| arm | t[s] | 2307.8 | 2352.5 | **2402.0** | 2451.0 | 2497.0 | 2541.0 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| G | `lnKE1` | 85.223 | 85.255 | **86.162** | 86.317 | 86.650 | 86.817 |
| G | `lnKEh` | 74.899 | 74.917 | 74.938 | 74.960 | 74.982 | 75.004 |
| G | `dt`[s] | 12.44 | 12.42 | 12.34 | 12.17 | 11.18 | 11.12 |
| H | `lnKE1` | 85.195 | 85.229 | **86.761** | 86.891 | 87.123 | 87.288 |
| H | `lnKEh` | 79.506 | 79.525 | 79.546 | 79.567 | 79.591 | 79.610 |

**The RADIAL kinetic energy jumps by a factor 2.5 (G) / 4.6 (H) in the one `hst` interval
that contains the restart, while the horizontal kinetic energy is perfectly continuous**
and `Gnr` at the same times is smooth (`lnKE1` 85.26 at 2357 -> 86.11 at 3533).  So the
restart injects a radial impulse.

**The mechanism is the MLT closure's relaxed profile.**  `problem/mlt_relax_time = 1.0e4`
s (2 turnovers) makes the applied subgrid flux follow the mean of the deficit rather than
its flicker: `fmlt1d` is relaxed toward the instantaneous target by `dt/relax_time` per
call (`src/pgen/red_giant.cpp:3400`), and it is **pgen state, not restart state** -- the
first source call of ANY process has `seed = !mlt_relax_seeded_` true
(`red_giant.cpp:3262,3409`) and sets `fmlt1d = ftar` outright.  The size of that reset is
in the dumps: `mltfaces_00.txt` (t = 0) against `mltfaces_01.txt` (written by the first
call after the restart at 0.5 turnover), `F_mlt/F_req`:

| r/R | t = 0 | after the restart |
| --- | --- | --- |
| 0.70 | 0.102 | **0.711** |
| 0.80 | 0.194 | **0.735** |
| 0.90 | 0.232 | **0.353** |

The applied subgrid convective flux jumps by a factor 1.5-7 in one step.  That is the
kick.  **Consequence for production: a `mlt_alpha > 0` run CANNOT be chained across
24-hour slots as it stands** -- every restart re-seeds the closure and hits the star with
a radial impulse.  Either `fmlt1d` has to go into the restart file (and
`mlt_relax_seeded_` with it), or `mlt_relax_time` has to be 0, or the legs have to be long
compared with the relaxation time.  Nothing else in this round depends on it: all the
physics numbers above are from single-piece runs.

## 4. THE GATES

| gate | result |
| --- | --- |
| alive to 5 turnovers | **FAIL** -- 1.077 (Gnr), 1.037 (Hnr) |
| `Gamma_rad -> 1` within 0.02 at 0.8 R | **not reached**: 1.007 at t = 0 (inside 0.02), 0.996 at 0.5 turnover -- but 0.90 R moves the WRONG way, 1.036 -> 1.084, and the total `Gamma` at the kappa peak rises 1.154 -> 1.247 over one turnover |
| `Mdot(0.9 R)` decreasing toward 0 | **FAIL**: +5.4e20 (0.25) -> +1.50e21 (0.5) -> +2.56e21 (0.75) -> **+2.83e21 g/s (1.00, the peak and the end value)**, monotone increasing |
| shell rho drift per turnover falling | **not settling**: FeCZ (0.68-0.92 R) -9.7 %/turnover, swept shell (0.92-0.96 R) +13.8 %/turnover, top atmosphere +5.4 %/turnover, all still growing when the run dies |
| `vceil` count | **0 through 1.00 turnover** in every arm (Gnr, Gv5, Hnr); it fires only inside the collapse.  The ceiling is not doing physics. |
| hand-over: resolved convection taking over? | **NO.** `F_res/F_req = 0.000` at 0.70/0.80/0.90 R at every dump; `F_mlt/F_req` = 0.10/0.19/0.23 at t = 0; `vr_rms/v_MLT` = 0.034 at 1 turnover (`mass_Gnr.txt` TABLE 3/4) |

So the honest reading of the whole round: **the floors are exactly the null-plus-safety-net
they were designed to be, and they are the right defaults to keep; but they do not make
this smoke model production-ready.**  The star does not settle -- `Gamma` at the kappa peak
RISES as the layer piles mass onto the rising side of the Fe bump, `Mdot(0.9 R)` grows
monotonically, and at t ~ 5.0e3 s a cell at r/R = 0.79 heats to 1-2e6 K and takes the
timestep with it.  `tests_r5` section 1.4 predicted exactly this ("with `mlt_alpha = 0`
the front keeps propagating"); what is new is that **`mlt_alpha = 1.5` only slows it**, and
that the end state is a local thermal runaway rather than a numerical velocity blow-up.

## 5. WHAT REMAINS BEFORE A PRODUCTION GRID CAN BE JUSTIFIED

1. **Diagnose the 0.79 R thermal runaway at t ~ 5.0e3 s.**  It is deterministic (four arms,
   two seeds, two ceilings, floors on and off, all within 0.04 turnover of each other, same
   cell index) and it is a HEATING event: `T` 1.7e5 -> 1.3-2.3e6 K in one cell.  The
   candidates are the subgrid flux deposition itself (`F_mlt` is applied as a divergence at
   that face), the two-stream/conduction blend at the same face, and the opacity table on
   the rising side of the Fe bump.  A single-column or thin-wedge reproducer at that radius
   is the cheap way in -- the event is not turbulent, so it should reproduce.
2. **Put `fmlt1d` (and `mlt_relax_seeded_`) into the restart file**, or the production run
   cannot be chained at all (section 3).  This is a blocker independent of (1).
3. Only then the settling question (5 turnovers, `Gamma_rad`, `Mdot -> 0`) and only then the
   resolution.  `vr_rms/v_MLT = 0.034` and `F_res/F_req = 0.000` say the resolved flow is
   nowhere near taking the 10-23 % of `L` the closure carries, so the hand-over budget
   cannot even be asked yet.
4. Keep `dfloor_keep_velocity = true` and `vceil = 1.0e8` in the input: they cost nothing
   (bitwise-quiet, `eos_vceil = 0`, structure unchanged to three digits) and they convert a
   212 c evacuated cell into a cell pinned at 7 v_MLT.

## Files

* `README.md` -- this file.
* `gate.py <armdir>`, `an6.py <armdir>` -- the analysers (see the header).
* `one.sh` -- one athena invocation, no restart chaining; `chain.sh` -- the 0.5-turnover
  leg driver that section 3 condemns (kept, because it is what produced the evidence).
* `Dnr/ Gnr/ Gv5/ Hnr/ G/ H/` -- `he4.hydro.hst`, `he4.log` (event log),
  `rt_profile.bin`, `rt_surface.bin`, `column_he4_*.txt`, `mltfaces_*.txt`.
* `table_<arm>.txt`, `mass_<arm>.txt` -- the tables as generated.
* `<arm>.log` -- the job logs (`<arm>.1.log` is a symlink `gate.py` uses to find them).
  The per-cycle `face budget:` lines have been STRIPPED to keep them small, so re-running
  `gate.py` now leaves its `L/Lstar` column empty; that column is preserved in the
  `table_<arm>.txt` files, which were generated from the full logs.  The startup block,
  the `dt COLLAPSE` blocks, the sponge-guard prints and the leg boundaries are all kept.
