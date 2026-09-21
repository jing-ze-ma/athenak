# `wb_phi_eff`: the well-balanced scheme against an EFFECTIVE gravity

Branch `wb-phieff`, 2026-09-21, viper, serial CPU Release, `-D PROBLEM=box_convection`,
gcc/14.

## What the option is

`<problem>/wb_phi_eff = true` gives the x1 well-balanced scheme its own potential

    Phi_eff(z) = g0*(z - zmin) - Integral a_rad(z') dz'          (g_eff = g0 - a_rad)

so that the background it subtracts before reconstruction, and the gravity source it
builds out of that background's own pressure drop, together deliver `-rho*g_eff`
instead of `-rho*g`.  `phicc0` / `phi0` -- the TRUE potential, which `<hydro>/etotgrav`
carries inside the conserved energy -- are untouched, and so is every other consumer of
them.  Off by default and then bitwise inert: `Hydro::phicc_wb` and
`Hydro::phi_wb_x1f` are shallow copies of `phicc0` and `phi0.x1f` unless
`Hydro::EnableWBEffectivePotential()` is called, so the default path runs the same
arithmetic on the same memory.

`a_rad(z)` comes from `<problem>/wb_arad_file`, a two-column ASCII file
(`z a_rad`, strictly increasing z, linear interpolation, constant extrapolation), read
on the host, resampled onto the pgen's fine grid and integrated with the trapezoid
rule; cell centres and x1 faces read the same running integral through the same linear
interpolant.  `Hydro::SetWBEffectivePotential(phicc, phix1f)` is there for the
milestone in which the M1 module, not a file, supplies the profile.

`<problem>/wb_arad_force = true` applies, as a Strang-split source in the place a
radiation module's force will act, the RESIDUAL
`rho*(a_rad_actual - a_rad_carried_by_Phi_eff)` -- for this file-driven test source the
whole `rho*a_rad` when `wb_phi_eff` is off, and identically zero when it is on.  Its
`v1` work term IS added to the total energy: `a_rad` is not in the etotgrav potential,
so a bare momentum source would take the work out of the internal energy.

## G1 -- default-off bitwise: **10 / 10 IDENTICAL**

The two `tests_gate_merge/g1_gate.sh` box commands (50 cycles, 16x16 transverse, modes
3 and 0) against `athenak/tests_gate_merge/g1_new_m{3,0}`:
`feczrt.hydro.hst`, `feczrt.user.hst`, `column_used.txt`, `rt_surface.bin`,
`rt_profile.bin`, both modes, all byte-identical.

## G2 -- the static 1-D column

`he_1d_wb.athinput`: box_w8's He column, 134 cells, `nx2 = nx3 = 1`, NO radiation at all
(no two-stream, no Conduction object, `eos_radiation = false`), seed and both sponges
off, `bc_mode = 3` on both walls.  `make_ic_and_arad.py` builds the gas-only IC (the
production `rho(z)`, `T(z)`, with `e = e_gas(rho,T)` found by a secant iteration against
the code's own column dump) and then DEFINES
`a_rad = g0 + (1/rho) dP_gas/dz` on the fine grid, so the column is hydrostatic under
`g_eff` by construction.  `ic_hse_retune` is 0 throughout: the column is used as it is.

`a_rad` versus `kappa_R F/c` (Rosseland table of the production run, F = 2.475202e15):

| region | max ratio | min | median |
| --- | --- | --- | --- |
| whole box | 1.62 | -0.26 | 1.000 |
| tau >= 2 | 1.06 | 0.87 | 1.000 |

`a_rad/g0` runs -0.044 .. 0.750, peaking at z = -9.33e7 (the Fe opacity bump), which is
the PLAN's 0.749 g.  The disagreement is confined to tau < 2, where the production EOS
taper had already switched the radiation pressure off.

Each arm: t = 300 s, i.e. 34.6 sound crossings of a scale height
(H_p = 4.97e7 cm, c_s = 5.75e6 cm/s), dt ~ 0.139 s.  v_MLT = 1.86e4 cm/s.

| arm | `wb_phi_eff` | `wb_arad_force` | cycles | max abs(v1) over the run | /c_s | /v_MLT | where | max abs(v1) at the end /v_MLT | max density drift |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A (the fix) | true | true | 2153 | 8.48e3 | 1.6e-3 | **0.46** | z = +2.94e7 | 0.43 | 2.8e-3 |
| B (option c) | false | true | 2297 | 1.98e5 | 3.8e-2 | **10.6** | z = +2.41e7 | 10.6 | 0.40 |
| C (control) | -- | -- | 2748 | 5.95e5 | 1.1e-1 | **32.0** | z = +3.20e7 | 9.3 | 0.27 |

Arm C is the present setup: the production column, the tapered EOS with the radiation
pressure in it, no `a_rad` anywhere.  **PASS**: arm A is 70x quieter than the control,
not merely within a factor 3 of it, and is well below v_MLT; arm B, the unfixed split,
sits at 10.6 v_MLT (the PLAN estimated up to 18).

One-cycle force balance, `|a|/g0` per cell after the first step from rest:
arm A **6.5e-3**, arm B **0.25** -- a factor 38.

## G3 -- restart

Arm A restarted from `wb1d.00001.rst` (t = 151 s) to t = 300 s: the overlapping history
lines agree to a max relative difference of **1.2e-10**, and the final density and
internal energy are BIT-IDENTICAL (only v1, whose magnitude is 1e-3 c_s, differs, by
7.7e-8 relative).  It is not bitwise -- but neither is the control: arm C's own restart
differs by 4.7e-13 and arm B's by 4.0e-11, so this is a pre-existing sub-roundoff
non-determinism of this configuration, not the new option.  A wrongly rebuilt Phi_eff
would show up as an O(1) difference in the first cycle, not 1e-10.

## G4 -- style

cpplint (repo config, `--filter=-build/include_subdir`) on the five changed files
against their `rt-integration` versions: `hydro.hpp` 114 -> 112, `hydro.cpp` 4 -> 4,
`hydro_wellbalance.cpp` 8 -> 8, `hydro_fluxes.cpp` 30 -> 30, `box_convection.cpp`
0 -> 0.  Zero new errors.  No tabs, no trailing whitespace, no `}}`, all files mode 644.

## Files here

* `he_1d_wb.athinput` -- the 1-D template (`IC_PROFILE_PLACEHOLDER`)
* `g2_arm{A,B,C}.athinput` -- the three arms
* `make_ic_and_arad.py` -- builds `ic/ic_s11.txt` (the IC) and `ic/wb_arad.txt`
* `run_g2.sh`, `analyze_g2.py` -- run and report
* `ic/ic_s11.txt`, `ic/wb_arad.txt` -- the products, committed so the arms can be rerun

Run outputs (`g2/`) are deliberately not committed; `run_g2.sh` regenerates them in
about three minutes on one core.
