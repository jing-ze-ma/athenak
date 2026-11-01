# TASK 2026-11-01 Caltech -> viper (USER): floor-mass tracer, please take it over (with your xe_tracer)

The USER asks viper to build the floor-mass tracer next to your x_e tracer, so the two passive scalars are designed together. Caltech has stopped its work on it.

## Goal
Remove the floor mass that remains with `<mhd>/valfven_noinfall_fac` (merged 1fb391a8).
- With the fix there is no radial rain (zero floored-to-unfloored infall).
- What remains is ~1.0e14 g/s exchanged SIDEWAYS from Alfven-floor cells into unfloored cells, mostly at night.
- Over 300 rot that would add ~5.4e21 g at fac 4, ~2.1x the mass of the 1e-5..1e-3 bar layers.

## Caltech WIP (not built, not gated, not merged)
Branch origin/floortrc-1101, commit 45010a2b, on fin51-r-fix b935e145. Feel free to redesign it.
1. **Key:** `<mhd>/floor_tracer` (read only when named). It is one passive scalar in the LAST slot, ifltr = nmhd + nscalars - 1, so the xe_tracer scalars can go before it. Optional keys `floor_tracer_margin` (2) and `floor_tracer_fmax` (0.01).
2. **Source:** in `GeneralMHD::ConsToPrim` (src/eos/general_mhd.cpp) and `IdealMHD::ConsToPrim` (src/eos/ideal_mhd.cpp), the `if (dfloor_used)` write-back adds u.d - d_in to the scalar. This catches dfloor and the Alfven floor (valfven_max_cgs).
3. **Drain:** new `MHD::FloorTracerDrain` (end of src/mhd/mhd_tasks.cpp), called in `MHD::ConToPrim` before the floor-mass tracking, and skipped on the frozen first conversion after a restart.
   - Where rho > margin * max(dfloor, B^2/v_A,max^2): dm = min(s, fmax*rho, rho - margin*rho_fl).
   - v and specific e_int are kept: rho, m and e_int are scaled by (1-f); E loses f(e_int + 1/2 m.v), plus dm*phi under etotgrav; the general-EOS pressure cache is scaled by (1-f).
   - The scalar is advected by the existing MHD scalar fluxes.
   - NOT checked: scalar handling in the outflow1 top ghosts and at the cs seams.
4. **Restart:** main.cpp records `global_variable::rst_mhd_nscalars`. NOT done: the MHD u0 reader in pgen.cpp has to read the file's slot count and zero the floor slot (copy the hydro `scalars_rst_missing` code, pgen.cpp ~270-295 and ~1120-1135).
5. **With the tracer on,** hst Mfloor (dhj.user.hst column 13) becomes the NET rate (added - drained).

## Gates we would use
- Key off: bitwise on CPU (r1m at C8x112, 20 cycles) and GPU vs b935e145.
- Key on: clean CPU run; restart from an rst file without the scalar; 0 <= s <= rho.
- **GPU smoke:** r1m production keys (fac 4, outflow1, rkg_dt_safety 0.5, rkg_lean), t=0 to 0.5 rot. Baseline p_f4 (1 H200):
  - Mfloor 1.64e14 g/s; dt 5.0 s at 0.4 rot; E_B 1.528e32 erg;
  - 1e-5..1e-4 bar band: 2.648e20 g / 2732 K / 2.53 km/s rms / 8.2 G.
  - p > 1e-5 bar must stay unchanged.

## Caltech plan
WASP-121b r1m stays HELD until your x_e tracer NOTE (and the floor tracer, if it is ready around the same time). We will launch with your recommended keys.
