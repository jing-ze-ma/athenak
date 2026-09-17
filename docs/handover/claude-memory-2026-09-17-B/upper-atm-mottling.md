---
name: upper-atm-mottling
description: The mottled T texture in the upper atmosphere is the H2 dissociation front and cold sinking plumes, NOT noise, floors or a meshblock artifact
metadata:
  type: project
---

Temperature maps of the upper atmosphere (both [[ck-limb-run]] and [[ck-grav-size-run]])
look speckled, worst in the terminator meridional slices. **Measured 2026-08-26: it is
physical structure, marginally resolved.** Do not re-diagnose it as noise.

**What it is.** Cold (T < 1800 K), dense, fast-sinking plumes of recombined H2 on the night
side of the terminator. Measured in the leading-terminator slab (|lon+90| < 20,
r = 1.25-1.40e10, ck_limb, three consecutive dumps):

| test | cold cells | random mask, same fraction |
|---|---|---|
| cold neighbours of 6 | 4.26 | 1.11 |
| isolated cells | 6.3 % | 29.2 % |
| overlap one rotation later | 80 % | -- |
| at dfloor or pfloor | **0.0 %** | -- |
| median v_r | **-5.2 km/s** | -0.9 km/s (hot cells) |

Clustered, persistent, sinking, and nowhere near a floor.

**Why it cannot look smooth, at any resolution.** At rho = 1e-11 the tabulated specific
internal energy runs 10^11.19 (1500 K) -> 10^12.07 (2000 K), a factor 7.7 across 500 K, then
only 1.7x more over all of 2000-5000 K. That step is the H2 dissociation reservoir. A cell
is on one branch or the other and the intermediate band is nearly empty, so a T map of the
transition has no middle colours to draw. **Pressure runs smoothly across the same region.**
This is the same stiffness as the "Gamma_1 varies 0.58 between adjacent cells" lead in
[[dhj-ck-eos-blowup]] -- which is CLOSED (it was the race), so do not reopen it.

**Three things it is NOT, each checked:** floors (0 % of cold cells floored, rho 55x dfloor);
a meshblock artifact (ring-to-ring jumps across the theta block edges j = 15/31/47 are the
same size as everywhere else); a checkerboard oscillation (see the table).

**"Especially near the poles" is the GRID, not the flow.** Cold-plume fraction peaks at LOW
latitude: 22.6 % for |lat| < 30, 16.0 % at 30-60, 9.3 % poleward of 60. With
f_stretch_theta = 3 a polar cell spans ~8 deg and an equatorial one 0.84 deg, so on a linear
latitude axis polar cells are fat blocks and equatorial ones hairlines, and the blocks read
as noise. Gouraud shading makes it worse -- re-render with `shading='nearest'` before
judging any texture in these slices.

**Separate and also real:** the outermost cells carry a floor-supported halo (36 % of polar
top-8 cells at dfloor in ck_limb, 15 % in ck_grav_size), whose pattern is 95-99 % identical
between dumps a rotation apart. It sits at log10 p ~ -8.5, i.e. **2.5 decades above the
1e-6 bar level**, and 0.00 % of cells ON the isobar are floored. Earlier notes saying "three
decades" were loose; 2.5 is the measured value.

**Open:** plume amplitude is not converged -- they are a few cells wide, and x1 cannot be
refined much further ([[nx1-ceiling-lds]]).
