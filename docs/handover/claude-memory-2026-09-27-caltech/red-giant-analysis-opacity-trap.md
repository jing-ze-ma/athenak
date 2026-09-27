---
name: red-giant-analysis-opacity-trap
description: "eoslib.get_kapr is solar_convection's ANALYTIC opacity and is 30-44x too high in a red giant interior; any F_rad computed with it in an analysis script is wrong by that factor. The red_giant RUN itself is clean -- it uses the tabulated stellar opacity."
metadata:
  type: feedback
---

Found 2026-09-08 when a flux profile said radiation carried only 3 % of the required flux
at the inner wall of `prod_topre` -- impossible, because that layer is subadiabatic
(grad - grad_ad = -0.077) and a stable layer must carry its flux radiatively.

**Cause: `tools/solar_convection/eoslib.py: get_kapr()` is the SOLAR_CONVECTION pgen's
analytic H- / Kramers / electron-scattering composite**, not the red giant's tabulated
stellar opacity. In a red giant interior it is wrong by a large factor, and the sign of
the error flips at the surface:

| | rho | T | kappa table | kappa eoslib | ratio |
|---|---|---|---|---|---|
| i=0 wall | 1.37e-1 | 5.49e6 | 0.778 | 33.98 | 43.7 |
| i=20 | 1.67e-2 | 1.50e6 | 9.04 | 380.8 | 42.1 |
| i=60 | 1.35e-3 | 2.82e5 | 380 | 10787 | 28.4 |
| i=280 | 2.38e-8 | 9.2e3 | 18.1 | 1.49 | **0.08** |

F_rad goes as 1/kappa, so the deep radiative flux was **understated 30-44x**. Corrected,
F_rad/F_req at the wall is **1.38, not 0.03** -- the bottom is a properly radiative,
slightly subadiabatic layer carrying the whole luminosity, as designed.

**Why:** the two pgens are different physics with the same-looking helper name, and the
analysis tools directory is shared. **How to apply:** in any red_giant analysis use the
run's own table, `data/stellar_opac/rosseland_gs98_x0.7_z0.014.txt`, read as log10 kappa
on (log10 T, log10 rho) with the header's `grid: nT nD lTmin dlT lDmin dlD`. The
scratchpad reader is `kaptab.py`. Never call `eoslib.get_kapr` on a red giant state.

## The RUN is not affected -- checked

`get_kapr` appears only in `cooling_convection.cpp` and `deep_hot_jupiter_rt_old.cpp`,
never in `red_giant.cpp`. The pgen loads `problem/opac_table` into the conduction module
(`rad_kappa_src = table_rho`; it FATALs on plain `table`), the grey two-stream takes
"opacity from the conduction module's table", and the startup log confirms
`every cell (ghosts included) is inside the opacity table's valid window logR = [-8, 1]`.

## What this retracted

An inner-boundary story built on the bogus 3 %: "the wall cannot shed the injected flux,
so it drives the deep motion". **Dead** -- radiation sheds it fine. What still stands and
is still unexplained: at t = 5e6 the layers i = 10-40 are neutral-to-SUBADIABATIC yet
carry v_rms 10-12x v_mlt; their mass flux rho|v_r|r^2 is ~5e26 against 2-4e23 at the
surface (1000x), so they are not the return flow of surface convection; and the growth
between t = 4 and 5e6 is 0.95-1.09 deep (saturated), 3.7-6.3 at i = 60-80, 1.00 at
i >= 160 -- i.e. spreading OUTWARD from the inside while the surface layer sits steady.
Next diagnostic: is it overturning or a trapped mode? See [[red-giant-flux-deficit-is-spinup]].
