# LHLLD benchmark suite: summary

Ad-hoc benchmarks for the low-dissipation HLLD solver (`mhd/rsolver = lhlld`, Minoshima
& Miyoshi 2021, JCP 446, 110639) against plain `hlld`.  These are NOT part of the
regression suite -- the low-Mach cases are far too slow -- they are run by hand from
this directory.  Details, command lines and tables are in the three companion files.

| benchmark | file | literature expectation | measured |
| --- | --- | --- | --- |
| Magnetized Gresho vortex, Mach scan | [results_gresho.md](results_gresho.md) | Leidi+22: LHLLD retention ~Mach-independent, HLLD degrading ~ Mach | NOT reproduced: both degrade then saturate; the LHLLD gain does rise monotonically as M falls (1.03, 1.10, 1.15 at M = 1e-1, 1e-2, 1e-3, beta 100).  The setup couples M to the Alfven Mach number, so it is not a clean test. |
| CPAW convergence | [results_cpaw_ot.md](results_cpaw_ot.md) | 2nd order must hold for both | PASS: ratio 0.291 (order 1.78) for both; solvers agree to 2e-5 |
| Orszag-Tang null test | [results_cpaw_ot.md](results_cpaw_ot.md) | the fix must be a null change on a supersonic solution | PASS: d1/d2 = 0.080, i.e. 8 per cent of hlld's own 128->256 discretization error |
| MHD Kelvin-Helmholtz growth rate | [results_kh.md](results_kh.md) | LHLLD resolution- and Mach-insensitive where HLLD is not | INCONCLUSIVE: the two solvers agree to 0.04 per cent (default Mach) and 0.35 per cent (10x lower), and BOTH are already resolution-insensitive.  Under rk3+wenoz the reconstruction, not the Riemann solver, sets the truncation error; a repeat at `plm` on a coarser grid is the way to see the effect. |

Identity check (not tabulated above): Brio-Wu boosted to a uniform Mach 5.7, where
chi = 1 everywhere, is BITWISE identical between `hlld` and `lhlld`.  Unboosted it
differs, i.e. the fix is live at rest.

NOT DONE: the optional Balsara advected MHD vortex of Leidi et al. 2022 sec. 4.2 (no
pgen exists for it), and the shock-detection factor theta of Minoshima & Miyoshi
eqns. (9)-(13), which needs transverse stencils the rsolver interface does not provide.
