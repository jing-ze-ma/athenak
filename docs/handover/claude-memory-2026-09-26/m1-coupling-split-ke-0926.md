---
name: m1-coupling-split-ke-0926
description: ROOT CAUSE of the He-box KE dt dependence (09-26) - the gas/radiation coupling split is 1st order when the energy exchange is stiff; real mode gamma0 1.5-2.2e-4/s; hesdirk2 overshoots, Strang multi-rate undershoots
metadata:
  type: project
---
/viper/ptmp2/jinma/kedt_0926/RESUME.md.
- **Mode.** Horizontal oscillation |k| 7-8, period ~118 s (plus radial acoustic modes). Driven by the radiation force (gas force off -> decays) and damped by the exchange (exchange off -> explodes).
- **Growth rate.** gamma0 + slope x coupling interval: hesdirk2 +2.0e-3 per s, be +5.5e-3, Strang multi-rate every 2 about -2e-3. All extrapolate to gamma0 = 1.5-2.2e-4 /s: real and dt-converged at this resolution.
- **Ruled out.** vimp, tolerances, one-pass, predictor, closure lag, opacity update, the transverse force.
- **Mechanism.** Adiabatic heating over the coupling interval, then re-equilibration in the implicit solve: first order in the stiff limit (driver.cpp:450-462 hesdirk2, 486-515 be; rad_m1_mr.cpp).
- **Fix needed.** A stiff-accurate coupling: implicit compression work in the radiation solve, or cancel the opposite-sign splittings. Until then, quote KE as the hesdirk2 / multi-rate bracket.
- The lower-bracket job 11990312 (multi-rate cfl 0.15/0.3 to 34000 s) was pending.
Related: [[m1-hesdirk2-cfl-recommendation]], [[m1-multirate-rejected-0926]].
