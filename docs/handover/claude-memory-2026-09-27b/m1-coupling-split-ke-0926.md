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
- **USER 09-26 goal for the coupling:** correct dt-independent growth AND 2nd order in time for all variables in the stiff regime; cfl 0.9 support does not matter; accuracy first, speed second. mix 1:2 (ke-dt-0926) = stop-gap only (growth right, v order 1.1-1.5, 1.54x cost). Next: per-stage IMEX-RK coupling (hydro tendency inside each implicit stage).
- **ROOT CAUSE UPDATE 09-27:** the dt dependence comes from the WB radiation-FORCE REFERENCE energy bookkeeping (rad_m1_coupling.cpp:38-60: the radiation pays the full work, the hydro WB source applies rho*a_ref with no energy term; they cancel only to O(dt^2) per step). The stage-consistent work is disabled when fref is on (rad_m1_implicit.cpp:7373).
- **Fix:** <rad_m1>/force_reference_work = split (branch ke-dt-0926, f66ea71b + d5f95e1c, default full). gamma 2.19/2.22/2.29/2.45/2.59e-4 at cfl 0.075..0.9, dt->0 limit 2.16e-4 (27x smaller slope); the lid is clean. rho, E, v1 are 2nd order; v_h is still order 0.88, due to the R-internal lag of the frozen vet_sc tensor + opacities (Eddington + converged opacity -> 1.2-1.4). The implicit_opac_update path itself breaks rho/E order and needs repair. Cost 46 vs 38.7 ms (to re-time).
- **Next:** stage-consistent closure + opacities; gate (d).
- **09-27: 2nd ORDER ACHIEVED** with force_reference_work=split (v2 c0389721/724b7ed4) + converged opacity update (implicit_opac_update, one_pass=0, tol 1e-12), vet_sc kept: v1 1.94/1.92, v_h 1.95/1.79, rho 1.94/1.87, E 1.69/2.07; errors ~100x below the default. Pending: growth-rate gate for this set, interleaved timing, recommended key set.
- **09-27 GATES (a)+(b) PASS:** split + implicit_opac_update=true + implicit_one_pass=0 + implicit_tol=1e-8. gamma 2.15e-4 at every cfl 0.075..0.9 (to 3 digits). 2nd order in all variables. tol 1e-8 = same accuracy as 1e-12 at 5.2 Picard passes; 1e-6 loses order; maxit 2 unstable. Cost/cycle ~2.3x (tol 1e-12: 3.0x) vs default; split alone costs ~0 but stays 1st order. At EQUAL ACCURACY ~15x cheaper than the default (cfl 0.9 vs 0.024). Pending: sph_wedge (guard + term 81e88e44), LE shocks, tst, interleaved tol 1e-8 timing at cfl 0.9.
