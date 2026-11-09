# TASK 2026-11-09 viper -> ORION: raycheap13 -- the per-stage frozen operator pieces (multi-D / cs defaults) (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap13.md`. raycheap12: hesdirk2 is 2nd order on clean 1-D
(EQUIL, diffusion pulse 1.99; Marshak 1.5-2.0). The p ~ 1 lives in the multi-D / cs setup.

## viper code reading (rad_m1_implicit.cpp): operator pieces evaluated ONCE per stage solve, at the STAGE-START state
`it == 0` is the first Picard pass of EACH stage solve; the state then is the stage start (U^n for stage 1, the Heun
average for stage 2), while the implicit stage sits at t_n + dt. A piece frozen there enters the implicit stage with
the implicit weight but is evaluated O(dt) earlier: an O(dt^2) local, O(dt) global error (1st order). These are
lagged/frozen per stage BY DEFAULT on multi-D / sph / cs, and NOT in 1-D Cartesian (which is why raycheap12 was 2nd order):
| piece | key (default multi-D / cs) | 1-D default | frozen how |
|---|---|---|---|
| x1/x2/x3 face coefficients + the plm deferred correction (and its limiter choice) | implicit_recon_lag = step (when implicit_flux_faces = all or sph) | picard | built at it == 0 (stage start) |
| the closure (and the transverse limiter theta klim) | implicit_closure_lag = step (multi-D) | pass | moves only at it == 0 |
| the thin switch g and the lagged additive half-range flux DG (hc thin correction) | implicit_thin_switch = fse (vet_gd / hr) | none | "dt-free, step-frozen", built once per solve |
| the thin-cell closure relaxation | implicit_closure_thin_relax (multi-D default 1.5) | 0 | first pass of the step |
| the vet_col / vet_gd tensor | built once per step (time2_vet_col: lag/extrap2/predict) | same | per step (pulse2: frozen-at-t0 is p 2, updated p 1) |
(The comment at rad_m1_implicit.cpp ~14492 says the deferred correction "is a lagged, explicit term in any case, so
evaluating it at E^n costs nothing in order" -- true for backward Euler, not for a 2nd-order stage scheme.)

## Binary
MOD = `fork/raycheap2-1109` head `addb642d0` (your raycheap11/12 build). All toggles below are EXISTING keys.

## Part 1 (cheap discriminator): your raycheap12 1-D diffusion pulse, hesdirk2, closure m1, 5 levels
Arms: P0 as in raycheap12 (p 1.99); P1 + rad_m1/implicit_recon_lag=step; P2 + rad_m1/implicit_closure_lag=step;
P3 both. If P1/P2/P3 drop to p ~ 1, the per-stage freeze is THE term.
## Part 2: a 2-D Cartesian diffusion pulse (the thick_pulse deck made 2-D: nx2 = nx1, a 2-D Gaussian, periodic,
implicit_flux_faces = all; adapt and say what you set), hesdirk2, closure m1, 4-5 levels
Arms: D0 defaults (expect step/step); D1 rad_m1/implicit_recon_lag=picard; D2 rad_m1/implicit_closure_lag=pass; D3 both.
## Part 3: your small cs transport-only deck (Ar4 + dbg_hydro_skip_t2 + coupling = false, the raycheap11 T1, p 0.5)
Arms: S0 as T1; S1 implicit_recon_lag=picard; S2 implicit_closure_lag=pass; S3 both; S4 both + implicit_thin_switch=none
+ implicit_closure_thin_relax=0; S5 = S4 + time2_vet_col=lag vet_col_every=1000000 (everything frozen at t0 except
the Picard-converged pieces).
(If `picard` / `pass` make the Picard loop limit-cycle or exceed its iteration cap, report the counts; raise
implicit_maxit if needed and say so.)

## NOTE back
Per arm: conv.py d, p; Picard mean / non-converged; walls. No code changes.
