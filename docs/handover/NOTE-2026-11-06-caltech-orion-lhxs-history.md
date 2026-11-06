# NOTE 2026-11-06 Caltech -> Orion (USER): AthenaK history for the overlay (reply to TASK-2026-11-06n)

`docs/handover/athenapp-bsg/orion3d/athenak_lhxs_history.txt`, from `prod_lhxs`. KE is on hst rows (~1000 s). dU, dEgrav and r_photo are `nan` except on 0.5 d dumps, because they need dumps; only t = 0 exists so far, and 0.5 d comes shortly. I'll refresh it with reports.
Definitions as in `orion3d/summary.txt`: pi/3 sr wedge, cgs; dU = change of sum(Pg/(g-1) + E_r) dV (M1 E, lab); dEgrav with GM = 2.62502e27; r_photo = tau 1 of shell-mean rho kappa_Ross from the outer face (T = (P/rho) mu/R on our float32 dumps, his table). Please check that against yours at t = 0 and say if it's not the same measure.
The run is AthenaK with M1 + vet_gd, kick-fixed v* mimic, RK2, and an atmosphere of g ~ 1/r^2 above 52.5 Rsun (USER 11-04), so it differs from his setup there. Please say so on the overlay.
