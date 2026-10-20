# NOTE 2026-10-20 raven -> caltech: answer to a381f5b6 (acc3d page energy residual) + L_z drift

## Energy residual (mk_prod.py:90, rE = dE - Ebnd - Erel + Emsk)
- Emsk is NOT inside Ebnd. In binaries up to 13206541 / 0fbfa31f, Emsk also books the kinetic energy that RyPerM1Post (sink_zero_v)
  removes from the sink's floor gas. That KE was put there by M1 coupling work the radiation never pays (outside the books), so
  creation and removal cancel, and booking the removal made the page formula wrong by exactly Emsk. Your no-Emsk number is right only
  because the two cancel.
- Fix (acclz-1020 8bd7982b, being integrated into prodmerge-1019b): that KE is no longer booked in Emsk -> mk_prod.py:90 is correct
  as written for runs on the new binary (rE 0.41 -> -5.3e-3 of |Ebnd| on a 39-cycle stream test; hydro.hst bitwise).
- For YOUR current run (13206541, hst written before the fix): show rE WITHOUT the +Emsk term (your 2.1e-3). Edit the page text to say so.
## L_z
- New lz.hst registers (Jgrv, Jcor, Jtr, Jrin, Jthw, Jspf, Jm1, Jres) close dJ to 1.4e-15 rel. Finding: at t = 0.01, dJ 8230 but
  Jgrv only 193; most of dJ is Jtr + Jcor (+8037): Coriolis source and frame part of the transport don't cancel on the grid.
  So don't read the drift as gravity torque. Registers arrive with the next gated fin51-r-fix (a NOTE will name the binary).
