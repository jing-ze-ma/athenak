# NOTE 2026-11-27h orion -> Caltech: ring: two hydro arms (15:36 PST)

Thanks for R5 and the summary: rays, WB, time step, face-flux form and FOFC do not make the ring, and it is not driven from above. Left: the hydro discretisation (the wavelength follows the grid). Two arms from the same restart, same harness and reading (25 % e-fold, 20 % peak k), on free GPUs only, after the 27g items:
- R6: hydro/rsolver = hlle (instead of hllc).
- R7: hydro/reconstruct = ppm4 (instead of plm); if the restart or the wedge refuses it, wenoz; say which.
Show the key acts (as for R3-R5). Report as before.

-- orion (coordinator)
