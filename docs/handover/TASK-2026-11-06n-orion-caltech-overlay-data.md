# TASK 2026-11-06 Orion -> Caltech (USER): AthenaK `prod_lhxs` history for an overlay on Linhao's Fig. 2

The USER wants the Orion Athena++ run and your AthenaK `doface_1105/prod_lhxs` overlaid on his Fig. 2, on one page, at about 1 d star time. It is built here; only your curves are missing.

**Please push `docs/handover/athenapp-bsg/orion3d/athenak_lhxs_history.txt` as soon as you have data, then with each of your reports.** Plain text, one header line, one row per ~1000 s hst or per 0.5 d dump, whichever you have:

`t[d] KE[erg] dU[erg] dEgrav[erg] r_photo[Rsun]`

- pi/3 sr wedge, cgs. dU = change of sum(Pg/(g-1) + E_r_lab) dV. dEgrav = change of -sum(GM rho/r) dV with r = x1v and your GM. Same as `orion3d/summary.txt`.
- Use `nan` for anything you can't compute cheaply, and say so in a one-line NOTE. KE alone is already useful.
