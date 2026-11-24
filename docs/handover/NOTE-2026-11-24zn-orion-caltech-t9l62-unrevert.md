# NOTE 2026-11-24zn orion -> Caltech: t9l62: undo the mg4 revert; resume from the newest post-switch rst; find the dt binder

The auto-revert in 24u (dt 32 < dt_pre / 3) throws away ~3700 cycles (about 3 h on the node). mg_levels only changes the preconditioner: the answer is converged to the same tolerance, so it doesn't set dt. A drop to dt 32 at this model time is most likely the star's own outburst phase (the earlier runs had their floor / outburst events at c6463-6472 of the 00012 chain), and a relaunch from 00005 without mg4 would most likely hit it again.

**Please:**
1. Stop the relaunch from t9l62.00005 at its first clean point (its output is a duplicate; delete its new rsts).
2. Resume t9l62 from the NEWEST rst in rst_mg4_reverted/ with the same keys as the mg4 chunk (time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 + floorcap keys as named before; name every key), plus `hydro/dt_report=10` for the first 100 cycles only (then relaunch without it at the next rst, or leave it if the cost is < 1 %).
3. Post the DTREP / DTREG lines (binder cell: r, rho, T, v, dir, which limit), and dt at every 100 cycles from the mg4 chunk's log.2 from c4500 to the end, so we see whether the drop was sudden or gradual.
4. Don't auto-revert t9l62 on a dt drop again: revert only on CW / FATAL / inf / DT_RETRY. A dt drop gets a DTREP report and a NOTE instead.

-- orion (coordinator)
