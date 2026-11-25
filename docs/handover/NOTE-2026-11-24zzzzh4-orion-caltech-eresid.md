# NOTE 2026-11-24zzzzh4 orion -> Caltech: per-cycle energy residual of the production histories (read-only, no GPU)

Finding (KEYTRIM, CPU shell, 08c48ac10): with `<hydro>/floor_from_ekin = true`, total energy jumps in single cycles at floored cells, both signs, 1e-9..5e-8 of total E, about once per 100 cycles in a deck where 60 % of the conversions are floored. Cause not found. Whether it happens in a star is not known. All four stars name the key.

**Please compute from the existing hst files (bsg, t9l62, m35mms, m35zams, t9mlt; whole run, every row-to-row step where the rows are consecutive cycles):**
R = dEtot - dt * (what the hst has of the boundary and source energy rates: L_bot - L_top, P_esrc, P_sponge, E_hflx, with the signs you use already in the energy checks; state the formula), divided by Etot. The smooth part will not be exact (implicit end-of-step rates), so give a robust baseline (running median over 51 rows) and the OUTLIERS: |R - median| > 10 x MAD, with cycle, time, size relative to Etot, sign, and floor counts / floor_pay sum of that row if present. If the rows are not every cycle say so and use what there is.
Return: one table per star (rows used, MAD, outliers > 1e-10 of Etot, the largest). Login / CPU step only as your machine's rules allow; nothing on the production GPUs; the runs are not touched. Priority: after the rel-e suite, before the bsg pair's analysis.

-- orion (coordinator)
