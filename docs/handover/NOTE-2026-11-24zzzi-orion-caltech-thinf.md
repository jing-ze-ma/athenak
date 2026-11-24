# NOTE 2026-11-24zzzi orion -> Caltech: USER: test the stars with clean_stag_thinf = true (P0 in 24zzzh)

The USER confirms the thinf arms of 24zzzh; make them P0 with Compton. Clarifications:
- The key is `rad_vet/clean_stag_thinf` (the stars' decks name false; the suite runs the fresh default true). On each of bsg (fresh, 300 cycles, against the T2 fresh c0-300 in hand), t9l62 (T2), m35mms (M3) and M35ZAMS (Z1): the true arm against false from the SAME rst, all else equal, every key named. Report true-vs-false: s/cyc, dt (min / mean), its and passes, CW / FATAL / inf / DT_RETRY, L_top / L_bot / L_in, KE_int and M_tot drift, and the CLEAN FWB / FCONE counts if the log has them.
- If a restart FATALs because the key can't change on a restart, say so and run that star's pair FRESH (300 cycles each) instead.
- Don't run any other thinf-related key (thinf_held, flim) unless the arm needs it; if it does, name it and say so.

-- orion (coordinator)
