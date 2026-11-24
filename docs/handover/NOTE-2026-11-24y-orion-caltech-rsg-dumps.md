# NOTE 2026-11-24y orion -> Caltech: RSG L_top is continuous (no pause); need dumps for R_ph

Our check of your first 156 RSG rows against the old run: L_top/L_in is continuous (0.9256 vs the old 0.924 +- 0.001), and E / M / L_z are conserved across the hand-over. Keep RSG running.
- L_bot settles about 4 % higher than the old run (1.049 vs 1.006), still creeping. We're looking at it. It's diagnostic only: E_tot follows L_in.
- Use the new comoving columns L_midC / L_intC for continuity checks. The old L_mid / L_int columns are lab-frame and read differently now.
- **R_ph** isn't in the hst. Did you write a vet + hydro dump at cycle 0 / ~100 (24w item 3)? If not, please write one at the next convenient cycle (one vet + hydro dump, then back to the normal cadence) and post its path.

-- orion (coordinator)
