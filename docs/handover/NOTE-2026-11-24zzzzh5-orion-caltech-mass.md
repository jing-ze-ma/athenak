# NOTE 2026-11-24zzzzh5 orion -> Caltech: add MASS to zzzzh4 (sharper detector)

KEYTRIM found the cause of the energy jumps at cell level on CPU: the FOFC-replaced flux of ONE face on the periodic phi seam differs between its two sides (mass 1e-5 of the flux, energy 1e-3 of it), so mass and energy are not conserved there. Rare; with floor_from_ekin on. Why the two sides differ is being printed. Please extend zzzzh4:
- per row: dM_tot vs the booked Mdot_top / Mdot_bot (and M_hflx), relative to M_tot, same outlier rule; m35zams's switch row had 3.7e-16 so mass should be round-off.
- if the rows show events: cycle numbers and whether FOFC counts / floor counts of that row are unusual.
- how many FOFC-flagged cells per cycle the logs print for t9l62, m35mms, m35zams, t9mlt, bsg (typical and max).
Read-only, no GPU, runs untouched.

-- orion (coordinator)
