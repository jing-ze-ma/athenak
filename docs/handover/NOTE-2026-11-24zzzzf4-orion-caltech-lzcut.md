# NOTE 2026-11-24zzzzf4 orion -> Caltech: exact L_z for the flux cut (fcone_lz + lz_spin): does it work? (USER: works -> default, else dropped)

Keys `rad_vet/clean_stag_fcone_lz=true rad_vet/clean_lz_spin=true` (on the cubed sphere also `clean_cs_lz_tel` if the binary asks: say what it asks). They have never run on a star. Binary: rel-a he_star_vet (5cfc3be5; if the keys are unknown there, say so and stop: they must exist at 95813cf85+). 
| arm | start | args |
|---|---|---|
| LZ0 (if HR's hst does not already serve) | HR's | HR word for word, +51 |
| LZ1 | hegiant.00205 | HR + the two keys, +51 |
| LZB0 / LZB1 | newest c1 bsg rst, a COPY | c1's keys, +100; LZB1 + the two keys (sp) |
Report: CW / FATAL / BE / fallback, dt, e_floor, floor_pay, L_top / L_mid / L_int vs the control at matched times, `CLEAN LZSPIN|FCONE_LZ|FCONE` lines (sums), gas + radiation L_z (Lz_tot, Lz_hflx, CLEAN LZB ledger) per cycle and the residual of bsg at c6050+ (the unexplained jumps: does LZB1 have fewer / smaller ones than LZB0 over the same 100 cycles? Take the 100 cycles from a rst where the residual is jumping (after c6050)).
**Works** = LZ1 and LZB1: CW / FATAL / BE / fallback 0, dt as control, L_top within 1.5 % of the control, gas + radiation L_z drift smaller than the control by >= 10x (bsg) or the cone's loss counted and returned; otherwise it is dropped. s/cyc with the keys vs without (cost). Priority: after HR150 / HRL; the bsg pair can share the node of the no-merge arms (same rst copy).

-- orion (coordinator)
