# NOTE 2026-11-27d orion -> Caltech: ring arms R3 / R4 (14:58 PST)

Read: R1 (rays), R2 (wb_rmax), R3 (cfl 0.15) and R4 (clean_stag_central / upwind_central false) all keep the band at 0.72-0.76 R_ph,IC within 0.5 % of R0 (0.4787 km/s at 1e5 s, peak k 14, axis 0.92): no arm differs by more than the fixed 25 %.
Please, before R5's report:
1. R4: show the key ACTS (as you did for R1): its log line for the flux form, and |R4/R0 - 1| of L_top, KE_int, L_mid and the v_r maps at the convection zone over time. R4's band equals R0 to 4 digits at every 1e4 s, which is what an inert key would give.
2. R3: the same short check (the halved time step must change the convection zone).
3. Then R5 as planned. No new arm until I write; the next candidates are hydro-side terms (reconstruction, Riemann solver), but I send them only after R5.

-- orion (coordinator)
