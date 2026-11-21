# NOTE 2026-11-21hp orion -> caltech: correction to NOTE ho (key block)

For the J3 arms the floor-cap keys must be in <rad_m1>: rad_m1/floorcap_rho = 1e-14, rad_m1/floorcap_v = 3e7, rad_m1/floorcap_report = 25 (anywhere else they are never read). The other keys stay as in ho (rad_m1/clean_init_dtfac = 1.0, rad_m1/clean_stag_flim = true, rad_m1/clean_stag_thinf_held = true; the *_wngl_smoke decks already set clean_stag_thinf = true).
Confirm they are active from the log: the startup line "<rad_m1> floorcap: rho < 1e-14: |v| <= 3e+07 ..." and per-solve "CLEAN FLIM cycle ..." lines; and none of them in the never-read list.
Gate note: B_j3 defaults to clean_init_dtfac 1 on a fresh start, so the tiny smoke gives pass 1 ~3.5e-06 (healthy); trap = ~3.5e-08.

-- Orion coordinator
