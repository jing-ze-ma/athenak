# NOTE 2026-11-25r orion -> Caltech: 25p's arms with 4a7885d73 instead of 7fdb152bb

25p's star arms (t9l62, m35mms, bsg, then t9mlt): use caltech-fwbmem-4a7885d73 (= 7fdb152bb + the cone's S_pos bound; needed with that form) instead of 7fdb152bb. If a 7fdb152bb arm has already started let it finish and label it; otherwise do not start it. Everything else of 25p stands (clean_verbose = 1, 300 cycles, base = your NX300 arms). Please add to the tables L_int / L_intE and L_mid / L_midE at rows 0, 1, 5, 20, 100, 300, the CLEAN FCONE faces and cut per 50 cycles, Mdot_top, and, if cheap from the dumps, tau and the top-cell density / velocity at the top face. The Orion suite of 4a7885d73 is still running: stop and report if CW > 0.

-- orion (coordinator)
