# NOTE 2026-11-21ev orion -> caltech (USER via coordinator, URGENT): start ALL launch smokes NOW, concurrently

Please start all five F_* smokes (NOTE eq, B_icbal, RT=1, 120 cycles) **now**, concurrently, one per free node: F_t9l62, F_t19l64, F_m35mms, F_t9l62_corona, F_t19l64_corona. MB_bsg62m is already running on 02-04.
- On any node still left: SQ_lq and SQ_gl (NOTE eu, seam A/B, 60 cycles), in parallel.
- **Hold** the v20.4 gate, CONVERR, ST/STD and every other non-launch arm until all of the above are running. If a non-launch arm holds a node an F_* arm needs, stop it; it is our own test arm.
- Keep R_LF running on 02-16.
Please commit each res the moment it lands, one commit per arm.

-- orion (for the USER)
