# NOTE 2026-11-21es orion -> caltech (coordinator, PRIORITY): F_* smokes get nodes FIRST; hold the v20.4 gate

The F_* final-recipe smokes (NOTE eq: F_t9l62, F_t19l64, F_m35mms, F_t9l62_corona, F_t19l64_corona, plus MB_bsg62m, all on B_icbal) gate the production launch. They must get nodes **before anything else**.
- **Do not start any new v20.4 gate arm (NOTE eo)** until all six F_* / MB arms have started on their own nodes.
- v20.4 arms that are already running: let them finish if they are nearly done. If one would block an F_* arm for more than a few minutes, stop it; it is our own test arm and can rerun later.
- After all F_* / MB arms have started, the order on free nodes is: the v20.4 gate, then the ST_* / STD_* stretched arms (NOTE er). Production preempts all of them on the go NOTE.
Please commit each F_* / MB result the moment it finishes.

-- orion (coordinator)
