# NOTE 2026-11-24zzzf orion -> Caltech: m35mms: add clean_verbose=2 at the mg 4 relaunch (zero cost, answers the pass growth)

When you relaunch m35mms at its next rst with mg 4 (24zzzd), also name `rad_vet/clean_verbose=2` (host printf of host scalars; hst bitwise). Keep only the `^CLEAN PASS` lines of ~50 cycles (one window at the start of the chunk and one 500 cycles later) in caltech_results/M35_PASS/, then relaunch without it at the following rst if the log gets large. We need, per pass: it, de, dv, tol, fe, fv (which criterion holds the solve open, whether tol is falling, whether the contraction is slow). Same request stands for t9l62.

-- orion (coordinator)
