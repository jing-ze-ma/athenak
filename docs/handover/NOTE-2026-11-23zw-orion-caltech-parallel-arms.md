# NOTE 2026-11-23zw orion -> Caltech: run the m35mms arms in parallel (USER: faster)

To speed things up, please use every free node at once instead of running in sequence:
- **N2 node A:** D1 (key on), as in 23zv.
- **hpc-sm-01-17 (bsg's node):** D2 (key + max_stage2). The M35ZAMS smoke there should be finished by now. If it is still running, stop it: that deck is on hold anyway (23zs) and its timing to date is enough. Please post its last s/cycle + memory line.
- **D1h** (key on, cfl 0.15, 1000 cycles): on whichever of the two frees first.
- N2 node B: the t9l62 arms S / L, as now.
Same build (caltech-m35fast-13679216a), args and returns as 23zv.

-- orion (coordinator)
