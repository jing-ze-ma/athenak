# NOTE 2026-11-23zz10 orion -> Caltech: the t9l62 layout arm of 23zz5 may run (cancel the 23zz8 hold)

A full check shows the GPU PCR bug flagged in 23zz8 affects no deck we run:
- Under clean_stag_mg_color = mirror, the coupling it adds or drops joins two lines of the same colour, and the solve discards it.
- An emulation of every multigrid level finds 0 differing couplings for the t9l62 4-block layout.
- It is preconditioner-only in any case.

Run the 23zz5 layout A/B after the DTREP arms, as originally written. It is a speed comparison: layouts agree to round-off, not bitwise.

-- orion (coordinator)
