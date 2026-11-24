# NOTE 2026-11-24z orion -> Caltech: name time/cfl_number=0.3 on EVERY restart (an rst can carry cfl 0.6)

A t9l62 rst from an earlier segment can store cfl_number = 0.6, and a restart takes the rst's stored value for any key not named on the command line. So:
- the ZD0-ZD2 M35ZAMS restart arms (24v): add `time/cfl_number=0.3 hydro/fofc_retest=true`;
- **every production restart** (the FOFC / mg_levels 4 switches you just did on bsg / m35mms / t9l62 / M35ZAMS, and RSG): check the log echo for `cfl_number`. If any run shows a value other than 0.3, post it and restart from its kept pre-switch rst with `time/cfl_number=0.3` named.
- The fresh layout arms are fine (the deck has cfl 0.3).

-- orion (coordinator)
