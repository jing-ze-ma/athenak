# NOTE 2026-11-23d orion -> Caltech: correction to 23c - put the IA9 keys INTO the deck, not on the command line

Found by our smoke of the same three fresh runs (Orion/viper, IA9 HIP build, original J4_{bsg,m35mms,t9l62} decks from PROD1121b):
on a FRESH start, a key given on the command line that is not already present in the input file is FATAL
("Parameter 'time2_tableau' in block 'rad_m1' on command line not found", parameter_input.cpp:435). Only restarts may add absent keys.
The J4 decks lack 12 of the keys listed in NOTE 23c.

So please:
- add every IA9 key from NOTE 23c as a line in each deck's <rad_m1> block (values as in 23c), and
- keep only time/cfl_number=0.3 (already a deck key) on the command line, or set it in <time>.

Our own fresh smokes of bsg / m35mms / t9l62 on IA9 are running now (55 cycles each); we post their numbers (CW, dt, s/cycle,
L_top/L_in) for comparison. Everything else in 23c stands (branch caltech-ia9-77551e823, CUDA 12.9 + g++13, smoke first).

-- orion (coordinator)
