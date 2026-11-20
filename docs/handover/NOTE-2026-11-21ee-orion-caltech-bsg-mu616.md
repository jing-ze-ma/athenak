# NOTE 2026-11-21ee orion -> caltech (USER): 2 BSG mu=0.616 arms at the VERY HEAD (they gate the BSG launch)

The USER asked for these two arms. The decks are in e6d62fe6c, under docs/handover/fresh_prod_1120/freshgo/decks/ (copy the new decks into fresh_prod_1120 first). Run them on B_v203, 1 node, NR=4, with your armg.sh copy and RT=1. Start them **before L300T_bsg (and L300_bsg) if those have not started yet**; let running arms finish.
1. `RT=1 armg.sh M150_bsg $B_v203 decks/bsg_mu616_smoke.athinput 150 output2/dcycle=150 output3/dcycle=150 output5/dcycle=75`
2. `RT=1 armg.sh M150T_bsg $B_v203 decks/bsg_mu616_smoke_thinf.athinput 150 output2/dcycle=150 output3/dcycle=150 output5/dcycle=75` (THINF)

Each arm is a straight 150-cycle run plus a restarted b leg (75 cycles, rst, 75 more): bins at cycle 150 only, rst at 75, about 4 min per arm. Please commit res.txt, the hst files and the log tails as before (the bitwise check of the restart against the straight run is part of the gate).

-- orion (for the USER)
