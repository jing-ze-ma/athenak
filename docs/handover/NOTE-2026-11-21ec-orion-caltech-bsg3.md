# NOTE 2026-11-21ec orion -> caltech (USER): 3 BSG arms at the VERY HEAD (needed by ~16:30 PST for the go decision)

S_bsg shows the inner-wall L_bot sign flipping; the go decision needs these three by ~16:30 PST. They run on B_v203, 1 node, NR=4, with armg.sh (Caltech copy). All keys already exist in the smoke decks.
1. `RT=0 armg.sh L300_bsg $B_v203 decks/bsg_smoke.athinput 300 output2/dcycle=100000 output3/dcycle=100000 output5/dcycle=100000`: key off, 300 cycles, hst every cycle, no bins or rst (~3 min).
2. `RT=0 armg.sh L300T_bsg $B_v203 decks/bsg_smoke_thinf.athinput 300 output2/dcycle=100000 output3/dcycle=100000 output5/dcycle=100000`: THINF.
3. `RT=0 KEEPBIN=1 armg.sh F0_bsg $B_v203 decks/bsg_smoke.athinput 2 output2/dcycle=1 output3/dcycle=1 output5/dcycle=100000 output7/dt=100`: the cycle-0/1/2 hydro_w + m1 + m1_face bins (double).
For (1) and (2) please commit the hst and res.txt as before.
For F0_bsg, please **do not commit 1-2 GB to git**. Instead:
- Preferred: slice the **inner 16 radial cells** of each variable (cycles 0-2) with numpy / your bin_convert.py into one small .npz (with the x1 / x2 / x3 coordinates) and commit that to caltech_results/F0_bsg/.
- If slicing isn't quick, keep the full bins on /resnick and commit only their path and sizes; we'll ask for a specific slice.

-- orion (for the USER)

**Addendum:** the F0_bsg .npz should also contain **x1v and x1f** and the **m1_face fx1 on the inner 16 (+1) radial faces** (cycles 0-2).
