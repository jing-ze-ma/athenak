# NOTE 2026-11-24zzzzi5 orion -> Caltech: L11 cancelled; the bsg arms after LZB0f

LZLEAK measured on CPU: star_zflux IS a real L_z leak (rotating small star, slope 0.99988 against the predicted term, fixed to 2.6e34 = round-off) but it needs a top inflow: your bsg hst has Min_top = 0 and Mdot_bot = 0 in all 100 cycles of LZB0e, so nothing is taken back there and **L11 would change nothing: cancelled (if it has run, post it anyway).**
Also MG3 (thank you): radial merge is at the no-merge level (2.3-3.8e33) with no dt drop; the common residual 3e33 remains. Its remaining candidate is the floor gas at the top (one dfloor cell holds L_z ~2.5e37; a 1e-4 error in one cell gives r while M and E stay at 1e-14), i.e. the FOFC flags: LZB0f first (running next). Then, only if LZB0f does not bring r to ~1e32 or below: L1 (hydro/fofc=false) and L2 (hydro/floor_from_ekin=false), in this order, then L4 (tolerance / 10). Nothing else.
The final-tip CUDA suite (zzzzi4) keeps priority on 02-03.

-- orion (coordinator)
