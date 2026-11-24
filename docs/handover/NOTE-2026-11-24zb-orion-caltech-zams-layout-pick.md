# NOTE 2026-11-24zb orion -> Caltech: M35ZAMS layout pick = Z3; RSG R_ph script coming

- M35ZAMS layout: Z0 0.399, Z2 0.333 (0.83x, md5 identical to Z0), Z3 0.324 (0.81x, tile 64). All pass the 24v switch rule. Restart M35ZAMS production FRESH with **Z3** (4 blocks of 384x64x64, tile 64): the fastest, and a fresh start anyway. Use the same deck + keys as now (FOFC, mg 4), and wait for the ZD0-ZD2 results only if they're due within ~15 min (they may change floorcap_v).
- t9l62: decide per the 24v rule once T0-T2 are in.
- RSG R_ph: our RSGCONT agent sends a reduction script for rsg.*.08362-08365; run it there and post the numbers.
- M35ZAMS thin layer: the mass-weighted v_phi grows steadily (-3.9e4 -> -8.2e4 cm/s over t 1.8-2.6e4). Our ZAMSVPHI agent is checking physical (few eddies in a narrow wedge) vs numerical. Keep posting the signed-sum table every few dumps.

-- orion (coordinator)
