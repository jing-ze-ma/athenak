# NOTE 2026-11-24t orion -> Caltech: USER: mg_levels 4 for all stars confirmed; R136 on GL 11x12 approved

1. **mg_levels 4:** USER confirmed for all stars (24s). Post s/cyc before / after per star. If any star is clearly slower (M35ZAMS is the one to watch), post it and we take it to the USER; don't switch back on your own.
2. **R136 (USER approved GL 11x12 + seam / wall keys + FOFC + mg 4):**
   - Binary: rc4, coming once its gates pass. R136's args target rc4 (the wall_sym extension for z-turning sets + the EMA ring fix are rc4-only).
   - rst: `/viper/ptmp2/jinma/r136fresh_1108/run/rst/r136.00149.rst` (11.3 GB, old implicit scheme, GL 10x12). It needs the same Globus route as RSG's rsg.01707 (USER), to `/resnick/groups/carnegie_poc/jingze/r136_rc4_1124/rst/`.
   - Args: `docs/handover/ia123rc3_decks/IA123_r136_restart.args` (95 keys, release names; 120 -> 132 rays on restart works). Layout: 4 MeshBlocks = one node, 4 H200.
   - Node: a freed test node once the PC / WD arms are done, or one of the RSG nodes while RSG's rst is still in transfer. Tell us which.

-- orion (coordinator)
