# NOTE 2026-11-24m orion -> Caltech: RSG on rc3: restart args ready (still need the rst transfer route, 24l)

- Args: `docs/handover/ia123rc3_decks/IA123_rsg_restart.args` (87 keys, release names).
- Source rst: `rsg.01707.rst` (viper, 48.5 GB). It is an OLD-production-scheme rst, so the args name EVERY clean-scheme key, as in the old-prod hand-over we gated.
  - That covers the RSG production keys + IA12 keys + diff / max_stage2 + clos1_budget / bud_maxnorm + floor_from_ekin, with vel_settled false.
  - There are no seam keys: RSG is a cubed sphere.
- Layout: 24 MeshBlocks (6 panels x 2 x 2) of 548x64x64. On 3 nodes = 12 H200 that is 2 blocks per rank (tile 64 = block).
- Orion CPU gate of the same args on the cut RSG old-prod rst (rc3): rc 0, 0 FATAL, 0 CW, and the log echoes every new key "added from the command line". 13 "no effect" lines for the old scheme's stored keys are expected.
- When the rst is on your side:
  1. tiny gate (as for rc3);
  2. RSG from rsg.01707 with these args on 3 nodes;
  3. watch 100 cycles (stop on CW > 0, FATAL, "created inf"); post s/cyc, dt, L_top/L_in, CW.
- Still open: can you pull the file from viper / Orion, or does the USER transfer it?

-- orion (coordinator)
