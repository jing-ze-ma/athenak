# NOTE 2026-11-21gw orion -> caltech: ray-intensity dump on the trapping H200 (2 one-minute runs, after NOTE gv)

Decides "rays wrong" vs "closure wrong". Decks: docs/handover/caltech_runs_1121/close/bsgs_n1.athinput (clean_init_n 1 = pass 0) and bsgs_n2.athinput
(clean_init_n 2 = pass 1); they are the BSG wngl smoke at 32^2, one MeshBlock, nlim 0 (keys are in the decks; do not
put clean_init_n on the command line, it is FATAL). Point the input-file paths at your fresh_prod_1120 BSG inputs if needed.
On a free test node, B_v203, 1 GPU / 1 rank:
  VGD_DUMP_I=$PWD/dumpI_h200_n1.bin athena -i bsgs_n1.athinput
  VGD_DUMP_I=$PWD/dumpI_h200_n2.bin athena -i bsgs_n2.athinput
Commit both .bin (gz; ~1x126x62x62x262 doubles each) to caltech_results/H200CHK/dumpI/ plus the clean_init lines.
If B_v203n (CUDA 12.9) is built by then, also run bsgs_n2 with it -> dumpI_h200new_n2.bin.
The newer-CUDA smoke (NOTE gt) and gv stay first.

-- Orion coordinator
