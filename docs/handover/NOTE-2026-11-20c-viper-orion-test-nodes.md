# NOTE 2026-11-20c viper -> orion (USER 11-19 ~17:30 PST "yes"): viper stops clean-cc; our 3 MI300A test nodes are yours

Times PST. The user decided viper stops its clean-cc (cell-layout) work and consolidates on your staggered line; our
findings follow in a separate handover NOTE (2026-11-20b) from our agent.

**Nodes for you now:** viper pHG 13255125, nodes **vipa1073, vipa1081, vipa1186** (3 nodes, 6 MI300A), free as soon as our
last test step stops (~17:35; check `squeue -s -j 13255125`: only steps on vipa1193/1203 (T19L6.4) and vipa1208
(dhj W178/W189) must remain). The allocation ends ~19:52.
**Tomorrow:** the pHG successor **13360518** (6 nodes, plain holder, PENDING): the same share, 3 of its 6 nodes; we tell you
the node names when it starts (T19L6.4 and the dhj runs take the other 3).

Rules as for the pRSG loan: `srun --jobid=<id> --overlap -w <your nodes>` steps only, never cancel the holder, do not touch
our steps on the other nodes or /viper/ptmp2/jinma/prod_1107 / t19l64_1119 / dhj_viper_1113. Please post a NOTE when you start.
