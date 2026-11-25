# NOTE 2026-11-25a orion -> Caltech: HG profile and old-vs-new time on H200 (SPEED50)


Idle H200 only, copies only, no rst / bin kept, HG only (RSG at Caltech is not allowed by the USER). Release binary: src b09a7b3c1 (A1+A2 in; build CUDA 12.9 + g++13, tiny BSG gate 3.686e-11), HG deck = the newest production HG keys Caltech runs, state = `hegiant.00205` (or the newest HG rst Caltech holds: say which and its cycle / dt). 4 H200, 6 MeshBlocks per GPU (24 blocks), or the layout Caltech ran before (say which).
1. **HGN**: 60 cycles plain (no outputs), ndiag 1: main-loop s/cycle over the last 50, dt per cycle (and the limiting cell if dt_report is free), passes and its per stage, FCT rounds, rays' per-build lines (vet_gd end lines), peak memory per GPU.
2. **HGO**: the same with the OLD HG production binary (1d1972805 line, prod-1118) and its own keys at the same rst, if a CUDA build of it exists at Caltech; if not, say so: no build is asked.
3. **PRF_hgn**: `profc.sh`, NCU_ON=0, NC=12 of HGN's state: res.txt, n/kname_r0..r3.txt, n/gaps_r*.txt, n/stats_r0.txt, CLEAN PASS / cycle lines with clean_verbose=2. Two nsys windows are better than one if the first cycle rotates the rays: give the window's rotation cycles.
4. **HG rotation cadence**: cycles of the last 400 with a rotation (`vet_rotate`), and the s/cycle of those vs the others (the 11-cycle 2 s spike of 11-21, H200PROF 5e: is it still there?).
Return ~1 node-h. -- orion (coordinator), drafted by SPEED50

Use idle nodes 02-03 / 02-09 / 02-16 (holders 7333072 etc.); the release binary is your rel-f 80d16d39 (b09a7b3c1). Coordinator: HG's dt reference is the OLD run's own dt at the same rst (equal to the new one, 117 s, at hegiant.00234 on viper; use the rst you have and report both dt).
