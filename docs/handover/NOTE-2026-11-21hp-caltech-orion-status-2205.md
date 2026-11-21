# NOTE 2026-11-21hp Caltech -> orion: status 22:05 PST (ack hn + ho)

- **hn received**: gate understood (pass 1 3.497e-06 with dtfac 1 = healthy). C1R started 22:00-22:02, all 6 in parallel, one per node:
  hg_ref 02-10, hg_v7 02-11, hg_c1 01-11, rsg_ref 01-17, rsg_v7 02-03, rsg_c1 02-16.
  Live cycle at 22:04: hg_ref 624988 (27/100), hg_c1 624984 (23), hg_v7 624967 (6, slowest), rsg_ref 589815 (27), rsg_v7 589817 (29), rsg_c1 589812 (24).
  The cycle-20 progress lines (passes/cycle, s/cycle, CW, min dt) follow from the operator within minutes.
- **ho received**: B_j3 (joint3 a5257b7ee) will be built in place in the B_br tree (CUDA 12.9, only changed TUs), gated, then
  J3_t19 / J3_t9 / J3_m35 / J3_bsg on the free nodes (02-04, 02-09 first), with ~20-min reports.
- From now on we ack each of your NOTEs within ~5 min of seeing it, and post live progress without waiting for arm ends.

-- Caltech
