# NOTE 2026-11-20dw Caltech -> orion (USER): node plan for fresh production ~17:00 PST (answer to 4c31e5bc0)

Checked 15:10 PST 11-20.

- **testH1 = job 6133593, node hpc-sm-01-11**, 1 node x 4 H200, ends **11-22 00:34 PST**.
- **The other holders end 11-21 17:46 PST** (~26.5 h from now):
  - hgGH3b 6133592: hpc-sm-01-17, 02-03, 02-04 (3 nodes);
  - rsgEH2 6133591: 02-09, 02-10;
  - r1H1 5840649: 02-11;
  - r1H2 5840650: 02-16 (R_LF is on it, started 14:54).
  
  They cannot be extended. A fresh production started at ~17:00 therefore gets about 24.5 h on them (31.5 h on testH1), unless it is
  restarted into the successor holders. Those are queued but still pending, with no start estimate: 7333072 (3 nodes), 7333073 (2 nodes),
  7333035/36/37 (1 node each), all 3-day. Plan for checkpoint/restart at the holder end: rst before 17:30 on 11-21.
- **Confirmed:** from ~17:00 every node except R_LF's (02-16) can go to fresh runs. The test arms running on them then (CONVERR / seeds /
  EW / resumed L2 / LX) are stopped as needed and logged for a later restart. RA_on / RA_off still go on the first free node once R_LF's
  592690 dump exists (~17:30). They are 30 cycles, so they hold that node only briefly.
- FRESHGO smoke arms: at the head of the queue on B_v203 as soon as the NOTE arrives, preempting CONVERR / seeds / EW.
- **Fresh production:** we prepare it from FRESHGO and launch **only on your explicit go** NOTE.
  - Every GPU job excludes hpc-sm-01-09.
  - Multi-node runs use OMPI_MCA_io=romio321.
  - Each rst write must fit the fixed 1 GiB piece logic (v20.2+). The rst write goes through romio321 multi-node, and the v20g2/g3 fix covered it.
- CET1_rsg: rerun once on a node checked idle with nvidia-smi (operator). The CONVERR script runs when the arms are done; the 4 text files go to caltech_results/CONVERR/.

-- Caltech
