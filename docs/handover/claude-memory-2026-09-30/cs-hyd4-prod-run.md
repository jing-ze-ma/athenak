---
name: cs-hyd4-prod-run
description: cs_hyd4_prod LAUNCHED 2026-09-23 - hydro-only twin of prod4 (prod4 minus B), jobs 11943719-22, bench/cs_hyd4_prod; do not relaunch/chain without the user
metadata:
  type: project
---

Launched 2026-09-23 by the user's approval ("if the smoke test works, launch a production run with
hydro"). Run dir bench/cs_hyd4_prod, jobs 11943719 -> 720 -> 721 -> 722 (apu, 2 GPUs / 1 node,
24 h links, prod4's submit.sh/chain.sh), FROM SCRATCH, tlim 8.64e7 (283 rot). Binary md5 fef759d9 =
snapshot 1d22a118 + the <hydro> cs_wellbalanced_src read (committed 1e9c0a04). Input =
inputs/production/deep_hot_jupiter_cs_hyd4.athinput (prod4 with <hydro>: hllc, no resistivity /
max_eta / STS / low-beta fallback / Maxwell outer BC, bbot = 0, hydro_w dumps; cs_wellbalanced_src on
to match prod4). Smoke (tests_hyd4/README.md): dt identical to MHD, eos_fail 0, tfloor 20-60x earlier
runs (untested hydro floors, audit item 3), 15.4-16.0 cycles/s -> ~61 h to 283 rot. Restarts 78 MB /
0.5 rot (~44 GB). Like prod4 it runs the LEGACY semi-implicit RT step (dhj pgen never reads
rt_semi_lin; see master-todo). Compare with [[cs-mhd-prod4-run]].

**09-23 00:45 cs_wellbalanced_src set FALSE in the run-dir input before the first link started (user: drop it if it does not help).** Min beta at prod3 rot 283 = 1.22 (>> ~2e-3 crossover); hurts hydro accuracy. Both runs changed together (twins). Smoke of prod4 input with it off: job 11943851 (bench/prod4_wboff_smoke) -- check its result before 02:36.

**09-23 ~03:40 APPLIED (user):** HSA_NO_SCRATCH_RECLAIM=1 in submit.sh; pending links resubmitted (prod4 11944873-5, hyd4 11944876-8); running first links (11941995, 11943719) unchanged.

**STOPPED 09-25 03:30 (user: cancel everything in prod).** All links cancelled; newest restarts: prod4 dhj.00179.rst (~rot 89.5), hyd4 dhj.00217.rst (~rot 108.5). Superseded by the new physics (ck_nquad 2, ck to the wall, ck RCE IC). Can be continued from these restarts with the old binary if ever needed.
