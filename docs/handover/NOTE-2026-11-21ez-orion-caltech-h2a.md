# NOTE 2026-11-21ez orion -> caltech: H200 speed arms on B_h2a (FREE nodes only, lowest priority)

These are short bitwise speed arms for H200SPD. Run them **only on free nodes**, after the launch smokes (NOTE ex) and the v20.4c gate (NOTE ey). Production has priority; never take its nodes.
- Build caltech-runs/h200spd = 105e3db5c as **B_h2a**. It is v20.3 d4b08dae4 + 2 commits: face-flux kernels with a compile-time direction (CleanDirCase), bitwise by construction. Use the same recipe as B_v203; incremental is fine. Please record the md5 and the provenance (HEAD, compiled TUs).
- Arm file: caltech_runs_1121/arms_h2a.txt.
  - H2A_rsg / H2A_hg: rt, 10 cycles. hst + final bins must be **bitwise** to V203L_rsg / V203L_hg (B_v203). Report s/cycle.
  - H2S_*_a / H2S_*_b: B_h2a vs B_v203, 20 straight cycles each, **same node, back to back**. Report steady s/cycle excluding the first and last cycles. The HG pair runs only if a node is free.
  - Optional, if quick: profc.sh ncu of CleanRightHandSide (cln_rhsaf) on B_h2a rsg, reporting IPC, the no_instruction stall and instructions per launch, against PRF_rsg c/.
- Results go to caltech_results/<TAG>/ as usual.

-- orion
