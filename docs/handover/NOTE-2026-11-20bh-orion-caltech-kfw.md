# NOTE 2026-11-20bh orion -> caltech (USER): KERN speed pairs (new build kern-fw) + 2 optional env-only arms

New build **B_kfw** = caltech-runs/kern-fw, commit 7ee3e80ddee422fb2ef00d92886fcc1964962994 (v20l plus 2 commits; it should be bitwise). Plain cmake line, as for B_v20l.
Arms are in `arms_kfw.txt`:
- RSG and HG pairs, B_v20l vs B_kfw, 20 cycles as in the BA arms. Round 1 together on comparable nodes, round 2 with the nodes swapped.
- Optional: KL8_rsg / KL2_rsg (env VGD_CS_LB, same binary as v20l).
Per arm: s/cycle, node, GPU peak, and hst md5. The md5 must be equal to the B_v20l arm of the same star (KL* vs KF_rsg_b_1).
Normal priority.

-- orion (for the USER)
