# NOTE 2026-11-20dm orion -> caltech (USER): build B_rsgdie now; 2 replay arms from R_LF's 592690 dump

- **Build now:** **B_rsgdie** = caltech-runs/rsgdie, commit f1855eaa8dcbd094168fcee6790dc2236b7e6d23, with the same he_star_m1 H200 recipe as B_v202 and provenance as usual. It is v20.2 code plus the key rad_m1/clean_retry_budget (default off = bitwise) and a diagnostic CLEAN WORST line, printed only after a failed backward-Euler solve.
- **Once R_LF has written its 592690 dump,** run two arms (in `arms_rsgdie.txt`) restarted **from that dump**, NR=4, 30 cycles (past the death at 592698), with the LF_rsg keys (floor_pay=rad_column) plus clean_nl_print=true:
  - RA_on: clean_retry_budget=true.
  - RA_off: clean_retry_budget=false. This is the control; it should die like LF_rsg.
  - They are short, so take any node, even ahead of the batch.
- Please return both run.log.gz files.

-- orion (for the USER)
