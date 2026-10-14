# NOTE Caltech -> Raven: queue cleared for bsgLIc (USER decision 10-14)

USER (10-14 ~21:00 PT): stop the other Caltech GPU jobs to make way for the Linhao BSG chain (bsgLIc 4487415, 4487610-13).
Stopped by scancel at 2026-10-13 21:01 PT; all data kept (bins, hst, logs, rst):
- **BSG fin3 link 2 (4434962)**: stopped at t 3.391e6 s (39.2 d), cycle 38757. Newest valid rst bsg3d.00040.rst t 3.3696e6 (39.0 d),
  cycle 38512 (rst_pick.py OK). Kept rst 00038-00040. No prune (link rc=x). The remaining 39 -> 57.4 d is NOT queued on Caltech.
- **He giant link 1 (4454171)**: stopped at t 6.79e5 s, cycle 36150 (heating issue, NOTE-2026-10-14-caltech-hegiant-heating.md).
  Newest valid rst hegiant.00015.rst t 6.48e5 s, cycle 34455. Kept rst 00013-00015.
- **E1 restart arms R0 4483862 / R1 4483863** (TASK e1div-1014 part 2): never started; removed. Dirs e1div_rst/R0, R1 (rst 00022 hard links) kept,
  can be resubmitted as they were.

bsgLIc link 1 still pending (Resources), estimate 2026-10-14 11:40.
