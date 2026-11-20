# NOTE 2026-11-20i orion -> caltech (USER): batch 4 = H200 layout / launch tuning (27 arms), after the queued arms

Files are in `docs/handover/caltech_runs_1121/batch4/`:
- `arms_batch4.txt`: arms, priorities and per-group rules are in its header.
- `armc2.sh`: armc.sh for 1-3 nodes (`NODES=<comma list>`, NR spread evenly, `NSYS=1` option, prints the hst md5s). Apply your
  armc_cal.sh changes (modules, res.py as a step) to it too.
- `nsyswrap.sh` (nsys per rank), `memwrap.sh`, `syncattr.py` + `r0lab.py` (nsys post-processing).

Builds (same cmake line as before; all on v20g 6b48f352c):
| build name | branch | commit |
|---|---|---|
| B_v20g | caltech-runs/v20g | 6b48f352cbdff12bd2c9d9fd708767ebf1a8aa5f (you have it) |
| B_p4 | caltech-runs/h200-p4 | dbf6a3053afd9cb1f8c034c3778b4464c070f016 |
| B_h200 | caltech-runs/h200 | a774c43b3c38d3f71e6a725fe03b5f073835a171 |
| B_h200mb2 | caltech-runs/h200 (same commit) | add `-D Athena_EXTRA_DEFS=ATHENA_CLEAN_FOR_MB=2` to the cmake line |

Priority: P1 gate (BA_* / BR_* / BP_*: the hst md5 of the four builds must be equal; rt and poison bitwise) > P2 layout scan
(LS_*, NR 2/3/6/8/12 on 1-3 nodes) > P3 nsys (PF_*) > P4 GDR on/off (TG_*, optional, last).
- Multi-node arms: `UCX_IB_GPU_DIRECT_RDMA=no`.
- HG arms drop `time/dt_retry_host=m1`, except BA_hg_v20g_m1.
- Timing arms: one arm per node at a time; not on hpc-sm-02-16.

Send back per arm: res.txt (hstmd5 lines, steady s/cycle, GPU memory peak) and info.txt; for PF_* also the syncattr/r0lab
output (not the .nsys-rep files unless small).

-- orion (for the USER)
