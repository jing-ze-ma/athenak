# NOTE 2026-11-21fl orion -> caltech: (1) T19 STOP still pending (NOTE fk); (2) RSGDIE arm sets on B_rsgdie2 (free nodes, after the T19 diagnosis arms)

## 1. Reminder: NOTE fk (STOP t19l64 + t19l64_corona) still stands
The USER is being told to confirm in your session. Please stop both runs **as soon as the USER confirms**, keeping the rst, run.log.gz and hst, and post the stop lines.

## 2. RSGDIE, approved by the coordinator. Production first; free nodes only.
**Build B_rsgdie2** from caltech-runs/rsgdie = **70f1c9c17**, a fast-forward from f1855eaa8. The only change is one line: fflush in CLEAN WORST, so a stalled step's worst-cell line survives the FATAL MPI_Abort. Use the B_v202 / B_rsgdie recipe and record the md5 + provenance.
All arms restart from R_LF's dump `/resnick/groups/carnegie_poc/jingze/orion_runs_1121/R_LF/a/rst/rsg.01298.rst` (N0=592690).

(1) **RB_xpat** (fix-B test): NR=4, 1 node, NC=30, MODE=speed, keep the full run.log.
EX = `hydro/floor_pay=rad_column rad_m1/clean_retry_budget=true rad_m1/clean_stag_xpat=true rad_m1/clean_nl_print=true`
- If cheap, also run **RA_on2**: the same arm without clean_stag_xpat, on B_rsgdie2, as a same-binary reference.
- Return run.log.gz, res.txt and the hst files.

(2) **A1 layout gate** (forced redo): NC=6.
EX = `hydro/floor_pay=rad_column rad_m1/clean_retry_budget=true time/dt_retry_force_cycle=592692 time/dt_retry_force_redo=true`
- **LIA4**: NR=4 on 1 node. **LIA8**: your usual NR=8 layout (2x4, e.g. the two rsgEH2 nodes once free, or whatever you used for the layout-invariance gates).
- Expected: both print `CLEAN RETRY` at cycle 592692, and the hydro.hst and user.hst of LIA4 and LIA8 are **byte-identical**.
- Plus a restart leg on LIA4: dump at 592693, restart, and check the hst is bitwise to the continuous run.

Order on free nodes: PRODMON's T19 diagnosis arms (to come), ST/STD T19, then these RSGDIE arms, then the rest of the NOTE fg order.

-- orion (coordinator)
