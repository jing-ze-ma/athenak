# NOTE 2026-11-20dn orion -> caltech (USER): CONVERR, 14 one-cycle arms on B_v202 (keep the bins)

The USER wants the per-cell solver error of ONE production step, at production tolerance vs tight. The arms are in `docs/handover/caltech_runs_1121/arms_converr.txt` (B_v202, rsg.01295 / hegiant.00139, 1 node each).
- Every arm is **KEEPBIN=1 FINALRST=0**. The final hydro_w and m1 bins are written in double precision, plus an EOS table dump in a/.
- Each arm is about 4-5 min, writing ~6 GB of bins (HG) or 5 GB (RSG).
- **USER:** slot them ahead of the long runs if that costs < 15 min, spread over the free nodes, in **file order** (CEP_rsg, CET3_rsg, CEP_hg and CET3_hg first).
- **KEEP the bins** until the analysis script (converr.py + converr.sh, which uses numpy and your bin_convert.py) has run on them. That script comes in a follow-up NOTE; only its small text output needs to come back.
- If CET3_hg or CET2_hg FATAL or dt_retry, just report it; CET1_hg is the fallback.
- Lines: CW / FATAL, nl_print summary, bins kept (paths).

-- orion (for the USER)
