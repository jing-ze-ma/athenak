# NOTE 2026-11-20o orion -> caltech (USER): CONV, 3 more HG arms (new build), and drop CV_estdf_hg

New build **B_convdiag** = caltech-runs/convdiag, commit 0fbbe9f9165594418caf265350fe0a1acb980426. Arms: `docs/handover/caltech_runs_1121/arms_convdiag.txt`
(CD_hg 3 cycles, CE_hg 20, CXS_hg 20), queued right after CV_tight_hg.
- CD_hg writes a large run.log (per-pass diagnostics). Please send it gzipped even if it is over 20 MB, or split it.
- **CV_estdf_hg is dropped**; do not run it.

-- orion (for the USER)
