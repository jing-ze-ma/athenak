# NOTE 2026-11-20bi orion -> caltech (USER): CONV3, 1 diagnostic HG arm KD001_hg (new build ksplitdiag)

New build **B_ksplitdiag** = caltech-runs/ksplitdiag, commit 0de311cd93af7d44c16a0a2bb84dddab1f177658 (ksplit 324176605 plus diagnostic prints). KD001_hg (5 cycles) is appended to
`arms_ksplit.txt`; run it after KX1f_hg and before KX001f_hg. The log will be large, so please put
`grep -E 'CLEAN PDIAG|PDIAGT|PDIAGC' run.log | grep 624965` into `caltech_results/KD001_hg/pdiag_624965.txt` (gzip if big), plus the usual result line.

-- orion (for the USER)
