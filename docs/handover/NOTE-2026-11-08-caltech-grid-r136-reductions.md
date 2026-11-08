# NOTE 2026-11-08 Caltech -> viper: RSG E vis3d grid + R136a1 B gap reductions sent (md5 OK)

Answers TASK-2026-11-08-caltech-rsg-vis3d-grid and TASK-2026-11-08-caltech-r136-gap-reductions. None of the requested dumps were missing.

| set | files | size | destination | md5sum -c on Raven |
|---|---|---|---|---|
| RSG E vis3d grid | r_05900 05920 05960 05980 06000 06020 06040 06060 06080 06100 06120 06140 06180 06200 06220 06240 06260 06280 06300 (19) | 138 MB | /raven/ptmp/jinma/pages_1107/cache/vis3d/rsgE_glob/red/ + MD5SUMS_grid_1108.txt | 19/19 OK |
| R136a1 B gap | red_00023..00030 + hg_00023..00030 (16) | 29 MB | /raven/ptmp/jinma/pages_1107/r136/work/cache/runB_prod/ + MD5SUMS_gap_1108.txt | 16/16 OK |

- Before sending, none of these files existed on Raven, and tar -k was used, so nothing was overwritten.
- Job 6467126 (expansion, CPU only) ran in 52 s:
  - RSG: gap_1108/red3d.py -> cs3d.load_reduced, REDUCE_VERSION c3, same cache key as the existing cache.
  - R136: gap_1108/r136red.py -> r136_page.reduce_dump + hgstyle_figs.reduce_hg, with the make_page_deep.sh environment. R136 times: 5.75-7.50 d.
- Scripts, log and md5 lists: $J/gap_1108/.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
