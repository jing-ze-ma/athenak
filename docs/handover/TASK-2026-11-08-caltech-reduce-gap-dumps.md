# TASK 2026-11-08 viper -> Caltech (USER: "ask caltech"): reduce the RSG E and HG gap dumps for the page movies

The RSG E and HG pages now continue your runs with our Raven segments. The history (hst) joins cleanly, but the
movies and slice/profile plots have a gap: these dumps were never reduced on Caltech, and they exist only on Caltech.

| run | dumps to reduce | times | where (Caltech) | restart we continued from |
|---|---|---|---|---|
| RSG E | 06304-06325 (22) | ~6304-6325 d | /resnick/groups/carnegie_poc/jingze/rsg_E_1020/prodH/bin | rsg.01295.rst, t = 546480000 s = 6325.000 d |
| HG global | 01730-01740 (11) | ~576.7-580 d | the HG global run's bin dir | hegiant.00139.rst, t = 50112000 s = 580.000 d |

Please reduce them with exactly the pipeline that made the existing caches (quickpage / reduce.py for red_NNNNN.npz,
and the vis3d reduction for r_NNNNN.npz), same settings, and rsync the outputs to Raven with md5 lists:

- RSG E red: `/raven/ptmp/jinma/pages_1107/cache/rsgE_page/out/data/red_NNNNN.npz`
- RSG E vis3d: `/raven/ptmp/jinma/pages_1107/cache/vis3d/rsgE_glob/red/r_NNNNN.npz` (r_06304 is already there; 06305-06325 needed)
- HG red: `/raven/ptmp/jinma/pages_1107/cache/hgglob_page/out/data/red_NNNNN.npz`
- HG vis3d: `/raven/ptmp/jinma/pages_1107/cache/vis3d/hggl_glob/red/r_NNNNN.npz`

Put `MD5SUMS_gap_1108.txt` next to each set. If any of these dumps no longer exists, list which, so the page can say so.
Do NOT send the raw dumps (too large); only the reductions. Reply with a NOTE on this branch when done.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
