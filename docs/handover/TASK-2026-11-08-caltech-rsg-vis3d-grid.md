# TASK 2026-11-08 viper -> Caltech (USER): 19 RSG E vis3d reductions for an exact 3-D movie grid

Thank you for the gap reductions (NOTE 943dba26a): md5 OK. All RSG E and HG movies are now continuous across the join.

The user wants every page movie to run at a constant pace in simulated days per video second, with no jitter. The 2-D
movies already do (one frame per dump; RSG 1 d, HG 1/3 d). The RSG E 3-D movie (`rsgE_3d*.mp4`, vis3d `cs3d.py`)
now uses exactly the dumps N = 0, 20, 40, ... (20 d per frame). The HG 3-D grid (every 6 dumps = 2 d) is already complete.

Your RSG E vis3d cache (`r_NNNNN.npz`) is sparse after ~5000 d (1-21 d between cached dumps). These 19 grid dumps
have no r_ file, and their raw dumps exist only on Caltech (`/resnick/groups/carnegie_poc/jingze/rsg_E_1020/prodH/bin`,
all before our restart at 6325.000 d):

```
05900 05920 05960 05980 06000 06020 06040 06060 06080 06100 06120 06140 06180 06200 06220 06240 06260 06280 06300
```

Please compute exactly these r_ files with your vis3d reduction. Use the same code and settings as the existing cache:
`cs3d.py` run `rsgE_glob`, `REDUCE_VERSION = 'c3'`, same opacity table and composition, i.e. `cs3d.load_reduced(run, n, fh,
key, cfg)` as in the reductions you already made. Then rsync them to

- `/raven/ptmp/jinma/pages_1107/cache/vis3d/rsgE_glob/red/r_NNNNN.npz`

with `MD5SUMS_grid_1108.txt` next to them (md5sum of the 19 files). Do not overwrite any r_ file that is already there.
If any of these dumps no longer exists, list which, and we will pick a coarser exact grid instead.
Reply with a NOTE on this branch when done. Our 3-D grid dumps from 06340 on are computed on Raven.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
