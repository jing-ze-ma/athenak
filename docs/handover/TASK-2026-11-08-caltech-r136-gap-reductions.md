# TASK 2026-11-08 viper -> Caltech (USER): reduce R136a1 B dumps 00023-00030 for the page movies

The R136a1 B page continues your run with our Raven segment (restart rst 8 at 7.51 d; our first dump is 00031 = 7.75 d).
The page sliders (eq, me) have a hole: your dumps 00023-00030 (5.75-7.5 d, every 0.25 d) were never reduced into the
page cache, and they exist only on Caltech. The user wants complete movies, with no held frames.

Please reduce exactly these 8 dumps with your page pipeline (`r136_page.py` / `hgstyle_figs.py` reductions, same
settings as the existing cache), giving for each N in 00023..00030:

- `red_000NN.npz` and `hg_000NN.npz`

and rsync them to `/raven/ptmp/jinma/pages_1107/r136/work/cache/runB_prod/` with `MD5SUMS_gap_1108.txt` (md5sum of
the 16 files) in the same directory. Do not overwrite any file that is already there (00000-00022 and 00031 are present).
If any of these dumps no longer exists, list which. Reply with a NOTE on this branch when done.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
