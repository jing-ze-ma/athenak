# Shared pipeline of the RSG E and HG global pages (cubed sphere)

Used by `../rsgE/` (page key `rsgE`) and `../hggl/` (page key `hggl`). Copies of `$J/<file>` unless noted.

| file | role |
|---|---|
| `pages.sh [status]` | one command for all live pages: `quickpage.sh rsgE`, `quickpage.sh hggl`, dhj_r1 `refresh.sh` |
| `quickpage.sh rsgE\|hggl [status\|text\|nocatch\|catchup]` | FAST text rebuild (~15 s): RSG `pageH/sync.sh`, `mk_figs.py` + `build_page.py` from hst + logs + already-reduced dumps -> `$P/<n>_glob/index.html`; then if dumps are unreduced, launches `catchup_job.sh` DETACHED in a held allocation from `pages_alloc.cfg` (`NSH=4 setsid nohup srun --overlap --jobid=<alloc> -w <node> -c 4 --mem=60G nice -n 19`), or sbatch if none runs |
| `catchup_job.sh rsgE\|hggl` | reduce new dumps (`reduce.py`, ONLY=<n>, up to 24 procs x RED_THREADS), 2-D movies (`mk_movies.py` limits pass -> NSH render shards -> encode/append parts), rsync movies to staging, `r3d_par.sh` (3-D), `quickpage.sh nocatch`, writes `$P/<n>_glob/READY` (line 1 = time + sim t; lines 2.. = files to publish) |
| `r3d_par.sh rsgE\|hggl` | 3-D movies: `vis3d/cs3d.py <run> --view single` then `--view layer` (24 workers, RED_CAP 16), copies `*_3d*.mp4` + posters to staging |
| `vis3d/cs3d.py` | 3-D renderer (vtk-osmesa 9.3.1 via `venv_vtk/pvpy.sh`). Per-dump reduction cache `vis3d/<run>/red/r_NNNNN.npz`, frames `vis3d/<run>/frames{,_layers}/f_NNNNN.{png,json}`, sticky colour limits `clim{,_layers}.json`, output `vis3d/<run>/out/<movie>{,_layers}.mp4/.json/_poster.png` -> copied to `<page>/out/publish/movies`. MAXFRAMES 300 (dumps subsampled evenly). |
| `vis3d/render_cs.sh`, `vis3d/README_cs3d.txt` | encode-only job and the renderer notes |
| `pagefast_lib/{figcache,hstcache,eosnb,ppcache}.py` | exact caches imported by `mk_figs.py` (0-d npz members cached in `out/data/.figcache_scal.pkl` keyed by name+size+mtime_ns; incremental hst parse with restart filter), numba EOS row bisection used by `eos_gen.py` (EOS_NUMBA=0 = numpy) |
| `vispy_4de8168d/bin_convert.py` | AthenaK `vis/python/bin_convert.py` at commit 4de8168d (the newer reader rejected bin format 1.2); `VISPY` points at this dir |
| `venv_vtk/pvpy.sh` | runs `$J/venv_vtk/bin/python` with a clean PYTHONPATH |
| `pagesync.sh`, `parred.sh`, `movies_par.sh`, `rereduce_job.sh`, `parpage.sh` | older/full-refresh helpers (pagesync = 3-D job -> page job; parred = fan-out reduce 8 dumps/job; movies_par = 2-D movies on one 32-core node; rereduce_job = full re-reduce after a reduce.py VER bump). Not used in the normal path any more. |
| `pages_alloc.cfg` | page -> held allocation/node for the in-alloc steps (Caltech job ids; replace with Raven ones) |

## Normal refresh (what Caltech did, ~1-30 min)

1. `bash quickpage.sh rsgE` (and `hggl`): text page now; catch-up detached for new dumps.
2. Publish `index.html` immediately (`url=` explicit).
3. When `quickpage.sh <n> status` says `READY`, publish `index.html` + the files in lines 2.. of `$P/<n>_glob/READY`
   (RSG: 2 calls to stay < 64 MB: index + 2-D parts, then index + 3-D files).

## Machine-specific bits to repoint

- `$J=/resnick/groups/carnegie_poc/jingze` hard-coded in: quickpage.sh, catchup_job.sh, r3d_par.sh, pages.sh,
  pagesync.sh, parred.sh, movies_par.sh, rereduce_job.sh, cs3d.py (`VIS`, `G`, `RUNS[...]['hydro'/'table'/'publish']`,
  `sys.path` to `rsgE_page/vispy_4de8168d` and `rsgE_page/src`), mk_figs.py / eos_gen.py (`sys.path` to pagefast_lib),
  reduce.py (`_HCF` = pagefast_lib/vet_hdr_t.json), mk_movies.py (`MATCH3D` default = vis3d/<run>/out/<movie>.json),
  build_page.py (footer text, rsgE `section_sg` reads `$J/rsgE_supergran_1101/fig*.png|jpg` -> committed in
  `../rsgE/static_inputs/`).
- Staging dir `$P=/resnick/home/jingze/ATHENAK/pages/<n>_glob` (quickpage.sh, catchup_job.sh, r3d_par.sh).
- Python: `$J/venv_vtk/bin/python` (python 3.11 venv: numpy, scipy, matplotlib, h5py, PIL, vtk-osmesa 9.3.1 for
  cs3d, numba 0.68 in `pagefast_lib/pylib` added to sys.path by eos_gen/eosnb). Any python with these works;
  vtk must be an OSMesa (offscreen) build. matplotlib backend Agg (`MPLBACKEND=Agg` exported); default fonts
  (DejaVu Sans), no custom font files.
- ffmpeg: `/central/software9/spack/.../ffmpeg-6.0-.../bin/ffmpeg` (update.sh `FFMPEG`, cs3d.py `FFMPEG`); libx264 needed.
- Slurm: `srun --overlap --jobid=<held alloc> -w <node> -N1 -n1 -c 4 --mem=40-60G nice -n 19` (Caltech held
  allocations from `pages_alloc.cfg`); `sbatch -A carnegie_poc -p expansion -c 32 --mem=200G -t 01:30:00` header in
  catchup_job.sh (fallback only). `squeue` is used for the "jobs" status line (`JOBNAME` env in update.sh).

## Required adaptation on Raven (the old dumps will not exist there)

- 2-D path is OK with the shipped caches: `reduce.py` only reduces dump numbers that have no `out/data/red_NNNNN.npz`
  (`new()` compares numbers), and `mk_figs.py` / `mk_movies.py` read all `red_*.npz` (Caltech history + new). Keep
  the mtimes when copying (figcache key). New Raven dumps must continue the dump numbering (they do, when resuming
  from the rst).
- 3-D path needs a small patch: `cs3d.py dump_list()` builds the frame list from the dump FILES
  (`RUNS[run]['hydro']` glob) and keys the reduction cache on realpath+size+mtime of the Caltech dumps. On Raven the
  old dumps are absent, so the 3-D movie would only show Raven dumps. Patch: add the numbers from
  `vis3d/<run>/red/r_NNNNN.npz` that have no dump, reusing the stored `key` (z['key'] minus the
  `|REDUCE_VERSION|table...` suffix), so `load_reduced` hits the cache. Frames re-render once (the cfg hash includes
  the repointed paths), from the cached reductions (cheap: ~2.6 s/frame).
- RSG `pageH/sync.sh` stitches prodG (< 1675 d) + prodH dumps/hst into `rsg_E_1020/pageH`; on Raven point
  `RUN/BINDIR/EXTRA_RUN/LOGDIR` at the Raven run and append its hst to the shipped stitched `pageH/rsg.user.hst` /
  `rsg.hydro.hst` (in pages_cache.tar), or extend sync_links.py with the Raven dir.
