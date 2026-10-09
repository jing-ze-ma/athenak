# BSG reproduction page scripts (viper `bsg_1001/figs/scripts_repro`, copied 10-09)

These scripts build the "BSG Reproduction Run" page (https://claude.ai/artifact/KkxYXphcf7mWaajATCd67V)
that compares a run with Ma, Bildsten & Jiang 2026. They are the copies viper ran on 10-08, with
only the edits listed under "Edits in this copy" (env vars whose defaults are the viper values, so
viper's results do not change).

## Files

| file | role |
| --- | --- |
| `update_page_repro.sh [pagedir]` | the refresh driver: bsg_figs -> cost_repro -> ke_compare -> gen_status -> compare57 -> movie/make_movie.sh -> movie_times; copies PNGs/frames/mp4 into `pagedir`, writes `pagedir/new_files.txt` (the list of files to publish) |
| `bsg_figs.py --arm repro` | Figs 2, 3, 4, 6 and the intensity map -> `out/ours_fig{2,3,4,6}.png`, `out/ours_int.png`; per-dump reductions cached in `cache/<basename of run dir>/red_NNNNN.npz`, `int2_NNNNN.npz`; `numbers.json` there; colour ranges in `movie/ranges.json` |
| `bin_convert_repro.py` | AthenaK `.bin` reader used by bsg_figs (a copy of `vis/python/bin_convert.py` that tolerates `=` inside header comments). The `vis/python` path that bsg_figs adds to `sys.path` is not used: it imports this local copy. |
| `cost_repro.py` | wall time, GPU-hours and projection from the run logs (`elapsed=` lines of `run.log`, `run.*.log`) + sacct -> `cache/cost.json` |
| `ke_compare.py` | KE_int(t) panel from `bsg3d.user.hst` of `repro` and `repro_full` -> `out/ours_ke.png`, `cache/ke.json` |
| `gen_status.py pagedir` | all page text and numbers -> `pagedir/status.json` (index.html fills itself from it) |
| `cmp57_dumps.py N [N ...]` | per-dump reductions for the 57.4-d comparison (`cache/cmp57/c_NNNNN.npz`). NOT called by update_page_repro.sh: run it once per new dump number in t = 30-57.41 d (viper ran it on dumps 61..118, 0.5 d apart) before compare57 |
| `compare57.py pagedir` | Figs 5 and 7 analogues + time-averaged L(r) -> `out/ours_cmp57.png`, `cache/cmp57.json`, keys cmp-body/w-cmp/n-cmp added to `pagedir/status.json` (needs the window 30-57.41 d to hold dumps, else it fails; the driver only prints "compare57 FAILED" and goes on) |
| `movie/make_movie.sh` | runs `movie/movie_frames.py --arm repro` (frames for every dump: f3 = phi=0 plane, f4 = r 40.64 Rsun shell, fint = intensity map) and ffmpeg -> `movie/fig3.mp4`, `fig4.mp4`, `figint.mp4` (+ gifs <= 8 MB) |
| `movie_times.py pagedir run` | `pagedir/movie/n.json` (frame count, dump time of each frame) |
| `index.html` | the page itself as published 10-08 (static; not generated; reads `status.json`, `movie/n.json`, `img/*.png`, `movie/*`). Its prose has viper/Raven-specific sentences (one "Raven"); the paper crops `img/paper_fig{2,3,4,6}.png` are copied from the artifact |
| `tess_s96_crop.png` | reference data: TESS sector 96 crop shown in Fig 6 (read by bsg_figs from its own dir) |

Not included: `paper_fig*.png` (on the artifact), `build_merged.py` (viper only: builds the merged
Raven+viper view of the viper run; skip it, see below), the one-off `patch*.py` that built index.html.

## Arguments and environment

All writes go under the script dir (`out/`, `cache/`, `movie/`) and `pagedir`. The run dir is read only.

| env | meaning | viper default | Caltech truerepro2 |
| --- | --- | --- | --- |
| `pagedir` ($1) | page directory | /viper/ptmp2/jinma/pages_1004/bsgrepro | your page dir |
| `BSG_REPRO_RUN` | the run shown | /viper/ptmp2/jinma/bsg_viper_1006/merged_repro | /resnick/groups/carnegie_poc/jingze/bsg_1009/prod |
| `BSG_REPRO_FULL_RUN` | second run in the KE panel only (viper: the unseeded force-work = full run) | /raven/ptmp/jinma/bsg_raven_1002/repro_4n | no twin: set it = `BSG_REPRO_RUN` (both curves coincide; edit the w-ke/n-ke text in gen_status.py) |
| `BSG_BUILD_MERGED` | merged-view builder | viper path | set to the empty string: `BSG_BUILD_MERGED=` |
| `BSG_INP` | input file name inside the run dir (grid stretch keys are checked against it) | bsg3d_repro.athinput | the run's own .athinput name |
| `BSG_IC` | the 1-D IC drawn in Fig 3 | $BSG_DATA_DIR/repro/ic3d/ic_ma2026.txt (md5 21672b25cd32f86bb75d3fd3cbd2506c) | the IC the run used (e.g. ../ic_ma2026_f0.txt from ../ic_ma2026_f0.txt.gz, md5 246a666e...) |
| `BSG_OPAC` | Ma's Rosseland table | $BSG_DATA_DIR/ma2026_files/athenak_tables/rosseland_ma2026_x0.7_z0.008_padded.txt (md5 450bc0c17a27c8e7e725b3ab7fb41587) | ../rosseland_ma2026_x0.7_z0.008_padded.txt (same md5) |
| `BSG_DATA_DIR` | base of the two defaults above | /viper/ptmp2/jinma/bsg_1001 | unused if BSG_IC and BSG_OPAC are set |
| `BSG_WORK` | where out/ cache/ live | the script dir | leave |
| `BSG_COST_PARTS` | JSON file: list of `{name, run, gpu, ngpu, layout, jobs:[...], host:"local"}` replacing the hard-coded viper/Raven parts in cost_repro.py | unset | e.g. `[{"name":"Caltech","run":"<run dir>","gpu":"H200","ngpu":N,"layout":"...","jobs":["4286643","4286644","4286645"],"host":"local"}]` |
| `BSG_RUN_LABEL`, `BSG_HW_LABEL` | run name and hardware text of the settings line in gen_status.py (non-merged runs) | run dir basename, "4 Raven nodes, 16 A100, 16 blocks" | your text |
| `BSG_RECHOOSE` | 1 = re-choose colour ranges from the newest dump and redraw all frames | 1 | 1 |
| `FFMPEG_MODULE`, `FPS` | ffmpeg module (load failure ignored; any ffmpeg on PATH works), movie fps | ffmpeg/7.1, 5 | |

Run-file names the scripts read (truerepro2 uses the same basename `bsg3d`):
`<run>/bsg3d.user.hst`, `<run>/bsg3d.x1grid.txt` (stretched-grid faces), `<run>/<BSG_INP>`,
`<run>/bin/bsg3d.hydro_w.NNNNN.bin`, `<run>/bin/bsg3d.m1.NNNNN.bin` (a dump is used only when both exist),
optional `<run>/bin/bsg3d.m1_vet.NNNNN.bin`, logs `<run>/run.log`, `<run>/run.*.log` (cost).
If the dumps are not under `bin/`, symlink a view dir. `gen_status.py` also reads `<run>/MERGED.txt`
only if it exists (viper merged view; absent at Caltech = the plain branch).

Text that stays viper-specific (edit by hand if wanted): the run description string `RUN` in
gen_status.py ("seeded reproduction run ... 1 % temperature seed, force reference work = split"),
the KE notes (unseeded twin), `TLIM = 4.96e6` s (57.4 d) in cost_repro.py, the KE legend in
ke_compare.py, and some sentences in index.html.

Example (Caltech):
```bash
export BSG_REPRO_RUN=/resnick/groups/carnegie_poc/jingze/bsg_1009/prod BSG_REPRO_FULL_RUN=$BSG_REPRO_RUN
export BSG_BUILD_MERGED= BSG_INP=<run input>.athinput BSG_IC=<IC txt> BSG_OPAC=<rosseland txt>
export BSG_COST_PARTS=$PWD/parts.json BSG_RUN_LABEL=truerepro2 BSG_HW_LABEL='...'
python3 cmp57_dumps.py <dump numbers with 30 <= t <= 57.41 d>    # once dumps exist there
bash update_page_repro.sh /path/to/page   # page dir holds index.html and img/paper_fig*.png
```

## Python dependencies

numpy, matplotlib (Agg), scipy (`scipy.optimize.curve_fit`, Fig 6 PSD fit) and **astropy**
(`astropy.timeseries.LombScargle`, Fig 6 periodogram; imported inside `psd_ls`, so only Fig 6 needs it,
but Fig 6 is on the page and bsg_figs aborts without it). ffmpeg for the movies. h5py is NOT needed:
bin_convert_repro.py imported it at top level only for its HDF5/XDMF writer, which the page never
calls; this copy imports it inside that writer. tqdm is optional (only in bin_convert's `__main__`).
cost_repro.py calls `grep` and `sacct` via subprocess (sacct missing = log wall only).

## Edits in this copy (relative to viper's files)

- bsg_figs.py: `BSG`, `FIGS`, `OPAC`, the repro arms' `inp`/`ic` from env (`BSG_DATA_DIR`, `BSG_WORK`,
  `BSG_OPAC`, `BSG_INP`, `BSG_IC`), defaults = viper values (FIGS default = the script dir, which on
  viper is scripts_repro).
- bin_convert_repro.py: `import h5py` moved into the HDF5 writer.
- gen_status.py: `S` = script dir (or `BSG_WORK`); `BSG_RUN_LABEL`/`BSG_HW_LABEL` in the non-merged branch.
- compare57.py: cache dir `cache/repro_sf_4n` -> `B.cache_dir(RUN)` (on viper `cache/merged_repro` is a
  symlink to `cache/repro_sf_4n`, so same files; at Caltech this is `cache/prod`).
- cost_repro.py: optional `BSG_COST_PARTS` JSON replaces `PARTS`.
- movie/movie_frames.py: `sys.path` = parent of its own dir; movie/make_movie.sh: `M` = its own dir,
  ffmpeg module optional.
- update_page_repro.sh: `S` = its own dir; build_merged.py only when `BSG_BUILD_MERGED` is non-empty.

Checked: `python3 -m py_compile` on every .py and `bash -n` on both .sh (no run of the pipeline here).
md5s of these files: MD5SUMS.
