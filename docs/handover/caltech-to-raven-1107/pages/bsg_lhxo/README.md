# BSG prod_lhxo page ("BSG Athena++-equivalent (his deck values)") — https://claude.ai/artifact/7mZFDbuJtsh96T8s64xMWF

THIS is the live page of the BSG production prod_lhxo (`$J/doface_1105/prod_lhxo`, Linhao deck values + vstar_nbe, GM 19.78, no
FOFC, CFL 0.3, finished 57.4 d at 11-07 06:14). Page v24 (last build 11-07 06:21: "OK | 57.40 d | wall 15.4 h", 116 frames,
575 files, 115.9 MB). The page dir `doface_1105/page` was cloned from `bsgli_1104/page_none` on 11-05 and repointed several times
(prod_doface -> prod_vsnsink -> prod_vsnmlat_ns -> thr1e7 -> cfl035 -> prod_lhxo on 11-06; see the make_page.sh header and the
`*.pre_*` backups not copied). The URL keeps its old title.

## Files (unchanged copies of `$J/doface_1105/page/...`)

- `refresh.sh [full]` — `srun --overlap --jobid=5719508 -w hpc-sm-01-05 -N1 -n1 -c 4 --mem=60G nice -n 19 bash -l make_page.sh $P`
  then writes `publish_files.txt` (index.html, status.json, movie/n.json, img/*.png, movie/*.mp4, movie/*/*.png). `full` sets
  BSG_RECHOOSE=1 (redraw all frames; colour limits are fixed constants read off Ma+2026 since 11-06, so it only redraws).
- `make_page.sh [pagedir]` — env + steps: modules (python 3.11.6, py-scipy 1.11.1, py-matplotlib 3.8.2) + venv
  `$J/bsg_1009/pagework/venv` (system site + astropy, for Fig 6), ffmpeg 6.0 on PATH; BSG_REPRO_RUN=BSG_REPRO_FULL_RUN=prod_lhxo,
  BSG_INP=`doface_1105/prep3d/bsgLI_lhx_exact_o.athinput`, BSG_IC=`ic/ic_linhao_7col.txt` (made by `tools/make_ic_linhao.py` from
  `bsgli_1014/files/linhao_profile_full.txt` if missing), BSG_OPAC=`bsgli_1014/files/rosseland_ma2026_x0.7_z0.008_padded.txt`,
  BSG_COST_PARTS=parts.json, labels. Steps: `pipe/cmp57_dumps.py` for new dumps with 30 <= t <= 57.41 d -> `pipe/cmp_arms.py`
  (img/ours_arms.png: this run vs the 'none' production bsgli_1104/prod and earlier arms) -> copy `pipe/index.html` ->
  `pipe/update_page_repro.sh $PAGE` -> `tools/dark_paper_fig4.py` -> `tools/inject_status.py`.
- `pipe/update_page_repro.sh` driver: bsg_figs.py --arm repro (Figs 2, 3, 4, 6, intensity map; per-dump caches) -> cost_repro.py ->
  ke_compare.py -> gen_status.py (all page text -> status.json) -> compare57.py (Figs 5/7 analogues, 30-57.4 d window) ->
  movie/make_movie.sh (movie_frames.py: f3 = phi=0 plane, f4 = r 40.64 Rsun shell, fint = intensity map; ffmpeg) -> movie_times.py.
  Full description of every script and env var: `pipe/README.md` (viper's original notes; still accurate apart from the paths).
- `pipe/index.html` (static page template that fills itself from status.json), `pipe/MD5SUMS`, `pipe/paper_fig6_src.png`,
  `pipe/tess_s96_crop.png` (Fig 6 reference data), `pipe/movie/{make_movie.sh,movie_frames.py,ranges.json,cadence.json}`,
  `tools/{build_template.py,dark_paper_fig4.py,inject_status.py,make_ic_linhao.py,refresh_text.sh}`, `lhxo_text_1106.py` (one-off text
  patch applied to gen_status/index on the lhxo repoint, kept for reference), `dtmap.py` (dt map helper), `parts.json`.
- Paper crops `img/paper_fig{2,3,4,5,6,7}.png` come from the artifact/staging dir (in pages_cache.tar under pages/bsg_doface/img).

## Inputs

- `prod_lhxo/bsg3d.user.hst` (KE_int, KEr_int, time; dt 1000 s), `bsg3d.x1grid.txt`, the deck, logs `run.log`/`run.*.log` (cost).
- dumps `bin/bsg3d.hydro_w.NNNNN.bin` + `bsg3d.m1.NNNNN.bin` every 0.5 d (116 dumps, 0-57.4 d), optional `bsg3d.m1_vet`.
- cmp_arms.py also reads the hst + per-dump reductions of the other arms (bsgli_1104/prod and the doface arms); their reductions are in
  `pipe/cache/arms/*.npz` (shipped), their hst files are NOT in the deps tar (only prod_lhxo's is): cmp_arms will fail on Raven
  unless those hst files are copied or the arms list is cut (it is optional: make_page.sh prints "cmp_arms FAILED" and continues).

## Caches (pages_cache.tar)

`doface_1105/page/pipe/cache/{prod_lhxo (730 MB: red_/int2_ per dump, numbers.json), cmp57 (44 MB), arms (2.7 MB), *.json, *.log,
cmp_arms.txt}`, `pipe/movie/frames` (80 MB, needed: the dumps are not on Raven), `pipe/out`, `doface_1105/page/ic`, prod_lhxo run
logs, staged `pages/bsg_doface`. Not shipped: caches of the abandoned arms (prod_doface, prod_vsnmlat_*, ~1 GB) and old frame dirs
`movie/frames_*_1106`.

## Publish

`index.html` with root `$P=/resnick/home/jingze/ATHENAK/pages/bsg_doface` and `files` = publish_files.txt minus index.html (575 files,
115.9 MB -> 2 publishes < 64 MB each; unchanged files are kept, so after the first full publish only new_files.txt is needed).

## Machine-specific bits

`W`, `P`, `BSG_*` paths, venv path, ffmpeg PATH, module names, `--jobid=5719508 -w hpc-sm-01-05`, `J` in cmp_arms.py,
`MPLCONFIGDIR=$W/pipe/cache/mpl`, BSG_HW_LABEL text (Caltech 4 H200).
