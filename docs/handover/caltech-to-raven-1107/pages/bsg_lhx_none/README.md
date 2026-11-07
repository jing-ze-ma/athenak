# BSG Linhao-exact (production) page — https://claude.ai/artifact/GyzAhm3DR7F81DyaBfvL2n

ARCHIVE page of the Linhao-exact 'none' run `$J/bsgli_1104/prod` (inner-face dissipation none, finished 57.4 d on 11-05, 115 dumps).
Last staged build 11-05 23:36, last publish 11-06. It is NOT the prod_lhxo page (that is 7mZFDbuJ, `../bsg_lhxo/`).
Sister test page (diss athena, `page_athena`, https://claude.ai/artifact/64r69zrYXGL5nB9aBmTAe7) not copied (same scripts, different run).

## Files (unchanged copies of `$J/bsgli_1104/page_none/...`)

Same pipeline as `../bsg_lhxo/` (that dir was cloned from this one), minus cmp_arms.py, plus:
- `refresh.sh [full]` also runs `bsgli_1104/run_cmp_mdot.sh` (srun in alloc 5719508: `cmp_mdot.py` -> `cmp_mdot.png`, copied to
  `img/ours_mdot_cmp.png`); copies in `bsgli_1104/`.
- `make_page.sh`: BSG_REPRO_RUN=`bsgli_1104/prod`, BSG_INP=`bsgli_1104/bsgLI_exact.athinput`, labels bsgLX-none; srun -c 8 --mem=64G.
- `merge_run.sh`, `upd_page.sh` (older drivers), `pipe/` (bsg_figs, cmp57_dumps, compare57, cost_repro, gen_status, ke_compare,
  movie_times, update_page_repro.sh, index.html, README.md, MD5SUMS, tess crop), `pipe/movie/`, `tools/`.

## Inputs / caches / publish

As bsg_lhxo, with run `bsgli_1104/prod` (hst/dumps every 0.5 d). Caches in pages_cache.tar: `bsgli_1104/page_none/pipe/cache/{prod
(730 MB), cmp57, *.json, *.log}`, `pipe/movie/frames`, `pipe/out`, `bsgli_1104/cmp_mdot.{png,txt}`, prod run logs, staged
`pages/bsg_lhx_none`. Publish index.html with root `$P/bsg_lhx_none` + publish_files.txt minus index.html.
Note: `bsgli_1104/prod/bsg3d.user.hst` and the deck are not in raven_deps_1107.tar; only needed if this page is ever rebuilt.
