# Caltech BSG true-repro page pipeline (bsg_1009/pagework), scripts only

Pushed 10-12 for TASK-2026-10-12-caltech-pagework (DeltaAI). Built the page https://claude.ai/artifact/BQ8LZRU8uV1Zw83GtZQQxK
(v1 10-09 at t 4.9 d, last refresh t 6.0 d) from the Caltech unseeded run bsg_1009/prod_noseed. Base = viper copy
`docs/handover/bsg-hrdet-1009/page` (see its README for every script's role). Layout here = Caltech dir layout:
`make_page.sh`, `parts.json` at top, the pipeline in `pipe/`.

## Files
- `make_page.sh` - driver (CPU node: `sbatch -A carnegie_poc -p expansion -c 8 --mem=64G -t 02:00:00 make_page.sh [pagedir]`).
  Sets the env (BSG_REPRO_RUN, BSG_INP, BSG_IC, BSG_OPAC, BSG_COST_PARTS, BSG_RUN_LABEL, BSG_HW_LABEL), runs cmp57_dumps.py on
  every new dump with 30 <= t <= 57.41 d (not cached), copies pipe/index.html into pagedir, then `pipe/update_page_repro.sh pagedir`.
  Caltech paths and modules (python 3.11 + scipy + matplotlib modules, a venv adding astropy, ffmpeg path) -> edit for DeltaAI.
- `parts.json` - cost_repro input (run dir, GPUs, job ids). NB it names bsg_1009/prod and jobs 4286643-45 (the seeded 10-09 chain);
  replace with the DeltaAI bsgprod run dir / job ids.
- `pipe/` - all scripts; unchanged vs base: update_page_repro.sh, bin_convert_repro.py, cost_repro.py, cmp57_dumps.py, compare57.py,
  movie_times.py, movie/make_movie.sh, movie/movie_frames.py (tess_s96_crop.png not repeated: take it from the base dir).

## Edits vs bsg-hrdet-1009/page (all marked "Caltech 10-09" in the code)
- `bsg_figs.py`: (1) skip a dump modified < 3 min ago (live run, half-written file); (2) hst reader drops a partially written last
  row; (3) Fig suptitle says "true reproduction run, no seed, no well-balancing, M1/VET".
- `gen_status.py`: default HW label = 1 Caltech node / 2 H200; RUN text describes the true reproduction deck (paper IC as is,
  F_r = 0 above r_photo, no seed, full M1 force, Ma zero-contact inner ghosts, Ma outflow top, no sponge, vet_gd); new "exactly
  laterally uniform" stage (all of vspread, rel_rho, T spread == 0); "grew from round-off" wording; a lateral-symmetry sentence from
  ke.json (lat_frac_max < 1e-6).
- `ke_compare.py`: second-run curve only if it is a different run dir; main label from BSG_KE_LABEL (default "true reproduction run
  (no seed)"); lateral KE (KE_int - KEr_int, user.hst) drawn when > 1e-10 KE_int, else a "laterally symmetric" text; writes
  lat_frac_max / ker_ratio_now / ke_05d / ke_min_after1 / t_min_after1_d into cache/ke.json (used by gen_status.py).
- `index.html`: the published Caltech page (92 changed lines vs viper's index.html: title, run/hardware prose, Caltech section,
  no Raven/viper sentences).
- **For bsgprod (seed 1e-6, he_seed_rad true, final method):** the "no seed"/"no well-balancing" wording in bsg_figs.py suptitle,
  gen_status.py RUN and BSG_KE_LABEL default is wrong for that run - set BSG_KE_LABEL / BSG_HW_LABEL and edit those two strings.
- Not pushed: venv/, pipe/cache, pipe/out, pipe/movie frames/mp4/gif, generated movie/ranges.json + cadence.json, logs.
