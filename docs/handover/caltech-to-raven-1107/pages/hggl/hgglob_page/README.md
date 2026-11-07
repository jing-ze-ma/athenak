# He Giant Global Caltech page (hgglob_page)

Results page for the He giant GLOBAL cubed-sphere run `/resnick/groups/carnegie_poc/jingze/hgglob_1021`
(deck `hgglob_fresh_clip.athinput` md5 06001321, binary 4de8168d md5 9a56b4a1, basename `hegiant`, squeue job name `hgG`).
The pipeline is a copy of `/resnick/groups/carnegie_poc/jingze/rsgE_page` (RSG E page) adapted to this run. It only READS the run dir.

## Rebuild

Run it as a CPU job on the expansion partition (not on the login node):

    cd /resnick/groups/carnegie_poc/jingze/hgglob_page
    sbatch -A carnegie_poc -p expansion -c 8 --mem=32G -t 01:30:00 -J hgGpage \
      -o /resnick/groups/carnegie_poc/jingze/hgglob_page/out/logs/page.%j.out \
      --wrap "bash /resnick/groups/carnegie_poc/jingze/hgglob_page/update.sh"

It is incremental: only new dumps are reduced (`out/data/red_NNNNN.npz`) and only new movie frames are rendered
(`out/frames/`). A dump younger than AGE (120 s) is skipped because it may still be being written. A lock file keeps
two builds from running at once.

## Output to publish

Publish the whole `out/publish/` directory: `index.html` (figures inlined as base64) and `movies/*.mp4`
(`hggl_eq`, `hggl_inner`, `hggl_phot`, `hggl_map`). The page loads the movies by the relative path `movies/<name>`.

Optional 3-D section: if `out/publish/movies/hggl_3d.mp4` exists when `build_page.py` runs, the page gets a "3-D view"
section with that movie (poster `movies/hggl_3d_poster.png` if present). The 3-D renders come from the separate
pipeline in `/resnick/groups/carnegie_poc/jingze/vis3d`: copy them into `out/publish/movies/` and rebuild (or just rerun
`build_page.py`, see below). Without the file the section is left out. The RSG E page does the same with `rsgE_3d.mp4` /
`rsgE_3d_poster.png`.

Only the page (seconds, after the movie files are in place):

    cd /resnick/groups/carnegie_poc/jingze/hgglob_page/out
    RUN=/resnick/groups/carnegie_poc/jingze/hgglob_1021/prod /resnick/groups/carnegie_poc/jingze/venv_vtk/bin/python \
      ../src/build_page.py img publish/index.html "page rebuilt $(date '+%Y-%m-%d %H:%M %Z')"

## Steps (src/)

- `reduce.py`: per dump (hydro_w + m1): shell means, L(r), T_gas from the EOS, R_ph / T_eff per column, native-cell
  equatorial and meridional rings, shells at 0.1, 0.2, 0.3, 0.6, 0.9, 0.98, 1.02 R_ph for the maps.
- `eos_gen.py`: numpy port of the code's general EOS with this deck's keys (X 0.0009496, Y 0.97752, a_metal 1.8047, H2 off,
  ionisation on, metal ionisation off, radiation off; table logd -21..-1, logT 3.3..6.8, dlog 0.02). Check:
  `python src/eos_gen.py <IC column>` gives mu 1.339625 and n_e/n_tot 0.655591 at the grid centre (run.log: 1.33962,
  0.655591) and reproduces the IC column T to max 1.5e-4 (median 3e-6).
- `mk_figs.py`: history (`hegiant.user.hst`), radial profiles (with the IC column), slices, photosphere band, Mollweide
  maps, status table (squeue `-n hgG`, run.log, link logs, `hegiant.log` event counters).
- `mk_movies.py`: the four movies, colour limits pooled over all dumps.
- `build_page.py`: the page text (star, grid, physics, seed, binary/deck, layout) and the assembly.

## Settings (update.sh, all overridable by env)

RUN / BINDIR / EXTRA_RUN = hgglob_1021/prod(/bin), LOGDIR = hgglob_1021 (link.*.out), ICF = the IC column
`files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt`, RPH0 = 7.94653e12 cm (IC header), RIN = 5.5656e11 cm (x1min, dashed
circle in the zooms), TFLOOR = 2400 K (deck tfloor_kelvin), BASENAME = hegiant, JOBNAME = hgG, PY = venv_vtk python,
VISPY = rsgE_page/vispy_4de8168d (bin reader, format 1.2), OUT = hgglob_page/out.
If the star, deck or binary change, edit the kv text in `src/build_page.py` and the EOS keys in `src/reduce.py`.
