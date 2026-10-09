# NOTE 2026-10-09 viper -> Caltech: BSG reproduction page scripts pushed

Answer to TASK-2026-10-09-viper-bsg-page-scripts.md (64a79371).

What: `docs/handover/bsg-hrdet-1009/page/` = viper's `bsg_1001/figs/scripts_repro` page pipeline:
update_page_repro.sh (driver), bsg_figs.py (Figs 2/3/4/6 + intensity map), bin_convert_repro.py (reader),
cost_repro.py, ke_compare.py, gen_status.py, cmp57_dumps.py + compare57.py (57.4-d comparison, Figs 5/7
analogues), movie/make_movie.sh + movie/movie_frames.py, movie_times.py, the published index.html
(static, fills from status.json) and tess_s96_crop.png (Fig 6 reference). README.md lists every
argument/env var, the run files read (bsg3d.user.hst, bsg3d.x1grid.txt, bin/bsg3d.{hydro_w,m1}.NNNNN.bin,
run*.log) and a Caltech example. MD5SUMS in the dir (bsg_figs.py 0af37d9f, gen_status.py 582f7ca7,
update_page_repro.sh a91019a5, index.html 9fca6260, README.md 554b8aca).

Edits vs viper: every viper path is an env var with the viper value as default (BSG_REPRO_RUN,
BSG_REPRO_FULL_RUN, BSG_IC, BSG_OPAC, BSG_INP, BSG_DATA_DIR, BSG_WORK, BSG_BUILD_MERGED, BSG_COST_PARTS,
BSG_RUN_LABEL, BSG_HW_LABEL); script dirs are found from `$0`/`__file__`; compare57 reads
cache/<run basename> instead of the hard-coded cache/repro_sf_4n; h5py import made lazy (only the
HDF5 writer, never used by the page). Checked: py_compile + bash -n only; the pipeline was not run.

Caveats:
- Fig 6 needs astropy (LombScargle); install it (pip --user) or Fig 6 / bsg_figs fail. h5py not needed.
- cmp57_dumps.py is not in the driver: run it on the dump numbers with 30 <= t <= 57.41 d first;
  until then compare57 fails (non-fatal in the driver).
- KE panel compares with a second run (viper: unseeded twin); set BSG_REPRO_FULL_RUN = BSG_REPRO_RUN and
  edit the w-ke/n-ke text in gen_status.py.
- Prose in gen_status.py (RUN string: seed, split, ...), cost_repro.py PARTS/TLIM and index.html still
  describe the viper/Raven run; set BSG_COST_PARTS / labels and edit text by hand.
- IC for Fig 3: point BSG_IC at the IC the run used (../ic_ma2026_f0.txt md5 246a666e differs from
  viper's ic_ma2026.txt md5 21672b25); BSG_OPAC = ../rosseland_ma2026_x0.7_z0.008_padded.txt (same md5
  450bc0c1 as viper's).
