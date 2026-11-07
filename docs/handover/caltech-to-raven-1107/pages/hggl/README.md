# He Giant Global Caltech page — https://claude.ai/artifact/LNtjtaL34oUhgcCRNp5V2M

Live version id 1794010989-9ef6 (11-07); last staged build 11-06 17:40, sim t = 4.981e7 s = 576.46 d (READY).
Same pipeline as RSG (adapted copy): read `../common_global/README.md` and `../rsgE/README.md`; differences only below.

## Files here (unchanged copies)

- `hgglob_page/update.sh` — env block: RUN=EXTRA_RUN=LOGDIR=`$J/hgglob_1021/prodG`, BINDIR=`prodG/bin`, BASENAME=hegiant,
  ICF=`hgglob_1021/files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt` (in raven_deps_1107.tar), RPH0=7.94653e12, RIN=5.5656e11 (inner
  wall, dashed in zooms), TFLOOR=2400, JOBNAME=hgGD,hgGS,hgGH3,hgGH1,hgGB, AGE=120, VISPY=`rsgE_page/vispy_4de8168d` (shared reader).
- `hgglob_page/src/`: `reduce.py`, `mk_figs.py`, `mk_movies.py`, `build_page.py`, `csplot.py`, `eos_gen.py`, `page.css` (no lomb.py,
  no add_shell.py, no supergran/modes sections).
- `hgglob_page/README.md` (Caltech README of this dir), `env_only.sh`.
- One-off tools for the rewind trap: `fix_stale_1106.sh`, `fix_stale_movies_1106.{py,sh}` (move caches of reused dump numbers aside
  after a restart from an earlier rst and rebuild the affected movie parts), `hg_reenc.sh` (full 2-D re-encode).

## Commands

`bash $J/quickpage.sh hggl` (text + detached catch-up) -> when READY publish `index.html` + READY files in ONE call (HG < 64 MB):
2-D parts `movies/hggl_{eq,inner,map,phot}_pNN.mp4` that changed + `movies/hggl_3d.mp4`, `hggl_3d_layers.mp4`, `hggl_3d_poster.png`,
`hggl_3d_layers_poster.png`.

## Inputs

- hst `hegiant.user.hst` (dt 25 s): same columns as RSG (time, dt, Etot, KE_int, KEr_int, L_bot, L_in, L_top, Mdot_top, M_tot,
  P_esrc, Picard, P_sponge); `hegiant.hydro.hst`.
- dumps `hegiant.hydro_w` + `hegiant.m1` (deck: 1 d), `hegiant.m1_vet` (5 d; matched by header time, "every 15th" hydro dump in
  practice); logs `run*.log`, `link.*.out`, `inalloc.out` in prodG.
- EOS composition X 0.0009496 / Y 0.97752 (eos_gen, cs3d), opacity `hgglob_1021/files/hg145-prod-1018/rosseland_mesakap_afe2dd8da2bd.txt` (cs3d).

## Caches

- bulk list: `hgglob_page/out/data/red_*.npz` (1730, 33 GB) + `.figcache_scal.pkl`; `vis3d/hggl_glob/red` (974, 8.4 GB).
- pages_cache.tar: `hgglob_page/out/frames` state (crf/stamp/txt/json/segs/limits.json), `out/img`, `out/publish`, `out_quick`,
  `vis3d/hggl_glob/{clim*.json,eos_fe.npy,out/}`, prodG run logs, staged `pages/hggl_glob`.
- NOT shipped: `out/stale_*`, `out/old_prodF_1022` (rewound/old runs, not shown on the page).
- Optional PNG frames 6 GB (bulk list).

HG note: prodG was rewound twice (11-02 strip at 315 d, 11-05/11-06 rewinds); the shipped caches are the cleaned ones. If Raven
resumes from rst 139 (5.0112e7 s, 580 d) dumps after 576.46 d may reuse numbers that exist in `out/data`: check `t` in the last
red npz against the first Raven dump and move overlapping numbers aside (see fix_stale_1106.sh).
