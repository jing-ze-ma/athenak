# R136a1 deep envelope run (B) page — https://claude.ai/artifact/J8TiLW4e42UMvjWCqQyvRh

Live version id 1794030610-d61a (55 files). Counter: v13 on 11-06 00:xx (2.83 d), several publishes since (last 11-07).
The first-look page https://claude.ai/artifact/AHgSPSmBA8TWhdLHNMtGbN (run128/prod, `refresh.sh`/`make_page.sh`, not copied)
is an archive.

## Files (unchanged copies of `$J/r136a1_1104/page/...`)

- `refresh_deep.sh [full]` — CPU step in the held allocation (`srun --overlap --jobid=$JID -w $ND -N1 -n1 -c 4 --mem=16G nice -n 19
  bash -l make_page_deep.sh $P [--full]`, JID/ND default 5544871 / hpc-sm-02-03), then writes `publish_files_deep.txt` (and a copy as
  `$P/publish_files.txt`).
- `make_page_deep.sh PAGEDIR [--full]` — env: R136_RUN=`$J/r136a1_1104/runB/prod`, R136_GM0=3.1017784815e28,
  R136_IC=`$J/r136a1_1104/ic/ic_r136a1_B.txt`, R136_LABEL, R136_TLIM=1.728e6, R136_BANNERHEAD; modules python/3.11.6, py-scipy 1.11.1,
  py-matplotlib 3.8.2 (Caltech spack names); copies `index_deep.html` -> `$PAGE/index.html`, runs `r136_page.py`, then `hgstyle_figs.py`.
- `r136_page.py` (definitions in its docstring: T from the run's general EOS, kappa_R, tau/R_ph/T_eff, comoving F0 (r136_f0.py thick-limit
  solver), Gamma_rad, Gamma_L, P_rad/P, L_rad/L_conv/L_MLT, Mdot_r) -> `img/{timeseries,profiles,mdot}.png`, movie frames
  `movie/eq/eq_NNNNN.png`, `movie/me/me_NNNNN.png`, `movie/n.json`, `status.json`.
- `r136_f0.py`, `hgstyle_figs.py` (11-05 USER: extra HG/RSG-style line plots `img/hgstyle_{profiles,hist,fluxbudget}.png`; uses
  `pagefast_lib/hstcache.py` from `../common_global/pagefast_lib/`), `index_deep.html` (static page; reads status.json, n.json, img/, movie/),
  `pipe/bin_convert_repro.py` (AthenaK .bin reader), `publish_files_deep.txt` (current list, 55 files).

## Inputs

- `r136.user.hst` columns time, Etot, KE_int, L_top, w_mlt (+ dt, KEr_int, L_in, Mdot_top, M_tot, P_esrc, Picard, P_sponge in
  hgstyle_figs); `run.log` / `run.*.log` (status, cost).
- dumps `bin/r136.hydro_w.NNNNN.bin` + `r136.m1.NNNNN.bin` every 2.16e4 s (0.25 d), `r136.m1_vet.*.bin` (P_rr, matched by time).
- EOS dump `$J/r136a1_1104/ic_cache/tables/eos_8f658626c61e/eos_dump.txt` (hard-coded `EOSF` in r136_page.py; in pages_cache.tar),
  IC `ic/ic_r136a1_B.txt` and opacity `opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt` (both in raven_deps_1107.tar).

## Publish

`index.html` with root `$P=/resnick/home/jingze/ATHENAK/pages/r136a1_deep` and `files` = publish_files_deep.txt minus index.html
(img/*.png, movie/eq/*.png, movie/me/*.png, movie/n.json, status.json). The "movies" are PNG frame sequences played by the page's JS.

## Caches

`cache/runB_prod/` (82 MB: red_00000-00022.npz, hg_00000-00022.npz, ranges.json) -> pages_cache.tar (path
`r136a1_1104/page/cache/runB_prod`), plus runB/prod run logs and the staged `pages/r136a1_deep`. Cache dir name is derived from
R136_RUN (`<parent>_<basename>` = `runB_prod`): on Raven either keep a run dir called `runB/prod` or rename the cache dir.
Dump numbers must continue from 00023 (the cache is keyed by dump number).

## Machine-specific bits

`W` and `PAGE` paths in refresh_deep.sh / make_page_deep.sh; `R` (= `$J/r136a1_1104`), `EOSF` in r136_page.py; `sys.path` to
`$J/pagefast_lib` in hgstyle_figs.py; `module load` names; `MPLCONFIGDIR=$W/cache/mpl`; JID/ND of the held allocation.
Run note: R136a1 B runs with `time2_adm_poslim=true` since 11-06 (name it at every restart).
