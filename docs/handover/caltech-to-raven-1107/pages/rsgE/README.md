# RSG E Caltech page — https://claude.ai/artifact/LBASg2iK2ArcLY1V2yBCAR

Live version id 1794016264-c620 (11-07); last staged build 11-06 17:48, sim t = 5.451e8 s = 6309.2 d (READY).
Shared scripts (quickpage/catchup/r3d_par/cs3d/pagefast_lib/vispy reader): `../common_global/` — read its README first.

## Files here (unchanged copies)

- `rsgE_page/update.sh` — env block (sourced by quickpage/catchup via `eval "$(grep -E '^(export|HERE=|R=)' update.sh)"`) and the
  slow full rebuild (reduce -> mk_figs -> mk_movies -> build_page). Defaults: RUN=`$J/rsg_E_1020/pageH`, BINDIR=`pageH/bin`,
  EXTRA_RUN=LOGDIR=`$J/rsg_E_1020/prodH`, ICF=`rsgE_page/ic/ic_rsg_g25_deep.txt`, RPH0=5.80019e13, JOBNAME=rsgED,rsgEB,rsgEH,rsgEH1,
  AGE=120, SITE, LAYOUT text, PY, VISPY, FFMPEG.
- `rsgE_page/src/`: `reduce.py` (VER 5), `mk_figs.py`, `mk_movies.py`, `build_page.py`, `csplot.py`, `eos_gen.py`, `lomb.py`,
  `add_shell.py` (adds a shell to existing red npz), `page.css`, `README.md` (original Raven RSGPAGE notes, partly outdated), `upd_page.sh` (old Raven driver).
- `rsg_E_1020/pageH/sync.sh`, `sync_links.py` — stitched run view (prodG < 1675 d + prodH from its 1675 d restart).
- `static_inputs/rsgE_supergran_1101/` — the 10 static supergranulation figures `build_page.section_sg()` embeds (reads them from
  `$J/rsgE_supergran_1101/`).
- The `modes/` figures (11 PNG, static snapshot to 6269 d from rsgmode_1106 + l-nu/nu-r diagrams) are in the staged page dir
  `pages/rsgE_glob/modes/` inside pages_cache.tar; `build_page.section_modes()` only links `modes/<file>`.

## Commands in order

1. `bash $J/quickpage.sh rsgE` -> `pageH/sync.sh` (stitch hst + dump links) -> `mk_figs.py data <out_quick>/img` -> `build_page.py`
   -> `$P/rsgE_glob/index.html`; if new dumps: detached `catchup_job.sh rsgE` (in-alloc).
2. catch-up: `reduce.py $RUN data` (ONLY=n per dump) -> `mk_movies.py` (LIMITS_ONLY pass, SHARD k/NS renders, final encode/append)
   -> rsync `out/publish/movies` -> `$P/rsgE_glob/movies` -> `r3d_par.sh rsgE` (cs3d single + layer) -> `quickpage.sh rsgE nocatch`
   -> `$P/rsgE_glob/READY`.
3. Publish (url explicit): call 1 `index.html` + the changed 2-D parts from READY; call 2 `index.html` + `movies/rsgE_3d.mp4`,
   `rsgE_3d_layers.mp4`, `rsgE_3d_poster.png`, `rsgE_3d_layers_poster.png` (+ `modes/*` only if changed).

## Inputs

- hst: `pageH/rsg.user.hst` columns time, dt, Etot, KE_int, KEr_int, L_bot, L_in, L_top, Mdot_top, M_tot, P_esrc, Picard,
  P_sponge (+ `rsg.hydro.hst` read as HH). hst dt 250 s.
- dumps: `rsg.hydro_w.NNNNN.bin` + `rsg.m1.NNNNN.bin` every 1 d (both needed per dump); `rsg.m1_vet.*.bin` every 5 d (own numbering,
  matched by header time: vet k = hydro 5k) for P_rr = D_rr E in the radiation enthalpy flux. Dumps younger than AGE=120 s skipped.
- logs: `run*.log` in RUN/EXTRA_RUN, `link.*.out` in RUN/LOGDIR (status block: cycles, s/cycle, FATAL/DIVERGED/NC/NaN counts, cost);
  `squeue` for job lines (JOBNAME).
- IC column `ic/ic_rsg_g25_deep.txt` (in pages_cache.tar under rsgE_page/ic), opacity table
  `rsg_E_1020/files/rsg-g25-prod/rosseland_mesakap_4c0e8b64d26a.txt` (cs3d T_eff; in raven_deps_1107.tar), general EOS composition
  X 0.653848 / Y 0.326395 (eos_gen, cs3d).

## What is computed

- reduce.py (per dump, ~16 MB npz): shell means/profiles (rho, T_gas via general-EOS inversion, v_r, KE, luminosity terms: Lr, EOS
  enthalpy flux, radiation enthalpy (E+P_rr)v_r, KE flux, comoving Lrcom, mean-flow Lbulk), photosphere R_ph/T_eff per column,
  equatorial slices, inner/photosphere shells, Mollweide maps (0.96 R shell), histograms, non-finite counts.
- mk_figs.py: time series (KE, energy/mass conservation, dt, Picard, Mdot_top, sponge/esrc power, L_top/L_in, Lomb-Scargle of the
  light curve), profiles, histograms, flux budget (VET dumps only), latest equatorial slice jpg, `table.json` status numbers.
- mk_movies.py: 4 2-D movies (eq slice +-1.5 R0, inner, photosphere, map), 160 dpi, one frame per dump, fps matched to the 3-D
  pace (MATCH3D json), parts <= 14 MB, append-encoded segments (frames/<movie>_segs/parts.json, crf in frames/<movie>.crf).
- cs3d.py: 3-D photosphere (T_eff, `hot`) + v_r wedge, single and layered views, 2000x1520, <= 300 frames.
- build_page.py: one HTML; figures base64-embedded, movies + modes linked as files.

## Caches needed for continuous history

- `rsgE_page/out/data/red_*.npz` (6304 files, 116 GB) + `.figcache_scal.pkl` -> bulk list (`$J/raven_copy_1107/pages_cache_bulk.lst`).
- `vis3d/rsgE_glob/red/r_*.npz` (2948, 22 GB) -> bulk list.
- In pages_cache.tar: `rsgE_page/out/frames/{*.crf,*.stamp,*.txt,*.json,moll_lut npy,<movie>_segs/,<movie>/limits.json}` (movie
  append state; closed parts are never re-encoded unless the limits are redone), `out/img`, `out/publish`, `out_quick`, `ic/`,
  `vis3d/rsgE_glob/{clim.json,clim_layers.json,eos_fe.npy,out/}`, `rsg_E_1020/pageH` (stitched hst + prefix caches; `bin/` holds
  symlinks to Caltech dumps, dangling on Raven), run logs of prodG/prodH, staged `pages/rsgE_glob` (incl. `modes/`).
- Optional PNG frames (29 GB) in the bulk list; without them the first Raven catch-up re-renders all 2-D frames from the npz
  (6304 x 4 frames; ~10-20 min on 24 cores) and may then fully re-encode.

## Machine-specific bits

See `../common_global/README.md`. Page-specific: update.sh `HERE`, `R`, `PY`, `VISPY`, `FFMPEG`, `SITE`/`LAYOUT` text (says Caltech H200),
`JOBNAME`; build_page.py footer and section_sg path; mk_movies.py MATCH3D default; cs3d.py `RUNS['rsgE_glob']` hydro glob, table, publish.
