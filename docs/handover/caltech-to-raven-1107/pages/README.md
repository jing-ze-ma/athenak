# Caltech page builders for Raven/viper (answer to TASK-2026-11-07-caltech-page-builders.md)

Everything here is an UNCHANGED copy of the builder that produced the live page (backups `*.pre_*`, `*.bak`, one-off
patches and `__pycache__` left out). Directory layout inside each page dir mirrors the Caltech layout relative to
`$J=/resnick/groups/carnegie_poc/jingze` (or `$P=/resnick/home/jingze/ATHENAK/pages` for staging dirs), so the
absolute paths in the scripts map one-to-one.

## URLs (checked with Artifact list 11-07 ~13:30 PST)

| dir | run | page title | URL | live? |
|---|---|---|---|---|
| `rsgE/` (+ `common_global/`) | RSG E, rsg_E_1020/prodH (stitched view pageH = prodG < 1675 d + prodH) | RSG E Caltech | https://claude.ai/artifact/LBASg2iK2ArcLY1V2yBCAR | yes, version id 1794016264-c620 |
| `hggl/` (+ `common_global/`) | He giant global, hgglob_1021/prodG | He Giant Global Caltech | https://claude.ai/artifact/LNtjtaL34oUhgcCRNp5V2M | yes, version id 1794010989-9ef6 |
| `r136a1_B/` | R136a1 deep run B, r136a1_1104/runB/prod | R136a1 deep envelope run (B) | https://claude.ai/artifact/J8TiLW4e42UMvjWCqQyvRh | yes, version id 1794030610-d61a |
| `bsg_lhxo/` | BSG prod_lhxo (doface_1105/prod_lhxo, Linhao deck values, finished 57.4 d 11-07 06:14) | BSG Athena++-equivalent (his deck values) | https://claude.ai/artifact/7mZFDbuJtsh96T8s64xMWF | **yes: this is the live page of prod_lhxo ("v24, finished at 57.4 d")** |
| `bsg_lhx_none/` | BSG Linhao-exact 'none' (bsgli_1104/prod, finished 57.4 d 11-05) | BSG Linhao-exact (production) | https://claude.ai/artifact/GyzAhm3DR7F81DyaBfvL2n | archive: run finished 11-05, last publish 11-06; it does NOT show prod_lhxo |
| `w121_dhj_r1/` | dhj round 1: p1 hydro runs of 5 planets + WASP-121b MHD r1m* (incl. r1m30) | Hot Jupiter Round 1 | https://claude.ai/artifact/A5omHrtMwF5roVtjNjDqq7 | yes (the only W121 r1m30 page; no separate r1m30 page exists) |

**Correction to the TASK table:** GyzAhm3D is the 'none' run (bsgli_1104/prod), not prod_lhxo. The page "v24 finished at
57.4 d" is 7mZFDbuJ (built by `doface_1105/page`, title "BSG Athena++-equivalent (his deck values)"; the page dir was
repointed from the doface/vstar arms to prod_lhxo on 11-06, see make_page.sh header). Both BSG runs are finished, so
both pages are final unless a new BSG run is pointed at the builder.

Page version counters: the Artifact service only shows opaque version ids (above). Last counters noted in the
handovers: RSG ~v125+ (v27 after the 11-02 account switch), HG ~v21+ (new account), R136a1 B v13 on 11-06 (more
publishes after), BSG 7mZF v24. For the "identical rebuild" check compare against `Artifact read`/`list scope=files`
of the live page, not the local staging: the live RSG index.html is 11231745 B vs staged 11231381 B and the live HG is
5704818 B vs staged 5645560 B (last publishes came from a slightly newer quickpage text build than the copy left in
staging), so expect small text/date differences in the status block.

## Caches

- Small caches (< 50 MB) are committed inside the page dirs (`w121_dhj_r1/cache`, `rsgE/static_inputs`, BSG
  `ranges.json`/`cadence.json`).
- `$J/raven_copy_1107/pages_cache.tar` (+ `.md5`): medium caches + run logs + the staged page dirs (see
  `pages_cache.tar.md` there for size/md5 and the content list). Paths relative to `$J`, except entries starting with
  `pages/` which are relative to `/resnick/home/jingze/ATHENAK`.
- `$J/raven_copy_1107/pages_cache_bulk.lst`: the per-dump reductions of the RSG/HG pages (~170 GB), NOT tarred
  (would double 170 GB on a 98 %-full filesystem). They are the only copy of the dump history (the dumps are not
  being copied to Raven), so ship them with rsync/fastpush directly if the pages must keep continuous history.

## User page conventions not in the scripts (from the Caltech memory)

1. No movie/media inlining (USER 10-24 "don't worry about phones, original resolution, frames and frequency"):
   publish `index.html` + media files through the Artifact `files` map at their relative paths. (Exception kept by
   the scripts: the RSG/HG static figures are base64-embedded by `build_page.py img()`; movies, posters and `modes/`
   are separate files.)
2. Per-file artifact cap 16 MB: 2-D movies are split into parts `<name>_<movie>_pNN.mp4` <= 14 MB; one publish
   <= 64 MB, so RSG (and any page > 64 MB) is published in 2 calls (index + changed 2-D parts, then index + 3-D files);
   unchanged files are kept by the service. After a skipped publish, publish all parts.
3. Always pass `url=` explicitly when publishing (the 11-02 account switch left old staging paths bound to dead URLs).
4. 3-D movies (vis3d style): black background, photosphere = true tau_R = 2/3 surface coloured by T_eff
   (RSG `hot`, HG `Blues`), 90x90 deg wedge cut showing v_r (coolwarm), layered view (outer 90x180, 90x135, 90x90;
   deeper shells in dT/<T>), relief ALWAYS true scale (never exaggerate; caption "photosphere at true scale"),
   T_eff smooth (RENDER_VERSION c3), 2000x1520. Videos do NOT loop. 2-D movie fps matched to the 3-D pace.
5. Results as plots with a one-sentence plain-language takeaway above each figure; never raw JSON/tables/key-name
   walls as content (USER 11-02).
6. RSG page priority over HG when both need a refresh (USER 10-26). Refreshes must be fast: quick text rebuild first,
   publish, then movies.
7. Report times as relative ("in ~X h"); reply to the user with one line per page (sim t + version).
8. No page work on the login node (USER 11-05): builds run as CPU steps in a held allocation (srun --overlap) or a
   small CPU job; never in a GPU production's CPU share beyond a few cores (catch-ups were 4 cores, nice 19).
9. After a restart from an EARLIER rst (dump numbers reused), move the page caches keyed by dump number for the
   rewound range aside (out/data red_*.npz, out/frames/*/NNNNN.png, vis3d/<run>/{red,frames,frames_layers}) or the
   page mixes old and new runs (TRAP 11-02, `hggl/hgglob_page/fix_stale_1106.sh` did this).
10. BSG pages: colour limits are fixed constants read off Ma+2026 (since 11-06); compare against Linhao's PLOTS (paper
   figs in img/paper_fig*.png), not his text.
11. R136a1: Gamma / L_rad from the COMOVING flux (lab flux is faint), the advected part goes into the enthalpy L
   (USER 11-05); the extra HG/RSG-style line plots (img/hgstyle_*.png) are additional, the original figures stay.
