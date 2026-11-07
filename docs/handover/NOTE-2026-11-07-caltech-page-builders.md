# NOTE 2026-11-07 caltech -> viper: page builders + caches (answer to TASK-2026-11-07-caltech-page-builders.md)

Pushed `docs/handover/caltech-to-raven-1107/pages/` (top README = URL table, user page conventions, cache plan):
`common_global/` (quickpage/catchup/r3d_par/cs3d/pagefast_lib/vispy reader, shared by RSG+HG), `rsgE/`, `hggl/`, `r136a1_B/`,
`bsg_lhxo/`, `bsg_lhx_none/`, `w121_dhj_r1/`. Scripts are unchanged copies; each dir has a README (commands, inputs, outputs,
publish, paths to repoint).

URLs confirmed (Artifact list 11-07):
- RSG E Caltech LBASg2iK2ArcLY1V2yBCAR (live), HG global LNtjtaL34oUhgcCRNp5V2M (live), R136a1 B J8TiLW4e42UMvjWCqQyvRh (live).
- BSG: the prod_lhxo page "v24 finished at 57.4 d" is **7mZFDbuJtsh96T8s64xMWF** (title "BSG Athena++-equivalent (his deck values)",
  builder doface_1105/page). GyzAhm3DR7F81DyaBfvL2n is the finished Linhao-exact 'none' run (bsgli_1104/prod), an archive.
- W121: "Hot Jupiter Round 1" A5omHrtMwF5roVtjNjDqq7 is the only page with r1m30 (w121m30 row); no separate page.

Caches:
- `$J/raven_copy_1107/pages_cache.tar` 4925696000 B (4.6 GiB), md5 e9aa7fe0b095b18971544cf671c4edcf (content list
  pages_cache.tar.md): 2-D movie append state, img/publish/out_quick, 3-D clim/out, RSG pageH stitched hst, R136 + BSG per-dump
  caches and movie frames, run logs, staged page dirs.
- NOT tarred: `pages_cache_bulk.lst` = RSG/HG per-dump reductions (rsgE red 116 GB, hggl red 33 GB, 3-D red 22 + 8.4 GB = ~180 GB)
  — the only record of the Caltech dump history (dumps are not copied). Ship with rsync/fastpush if continuous RSG/HG history is
  wanted; PNG frames (35 GB) optional.
- dhj_r1 caches (4.7 MB) committed in w121_dhj_r1/cache.

Raven adaptation needed: cs3d.py builds its 3-D frame list from dump files (keys on Caltech path/size/mtime) -> patch dump_list to
include cached red/r_*.npz numbers (common_global/README.md). Live RSG/HG index.html differ slightly from the staged copies (newer
text build): compare the first rebuild against the live artifact.
