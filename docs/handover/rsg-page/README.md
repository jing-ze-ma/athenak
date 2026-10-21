# RSG cubed-sphere page pipeline (Raven RSGPAGE 10-19; "RSG g25 Production" page)

upd_page.sh -> reduce.py (per-dump reduction, panel-major cs blocks: nmb/6 blocks per panel, gnomonic equiangular, handles any
2x2 / 1x1 blocks-per-panel layout from the block logical positions) -> mk_figs.py (img/) -> build_page.py (single-file
publish/index.html, figures inlined). Movies are NOT part of this pipeline (Raven page has none yet).
Run elsewhere with INLINE=1 (or NO_SBATCH=1) and set the env: RUN (dir with rsg.user.hst, run.log, link.*.out), BINDIR
(rsg.{hydro_w,m1}.NNNNN.bin), OUT (work dir), ICF (IC column, e.g. docs/handover/rsg-E-1020/ic_rsg_g25_deep.txt),
OLDRUN (optional old-seed run for the KE comparison; unset/empty skips it), VISPY (AthenaK vis/python), PY (python with numpy +
matplotlib), SITE, LAYOUT (text for the page). Edit the footer string in build_page.py (Raven paths) for your site.

## T_gas (USER: pages show gas temperature)
tgas.py inverts the code's own general-EOS table dump. Make the dump on your site with any 4de8168d CPU or GPU binary built
WITHOUT a user pgen (shock_tube is built in): `athena -i eos_dump_rsgE.athinput` (1 cell, writes eos_dump_rsgE.txt at setup;
output-only, no production deck change). Raven reference md5 of that file: see the NOTE. Then `python tgas.py eos_dump_rsgE.txt`
prints 4 round-trip checks (T -> eint -> T to 0.1 K).
