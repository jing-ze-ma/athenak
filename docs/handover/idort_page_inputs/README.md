# idort_page_inputs (viper -> Orion, 2026-11-16; answers QUESTION-2026-11-16-orion-viper-idort-page-inputs)

Raven-only inputs of the IDORT suite page (idort_page_pipeline/figsrc/extract.py + plot.py). Text, small npz and PNG
only (2.8 MB). Raven source path of every item is given so extract.py's paths can be remapped.

## decks/ (one directory per run: in.athinput as run, plus whatever metadata the run dir had)
- in.athinput = the input the run used. idort1020 runs (ddf, ddn, ddiff, sdifk, eq*, hcyl, cloud*) were run as
  `athena -i in.athinput` (all FIN3 keys in the deck; bin.txt = CPU binary + md5, sum.txt = rc/fatal/dt_retry/NC line).
- fin3_1012 bub/bubu runs: in.athinput + cmd.txt = the full command line with the FIN3 keys as command-line overrides
  (args.txt, fixkeys.txt as run). Apply cmd.txt's overrides on top of in.athinput.

| dir | Raven source | page function |
| --- | --- | --- |
| ddf_FIN3_d | idort1020/runs/cart/ddf/FIN3_d/p (= p19runs/cart/ddf/FIN_d) | extract.dif, plot.f10 (4.2.1) |
| ddn_FIN3_d | idort1020/runs/cart/ddn/FIN3_d/p | extract.dif, plot.f11 (4.2.1) |
| ddiff_FIN3_d | idort1020/runs/cart/ddiff/FIN3_d/p | extract.dif (profiles) |
| sdifk_FIN3_d | idort1020/runs/cart/sdifk/FIN3_d/p | extract.dif, plot.f09 (static diffusion; not in task_* packages) |
| hcyl_FIN3_n512 | idort1020/runs/cart/hcyl/FIN3_n512/p | extract.hcyl, plot.f05 (512^2 arm; g8x32 deck is in task_1115_cc_idort_head) |
| eq1_FIN3_fine, eq2_FIN3_fine | idort1020/runs/cart/eq{1,2}/FIN3_fine/p | extract.eq, plot.f13/f14 (4.3.1) |
| eq3t_FIN3_t | idort1020/runs/cart/eq3t/FIN3_t/p | extract.eq, plot.f15 (moving equilibrium) |
| cloud_FIN3_r12g2 | idort1020/runs/cart/cloud/FIN3_r12g2/p | extract.cloud arm FIN, plot.f19 |
| cloud_b540 | idort1020/runs/xtra/cloud_b540/p | extract.cloud arm FINb540, plot.f19 (4.3.4) |
| bub_FIN3_r100k5, bub_FIN3_r100k3, bub_FIN3_g8x32k3 | fin3_1012/cart281/bub/<run> | extract.bub (FIN_r100k5 etc.), plot.f12 row 1 |
| bubu_FIN3_r100, bubu_FIN3_r200 | fin3_1012/cart281/bubu/<run> | plot.f12 rows 2-3 (4.2.2 outflow) |

Run name mapping: page names FIN_* = Raven run dirs FIN3_*. plot.f12 reads P/fin3runs/cart/{bub,bubu}/FIN_*; that link
tree does not exist on Raven any more; the runs are fin3_1012/cart281/{bub,bubu}/FIN3_*.

## ana/ (modules extract.py imports that were not shipped before)
| file | Raven source | used by |
| --- | --- | --- |
| ddts.py | suiteF_1012/ana/mq/ddts.py | extract.dif (ana, rd, V, L, BG) |
| eqts.py | suiteF_1012/ana/mq/eqts.py | extract.eq (series) |
| shockcmp.py | suiteF_1012/ana/shockcmp.py | extract.shock (ref, run, jump); REF path -> ref/idort_rad_shock_sol/M%s |
| bub.py + bub_exact_128_2.npz, bub_exact_512_8.npz, bub_sources.json | suiteF_1012/ana | extract.bub (metrics); bub.py loads the npz from its own dir; it imports bin_convert from a Raven path (use the shipped task_1115_cc_idort_head/ana/bin_convert.py) |
| patm_eval.py | suiteF_1012/sp/ana/patm_eval.py | extract (sp atmosphere section; imports spl, shipped in beams-1110/ana/sp) |
Already shipped elsewhere: cb, lcommon, shd, hcyl, s2, bin_convert (task_1115_cc_idort_head/ana), spl, spt, s3_eval
(beams-1110/ana/sp), pipe (suite-pipe-1011).

## ref/
- idort_rad_shock_sol/M1p2, M2, M3: IDORT Zenodo radiating-shock transport reference ({x,rho,p,E,u,rF}_4096.dat +
  radshock.cpp), from suite_1011/ext/idort_mini/x/IDORT_mini/shock_1d/Rad_shock_sol; shockcmp.ref(m) reads
  {x,rho,p,E,u}_4096.dat (plot.f18).
- idort_shadow_2d/irradiation_003/005/008.png: Zenodo shadow_2d/ds_selected snapshots (extract.cloud idort_png,
  plot.f19 right column; extract.cloud reads the colour bar off these PNGs, pixel calibration in its comments).

## cache/
- bubeval_2.json: idort_page/cache/bubeval_2.json (extract.bub 'eval'; bubeval.py output).

## Skipped (too large; regenerable)
- idort_page/cache/bubref_2.npz (25 MB): the retarded periodic/open free-streaming bubble reference; regenerate with
  idort_page_pipeline/figsrc/bubref.py (every=2). extract.bub and plot.f12 load it.
- run outputs (bin/, tab/) of all the runs above: not shipped (rerun on the clean scheme).
