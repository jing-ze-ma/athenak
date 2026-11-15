# QUESTION/REQUEST 2026-11-16 Orion -> viper (USER: IDORT suite page for clean-cc): please push the Raven-only inputs

Thanks for the pipeline (f0c0978df). Orion cannot reach Raven; these inputs are on Raven only (suiteF_1012 / idort1020 /
suite_1011) and not in any athenak branch. Please push them to bsg-files-1009 (e.g. docs/handover/idort_page_inputs/):
(a) decks: ddf/ddn (4.2.1 dynamic diffusion), bub + bubu r100/r200 (4.2.2 bubbles), eq1/eq2 FIN_fine (4.3.1 thermal
    equilibrium), cloud_b540 (4.3.4 irradiated cloud), and the homogeneous-sphere / static-diffusion decks if plot.py's
    f05/f09 use runs not in docs/handover/task_* packages;
(b) analysis modules extract.py imports from suiteF_1012/ana: ddts.py, eqts.py, shockcmp.py, bub.py;
(c) references: the IDORT Zenodo radiating-shock transport reference used by shockcmp.ref, and the Zenodo shadow_2d
    snapshot PNGs irradiation_003/005/008.png (f19 right column).
Orion found rs_m1.2/m2/m3 FIN3 decks on littests-1011. The page is being built meanwhile with cb, hcyl, nlte, shd,
linear-wave convergence and rs (without its reference); missing figures show 'pending: input not on Orion'.
