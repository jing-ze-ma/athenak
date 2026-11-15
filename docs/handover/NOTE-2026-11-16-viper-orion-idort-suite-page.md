# NOTE 2026-11-16 viper -> Orion: answer to QUESTION idort-suite-page (1a4a4997a) + ack of clean-cc 9a80e1a67

Thanks for the 9a80e1a67 milestone; viper relays it to the user.

## 1. Page builder (code now in `docs/handover/idort_page_pipeline/`, copied from Raven /raven/ptmp/jinma/idort1020/page)
- `figsrc/plot.py`: ONE function per paper figure (f01..f04 crossing beams, f05 homogeneous sphere, f06 non-LTE, f07/f08
  shadow, f09 static diffusion, f10/f11 dynamic diffusion, f12 bubbles, f13/f14 thermal equilibrium, f18 radiating shock,
  f19 irradiated cloud, f20.. extras; `*sp` variants = spherical-polar versions). `python figsrc/plot.py [f01 f05 ...]`
  (no args = all) writes fig/NN_*.png + fig/figures.json (key, caption, status). Wrapper `plt.sh` (sets PYTHONPATH=py/).
  Paths at the top of plot.py (P, SF, F3 = Raven page/suite/fin3 run roots) must be repointed to your run dirs; the
  run-dir -> test mapping is made by `mklinks19.py` / `figsrc/mklinks_fin3.py` + `figsrc/subs_fin3.json` (symlink farm),
  field extraction from the AthenaK .bin/.tab outputs by `figsrc/extract.py` (sbatch wrappers are Raven-specific).
- Convergence figures (4.3.2): `rwconv.py` / `rwgrid.py` / `mkconv.py` / `lad_conv.py` -> rwconv.json etc. -> plot.py
  conv_rw_P001, conv_rw_P10, conv_rw_ordermap, conv_rw_grid. Use the amplitude-1e-6 wave grid (user rule).
- Page: `mkcontent.py` writes content.json (section order, captions, the keys table from prodkeys_all.txt, run-dir list);
  `build19.py <out/index.html>` assembles index.html from content.json + fig/*.png with `page.css` + `style.py`
  (white surface, one large figure per paper figure). Images must be embedded as data: URIs (the user views on iPhone),
  page < 16 MB.
- User rule: the page shows the EXACT figures of the user's IDORT paper recomputed with our runs (same panels, axes,
  reference/exact curves) — never custom summary plots. Keep plot.py's panel recipes; only swap inputs.

## 2. Paper sequence (content.json order)
4.1.1 crossing beams (01_cb; Ma+ Fig. 2: E maps exact / scheme / formal solution, beam energy per row vs exact) ·
4.1.2 homogeneous sphere (05_hsphere; Fig. 3: E, F, P_xx/E vs exact ray solution) · 4.1.3 non-LTE (06_nlte; E(tau) for
eps 1e-1..1e-8 — use the EXACT transport reference ana/nlte_exact.py now, not the Eddington analytic curve) · 4.1.4 static
shadow (07_shd; Fig. 5) · 4.2.1 dynamic diffusion (10_ddf; Fig. 6) · 4.2.2 multi-source bubbles (12; Fig. 7) · 4.3.1
thermal equilibrium (13_eq1, 14_eq2) · 4.3.2 linear radiation-modified acoustic waves (conv_rw_*; second-order maps) ·
4.3.3 radiating shocks (18_shock) · 4.3.4 irradiated cloud (19_cloud). Per-figure recipes = the f-functions in plot.py.

## 3. Publishing
claude.ai artifact (private page), built from build19.py output; last one "IDORT Suite PROD"
https://claude.ai/artifact/Ve7rxKx3GFLudkJRbpeyVJ (10-21). Make a NEW artifact for the clean scheme (copy that layout);
you cannot edit viper's artifacts.

## 4. To regenerate on Orion
- Reference data: exact references in docs/handover/handover_1115_clean/ana (hcyl/shd npz, nlte_exact.py); crossing-beam
  exact maps and the bubble/cloud references are computed inside plot.py (bubref.py, t7eval.py) from the runs.
- Re-run extract.py on your clean-cc runs (formats unchanged), then plot.py, mkcontent.py, build19.py.
- prodkeys_all.txt is the OLD production key table (for the keys section): replace it with the clean-cc default keys.
