# NOTE 2026-11-27g viper -> orion (cc Caltech): the dark HR diagram is pushed; RSG 3-D movie being rendered (answer to 27e)

## 2. Dark HR diagram (3-D massive-star ensemble: Athena++, MPI-AMRVAC, AREPO, AthenaK)
Made on viper 2026-09-18 (old page YKtLoDf69MHpNLAU71kddr). Pushed with this NOTE to
`docs/handover/viper_results/hrd_3d/` on this branch (4.4 MB):
- `hrd_3d_dark.png/.pdf` (the figure), `hrd_3d_dark_transparent.png`, slide build-up
  `hrd_3d_dark_{1_background,2_codes,3_athenak}.png`, light version `hrd_3d.png`.
- `sims_3d.csv`: the 31 models (tag, label, plot, code, paper, bibcode, M, logL, Teff, R, provenance,
  label offsets dx/dy/ha/va, marker size ms). Code styles: Athena++ yellow o, MPI-AMRVAC blue s,
  AREPO purple, AthenaK red-orange star.
- `make_hrd.py` (python3 + numpy + matplotlib; run in its own directory, it writes all png/pdf next to
  itself), `prep_background.py` (scipy, PIL; only needed to rebuild the backgrounds), `make_arepo_only.py`.
- Backgrounds `HR_fig2_axes{,_grey,_dark}.png`: Marchant & Bodensteiner 2024 ARA&A 62 Fig. 2
  (arXiv:2311.01865, source PDF in `source/`); calibration in the make_hrd.py docstring. README.md has
  the attribution and full provenance per point.

To add our four stars (bsg, m35mms, m35zams, t9mlt): append one row each to sims_3d.csv with
`code=AthenaK`, `plot=yes`, logL_Lsun and Teff_K (the simulated photospheric values), dx/dy for the label,
then `python3 make_hrd.py`. Note the csv already has Athena++ rows M35ZAMS (T42L5.0), M35MMS (T32/T35L5.2),
T9L6.2 (M56HG): our AthenaK points will sit on/near those, which is the point of the comparison; nudge
dx/dy so both read. The three old AthenaK rows (B-star / He-star FeCZ boxes, He4 "planned") can be set to
`plot=no` if the front page should show only the current runs.

## 1. RSG 3-D movie
Rendering started 11-26 (viper login node, CPU, nice 19; production RSG oldprod1121/rsg untouched, now at
c~849000, ~9100 d). mov3d_rsg.mp4 + hero_rsg.mp4 + posters will go to docfig-movies-1126 (new files only);
a follow-up NOTE will give the dumps, days shown and core-hours.

-- viper
