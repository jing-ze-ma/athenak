# NOTE 2026-11-27i viper -> orion (cc Caltech): use the dark HR diagram WITHOUT points (correction to viper 27g)

The USER: the dark HR diagram they meant is the background with no dots. It is in
`docs/handover/viper_results/hrd_3d/` on this branch:
- `hrd_3d_dark_1_background.png` / `.pdf` / `_transparent.png`: the Marchant & Bodensteiner 2024 Fig. 2
  regions + tracks on black, axes in white, no simulation points. Same framing and calibration as the
  other hrd_3d figures, so points can be overlaid exactly.

To put only our four stars (bsg, m35mms, m35zams, t9mlt) on it: in `sims_3d.csv` set `plot=no` on the three old
AthenaK rows, append the four new rows with `code=AthenaK`, then in make_hrd.py's directory:

    python3 -c "import make_hrd as m; r=m.load('sims_3d.csv'); m.draw(r, True, 'hrd_ours_dark', dark=True, codes=['AthenaK'])"

(`codes=[...]` limits the drawn points to those codes; `codes=[]` reproduces the empty background.)
Ignore the 27g advice about nudging labels next to the Athena++ points; those points are not drawn this way.

-- viper
