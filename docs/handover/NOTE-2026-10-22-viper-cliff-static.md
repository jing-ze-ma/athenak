# NOTE 2026-10-22 viper -> caltech: ADDENDUM to 8d78b5c5: the bare cliff layer breaks symmetry from round-off on BOTH cs and sp (seed 0, no MLT scaffold, arith faces)

Same static runs (he_star_m1 RSG E, he_seed 0, MLT scaffold off, radiation coupled, f55be371 = arith face opacity): the interior
0.3-0.95 R_ph stays symmetric to round-off, but the cliff layer amplifies round-off into column-local structure within 1-2 d:
| r / R_ph | lateral v_r rms, cs N=32 (1.0 d -> to 7.4 d) | sp wedge 2.5 deg (2.0 d -> to 15 d) |
|---|---|---|
| 0.99  | 2.4e3 -> 1.3-1.6e3 cm/s | 1.9e3 -> 1.4-3.5e3 cm/s |
| 0.995 | 1.3e2 -> 3e2 cm/s       | 2.8e2 -> 2-5e2 cm/s     |
| 1.0   | 42 -> 1-2e2 cm/s        | 63 -> 1.5-4.7e2 cm/s    |
on top of an m=0 infall of the top layers (-1.5e5..-5e5 cm/s; no scaffold), identical in cs and sp to 2-3 digits. corr with cell
geometry is small in both (cs dOmega +0.03..+0.10; sp sin(theta) -0.09..+0.14).
=> The cliff instability is NOT a cubed-sphere effect (the sp grid has it too); it does not need the MLT scaffold deposit (the bare
cliff with arith faces is unstable at this resolution), consistent with the face-flux monotonicity criterion that harmonic /
kirchhoff are built to satisfy. cs only maps the cell-size-dependent growth onto the panel pattern. Implication: sp production decks
with arith faces (HG wedge hg145, AG Car / BSG-type) are exposed as well.
Next on viper (USER GO): the same seed-free static arms with implicit_face_opac = harmonic and = kirchhoff on cs N=32 and the sp wedge:
does the face option alone keep the bare cliff symmetric on both grids? Files: /viper/ptmp2/jinma/photinst_1022 (scr/static.py,
runs/hse32, runs/hse_sp).
