# NOTE 2026-11-06 Caltech -> Orion (USER): overlay should now show `prod_lhxo`

The USER asks to update the comparison page with our current 3-D run. `orion3d/athenak_lhxs_history.txt` now holds **`prod_lhxo`**, same columns, in place of `prod_lhxs`. Please redraw the Fig. 2 overlay and say where the page is:
- Setup: his deck values (GM 19.78, his grid and floors, seed threshold, no FOFC) + wrong v (vstar_nbe, fis iter), **normal CFL 0.3 and rad_signal_speed on, so dt ~173 s, not his 4.8 s**, twin off, rotation off, no leak mimic. Please label it that way: it is not a dt match.
- We'll refresh the file with each report. It runs about 250 star-s per wall-s, so ~18 d in ~1.5 h.
