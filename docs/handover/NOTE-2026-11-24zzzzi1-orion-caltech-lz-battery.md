# NOTE 2026-11-24zzzzi1 orion -> Caltech: bsg L_z leak: single-toggle battery on 02-09 (USER: find the leak)

USER order. Base arm = LZB0e (rel-e 08c48ac10 star_vet 8737e5b0, bsg3d.00004.rst c14605, +100, 4 ranks, its exact args); r RMS 3.04e33. After LZB0f (fix build), run one after another on 02-09, ~95 s each, each with ONE thing changed from LZB0e and everything else identical; stop the battery if one arm FATALs, post it, go to the next:
| arm | change |
|---|---|
| L1 | hydro/fofc=false |
| L2 | hydro/floor_from_ekin=false |
| L3 | hydro/merge_lat_f=1,1,1 (no merge; LZB0e's dt) |
| L4 | rad_vet/clean_conv_max at its value / 10 (name the stored value) |
| L5 | time/cfl_number=0.15, 200 cycles |
| L6 | rad_vet/clean_wb_exch_local=false |
| L7 | sponge off, if the deck has one (name the key) |
| L8 | hydro/vceil=1e30 (ceiling never acts) |
| L9 | rad_vet/floorcap_* off (all floorcap keys named as none / large), if the deck has them |
| L10 | hydro/floor_pay=none if it does not FATAL |
If a key is refused on restart, say so and skip that arm (no reshaping). Per arm: r RMS, max, sum, sum/RMS, dLz_tot, L_top / KE vs LZB0e, CW / FATAL, the arm's dt. Diagnostics only: no setting is being proposed. Reading: the arm(s) that bring r to about 3e31 (round-off) name the term; several reductions = several terms.
In between, LZB0f (zzzzh8) stays first. The final-sha suite (not before 22:10 PST) takes priority over the battery.

-- orion (coordinator)
