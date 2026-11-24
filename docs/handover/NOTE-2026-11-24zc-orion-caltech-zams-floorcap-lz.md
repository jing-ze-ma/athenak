# NOTE 2026-11-24zc orion -> Caltech: M35ZAMS fresh Z3 restart: add rad_vet/floorcap_lz=true (L_z exact)

For the M35ZAMS fresh start on the Z3 layout (24zb), add **`rad_vet/floorcap_lz = true`** to the deck. Without it, the floor velocity cap scales v_phi and removes L_z in floor gas, which breaks exact L_z (absolute rule). m35mms: add the same key at its next rst.

Also (no action from you, FYI): the M35ZAMS deck names both GL 7x18 and vet_gd_ls = 10, and the code silently uses GL 7x18. The USER decides which set when back. Keep the deck's sets as they are for this restart (the run has been GL 7x18 all along).

ZAMSFAST's dt finding: late-time dt (2.49 s) is set by floor gas held at floorcap_v = 3e7 near the photosphere; the physical ceiling is ~3.5 s at the inner wall. A dt lever for floor gas waits for the USER. Run ZD0 / ZD1 as planned (they show it).

-- orion (coordinator)
