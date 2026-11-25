# NOTE 2026-11-24zzzzd4 orion -> Caltech: cs pattern: the cause found on Orion, arm L

THINFCS measured on the small star3d_cs: a constant panel-locked lateral radiation force (6.5-7.5e-4 g at every radius, from the first step) at the THICK lateral faces (equiangular cells are not mirror images about the face). It was solved on 11-20 (rad_m1 key `clean_cs_lat_wb`) and removed as unused on 11-24; restored as a key, default off. Test it on the hero deck:
| arm | binary | deck |
|---|---|---|
| L | branch `caltech-latwb-83663792f` (= 53290b5af + the key; own build dir, CUDA 12.9 + g++ 13, he_star_vet, gate 3.686e-11, md5) | B's deck, nlim 600, and `clean_cs_lat_wb = true` in <rad_vet> (a fresh deck must name it) |
Score as E / H: table with m4 / phase / sym_rms / rms_all at c150 / 300 / 450 / 600, hst L_top / KE vs Bp, CW / FATAL / DT_RETRY; keep the c600 dump. **Rule:** it passes if m4(200-400) at c600 falls by at least a factor 3 against Bp (1.391e-4), L_top and KE stay within 1 % of Bp and CW / FATAL / DT_RETRY are 0. A residual radial panel pattern (2.9e-5 g on Orion) is expected: report it, do not score it.
Then 60 cycles restart check from a rst of L (chain 30 + rst + 30 vs 60: hst printed bitwise) if a node is free. After H; production first. Not for production.

-- orion (coordinator)
