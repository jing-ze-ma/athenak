# NOTE 2026-10-09 (viper -> Caltech, DeltaAI, Orion): sp theta-periodic seam mass leak + the new thin-atmosphere settings

## 1. Mass leak at the theta-periodic seam (affects the He giant N897 wedge)
Every sp run with `ix2_bc = ox2_bc = periodic` (He giant, AG Car, BSG Ma26 wedges) leaks mass/momentum/energy through
the theta seam once gas crosses it: the angular PLM uses sin-weighted centroids x2v, and the ghost centroids are not the
periodic images of the interior ones, so F2(js) != F2(je+1). Measured (AG Car A 480x8x8 seeded): dM/M ~3-4e-10 per
cycle, energy residual 1.8-2.8e-4 L_in dt per cycle; it is a persistent seam source (scales ~ dtheta^2 at higher
resolution). Fix: `<hydro>/sp_x2_periodic_image = true` (and `<mhd>/...`), default false = bitwise; with it mass and
energy close to round-off (1e-14, 1e-8). Code: branch rc-1009 (commit 36f18f1d on the fork once pushed; until then
patch in viper /viper/ptmp2/jinma/mass_1009/fix_sp_x2_periodic_image.patch). Also: on sp, `hydro/reconstruct` is
ignored in the angular sweeps (always PLM).
Caltech He giant: the running N897 chain has this leak. Suggestion (user decides): at a link boundary, restart with a
binary that has the key and `hydro/sp_x2_periodic_image = true` in the restart's -i input (the key must exist in the
input file, command-line keys absent from the input are FATAL on fresh starts — on restarts check first); smoke 10
cycles first; the change is confined to the seam cells.

## 2. Thin-atmosphere radiation (for any he_star_m1 run with vet_col/vet_gd)
`vet_source_noesrc = true` + `vet_gd_twin = true` (both needed together) fix the formal-solution over-emission and the
vet_gd ray effect; on AG Car A (3-D, Raven) L = 1.00 above R_ph with no inward cells, 53 % cheaper than blend alone.
For the He giant only `vet_gd_twin` matters (its MLT scaffold is already zero, so noesrc is a no-op); it runs central,
and central->blend restarts diverge (earlier NOTE), so no change is suggested mid-run without a viper test.
