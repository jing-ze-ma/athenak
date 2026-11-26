# NOTE 2026-11-27zg orion -> Caltech: scspeed, decks fixed (your 27ze smoke FATAL)

Thanks for stopping on the FATAL. Its cause: the geodesic set has no azimuthal period, so the 0.4 rad
phi wedge of sp.athinput cannot have exact phi-wall maps. Both decks now use the same **product
set GL 8 x 20, unstaggered (160 rays)** (`vet_gd_gl_nmu = 8`, `vet_gd_gl_nphi = 20`,
`vet_gd_gl_stagger = false`). On Cartesian (`vet_sc_set = gd`) its pole is on x1, with no ray in the
layers. The sp phi wedge is now one azimuthal period, 2 pi/20 = 0.31416, so the phi-wall maps are exact.
The updated files are on this branch: `docs/handover/scspeed-1127/cart.athinput` and `sp.athinput`.
Please `git pull` before you run.

Run everything in NOTE 27ze unchanged (steps 0-5, the same overrides, the same report lines),
using the new decks. The vet_gd log line now names the GL set, not geodesic 4.
If step 0 FATALs again, post the lines and stop.

-- orion (SCSPEED, for the coordinator)
