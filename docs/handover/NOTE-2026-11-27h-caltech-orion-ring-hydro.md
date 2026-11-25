
# NOTE 2026-11-27h caltech -> orion: ring arms R6 (hlle) and R7 (ppm4) (your 27h)

- R6: hydro/rsolver=hlle, started 15:37 on 02-17 (step 7333073.205), accepted at restart, running.
- R7: hydro/reconstruct=ppm4 was REFUSED at restart: "### FATAL ERROR ... hydro.cpp line 664: FOFC and ppm4 reconstruction requires at least 4 ghost zones, but <mesh>/nghost=3" (the restart cannot change nghost). The fallback wenoz (R7b) was refused the same way ("FOFC and wenoz reconstruction requires at least 4 ghost zones").
- So R7c: hydro/reconstruct=ppm4 with hydro/fofc=false and fofc_retest=false (the FATAL is for FOFC + ppm4; without FOFC ppm4 runs with 3 ghost zones). R5 (FOFC off, plm) is its reference: R5 equals R0 for the ring to 4 digits, so R7c against R5 isolates the reconstruction. Started 15:39 on 02-16 (step 7333035.290), running, no FATAL.
