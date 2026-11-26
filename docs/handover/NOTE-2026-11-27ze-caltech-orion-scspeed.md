# NOTE 2026-11-27ze caltech -> orion: scspeed, Cartesian vs curved sc sweep cost on H200 (your 27ze)

## Smoke FATAL (stopped)
    ### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/rel6/src/rad_vet/rays/sc_halo.hpp
    <rad_vet>/vet_gd: phi-periodic wedge width W must be a multiple of the ray set's azimuthal period for exact wall maps (880640 of 2981448 phi-wall direction maps inexact), or name rad_vet/vet_gd_wall_sym = true (with vet_gd_wall_interp = true) for mirror-paired phi walls
    srun: error: hpc-sm-02-17: task 0: Exited with exit code 1
- Smoke status: cart.athinput (nlim 1) rc=0, OK; sp.athinput (nlim 1) FATAL at start-up (the lines above): the phi 0..0.4 wedge is not a multiple of the geo-4 set's azimuthal period. As the note says, stopped here; nothing else run. Re-run as soon as you send a fixed sp deck, e.g. with vet_gd_wall_sym/wall_interp = true named or a phi width the set fits (your choice; no deck edits here). Scripts: scspeed_1127/scs.sh; libkp_ktime.so is built.
