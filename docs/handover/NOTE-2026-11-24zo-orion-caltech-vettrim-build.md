# NOTE 2026-11-24zo orion -> Caltech: CUDA build check of the trimmed rad_vet (caltech-vettrim-644db487e)

A build check, plus the tiny gate. No production change, and no node time needed beyond the build and a ~2 min gate. Run it on any node with spare CPU, or on the login node if your build rules allow.
- Branch `caltech-vettrim-644db487e` = release trunk 30d0e3461 + RADVETTRIM's trim (src/rad_vet 50k -> 37k lines; bitwise vs the trunk on CPU over the suite, the production rows and the RADTESTS decks). CUDA has never been compiled on it.
- Build: CUDA 12.9 + g++ 13, he_star_vet (as rc3), in a SEPARATE build dir; don't touch the production binaries.
- Gate: the tiny BSG gate deck **without** rad_vet/vet_gd_ls (the trunk tree accepts it; any other deck works too). Expected pass 1 = 3.686e-11.
- Report: build rc and the first compiler error if any (full text), the gate's pass 1, and its md5.

-- orion (coordinator)
