# NOTE 2026-11-26t orion -> Caltech: CUDA build + suite of the corrected candidate (next5)

The thin-flux form c913e1dcc and the weight fix 3cf6d231c are OUT: Orion's refinement test of a closed rotating scattering shell at large c dt / dx showed an artifact of the form (radial KE 5.8x too large on the coarse grid, the gas thinned to rho 3.5 of 100; it converges away with grid and step). Your next4 result (58/60) therefore stands for the dropped version only; please keep it, no more work on 8651496c3.

Branch **caltech-next5-b5ade8753** = release c3f133b44 + dead column code removed + key removal and defaults (incl. the kappa deck fix) + star decks, test counts. Nothing else. Please on ONE idle test node (not M35C's, not a production node):
1. Build (star_vet + built_in_pgens, CUDA 12.9 + g++ 13, .prov b5ade8753), tiny gate (3.686e-11).
2. Full suite (59 tests) with the same --only split. **Expected: no FAIL: order_time and star_smoke use the recorded references, which this candidate must reproduce** (star_smoke dev about 1.1e-4 of 2e-4 is the known move of the two budget defaults). Report any FAIL at once.
3. Values against next3 (9b0c24989), all differences, per test.
If the Orion quick gate fails I cancel by a follow-up (NOTE 26u). Do not start the suite before your tiny gate passes.

-- orion (coordinator)
