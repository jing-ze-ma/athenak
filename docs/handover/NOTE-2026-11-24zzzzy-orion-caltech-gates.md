# NOTE 2026-11-24zzzzy orion -> Caltech: release gates on idle nodes (Orion is full; USER)

Branch `caltech-starname-4a26282a6` = vet-inputs (decks and tests moved) + vet-keytrim-decks + vet-rmrun + the he_* -> star_* rename (files, keys, PROBLEM). Never run on a GPU and never built with CUDA; its CPU gates pass on Orion (quick 16/16, full 37/37).
1. **CUDA build** (12.9 + g++ 13, HOPPER90) in its own dir with `-DPROBLEM=star_vet` AND once with `-DPROBLEM=he_star_vet` (must map to star_vet with one message); plus the built_in pgens build for the suite.
2. **Suite on a GPU** (tests/vet_clean/run_tests.py, the full list, 1-2 H200 of an idle node): PASS/FAIL list, md5.
3. **star3d_sp and star3d_cs templates** (now inputs/stars/): 3 cycles on 1 GPU each, CW / FATAL 0; and the old-key test: one of your existing star decks (he_* keys) starts with the NOTE lines and runs 3 cycles.
Idle nodes only; after the candidate timings (24zzzzu) and the sdirk2 fix-2 build. If a build fails, post the first error and go on with the others. Nothing on production decks.

-- orion (coordinator)
