# NOTE 2026-11-26m orion -> Caltech: CUDA build + suite of the combined candidate (next4)

Branch caltech-next4-8651496c3 = release c3f133b44 + dead column code removed + key removal (defaults and removed keys) + thin-flux form c913e1dcc + flux-weight fix 3cf6d231c + star decks and test counts (28 commits; author names still to be unified, test branch only). Please, on ONE idle test node (02-09 or 02-16; never the M35C node or a production node):

1. Build (star_vet + built_in_pgens, CUDA 12.9 + g++ 13, .prov of 8651496c3) and tiny gate (pass 1 = 3.686e-11). About 37 min at 2 builds in parallel as for next3.
2. Full suite (60 tests), the same --only list split as next3. **Known FAIL: order_time and star_smoke (old recorded references; star_smoke's stored L_int(t=0) is the defective F; re-recording is approved by the USER on the final tip). Anything else that FAILs: print it as it comes.** steady_flux is new: expect thick ~1e-15, total about 1.1.
3. Print, against next3 (9b0c24989), all values that differ, per test.

If the Orion CPU quick gate fails I will cancel this NOTE by a follow-up (please check for "NOTE 26n" before starting the suite). Do not start the suite before the tiny gate passes.

-- orion (coordinator)
