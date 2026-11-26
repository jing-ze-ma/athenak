# NOTE 2026-11-27u orion -> Caltech: the rel82 smoke gate fails on a removed key (expected), how to gate (17:06 PST)

Likely cause of 'GATE FAIL rel82' with an empty pass1: the key removal in this release turns removed keys into rules. `clean_face_hctau` is now fixed at true (src/rad_vet/rad_vet_params.cpp:144: "+true"); a deck that names it with another value stops at start-up with a FATAL naming the key (by design). Your tiny gate deck sets clean_face_hctau = false, so the run stops before pass 1. Please check the log for that FATAL. If so:
- run the gate deck WITHOUT that line on BOTH binaries (b5ade8753 and rel82): the two pass-1 values must be equal (the key removal is bitwise to b5ade8753 for the same physics); the old 4.020e-06 reference belongs to hctau = false and is not comparable.
- any other key in your decks that the removal touched gives the same FATAL with its name: drop the line if its value equals the fixed one, else tell me.
Then go on with the suite (the suite's own decks were updated with the removal).

-- orion (coordinator)
