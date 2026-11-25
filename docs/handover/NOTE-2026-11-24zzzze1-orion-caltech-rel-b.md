# NOTE 2026-11-24zzzze1 orion -> Caltech: GPU gate of batch B (names, layout, tests, keys on top of rel-a)

Branch `caltech-rel-b-c8a0fc099` = rel-a c880e2703 + 30 commits (tests 40; CPU 40/40; src not changed by the last three commits). After the rel-a suite and zzzzd7:
1. CUDA 12.9 + g++ 13, own build dir: built_in and star_vet (PROBLEM = star_vet; he_star_vet is the same source).
2. **The tiny gate deck needs `clean_face_hctau = false` in its radiation block** (B makes it default true; without the key the value is 3.720e-11, with it 3.686e-11): add it to your copy of the gate deck, say so.
3. Full GPU suite, 40 tests incl. star_smoke (its reference is recorded on CPU; the solver-tolerance difference CPU / GPU is not known: report its dev, bound 2e-4, don't change anything) and radwave N80.
4. Old-rst check as for rel-a (8 + 8, old rst from the release or 95813cf85, keys absent).
Priority: after zzzzd7, before the speed arms only if the node is free; production first.

-- orion (coordinator)
