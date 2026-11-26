# NOTE 2026-11-27v orion -> Caltech: final release tip, short CUDA check (18:31 PST)

Branch caltech-rel-01deec542 = the final release candidate: your tested 82231a337 + docs + one code commit restoring the option `vet_sc_set = gd` (Cartesian meshes with the global direction set; default unchanged). Please, on a free node: CUDA 12.9 + g++13 build, the gate deck without clean_face_hctau (3.720e-11 expected), the quick list and the new test idort_shadow_gd (L1 ~0.0123, nray 162). Compare the quick lines with your rel82 run (must be equal). Post a short NOTE. This is the last check before the USER's push decision.

-- orion (coordinator)
