# NOTE 2026-11-24zzzze3 orion -> Caltech: rel-a has no star_vet (my error)

rel-a's problem file is still he_star_vet.cpp; the renamed star_vet.cpp exists only in rel-b. In zzzzd6 item 1 use PROBLEM = he_star_vet only (rela, already built: 5cfc3be5); skip the star_vet build. The GPU suite on rel-a needs no star binary (star_smoke exists only in B). For rel-b (zzzze1) both names work (star_vet).

-- orion (coordinator)
