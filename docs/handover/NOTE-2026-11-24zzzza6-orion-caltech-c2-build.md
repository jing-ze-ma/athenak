# NOTE 2026-11-24zzzza6 orion -> Caltech: build c2 (= 95813cf85 + the ray-hole cure), arms on hold

T9LATE found and cured the t9l62 clean_init stall (thinf in ray holes; 3 passes on the Orion reproducer). USER: start t9l62 with it, and restart bsg / M35ZAMS / m35mms fresh with it. **Build now** branch `caltech-c2-d4ae795` (95813cf85 + cherry-pick of T9LATE's 6211aaa56; src diff 10/-5 lines in rad_vet_clean.cpp; CUDA 12.9 + g++ 13, own dir, tiny gate pass1, md5, .prov). Do not run arms until orion says the gates pass. The running c1 stars (95813cf85) go on untouched until then; you get a separate note for the stop + fresh restart (c2).
When orion releases it: (1) t9l62 fresh, 1 cycle then 100 cycles (decks candidate_1124/IA123_t9l62): clean_init converged, CW 0; (2) preflights of the other three on c2 (100 cycles each, as before).

-- orion (coordinator)
