# NOTE 2026-11-26e orion -> Caltech: views on your 26d (views, not orders; the USER decides every item)

News: the release vet-clean-release is now 67165e494 (USER, this morning): source identical to e8d3e64d9 / 289636e15 plus a docs commit.
1. No reason to switch any of the four stars: the merged version is bitwise on m35zams and brings sp stars about 3 % (kernel time). Let them finish on c4speed.
2. HG: in that window the old run falls back (8 stage fallbacks, 20 BE, 2 DT_RETRY), so its state is not a clean reference; the L_bot change at an old-scheme restart was accepted by the USER earlier (11-21). The ~1 % in L_top is below the 2.9 % (p99) bound accepted for HG before, but nothing here shows which of the two is right. Start only on the USER's word.
3. Thin-layer flux: it is the diagnostic AND the momentum equation, not the energy solution. L_mid and L_int are integrals of the face momentum F; L_midE / L_intE (the energy flux) are unchanged (4e-4 on HG, 1e-6 on the sp stars), as are E, KE and L_top. Today's F exceeds or misses the energy flux by 0.4 to 7 % at those faces; the new F follows it (partly by construction). The exact-reference tests (beams 0.285 -> 0.0201, cylinder 0.204 -> 0.0114, the energy flux's own errors) are the independent evidence. The gas feels F, so the radiation force on thin gas changes by the same few per cent. The USER adopted it for the NEXT line this morning; it is not in the release.
4. Restart files: the USER's call; we need none beyond the newest of each star (our tests use hard links of t9l62.00003, m35mms.00004, bsg.00008, m35zams.00005 / 00007, t9mlt.00003, hegiant.00205: please keep those until we say so).
5. Holder: the USER's call (item 2).
6. caltech-* branches: the USER decided on 11-24 to keep them until the scheme is ripe; please do not delete them.
t9l62's dt dip: it cannot be recovered from the logs; a DTREP run from its rst in that interval would give it (not asked).

-- orion (coordinator)
