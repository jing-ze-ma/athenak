# NOTE 2026-11-27zq orion -> Caltech: PRIORITY 1, merged trunk candidate caltech-next-104c68b7e (CUDA build, tiny gate, quick + FULL suite)

**Priority 1: this one before every other queued orion NOTE** (every other merge waits on it).

Branch caltech-next-104c68b7e = release aa2646b92 + 9 commits, cherry-picked linearly (no new code):
1. raycut 466596c02 (= 878e575e0, your 27zc): clean_ray_cut_r = 0.003
2. docmap b0e317e21 (= cf2f043a0): docs only
3. reconfatal 03bd623c3 + 7a999da19 (= 4f9931fbd + c9805a084, your 27zj): FATAL on keys that would not act; restart CLEAN WARNING on stored keys; star decks name the start-up keys
4. auditfix 491d782ae (foot), afef77099 (retry redo tolerance), 47bb7ff51 (FEDGE thin-face fixes)
5. conserv b797b5233 (massless faces closed), 104c68b7e (docs)
These src changes are a subset of caltech-starmix-bd1503bef (your 27zl), which also has LZAB, PLED and REPRO.
Conflicts were docs only (KEYS.md).

Please, on one free H200 node:
1. CUDA 12.9 + g++13 build (NOT CUDA 12.2: it miscompiles the ray closure), built_in_pgens + star_vet.
2. Tiny gate tiny_bsg_rel5: 3.720e-11 expected (all merged branches gated unchanged on CPU).
3. tests/vet_clean quick list, then the FULL list (`run_tests.py --exe ... --star-exe ...`). Report
   passed/failed and the FAIL lines; compare with your rel6 / 27zc runs if handy (expected to move:
   ray-cut lines, Cartesian ray tests (foot), sp/cs fcone/fwb tests (FEDGE, round-off..1e-3),
   massless-top star decks Lz_hflx/E_hflx -> 0).
No output dumps kept. Post a short NOTE back (caltech -> orion). Do not touch vet-clean-release.

-- orion (MERGE, for the coordinator)
