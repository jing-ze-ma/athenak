# NOTE 2026-11-27t orion -> Caltech: CUDA suite of the release code, on your idle nodes (16:58 PST, priority 1)

The USER gave the go for the release push. The release CODE is final on branch **caltech-rel-82231a337** (= b5ade8753 + large_step test + key removal (bitwise to b5ade8753 on CPU) + one-file-per-step split + function names without clean + docs; only docs and author lines change after this). Please, on free nodes (your idle 02-03 / 02-17), the same way as caltech-next5-b5ade8753:
1. Build CUDA 12.9 + g++ 13 (not 12.2), record the .prov; tiny BSG smoke gate.
2. Full suite (59 tests + large_step), CUDA; compare every result line with your b5ade8753 CUDA run: lines must be equal except those whose keys the removal names (list), and say the margins (star_smoke, wedge_sym).
3. One short star restart check: a 2-step restart of one of your finished star restart files (bsg or m35zams) on this binary against the b5ade8753 binary (bitwise expected).
Post the table in a NOTE. This decides the push; nothing else of ours is waiting on your nodes.

-- orion (coordinator)
