# NOTE 2026-11-23zc orion -> Caltech: use the idle t9l62 node now: build IA12.1 candidate, tiny BSG gate, short M35ZAMS 128² GL smoke

The USER asked to use the idle t9l62 node (holder 7333072, hpc-sm-02-09) now.
- This is the IA12.1 CANDIDATE: the same binary that is in its final gates on our side tonight (~21:30 PST).
- The t9l62 restart on it follows only after the table passes (separate NOTE).
- bsg and m35mms keep running on IA11. Don't touch them.

## 1. Build
- Branch `caltech-ia121-21e176c67` on GitHub jing-ze-ma/athenak (= Orion clean-stag2-interim-e 21e176c67).
- CUDA 12.9 + g++ 13, same CMake line as IA11. NOT CUDA 12.2.1 + g++ 11.5: that miscompiles on H200.

## 2. Tiny BSG gate (as before)
- With the gate deck unchanged (it names no new keys, so every IA12.1 addition is off), pass 1 should equal the IA11 / IA12 value 3.686e-11.
- If it differs, stop and post it.

## 3. Short M35ZAMS 128² smoke on that node (4 H200)
- Deck: `docs/handover/caltech_ia121/IA121_m35zams128_GL.athinput` (md5 in MD5SUMS there).
  - Paper S22 grid: 384 x 128 x 128, GL 7x18 staggered rays, rotation off, IA12.1 keys, m35mms outputs.
  - The IC and tables are your fresh_prod_1120 paths (/resnick/groups/carnegie_poc/jingze/fresh_prod_1120/m35zams/ic_m35zams_ideal.txt and .../tables/). Check they exist.
- Layout: the deck's 2x2 MeshBlocks of 384x64x64 on 4 GPUs (tile 64). We gated this layout bitwise to 1 block.
- Run from t=0, `-t 01:30:00` (or ~1000 cycles), output as in the deck.
- Report:
  - build OK + gate value;
  - s/cycle (cycles 20-end);
  - dt;
  - CLEAN WARNING / FATAL / stage fallbacks / dt_retry counts;
  - Ltop/Lin and Lbot/Lin every ~100 cycles;
  - GPU memory per GPU;
  - the first 3 CLEAN WARNING context blocks if any.
- Keep its rst; it may become the production start.

## 4. Then
When our IA12.1 table passes (~21:30 PST), we post the t9l62 restart NOTE (from t9l62.00012.rst, every key named). Keep the node for it.

-- orion (coordinator)
