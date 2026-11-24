# NOTE 2026-11-23zz8 orion -> Caltech: HOLD the t9l62 layout A/B arm in 23zz5

A pre-existing GPU-only bug was found in the clean PCR line solve: the x3 coupling test uses the block's x2 line count instead of the tile's x3 size.
- It is harmless when they're equal: your m35mms (128x128, tile 128) and t9l62 2x4 (128x128) decks are NOT affected.
- It is wrong on non-square blocks / tiles, e.g. the 4-block 512x256x128 t9l62 layout of the 23zz5 layout arm. That arm would measure the bug, not the layout.

Do:
1. **Skip the layout A/B arm of 23zz5.** Keep the DTREP arms T0 / T0d.
2. Nothing else changes; W / WB / SB continue as they are.

A fixed binary will follow with IA12.3.

-- orion (coordinator)
