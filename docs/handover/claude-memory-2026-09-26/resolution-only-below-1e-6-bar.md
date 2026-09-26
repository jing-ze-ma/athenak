---
name: resolution-only-below-1e-6-bar
description: dhj grid design rule (user, restated 09-25) -- cells per scale height count ONLY for p > 1e-6 bar; above 1e-6 bar no constraint, let dr grow
metadata:
  type: feedback
---
Every radial-resolution statement for the dhj (cells per H, "3/H", "10/H", SPARC-like / UM-like grid design, e-fold
counts, grid comparisons) considers ONLY the layers below the 1e-6 bar level (p > 1e-6 bar). Above 1e-6 bar there is
no requirement: dr may grow freely (keep the stretch monotonic and smooth).

**Why:** the user's requirement since 08-26 ([[dhj-grid-resolution-design]]: "10 cells per scale height BELOW the
1e-6 bar level"); on 09-25 I designed a "SPARC-like 3/H" grid over the whole domain (22.9 e-folds incl. 1e-9..1e-6 bar,
nx1 72 instead of ~64) and the user had to correct it.
**How to apply:** when computing nx1, e-folds, cells-per-H tables or dt for a grid, integrate/evaluate only over
p > 1e-6 bar; mention the top band only as information, never as a target.
