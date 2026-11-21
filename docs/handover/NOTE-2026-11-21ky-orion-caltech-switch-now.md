# NOTE 2026-11-21ky orion -> caltech: CW rising on all four runs — apply NOTE kw NOW on the current holders (don't wait for 17:46)

PROD 15:55: t9l62 CW 4->12, t9r CW 0->3 with dt 126 -> 25, m35mms 4->5, f1k 0->2 — the stage-stall class the gas Newton switch cures
(V214 D_c). Please build B_p1 (branch caltech-prod-bf4f917d4, CUDA 12.9 + g++13, BSG smoke gate) now and restart all four runs from
their newest rst on the CURRENT holders with the NOTE kw keys (kp keys + clean_gas_newton_switch=4 + clean_gas_newton_switch_min=6
named; export VGD_XFE=0; t9r cfl 0.6 + tfloor 1000 K, others cfl 0.9; no xpat, no damp). The successors then just continue.
Keep the current rst backed up (hard link) before restarting. Post the first-50-cycle tables.

-- orion (for the USER)
