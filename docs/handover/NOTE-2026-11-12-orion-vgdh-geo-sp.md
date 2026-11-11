# NOTE 2026-11-12 ORION -> viper: vgdh-geo sp (item 2, last row) -- on sp more directions DO lift the finer pulse_u pair:
# CAND 0.45 -> geo3 1.90 / geo4 1.36 (F1 0.56 -> 1.50 / 1.15), unlike cs (geo3/geo4 = CAND there)

Last sp rows of 99e55db19 / 242985932. fork/vgdlay-1110 453f879bd, pulse_u SP_FULL (angular 4 8 16 32 for nx1 64..512), one
angular MeshBlock per level, CAND + vet_gd_geo, rotation off. 4/4 rc 0 each.

| arm (pulse_u full C) | E p (64-128-256, 128-256-512) | F1 p | d at the finest pair (E) |
|---|---|---|---|
| CAND (geo 2, 42 dirs) | 1.14 0.45 | 0.87 0.56 | see 242985932 |
| CAND-geo3 (92) | 0.77 **1.90** | 0.89 **1.50** | 4.10e-3 |
| CAND-geo4 (162) | 1.15 **1.36** | 1.07 1.15 | 3.98e-3 |
| CAND-sw5 | 1.30 1.53 | 1.09 1.86 | 6.91e-3 |

## Reading
On cs the 42 -> 162 directions left pulse_cs R3 at 1.07 1.72 (99e55db19); on sp the same change takes the finer pulse_u
pair from 0.45 to 1.4-1.9, with ~2x smaller absolute differences than sw5 at the finest pair. So on sp part of the cap IS the
angular resolution of the vet_gd set (the sp wedge's ray sweep sees the fixed set at a much coarser angular spacing per cell
than the cs panels). The obvious next sp arm is **CAND-sw5 + geo3** (combines the w-variation fix with the angular one);
I can run it on the same layout (~5-6 h for the 512 level on one rank) on your GO.
