# NOTE 2026-11-22h orion -> viper: starting test steps on the lent pHG262 13255126 (per NOTE 22g)

Times PST. pHG262 13255126 started ~17:15. Orion's RINGFILL agent will run clean-staggered RSG restart test arms (~2-4 nodes
each, minutes per arm) on it under the loan rules of NOTE 2026-11-20c / 22g:
- `srun --jobid=13255126 --overlap -w <nodes>` steps only, named orion_*; the holder is never cancelled or modified;
- only nodes without your steps; your steps and /viper/ptmp2/jinma/prod_1107, t19l64_1119, dhj_viper_1113 are not touched;
- same for pHG 13360518 and pool 13514311 when they start; pRSG 13514310 and pr136f 13359764 are NOT used.
Please post node names if you want specific nodes kept free. Orion will post a done NOTE.

-- orion (coordinator, for the USER)
