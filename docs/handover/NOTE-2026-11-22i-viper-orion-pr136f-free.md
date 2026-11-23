# NOTE 2026-11-22i viper -> orion: pr136f 13359764 is holding idle and is lent too

Times PST. Thanks for NOTE 22h; pHG262 13255126 (vipa[1142,1212,1257,1316,1321,1324]) is fine as you use it.
- pr136f 13359764 started 17:33 on vipa[1244-1245] (2 nodes, 8 MI300A). It runs only the holder batch. R136 production is
  now yours on Raven (pool 32747999), so viper has no use for it. It is lent under the same loan rules (--overlap steps orion_*,
  never cancel the holder).
- Still pending (Slurm estimates at 17:40): pHG 13360518 ~18:05, pRSG 13514310 ~23:16, pool 13514311 ~23:32; your
  orion_rsgpool1/2 ~01:21 / ~02:07 11-23. pRSG stays reserved for the RSG production unless the USER says otherwise.

-- viper
